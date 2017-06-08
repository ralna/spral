#include <stdlib.h>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <cublas_v2.h>   /* CUBLAS public header file  */

typedef size_t ptr_t;

// basic node matrix data
template< typename ELEMENT_TYPE >
struct node_data {
  ELEMENT_TYPE* ptr_v; // matrix elements
  int ld; // leading dimension
  int nrows;
  int ncols;
};

// tile data needed for tile presolve step
template< typename ELEMENT_TYPE >
struct tile_presolve_data {
  ELEMENT_TYPE* ptr_diag; // diagonal tile elements' array pointer
  ELEMENT_TYPE* ptr_offd; // off-diagonal tile elements' array pointer
  int ldd; // leading dimension of the diagonal tile elements' array
  int ldo; // leading dimension of the off-diagonal tile elements' array
  int nrows;
  int ncols;
};

// input data for cu_multinode_dgemm_n and cu_multinode_solve_n/t
template< typename ELEMENT_TYPE >
struct node_solve_data {
  // pointers are used by cu_multinode_dgemm_n
  ELEMENT_TYPE* ptr_a;
  ELEMENT_TYPE* ptr_b;
  ELEMENT_TYPE* ptr_u;
  ELEMENT_TYPE* ptr_v;
  // leading dimensions
  int lda;
  int ldb;
  int ldu;
  int ldv;
  // sizes
  int nrows;
  int ncols;
  int nrhs;
  // this CUDA block data offset
  int offb;
  // array offsets are used by cu_node_solve_n/t
  long off_a;
  int off_b;
  int off_u;
  int off_v;
};

///////////////////////////////////
namespace spral { namespace ssids {
///////////////////////////////////

// data for inverting the node matrix' diagonal block
template< typename ELEMENT_TYPE >
struct l_inv_data {
  ELEMENT_TYPE* ptr_l; // elements of the diagonal block
  ELEMENT_TYPE* ptr_i; // elements of its inverse
  // leading dimensions
  int ldl;
  int ldi;
  // size
  int n;
  // CUDA block's local number
  int block;
  // # CUDA blocks for this node matrix
  int nb;
};

extern __shared__ char SharedMemory[];

/*
 *
 Auxiliary kernels used on solve phase.
 *
 */

// scales indexed rows of a matrix
template< typename ELEMENT_TYPE >
__global__ void
cu_scale( 
  const int nrows, const int ncols, 
  ELEMENT_TYPE *const a, const int lda, 
  ELEMENT_TYPE *const s, 
  int *const ind
){
  for ( int y = threadIdx.y + blockIdx.y*blockDim.y; y < ncols; 
        y += blockDim.y*gridDim.y )
    for ( int x = threadIdx.x + blockIdx.x*blockDim.x; x < nrows; 
          x += blockDim.x*gridDim.x )
      a[ind[x] - 1 + y*lda] *= s[x];
}

// gathers D**(-1) from nodes' data into a coniguous array
template< typename ELEMENT_TYPE >
__global__ void
cu_gather_diag( const int n, ELEMENT_TYPE *const src, ELEMENT_TYPE *const dst, long *const ind )
{
  for ( int x = threadIdx.x + blockIdx.x*blockDim.x; x < n;
        x += gridDim.x*blockDim.x ) {
    const long i = ind[x] - 1;
    dst[2*x] = src[i];
    dst[2*x + 1] = src[i + 1];
  }
  // We add a dummy 0.0 value to left of d_11 that allows simpler application
  // later: memory is offset on call so this is safe!
  if ( blockIdx.x == 0 && threadIdx.x == 0 )
    dst[-1] = 0.0;
}

// gathers backward solve rhs/solution vector; rhs part is multiplied by D**(-1)
// v = D * u
template< typename ELEMENT_TYPE >
__global__ void
cu_gather_dx( 
  const int nrows, const int ncols,
  ELEMENT_TYPE *const d,
  ELEMENT_TYPE *const u, const int ldu,
  ELEMENT_TYPE *const v, const int ldv,
  const int *const drow,
  const int *const urow
){
  for ( int x = threadIdx.x + blockIdx.x*blockDim.x; x < nrows; 
        x += blockDim.x*gridDim.x ) {
    const int dr = drow[x];
    const int ur = urow[x];
    if ( dr ) {
      const int i = 2*(dr - 1) + 1;
      for ( int y = threadIdx.y + blockIdx.y*blockDim.y; y < ncols; 
            y += blockDim.y*gridDim.y ) {
        ELEMENT_TYPE s = d[i]*u[ur - 1 + ldu*y];
        if ( d[i - 1] )
          s += d[i - 1]*u[urow[x - 1] - 1 + ldu*y];
        if ( d[i + 1] )
          s += d[i + 1]*u[urow[x + 1] - 1 + ldu*y];
        v[x + ldv*y] = s;
      }
    }
    else {
      for ( int y = threadIdx.y + blockIdx.y*blockDim.y; y < ncols; 
            y += blockDim.y*gridDim.y )
        v[x + ldv*y] = u[ur - 1 + ldu*y];
    }
  }
}

// applies D**(-1) for a partial solve
template< typename ELEMENT_TYPE >
__global__ void
cu_apply_d( 
  const int nrows, const int ncols,
  ELEMENT_TYPE *const d,
  ELEMENT_TYPE *const u, const int ldu,
  ELEMENT_TYPE *const v, const int ldv,
  const int *const urow
){
  for ( int x = threadIdx.x + blockIdx.x*blockDim.x; x < nrows; 
        x += blockDim.x*gridDim.x ) {
    const int ur = urow[x];
    const int i = 2*x + 1;
    for ( int y = threadIdx.y + blockIdx.y*blockDim.y; y < ncols; 
          y += blockDim.y*gridDim.y ) {
      ELEMENT_TYPE s = d[i]*u[ur - 1 + ldu*y];
      if ( d[i - 1] )
        s += d[i - 1]*u[urow[x - 1] - 1 + ldu*y];
      if ( d[i + 1] )
        s += d[i + 1]*u[urow[x + 1] - 1 + ldu*y];
      v[x + ldv*y] = s;
    }
  }
}

// gathers rows of a matrix
template< typename ELEMENT_TYPE >
__global__ void
cu_gather(
  const int nrows,
  const int ncols,
  ELEMENT_TYPE *const src, const int lds,
  ELEMENT_TYPE *const dst, const int ldd,
  int *const ind
){
  for ( int x = threadIdx.x + blockIdx.x*blockDim.x; x < nrows;
        x += gridDim.x*blockDim.x ) {
    const int i = ind[x] - 1;
    for ( int y = threadIdx.y + blockIdx.y*blockDim.y; y < ncols;
          y += gridDim.y*blockDim.y )
      dst[x + y*ldd] = src[i + y*lds];
  }
}

// scatters rows of the sum of two matrices
template< typename ELEMENT_TYPE >
__global__ void
cu_scatter_sum(
  const int nrows,
  const int ncols, 
  ELEMENT_TYPE *const u, const int ldu,
  ELEMENT_TYPE *const v, const int ldv,
  ELEMENT_TYPE *const dst, const int ldd,
  const int *const ind
){
  for ( int x = threadIdx.x + blockIdx.x*blockDim.x; x < nrows;
        x += gridDim.x*blockDim.x ) {
    const int i = ind[x] - 1;
    for ( int y = threadIdx.y + blockIdx.y*blockDim.y; y < ncols;
          y += gridDim.y*blockDim.y )
      dst[i + y*ldd] = u[x + y*ldu] + v[x + y*ldv];
  }
}

// scatters rows of a matrix
template< typename ELEMENT_TYPE >
__global__ void
cu_scatter(
  const int nrows,
  const int ncols,
  ELEMENT_TYPE *const src, const int lds,
  ELEMENT_TYPE *const dst, const int ldd,
  const int *const ind
){
  for ( int x = threadIdx.x + blockIdx.x*blockDim.x; x < nrows;
        x += gridDim.x*blockDim.x ) {
    const int i = ind[x] - 1;
    for ( int y = threadIdx.y + blockIdx.y*blockDim.y; y < ncols;
          y += gridDim.y*blockDim.y )
      dst[i + y*ldd] = src[x + y*lds];
  }
}

/*
 *
 Pre-solve kernels below do some post-processing with the L-factor
 in order to accelerate the subsequent solve(s) with it.
 Each node matrix

   | L_d |
   | L_o |,

 where L_d is the square dense diagonal block and L_o is compactly stored
 sparse off-diagonal block, is replaced by

   |    L_d**(-1)   |
   | -L_o L_d**(-1) |

 so that the triangular solves are performed by gemm-like CUDA kernels
 rather than trsv-like kernels.
 
 Below I always stands for an identity of a size determined by the context.

 In order to invert L_d's, we split them into square tiles of size <tile_size>
 and perform block-backward solve for the systems

   L_d**T X = I,
   
 by operations on tile-rows of the extended matrix L_e**T = | L_d**T I | or,
 equivalently, on tile-columns of L_e.

 At the first pre-solve step, we simultaneously apply the inverses of the
 diagonal tiles to the respective tile-columns of all extended matrices L_e.
 Note that this replaces the diagonal tiles T_d of the upper half of L_e with
 identity tiles and the diagonal tiles of the lower half of L_e with T_d**(-1).

 Next, we split L_e into same-tile-width groups (cf. presolve_lwork and
 presolve_first in factor_gpu.f90) and, within each group, simultaneously apply
 to each L_e tile-column operations that eliminate the off-diagonal tiles of 
 the upper half of L_e, thus transforming it into the identity and the lower
 half into L_d**(-1).

 Once L_d**(-1) are computed, we simultaneously apply them to respective L_o.
 *
 */

/*
 *
 The next two kernels are involved in 'tile pre-solve step', whereby the 
 inverses of the diagonal tiles are computed and applied to the respective
 tile-columns of L_e by solving simultaneously triangular systems
 
   T_d**T X**T = | T_o**T I |,
   
 where each T_d is a diagonal tile of L_d and T_o is formed by respective
 sub-diagonal ones. The solve is performed by operations on columns of
 the extended matrix T_e, where T_e**T = | T_d**T T_o**T I |.

 Each CUDA block is provided with T_d and either a tile from T_o or
 an identity tile and computes the respective tile of X.
 *
 */

// This kernel forms identities for the tile presolve step.
// blockDim.x = blockDim.y = tile_size
template< typename ELEMENT_TYPE >
__global__ void
cu_init_presolve( node_data< ELEMENT_TYPE > *data )
{
  data += blockIdx.x;
  const int n = data->ncols;  
  const int ld = data->ld;
  ELEMENT_TYPE *const v = data->ptr_v;
  for ( int x = threadIdx.x; x < n; x += blockDim.x )
    for ( int y = threadIdx.y; y < n; y += blockDim.y )
      v[x + y*ld] = (x == y) ? 1.0 : 0.0;
}

// This kernel performs solve T_d**T Xi = Yi**T, where Yi is either
// a tile from T_o or an identity tile, by operations on columns of
// T_ei, where T_ei**T = | T_d**T Yi**T |
// blockDim.x = tile_size*tile_size
// shared memory size 2*tile_size*tile_size
template< typename ELEMENT_TYPE, bool NONUNIT_DIAG > 
__global__ void
cu_tile_presolve( 
  const int tile_size,
  tile_presolve_data< ELEMENT_TYPE > *data
){
  const int LDW = 2*tile_size;

  volatile ELEMENT_TYPE *const work = (volatile ELEMENT_TYPE*)SharedMemory;

  data += blockIdx.x;

  const int ldo   = data->ldo;
  const int nrows = data->nrows;
  const int ncols = data->ncols;

  ELEMENT_TYPE *const offd = data->ptr_offd;
  
  const int x = threadIdx.x % tile_size;
  const int y = threadIdx.x / tile_size;
  
  {
    const int ldd = data->ldd;
    work[x + LDW*y] = 
      (x < ncols && y <= x) ? data->ptr_diag[x + ldd*y] : 0.0;
  }
  work[x + tile_size + LDW*y] = 
    (x < nrows && y < ncols) ? offd[x + ldo*y] : 0.0;
  __syncthreads();
  
  if ( NONUNIT_DIAG ) {
    if ( x != y && y < ncols )
      work[x + LDW*y] /= work[y + LDW*y];
    work[x + tile_size + LDW*y] /= work[y + LDW*y];
    __syncthreads();
  }
  
  for ( int m = ncols - 1; m > 0; m-- ) {
    if ( y < m )
      work[tile_size + x + LDW*y] -= 
        work[m + LDW*y] * work[tile_size + x + LDW*m];
    __syncthreads();
  }
  
  if ( x < nrows && y < ncols )
    offd[x + ldo*y] = work[x + tile_size + LDW*y];
}

// This kernel replaces the diagonal block of each node matrix
// with its inverse
template< typename ELEMENT_TYPE >
__global__ void
cu_multi_l_inv_copy( l_inv_data< ELEMENT_TYPE > *data, const int tile_size )
{
  data += blockIdx.x;
  const int n     = data->n;
  const int ldl   = data->ldl;
  const int ldi   = data->ldi;
  const int block = data->block;
  const int nb    = data->nb;
  
  for ( int i = threadIdx.x + blockDim.x*(threadIdx.y + block*blockDim.y);
        i < n*n;
        i += nb*blockDim.x*blockDim.y ) {
    const int x = i % n;
    const int y = i / n;
    data->ptr_l[x + ldl*y] = data->ptr_i[x + ldi*y];
  }
}

// This kernel computes (w - step + 1)-th tile-row of the inverse of L_d,
// where w is L_d's width in tiles, by eliminating (w - step + 1)-th
// tile row of the extended matrix L_e. This is done by subtracting from 
// the first (w - step) tile columns of L_e the same columns of
// the product of its (step)-th tile column by its (w - step + 1)-th
// tile row, since the upper half of L_e has identity tiles on the
// diagonal after the tile pre-solve step. Only non-zero tiles are
// computed and subtracted.
template< typename ELEMENT_TYPE >
__global__ void
cu_multi_l_inv( l_inv_data< ELEMENT_TYPE > *data, const int tile_size, const int step )
{
  data += blockIdx.x;

  const int tid = threadIdx.x + blockDim.x*threadIdx.y;
  ELEMENT_TYPE s[16];

  __shared__ volatile ELEMENT_TYPE as[128], bs[128];
  __shared__ ELEMENT_TYPE *volatile ptr_l;
  __shared__ ELEMENT_TYPE *volatile ptr_i;
  __shared__ volatile int n, m, k, ldl, ldi, block, bx, by;

  if ( tid == 0 ) {
    n = data->n;
    const int tiles = (n + (tile_size - 1))/tile_size;
    m = (tiles - step)*tile_size;
    n -= m;
    k = min(n, tile_size);
    ldl = data->ldl;
    ldi = data->ldi;
    ptr_l = data->ptr_l;
    ptr_i = data->ptr_i;
    block = data->block;
    bx = (n - 1)/32 + 1;
    by = block / bx;
    bx = block % bx;
  }
  __syncthreads();
  if ( m < 1 || by > (m + 31)/32 )
    return;

  for ( int i = 0; i < 16; i++ )
    s[i] = 0;

  for ( int i = 0; i < k; i += 4 ) {
    {
      int x, y;
      x = tid % 16;
      y = tid / 16;
      if ( 32*bx + x < n && i + y < k )
        as[x + 32*y] = ptr_i[m + 32*bx + x + ldi*(m + i + y)];
      else
        as[x + 32*y] = 0.0;
      x += 16;
      if ( 32*bx + x < n && i + y < k )
        as[x + 32*y] = ptr_i[m + 32*bx + x + ldi*(m + i + y)];
      else
        as[x + 32*y] = 0.0;
      x = tid % 4;
      y = tid / 4;
      if ( i + x < k && 32*by + y < m )
        bs[x + 4*y] = ptr_l[m + i + x + ldl*(32*by + y)];
      else
        bs[x + 4*y] = 0.0;
      y += 16;
      if ( i + x < k && 32*by + y < m )
        bs[x + 4*y] = ptr_l[m + i + x + ldl*(32*by + y)];
      else
        bs[x + 4*y] = 0.0;
    }
    __syncthreads();
    
    for ( int jy = 0; jy < 4; jy++ ) {
      for ( int iy = 0; iy < 4; iy++ ) {
        s[jy*4    ] += as[threadIdx.x + 32*iy     ]*bs[4*threadIdx.y + 32*jy + iy];
        s[jy*4 + 1] += as[threadIdx.x + 32*iy +  8]*bs[4*threadIdx.y + 32*jy + iy];
        s[jy*4 + 2] += as[threadIdx.x + 32*iy + 16]*bs[4*threadIdx.y + 32*jy + iy];
        s[jy*4 + 3] += as[threadIdx.x + 32*iy + 24]*bs[4*threadIdx.y + 32*jy + iy];
      }
    }
    __syncthreads();
  }

  for ( int iy = 0; iy < 4; iy++ )
    for ( int ix = 0; ix < 4; ix++ ) {
      const int x = threadIdx.x + (ix + bx*4)*8;
      const int y = threadIdx.y + (iy + by*4)*8;
      if ( y < m && x < n )
        ptr_i[m + x + y*ldi] -= s[ix + iy*4];
    }
}

// This kernel simultaneously computes the products of several pairs of
// non-transposed matrices multiplied by a scalar alpha. It is used for
// computing -L_o L_d**(-1) on the last pre-solve step.
template< typename ELEMENT_TYPE >
__global__ void
cu_multinode_dgemm_n( 
  node_solve_data< ELEMENT_TYPE > *data, 
  const ELEMENT_TYPE alpha,
  const int off )
{
  data += blockIdx.x;

  const int tid = threadIdx.x + blockDim.x*threadIdx.y;
  ELEMENT_TYPE s[16];

  __shared__ volatile ELEMENT_TYPE as[128], bs[128];
  __shared__ volatile int n, m, k, lda, ldb, ldu, offb, bx, by;

  if ( tid == 0 ) {
    n = data->nrows;
    k = data->ncols;
    m = data->nrhs;
    lda = data->lda;
    ldb = data->ldb;
    ldu = data->ldu;
    offb = data->offb;
    bx = (n - 1)/32 + 1;
    by = (off + blockIdx.x - offb) / bx;
    bx = (off + blockIdx.x - offb) % bx;
  }
  __syncthreads();
  if ( by > (m + 31)/32 )
    return;

  for ( int i = 0; i < 16; i++ )
    s[i] = 0;

  for ( int i = 0; i < k; i += 4 ) {
    {
      int x, y;
      x = tid % 16;
      y = tid / 16;
      if ( 32*bx + x < n && i + y < k )
        as[x + 32*y] = data->ptr_a[32*bx + x + lda*(i + y)];
      else
        as[x + 32*y] = 0.0;
      x += 16;
      if ( 32*bx + x < n && i + y < k )
        as[x + 32*y] = data->ptr_a[32*bx + x + lda*(i + y)];
      else
        as[x + 32*y] = 0.0;
      x = tid % 4;
      y = tid / 4;
      if ( i + x < k && 32*by + y < m )
        bs[x + 4*y] = data->ptr_b[i + x + ldb*(32*by + y)];
      else
        bs[x + 4*y] = 0.0;
      y += 16;
      if ( i + x < k && 32*by + y < m )
        bs[x + 4*y] = data->ptr_b[i + x + ldb*(32*by + y)];
      else
        bs[x + 4*y] = 0.0;
    }
    __syncthreads();
    
    for ( int jy = 0; jy < 4; jy++ ) {
      for ( int iy = 0; iy < 4; iy++ ) {
        s[jy*4]     += as[threadIdx.x + 32*iy     ]*bs[4*threadIdx.y + 32*jy + iy];
        s[jy*4 + 1] += as[threadIdx.x + 32*iy + 8 ]*bs[4*threadIdx.y + 32*jy + iy];
        s[jy*4 + 2] += as[threadIdx.x + 32*iy + 16]*bs[4*threadIdx.y + 32*jy + iy];
        s[jy*4 + 3] += as[threadIdx.x + 32*iy + 24]*bs[4*threadIdx.y + 32*jy + iy];
      }
    }
    __syncthreads();
  }

  for ( int iy = 0; iy < 4; iy++ )
    for ( int ix = 0; ix < 4; ix++ ) {
      const int x = threadIdx.x + (ix + bx*4)*8;
      const int y = threadIdx.y + (iy + by*4)*8;
      if ( x < n && y < m )
        data->ptr_u[x + y*ldu] = alpha*s[ix + iy*4];
    }
}

/*
 *
 This kernel simultaneously computes the products of several pairs of
 non-transposed matrices. It is used for applying pre-processed node
 matrices

   |    L_d**(-1)   |
   | -L_o L_d**(-1) |

 to corresponding gathered portions of the rhs matrix during the forward
 solve for a sufficiently large number of rhs.
 
 Each CUDA block computes a 64x8 tile of one product using 16x2 threads,
 each one computing 4x4 elements of the product.
 
 The upper and lower parts of the product, corresponding to the respective 
 parts of the preprocessed node matrix, are stored separately, as they
 are subsequently used in a different manner by the multi-frontal forward
 solve.
 *
 */

template< typename ELEMENT_TYPE >
__global__ void
cu_multinode_solve_n_16x2( 
  const int m,
  double *const a, double *const b, double *const u, double *const v,
  node_solve_data< ELEMENT_TYPE > *data,
  const int off
)
{
  data += blockIdx.x;

  const int tid = threadIdx.x + blockDim.x*threadIdx.y;
  ELEMENT_TYPE s[16];

  __shared__ volatile ELEMENT_TYPE as[256], bs[32];
  __shared__ volatile int n, k, lda, ldb, ldu, ldv, offb, bx, by;
  __shared__ volatile int off_a, off_b, off_u, off_v;

  if ( tid == 0 ) {
    n = data->nrows;
    k = data->ncols;
    off_a = data->off_a;
    off_b = data->off_b;
    off_u = data->off_u;
    off_v = data->off_v;
    lda = data->lda;
    ldb = data->ldb;
    ldu = data->ldu;
    ldv = data->ldv;
    offb = data->offb;
    bx = off + blockIdx.x - offb;
    by = blockIdx.y;
  }
  __syncthreads();
  if ( by > (m + 7)/8 )
    return;

  for ( int i = 0; i < 16; i++ )
    s[i] = 0;

  for ( int i = 0; i < k; i += 4 ) {
    {
      int x, y;
      x = tid % 8;
      y = tid / 8;
      for ( int j = 0; j < 8; j++, x += 8 ) {
        if ( 64*bx + x < n && i + y < k )
          as[x + 64*y] = a[off_a + 64*bx + x + lda*(i + y)];
        else
          as[x + 64*y] = 0.0;
      }
      x = tid % 4;
      y = tid / 4;
      if ( i + x < k && 8*by + y < m )
        bs[x + 4*y] = b[off_b + i + x + ldb*(8*by + y)];
      else
        bs[x + 4*y] = 0.0;
    }
    __syncthreads();
    
    for ( int jy = 0; jy < 4; jy++ ) {
      for ( int iy = 0; iy < 4; iy++ ) {
        s[jy*4]     += as[threadIdx.x + 64*iy     ]*bs[4*(threadIdx.y + 2*jy) + iy];
        s[jy*4 + 1] += as[threadIdx.x + 64*iy + 16]*bs[4*(threadIdx.y + 2*jy) + iy];
        s[jy*4 + 2] += as[threadIdx.x + 64*iy + 32]*bs[4*(threadIdx.y + 2*jy) + iy];
        s[jy*4 + 3] += as[threadIdx.x + 64*iy + 48]*bs[4*(threadIdx.y + 2*jy) + iy];
      }
    }
    __syncthreads();
  }

  for ( int iy = 0; iy < 4; iy++ )
    for ( int ix = 0; ix < 4; ix++ ) {
      const int x = threadIdx.x + (ix + bx*4)*16;
      const int y = threadIdx.y + (iy + by*4)*2;
      if ( y < m )
        if ( x < k )
          u[off_u + x + y*ldu] = s[ix + iy*4];
        else if ( x < n )
          v[off_v + x - k + y*ldv] = s[ix + iy*4];
    }
}

/*
 *
 This kernel simultaneously computes the products of several matrix-vector pairs.
 It is used for applying pre-processed node matrices

   |    L_d**(-1)   |
   | -L_o L_d**(-1) |

 to corresponding gathered portions of the rhs matrix during the forward
 solve for a small number of rhs.

 The upper and lower parts of the product, corresponding to the respective
 parts of the preprocessed node matrix, are stored separately.
 *
 */

template< typename ELEMENT_TYPE, unsigned int NRHS >
__global__ void
cu_multinode_solve_n_few_64( 
  double *const a, double *const b, double *const u, double *const v,
  node_solve_data< ELEMENT_TYPE > *data,
  const int off
){
  ELEMENT_TYPE s[NRHS];
  int n, k, lda, ldb, ldu, ldv, offb;
  int off_a, off_b, off_u, off_v;
  
  data += blockIdx.x;
  n = data->nrows;
  k = data->ncols;
  off_a = data->off_a;
  off_b = data->off_b;
  off_u = data->off_u;
  off_v = data->off_v;
  lda = data->lda;
  ldb = data->ldb;
  ldu = data->ldu;
  ldv = data->ldv;
  offb = data->offb;

  int x = threadIdx.x + 64*(off + blockIdx.x - offb);
  if ( x >= n )
    return;

  for ( int i = 0; i < NRHS; i++ )
    s[i] = 0.0;
  for ( int y = 0; y < k; y++ ) {
    for ( int i = 0; i < NRHS; i++ )
      s[i] += a[off_a + x + y*lda]*b[off_b + y + ldb*i];
  }

  if ( x < k ) {
    for ( int i = 0; i < NRHS; i++ )
      u[off_u + x + ldu*i] = s[i];
  }
  else {
    for ( int i = 0; i < NRHS; i++ )
      v[off_v + x - k + ldv*i] = s[i];
  }
}

/*
 *
 This kernel simultaneously computes the products of several 
 transposed-matrix-vector pairs. It is used for applying transposes of
 pre-processed node matrices

   |    L_d**(-1)   |
   | -L_o L_d**(-1) |

 during the backward solve for a sufficiently large number of rhs.

 The kernel computes two parts u and v of each product, so that the
 actual product is u + v.
 *
 */

template< typename ELEMENT_TYPE >
__global__ void
cu_multinode_solve_t_mp_2x2x8( 
  int m,
  double *const a, double *const b, double *const u, double *const v,
  node_solve_data< ELEMENT_TYPE > *data,
  const int off
){
  data += blockIdx.x;

  const int tid = threadIdx.x 
                 + threadIdx.y*blockDim.x
                 + threadIdx.z*blockDim.x*blockDim.y;
  ELEMENT_TYPE s[16];

  __shared__ volatile ELEMENT_TYPE ws[512];
  __shared__ volatile int n, k, lda, ldb, ldu, ldv, offb, bx, by;
  __shared__ volatile int off_a, off_b, off_u, off_v;

  if ( tid == 0 ) {
    k = data->nrows;
    n = data->ncols;
    off_a = data->off_a;
    off_b = data->off_b;
    off_u = data->off_u;
    off_v = data->off_v;
    lda = data->lda;
    ldb = data->ldb;
    ldu = data->ldu;
    ldv = data->ldv;
    offb = data->offb;
    bx = off + blockIdx.x - offb;
    by = blockIdx.y;
  }
  __syncthreads();
  if ( by > (m + 7)/8 )
    return;

  for ( int i = 0; i < 16; i++ )
    s[i] = 0;

  for ( int row = 0; row < k; row += 32*gridDim.z ) {
    {
      int x = tid % 4;
      const int y = tid / 4;
      const int i = row + 32*blockIdx.z;
      for ( int j = 0; j < 8; j++, x += 4 ) {
        if ( i + x < k && 8*bx + y < n )
          ws[x + 32*y] = a[off_a + i + x + lda*(8*bx + y)];
        else
          ws[x + 32*y] = 0.0;
        if ( i + x < k && 8*by + y < m )
          ws[256 + x + 32*y] = b[off_b + i + x + ldb*(8*by + y)];
        else
          ws[256 + x + 32*y] = 0.0;
      }
    }
    __syncthreads();

    for ( int j = 0; j < 4; j++ ) {
      for ( int i = 0; i < 4; i++ ) {
        const int l = i + 4*threadIdx.z;
        s[j*4    ] += ws[32*(threadIdx.x    ) + l]*ws[256 + 32*(threadIdx.y + 2*j) + l];
        s[j*4 + 1] += ws[32*(threadIdx.x + 2) + l]*ws[256 + 32*(threadIdx.y + 2*j) + l];
        s[j*4 + 2] += ws[32*(threadIdx.x + 4) + l]*ws[256 + 32*(threadIdx.y + 2*j) + l];
        s[j*4 + 3] += ws[32*(threadIdx.x + 6) + l]*ws[256 + 32*(threadIdx.y + 2*j) + l];
      }
    }
    __syncthreads();
  }

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < 4; i++ )
      ws[threadIdx.x + 2*i + (threadIdx.y + 2*j)*8 + threadIdx.z*64] = 
        s[i + 4*j];
  __syncthreads();
  
  s[0] = 0.0;
  s[1] = 0.0;
  for ( int i = 0; i < 8; i++ ) {
    s[0] += ws[tid + i*64];
    s[1] += ws[tid + i*64 + 32];
  }
  const int x = tid % 8 + bx*8;
  int y = tid / 8 + by*8;
  if ( blockIdx.z ) {
    if ( x < n && y < m )
      v[off_v + x + y*ldv] = s[0];
    y += 4;
    if ( x < n && y < m )
      v[off_v + x + y*ldv] = s[1];
  }
  else {
    if ( x < n && y < m )
      u[off_u + x + y*ldu] = s[0];
    y += 4;
    if ( x < n && y < m )
      u[off_u + x + y*ldu] = s[1];
  }
}

/*
 *
 This kernel simultaneously computes the products of several 
 transposed-matrix-vector pairs. It is used for applying transposes of
 pre-processed node matrices

   |    L_d**(-1)   |
   | -L_o L_d**(-1) |

 during the backward solve for one rhs.

 The kernel computes two parts u and v of each product, so that the
 actual product is u + v.
 *
 */

template< typename ELEMENT_TYPE >
__global__ void
cu_multinode_solve_t_one_8x8( 
  double *const a, double *const b, double *const u, double *const v,
  node_solve_data< ELEMENT_TYPE > *data,
  const int off
){
  ELEMENT_TYPE s;
  __shared__ volatile ELEMENT_TYPE ws[64];
  __shared__ volatile int n, k, lda, offb; //, bx;
  __shared__ volatile int off_a, off_b, off_u, off_v;

  data += blockIdx.x;
  if ( threadIdx.x == 0 && threadIdx.y == 0 ) {
    k = data->nrows;
    n = data->ncols;
    off_a = data->off_a;
    off_b = data->off_b;
    off_u = data->off_u;
    off_v = data->off_v;
    lda = data->lda;
    offb = data->offb;
  }
  __syncthreads();

  const int y = threadIdx.y + 8*(off + blockIdx.x - offb);

  s = 0.0;
  if ( y < n )
    for ( int x = threadIdx.x + 8*blockIdx.y; x < k; x += 16 )
      s += a[off_a + x + y*lda] * b[off_b + x];

  ws[threadIdx.x + 8*threadIdx.y] = s;
  __syncthreads();
  
  if ( threadIdx.x == 0 && y < n ) {
    s = 0.0;
    for ( int x = 0; x < 8; x++ )
      s += ws[x + 8*threadIdx.y];
    if ( blockIdx.y )  {
      v[off_v + y] = s;
    }
    else {
      u[off_u + y] = s;
    }
  }
  __syncthreads();
}

/*
 *
 This kernel simultaneously computes the products of several 
 transposed-matrix-vector pairs. It is used for applying transposes of
 pre-processed node matrices

   |    L_d**(-1)   |
   | -L_o L_d**(-1) |

 during the backward solve for a small number of rhs.

 The kernel computes two parts u and v of each product, so that the
 actual product is u + v.
 *
 */

template< typename ELEMENT_TYPE, unsigned int NRHS >
__global__ void
cu_multinode_solve_t_few_8x8( 
  double *const a, double *const b, double *const u, double *const v,
  node_solve_data< ELEMENT_TYPE > *data,
  const int off
){
  ELEMENT_TYPE s[NRHS];
  __shared__ volatile ELEMENT_TYPE ws[64*NRHS];
  int n, k, lda, ldb, ldu, ldv, offb;
  int off_a, off_b, off_u, off_v;

  data += blockIdx.x;
  k = data->nrows;
  n = data->ncols;
  off_a = data->off_a;
  off_b = data->off_b;
  off_u = data->off_u;
  off_v = data->off_v;
  lda = data->lda;
  ldb = data->ldb;
  ldu = data->ldu;
  ldv = data->ldv;
  offb = data->offb;

  const int y = threadIdx.y + 8*(off + blockIdx.x - offb);

  for ( int i = 0; i < NRHS; i++ )
    s[i] = 0.0;
  if ( y < n )
    for ( int x = threadIdx.x + 8*blockIdx.y; x < k; x += 16 ) {
      for ( int i = 0; i < NRHS; i++ )
        s[i] += a[off_a + x + y*lda] * b[off_b + x + ldb*i];
    }

  for ( int i = 0; i < NRHS; i++ )
    ws[threadIdx.x + 8*threadIdx.y + 64*i] = s[i];
  __syncthreads();
  
  if ( threadIdx.x == 0 && y < n ) {
    for ( int i = 0; i < NRHS; i++ )
      s[i] = 0.0;
    for ( int x = 0; x < 8; x++ ) {
      for ( int i = 0; i < NRHS; i++ )
        s[i] += ws[x + 8*threadIdx.y + 64*i];
    }
    if ( blockIdx.y )  {
      for ( int i = 0; i < NRHS; i++ )
        v[off_v + y + ldv*i] = s[i];
    }
    else {
      for ( int i = 0; i < NRHS; i++ )
        u[off_u + y + ldu*i] = s[i];
    }
  }
  __syncthreads();
}

// C interface for cu_init_presolve
void
cuda_init_presolve( 
   const cudaStream_t stream,
   const int nblocks, 
   node_data< double > *const data 
){
  const dim3 threads(8, 8);
  for ( int i = 0; i < nblocks; i += 65535 ) {
    const int blocks = min(65535, nblocks - i);
    cu_init_presolve< double ><<< blocks, threads, 0, stream >>>( data + i );
  }
}

// C interface for cu_tile_presolve
void
cuda_tile_presolve( 
  const cudaStream_t stream,
  const int nblocks,
  const int tile_size, 
  tile_presolve_data< double > *const data,
  const int nud
){
  const int nt = tile_size*tile_size;
  const int sms = 2*nt*sizeof(double);
  for ( int i = 0; i < nblocks; i += 65535 ) {
    const int blocks = min(65535, nblocks - i);
    const dim3 threads(24,24);
    if ( nud )
      cu_tile_presolve< double, true >
        <<< blocks, nt, sms, stream >>>( tile_size, data + i );
    else
      cu_tile_presolve< double, false >
        <<< blocks, nt, sms, stream >>>( tile_size, data + i );
  }
}

// C interface for cu_multi_l_inv
void
cuda_multi_l_inv( 
  const cudaStream_t stream, 
  const int nblocks,
  l_inv_data< double > *const data,
  const int tile_size,
  const int step
){
  const dim3 threads(8, 8);
  
  for ( int i = 0; i < nblocks; i += 65535 ) {
    const int blocks = min(65535, nblocks - i);
    cu_multi_l_inv< double >
      <<< blocks, threads, 0, stream >>>( data + i, tile_size, step );
  }
}

// C interface for cu_multi_l_inv_copy
void
cuda_multi_l_inv_copy( 
  const cudaStream_t stream, 
  const int nblocks,
  l_inv_data< double > *const data,
  const int tile_size
){
  const dim3 threads(8, 8);
  
  for ( int i = 0; i < nblocks; i += 65535 ) {
    const int blocks = min(65535, nblocks - i);
    cu_multi_l_inv_copy< double >
      <<< blocks, threads, 0, stream >>>( data + i, tile_size );
  }
}

// computes the number of tiles needed for performing the tile pre-solve step
// on a nrows-by-ncols node matrix
int
tile_presolve_ntiles( const int tile_size, const int nrows, const int ncols )
{
  const int ntcols = (ncols + (tile_size - 1))/tile_size;
  return ntcols + (ntcols*(ntcols - 1))/2;
}

// computes the number of CUDA blocks needed to invert the diagonal block
// L_d of size n
int
l_inv_nblocks( const int tile_size, const int n )
{
  const int nt = (n + (tile_size - 1))/tile_size;
  int nb = 1;
  for ( int t = 1; t < nt; t++ ) {
    const int nx = t*tile_size;
    const int ny = (nt - t)*tile_size;
    nb = max(nb, ((nx + 31)/32)*((ny + 31)/32));
  }
  return nb;
}

// fills the input data structure for the tile pre-solve step
int 
tile_presolve_setup( const int tile_size, 
                    const int nrows, const int ncols, 
                    double *const node, const int ldn,
                    double *const invd, const int ldi,
                    tile_presolve_data< double > *tile_data,
                    const int off
)
{
  const int ntcols = (ncols + (tile_size - 1))/tile_size;
  int k = 0;
  
  tile_data += off;
  for ( int j = 0; j < ntcols; j++ ) {
    for ( int i = j + 1; i < ntcols; i++, k++ ) {
      tile_data[k].ptr_diag = node + j*tile_size*nrows + j*tile_size;
      tile_data[k].ptr_offd = node + j*tile_size*nrows + i*tile_size;
      tile_data[k].ldd = ldn; //nrows;
      tile_data[k].ldo = ldn; //nrows;
      tile_data[k].nrows = min(tile_size, ncols - i*tile_size);
      tile_data[k].ncols = min(tile_size, ncols - j*tile_size);
    }
  }
  for ( int j = 0; j < ntcols; j++, k++ ) {
    tile_data[k].ptr_diag = node + j*tile_size*nrows + j*tile_size;
    tile_data[k].ptr_offd = invd + j*tile_size*ncols + j*tile_size;
    tile_data[k].ldd = ldn; //nrows;
    tile_data[k].ldo = ldi; //ncols;
    tile_data[k].nrows = min(tile_size, ncols - j*tile_size);
    tile_data[k].ncols = min(tile_size, ncols - j*tile_size);
  }
  
  return off + k;
}

///////////////////////////////////
}} // namespace spral::ssids
//////////////////////////////////

using namespace spral::ssids;

extern "C" {

// gathers <nrows> rows of a sparse matrix <ncols> columns wide
// into a dense <nrows>-by-<ncols> matrix
void spral_ssids_gather(const cudaStream_t stream, const int nrows, const int ncols,
      double *const src, const int lds, double *const dst, const int ldd, int *const ind) {
  int nt = ((nrows*ncols + 31)/32)*32;
  if ( nt > 1024 )
    nt = 1024;
  int ty = ((ncols + 3)/4)*4;
  if ( ty > nt/4 )
    ty = nt/4;
  const int tx = nt/ty;
  const dim3 threads(tx, ty);
  const int nx = (nrows + (tx - 1))/tx;
  const int ny = (ncols + (ty - 1))/ty;
  const dim3 grid(nx, ny);
  cu_gather< double ><<< grid, threads, 0, stream >>>
    ( nrows, ncols, src, lds, dst, ldd, ind );
}

// gathers nodes' D factor pieces from the level data array into
// a coniguous array of diagonal and offdiagonal values' pairs
void spral_ssids_gather_diag(const cudaStream_t stream, int n, double *const src,
      double *const dst, long *const ind) {
   int nt = ((n + 31)/32)*32;
   if ( nt > 1024 )
      nt = 1024;
   const int nb = (n + (nt - 1))/nt;
   cu_gather_diag< double ><<< nb, nt, 0, stream >>>( n, src, dst, ind );
}

// gathers <nrows> rows from rhs/solution matrix of the backward solve
// into a dense <nrows>-by-<ncols> matrix;
// rhs part is simultaneously multiplied by D
void spral_ssids_gather_dx(const cudaStream_t stream, const int nrows, const int ncols,
      double *const d, double *const u, const int ldu, double *const v, const int ldv, int *const indd, int *const indx) {
  int nt = ((nrows*ncols + 31)/32)*32;
  if ( nt > 1024 )
    nt = 1024;
  int ty = ((ncols + 3)/4)*4;
  if ( ty > nt/4 )
    ty = nt/4;
  const int tx = nt/ty;
  const dim3 threads(tx, ty);
  const int nx = (nrows + (tx - 1))/tx;
  const int ny = (ncols + (ty - 1))/ty;
  const dim3 grid(nx, ny);
  cu_gather_dx< double ><<< grid, threads, 0, stream >>>
    ( nrows, ncols, d, u, ldu, v, ldv, indd, indx );
}

// gathers <nrows> rows from rhs/solution matrix of the backward solve
// multiplies by D and puts into a dense <nrows>-by-<ncols> matrix
void spral_ssids_apply_d(const cudaStream_t stream, const int nrows, const int ncols,
      double *const d, double *const u, const int ldu, double *const v, const int ldv, int *const indx) {
  int nt = ((nrows*ncols + 31)/32)*32;
  if ( nt > 1024 )
    nt = 1024;
  int ty = ((ncols + 3)/4)*4;
  if ( ty > nt/4 )
    ty = nt/4;
  const int tx = nt/ty;
  const dim3 threads(tx, ty);
  const int nx = (nrows + (tx - 1))/tx;
  const int ny = (ncols + (ty - 1))/ty;
  const dim3 grid(nx, ny);
  cu_apply_d< double ><<< grid, threads, 0, stream >>>
    ( nrows, ncols, d, u, ldu, v, ldv, indx );
}

// multiplies several matrices simultaneously
void spral_ssids_multinode_dgemm_n(const cudaStream_t stream, const int nblocks,
      node_solve_data< double > *const data, const double a) {
  const dim3 threads(8, 8);
  
  for ( int i = 0; i < nblocks; i += 65535 ) {
    const int blocks = min(65535, nblocks - i);
    cu_multinode_dgemm_n< double ><<< blocks, threads >>>( data + i, a, i );
  }
}

// performs forward nodal solve for several nodes simultaneously
// using the modified nodal matrices from presolve;
// the upper part of the solution that corresponds to the diagonal block
// of the nodal matrix is put into an array poined by ptr_u
// member of the node_solve_data structure,
// the rest goes into one pointed by ptr_v
void spral_ssids_multinode_solve_n(const cudaStream_t stream, const int nblocks,
      const int nrhs, double *const a, double *const b, double *const u, double *const v,
      node_solve_data< double > *const data) {
  for ( int i = 0; i < nblocks; i += 65535 ) {
    const int blocks = min(65535, nblocks - i);
    if ( nrhs > 14 ) {
      const dim3 threads(16, 2);
      const int ny = (nrhs + 7)/8;
      const dim3 grid(blocks, ny);
      cu_multinode_solve_n_16x2< double ><<< grid, threads, 0, stream >>>
          ( nrhs, a, b, u, v, data + i, i );
    }
    else if ( nrhs == 14 ) {
      cu_multinode_solve_n_few_64< double, 14 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 13 ) {
      cu_multinode_solve_n_few_64< double, 13 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 12 ) {
      cu_multinode_solve_n_few_64< double, 12 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 11 ) {
      cu_multinode_solve_n_few_64< double, 11 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 10 ) {
      cu_multinode_solve_n_few_64< double, 10 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 9 ) {
      cu_multinode_solve_n_few_64< double, 9 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 8 ) {
      cu_multinode_solve_n_few_64< double, 8 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 7 ) {
      cu_multinode_solve_n_few_64< double, 7 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 6 ) {
      cu_multinode_solve_n_few_64< double, 6 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 5 ) {
      cu_multinode_solve_n_few_64< double, 5 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 4 ) {
      cu_multinode_solve_n_few_64< double, 4 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 3 ) {
      cu_multinode_solve_n_few_64< double, 3 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 2 ) {
      cu_multinode_solve_n_few_64< double, 2 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
    else {
      const dim3 threads(64, 2);
      cu_multinode_solve_n_few_64< double, 1 ><<< blocks, 64, 0, stream >>>
         ( a, b, u, v, data + i, i );
    }
  }
}

// performs backward nodal solve for several nodes simultaneously
// using the modified nodal matrices from presolve;
// the solution is the sum of two parts placed into arrays pointed to
// by ptr_u and ptr_v members of node_solve_data structure
void spral_ssids_multinode_solve_t(const cudaStream_t stream, const int nblocks,
      const int nrhs, double *const a, double *const b, double *const u, double *const v,
      node_solve_data< double > *const data) {
  for ( int i = 0; i < nblocks; i += 65535 ) {
    const int blocks = min(65535, nblocks - i);
    if ( nrhs > 14 ) {
      const dim3 threads(2, 2, 8);
      const int ny = (nrhs + 7)/8;
      const dim3 grid(blocks, ny, 2);
      cu_multinode_solve_t_mp_2x2x8< double ><<< grid, threads, 0, stream >>>
        ( nrhs, a, b, u, v, data + i, i );
    }
    else if ( nrhs == 14 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 14 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 13 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 13 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 12 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 12 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 11 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 11 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 10 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 10 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 9 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 9 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 8 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 8 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 7 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 7 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 6 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 6 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 5 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 5 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 4 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 4 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 3 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 3 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else if ( nrhs == 2 ) {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_few_8x8< double, 2 ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
    else {
      const dim3 threads(8, 8);
      const dim3 grid(blocks, 2);
      cu_multinode_solve_t_one_8x8< double ><<< grid, threads, 0, stream >>>
        ( a, b, u, v, data + i, i );
    }
  }
}

// performs indexed scaling of <nrows>-by-<ncols> matrix
void spral_ssids_scale( const int nrows, const int ncols, double *const a, const int lda, double *const s,
      int *const ind ) {
  const dim3 threads(8,8);
  const int nx = min(64, (nrows + 7)/8);
  const int ny = min(8, (ncols + 7)/8);
  const dim3 grid(nx,ny);
  cu_scale< double ><<< grid, threads >>>( nrows, ncols, a, lda, s, ind );
}

// the opposite of spral_ssids_gather
void spral_ssids_scatter(const cudaStream_t stream, const int nrows, const int ncols,
      double *const src, const int lds, double *const dst, const int ldd, int *const ind) {
  int nt = ((nrows*ncols + 31)/32)*32;
  if ( nt > 1024 )
    nt = 1024;
  int ty = ((ncols + 3)/4)*4;
  if ( ty > nt/4 )
    ty = nt/4;
  const int tx = nt/ty;
  const dim3 threads(tx, ty);
  const int nx = (nrows + (tx - 1))/tx;
  const int ny = (ncols + (ty - 1))/ty;
  const dim3 grid(nx, ny);
  cu_scatter< double ><<< grid, threads, 0, stream >>>
    ( nrows, ncols, src, lds, dst, ldd, ind );
}

// scatters the sum of the two backward solution parts
// produced by spral_ssids_multinode_solve_t
void spral_ssids_scatter_sum(const cudaStream_t stream, const int nrows, const int ncols,
      double *const u, const int ldu, double *const v, const int ldv, double *const dst, const int ldd, int *const ind) {
  int nt = ((nrows*ncols + 31)/32)*32;
  if ( nt > 1024 )
    nt = 1024;
  int ty = ((ncols + 3)/4)*4;
  if ( ty > nt/4 )
    ty = nt/4;
  const int tx = nt/ty;
  const dim3 threads(tx, ty);
  const int nx = (nrows + (tx - 1))/tx;
  const int ny = (ncols + (ty - 1))/ty;
  const dim3 grid(nx, ny);
  cu_scatter_sum< double ><<< grid, threads, 0, stream >>>
    ( nrows, ncols, u, ldu, v, ldv, dst, ldd, ind );
}

// simultaneously computes inverses of the diagonal blocks of several
// tiled nodal matrices of the same tile-width
int spral_ssids_multi_Ld_inv(const cudaStream_t stream, const int nnodes,
      node_data< double > *const data, const int tile_size, double *const d_work) {
  int nblocks;
  int nrows, ncols;
  int maxcols;
  int k, n;
  int status;

  l_inv_data< double > *inv_data = 0;
  l_inv_data< double > *d_inv_data = 0;

  double *d_node_ptr = 0;
  double *d_invd_ptr = 0;
  
  status = 0;
  
  for (;;) {

    maxcols = 0;
    for ( n = 0; n < nnodes; n++ ) {
      ncols = data[n].ncols;
      maxcols = max(maxcols, ncols);
    }

    nblocks = 0;
    for ( n = 0; n < nnodes; n++ ) {
      ncols = data[n].ncols;
      nblocks += l_inv_nblocks(tile_size, ncols);
    }

    inv_data = new l_inv_data< double > [nblocks];
    k = 0;
    d_invd_ptr = d_work;
    for ( n = 0; n < nnodes; n++ ) {
      d_node_ptr = data[n].ptr_v;
      ncols = data[n].ncols;
      nrows = data[n].nrows;
      const int nb = l_inv_nblocks(tile_size, ncols);
      for ( int i = 0; i < nb; i++, k++ ) {
        inv_data[k].ptr_l = d_node_ptr;
        inv_data[k].ptr_i = d_invd_ptr;
        inv_data[k].ldl = nrows;
        inv_data[k].ldi = ncols;
        inv_data[k].n = ncols;
        inv_data[k].block = i;
        inv_data[k].nb = nb;
      }
      d_invd_ptr += ncols*ncols;
    }

    n = nblocks*sizeof(l_inv_data< double >);
    status = cudaMalloc((void**) &d_inv_data, n);
    if ( status )
      break;
    status = cudaMemcpyAsync(d_inv_data, inv_data, n, cudaMemcpyHostToDevice,
      stream);
    if ( status )
      break;

    maxcols = (maxcols - 1)/tile_size + 1;

    int step;
    for ( step = 1; step < maxcols; step++ ) {
      cuda_multi_l_inv( stream, nblocks, d_inv_data, tile_size, step );
    }

    cuda_multi_l_inv_copy( stream, nblocks, d_inv_data, tile_size );

    break;
  }

  if ( d_inv_data )
    cudaFree( d_inv_data );
  if ( inv_data )
    delete [] inv_data;

  return status;
}

// initializes the inversion of the previous kernel by inverting
// the diagonal tiles of blocks to be inverted and applying the inverses
// to respective tile columns
int spral_ssids_multi_Ld_inv_init(const cudaStream_t stream, const int nnodes,
      node_data< double > *const data, const int tile_size, const int nud, double *const d_work) {
  int ntiles;
  int nrows, ncols;
  int maxcols;
  int k, n;
  int status;

  node_data< double > *init_data = 0;
  node_data< double > *d_init_data = 0;

  tile_presolve_data< double > *tile_data = 0;
  tile_presolve_data< double > *d_tile_data = 0;

  double *d_node_ptr = 0;
  double *d_invd_ptr = 0;
  
  status = 0;
  
  for (;;) {

    init_data = new node_data< double > [nnodes];
    maxcols = 0;
    d_invd_ptr = d_work;
    for ( n = 0; n < nnodes; n++ ) {
      ncols = data[n].ncols;
      init_data[n].ptr_v = d_invd_ptr;
      init_data[n].ncols = ncols;
      init_data[n].ld = ncols;
      maxcols = max(maxcols, ncols);
      d_invd_ptr += ncols*ncols;
    }
    n = nnodes*sizeof(node_data< double >);
    status = cudaMalloc((void**) &d_init_data, n);
    if ( status )
      break;
    status = cudaMemcpyAsync(d_init_data, init_data, n, cudaMemcpyHostToDevice,
      stream);
    if ( status )
      break;
    cuda_init_presolve( stream, nnodes, d_init_data );

    ntiles = 0;
    for ( n = 0; n < nnodes; n++ ) {
      nrows = data[n].nrows;
      ncols = data[n].ncols;
      ntiles += tile_presolve_ntiles(tile_size, nrows, ncols);
    }

    tile_data = new tile_presolve_data< double > [ntiles];    
    k = 0;
    d_invd_ptr = d_work;
    for ( n = 0; n < nnodes; n++ ) {
      d_node_ptr = data[n].ptr_v;
      nrows = data[n].nrows;
      ncols = data[n].ncols;
      k = tile_presolve_setup(tile_size, nrows, ncols, 
                              d_node_ptr, nrows, d_invd_ptr, ncols,
                              tile_data, k);
      d_invd_ptr += ncols*ncols;
    }

    n = ntiles*sizeof(tile_presolve_data< double >);
    status = cudaMalloc((void**) &d_tile_data, n);
    if ( status )
      break;
    status = cudaMemcpyAsync(d_tile_data, tile_data, n, cudaMemcpyHostToDevice,
      stream);
    if ( status )
      break;
    
    cuda_tile_presolve( stream, ntiles, tile_size, d_tile_data, nud );

    break;
  }

  if ( d_init_data )
    cudaFree( d_init_data );
  if ( d_tile_data )
    cudaFree( d_tile_data );
  if ( init_data )
    delete [] init_data;
  if ( tile_data )
    delete [] tile_data;

  return status;
}

// sets up data for spral_ssids_multinode_dgemm_n
int spral_ssids_multinode_dgemm_setup(const int nrows, const int ncols, const int nrhs, double *const a,
      const int lda, double *const b, const int ldb, double *const u, const int ldu, double *const v, const int ldv,
      node_solve_data< double > *data, const int off) {
  const int nr = (nrows + 31)/32;
  const int nc = (ncols + 31)/32;

  data += off;
  for ( int i = 0; i < nr*nc; i++ ) {
    data[i].ptr_a = a;
    data[i].ptr_b = b;
    data[i].ptr_u = u;
    data[i].ptr_v = v;
    data[i].lda = lda;
    data[i].ldb = ldb;
    data[i].ldu = ldu;
    data[i].ldv = ldv;
    data[i].nrows = nrows;
    data[i].ncols = ncols;
    data[i].nrhs = nrhs;
    data[i].offb = off;
  }
  return off + nr*nc;
}

} // end extern "C"
