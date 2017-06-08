#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>

#include "../cuda/cuda_check.h"
#include "kernel_datatypes.h"

#define HOGG_ASSEMBLE_TX 128  // Block height
#define HOGG_ASSEMBLE_TY 8    // Block width
#define HOGG_ASSEMBLE_NTX 32  // Number of threads x
#define HOGG_ASSEMBLE_NTY 4   // Number of threads y
#define ADD_DELAYS_TX 32
#define ADD_DELAYS_TY 4

namespace spral { namespace ssids {

struct load_nodes_type {
  int nnz;    // Number of entries to map
  int lda;    // Leading dimension of A
  int ldl;    // Leading dimension of L
  double *lcol; // Pointer to non-delay part of L
  int offn;   // Offset into nlist
  long offr;  // Offset into rlist
};

/*
 * Perform assembly according to nlist:
 * lval( nlist(2,i) ) = val( nlist(1,i) )    (in Fortran)
 *
 * Each block handles one node (regardless of size!!!)
 * Note: modified value lval is passed in via pointer in lndata, not as argument
 */
__global__ void
cu_load_nodes(
      const struct load_nodes_type *lndata,
      const int *nlist,
      const double *const aval) {

   lndata += blockIdx.x;
   const int nnz = lndata->nnz;
   const int lda = lndata->lda;
   const int ldl = lndata->ldl;

   nlist += 2*lndata->offn;
   double *const lval = lndata->lcol;
  
   for ( int i = threadIdx.x; i < nnz; i += blockDim.x ) {
      // Note: nlist is 1-indexed, not 0 indexed, so we have to adjust
      const int r = (nlist[2*i+1] - 1) % lda; // row index
      const int c = (nlist[2*i+1] - 1) / lda; // col index
      const int sidx = nlist[2*i+0] - 1; // source index
      lval[r + c*ldl] = aval[sidx];
   }
}

/*
 * Perform assembly according to nlist:
 * lval( nlist(2,i) ) = val( nlist(1,i) )    (in Fortran)
 * with the added twist of needing to perform a scaling at the same time
 *
 * Each block handles one node (regardless of size!!!)
 * Note: modified value lval is passed in via pointer in lndata, not as argument
 */
__global__ void
cu_load_nodes_sc(
      const struct load_nodes_type *lndata,
      const int *nlist,
      const int *rlist,
      const double *const scale,
      const double *const aval) {

   lndata += blockIdx.x;
   const int nnz = lndata->nnz;
   const int lda = lndata->lda;
   const int ldl = lndata->ldl;

   nlist += 2*lndata->offn;
   double *const lval = lndata->lcol;
   rlist += lndata->offr;
  
   for ( int i = threadIdx.x; i < nnz; i += blockDim.x ) {
      // Note: nlist and rlist are 1-indexed, not 0 indexed, so we adjust
      const int r = (nlist[2*i+1] - 1) % lda; // row index
      const int c = (nlist[2*i+1] - 1) / lda; // col index
      const int sidx = nlist[2*i+0] - 1; // source index
      const double rs = scale[rlist[r] - 1]; // row scaling
      const double cs = scale[rlist[c] - 1]; // col scaling
      lval[r + c*ldl] = rs * aval[sidx] * cs;
   }
}

// BLOCK_SIZE = blockDim.x
// maxabs must be initialized to zeros
template< typename ELEMENT_TYPE, unsigned int BLOCK_SIZE >
__global__ void
cu_max_abs( const long n, ELEMENT_TYPE *const u, ELEMENT_TYPE *const maxabs )
{
  __shared__ volatile ELEMENT_TYPE tmax[BLOCK_SIZE];
  
  tmax[threadIdx.x] = 0.0;
  for ( long i = threadIdx.x + blockDim.x*blockIdx.x; i < n; 
        i += blockDim.x*gridDim.x ) {
    const ELEMENT_TYPE v = fabs(u[i]);
    if ( v > tmax[threadIdx.x] )
      tmax[threadIdx.x] = v;
  }
  __syncthreads();
  
  for ( int inc = 1; inc < BLOCK_SIZE; inc *= 2 ) {
    if ( 2*inc*threadIdx.x + inc < BLOCK_SIZE
        && tmax[2*inc*threadIdx.x + inc] > tmax[2*inc*threadIdx.x] )
      tmax[2*inc*threadIdx.x] = tmax[2*inc*threadIdx.x + inc];
    __syncthreads();
  }
  if ( threadIdx.x == 0 && tmax[0] > 0.0 )
    maxabs[blockIdx.x] = tmax[0];
}


/* Following data type describes a single child-parent assembly */
struct assemble_cp_type {
  // Parent data
  int pvoffset; // Offset to start of parent node values
  double *pval; // Pointer to non-delay part of parent L
  int ldp; // Leading dimension of parent

  // Child data
  int cm; // Number of rows in child
  int cn; // Number of columns in child
  int ldc; // Leading dimension of child
  long cvoffset; // Offset to start of child node values
  double *cv; // Pointer to start of child node values

  // Alignment data
  int *rlist_direct; // Pointer to start of child's rlist
  int *ind; // Pointer to start of child's contribution index

  // Sync data
  int sync_offset; // we watch sync[sync_offset]
  int sync_wait_for; // and wait for it to have value >= sync_wait_for
};

/* Following data type describes actions of single CUDA block */
struct assemble_blk_type {
  int cp; // node we're assembling into
  int blk; // block number of that node
};

/* Used to force volatile load of a declared non-volatile variable */
template <typename T_ELEM>
__inline__ __device__ T_ELEM loadVolatile(volatile T_ELEM *const vptr) {
  return *vptr;
}

/* Performs sparse assembly of a m x n child into a parent as dictated by
 * rlist_direct (supplied as part of cpdata).
 *
 * A lookup is performed in blkdata to determine which child-parent assembly
 * is to be performed next, and which block of that assembly this is.
 *
 * next_blk is used to ensure all blocks run in exact desired order.
 * sync[] is used to ensure dependencies are completed in the correct order.
 */
template <unsigned int blk_sz_x, unsigned int blk_sz_y,
          unsigned int ntx, unsigned nty>
void __global__ assemble(
      const struct assemble_blk_type *blkdata, // block mapping
      const struct assemble_cp_type *cpdata, // child-parent data
      const double *const children, // pointer to array containing children
      double *const parents, // pointer to array containing parents
      unsigned int *const next_blk, // gmem location used to determine next block
      unsigned int *const sync // sync[cp] is #blocks completed so far for cp
      ) {
   // Get block number
   __shared__ volatile unsigned int mynext_blk;
   if (threadIdx.x==0 && threadIdx.y==0)
      mynext_blk = atomicAdd(next_blk, 1);
   __syncthreads();

   // Determine global information
   blkdata += mynext_blk;
   cpdata += blkdata->cp;
   const int blk = blkdata->blk;
   const int nx = (cpdata->cm + (blk_sz_x - 1)) / blk_sz_x; // number of blocks high child is
   const int bx = blk % nx; // coordinate of block in x direction
   const int by = blk / nx; // coordinate of block in y direction
   const int ldc = cpdata->ldc;
   const int ldp = cpdata->ldp;

   // Initialize local information
   const int m = min(blk_sz_x, cpdata->cm - bx*blk_sz_x);
   const int n = min(blk_sz_y, cpdata->cn - by*blk_sz_y);
   const double *const src = 
      cpdata->cv + ldc*by*blk_sz_y + bx*blk_sz_x;
   double *const dest = cpdata->pval;
   int *const rows = cpdata->rlist_direct + bx*blk_sz_x;
   int *const cols = cpdata->rlist_direct + by*blk_sz_y;

   // Wait for previous child of this parent to complete
   if (threadIdx.x==0 && threadIdx.y==0) {
     while (loadVolatile(&sync[cpdata->sync_offset]) < cpdata->sync_wait_for) /**/;
   }
   __syncthreads();

   // Perform assembly
   for (int j=0; j<blk_sz_y/nty; j++) {
      if ( threadIdx.y+j*nty < n ) {
         const int col = cols[threadIdx.y+j*nty]-1;
         for (int i=0; i<blk_sz_x/ntx; i++) {
            if ( threadIdx.x+i*ntx < m ) {
               const int row = rows[threadIdx.x+i*ntx]-1;
               dest[row + col*ldp] += 
                  src[threadIdx.x+i*ntx + (threadIdx.y+j*nty)*ldc];
            }
         }
      }
   }

   // Record that we're done
   __syncthreads();
   if (threadIdx.x==0 && threadIdx.y==0) {
      atomicAdd(&sync[blkdata->cp], 1);
   }
}

/* Same for the solve phase, multiple rhs */
template <unsigned int SIZE_X, unsigned int SIZE_Y,
          unsigned int STEP_X, unsigned STEP_Y>
void __global__ cu_assemble4solve(
      const int nrhs,
      const struct assemble_blk_type *blkdata, // block mapping
      const struct assemble_cp_type *cpdata, // child-parent data
      double *const pval,
      const double *const cval,
      unsigned int *const next_blk, // gmem location used to determine next block
      unsigned int *const sync // sync[cp] is #blocks completed so far for cp
      ) {
   // Get block number
   __shared__ volatile unsigned int mynext_blk;

   if ( threadIdx.x == 0 && threadIdx.y == 0 )
      mynext_blk = atomicAdd(next_blk, 1);
   __syncthreads();

   // Determine global information
   blkdata += mynext_blk;
   const int cp = blkdata->cp;
   cpdata += cp;

   // Wait for previous child of this parent to complete
   if ( threadIdx.x == 0 && threadIdx.y == 0 ) {
     while (loadVolatile(&sync[cpdata->sync_offset]) < cpdata->sync_wait_for) /**/;
   }
   __syncthreads();

   const int cm = cpdata->cm;
   if ( cm < 1 || nrhs < 1 ) {
     if ( threadIdx.x == 0 && threadIdx.y == 0 )
        atomicAdd(&sync[cp], 1);
     return;
   }

   const int blk = blkdata->blk;
   const int ldc = cpdata->ldc;
   const int ldp = cpdata->ldp;

   // Initialize local information
   const int m = min(SIZE_X, cm - blk*SIZE_X);
   double *const dest = pval + cpdata->pvoffset;
   const double *const src = cval + cpdata->cvoffset;
   int *const rows = cpdata->rlist_direct;
   int *const ind = cpdata->ind + blk*SIZE_X;

   for ( int i = threadIdx.x; i < m; i += STEP_X ) {
      const int indx = ind[i] - 1;
      const int row = rows[indx] - 1;
      for ( int j = threadIdx.y; j < nrhs; j += STEP_Y )
         dest[row + j*ldp] += src[indx + j*ldc];
   }

   // Record that we're done
   __syncthreads();
   if ( threadIdx.x == 0 && threadIdx.y == 0 )
      atomicAdd(&sync[cp], 1);
}

/* Same for the solve phase, one rhs */
template< unsigned int SIZE_X, unsigned int STEP_X >
void __global__ cu_assemble4solve_onerhs(
      const struct assemble_blk_type *blkdata, // block mapping
      const struct assemble_cp_type *cpdata, // child-parent data
      double *const pval,
      const double *const cval,
      unsigned int *const next_blk, // gmem location used to determine next block
      unsigned int *const sync // sync[cp] is #blocks completed so far for cp
      ) {
   // Get block number
   __shared__ volatile unsigned int mynext_blk;

   if ( threadIdx.x == 0 )
      mynext_blk = atomicAdd(next_blk, 1);
   __syncthreads();

   // Determine global information
   blkdata += mynext_blk;
   const int cp = blkdata->cp;
   cpdata += cp;

   // Wait for previous child of this parent to complete
   if ( threadIdx.x == 0 ) {
     while (loadVolatile(&sync[cpdata->sync_offset]) < cpdata->sync_wait_for) /**/;
   }
   __syncthreads();

   const int cm = cpdata->cm;
   if ( cm < 1 ) {
     if( threadIdx.x == 0 )
        atomicAdd(&sync[cp], 1);
     return;
   }

   const int blk = blkdata->blk;

   // Initialize local information
   const int m = min(SIZE_X, cm - blk*SIZE_X);
   double *const dest = pval + cpdata->pvoffset;
   const double *const src = cval + cpdata->cvoffset;
   int *const rows = cpdata->rlist_direct;
   int *const ind = cpdata->ind + blk*SIZE_X;

   for ( int i = threadIdx.x; i < m; i += STEP_X ) {
      const int indx = ind[i] - 1;
      const int row = rows[indx] - 1;
      dest[row] += src[indx];
   }

   // Record that we're done
   __syncthreads();
   if ( threadIdx.x == 0 )
      atomicAdd(&sync[cp], 1);
}

struct assemble_delay_type {
  int dskip; // Number of rows to skip for delays from later children
  int m; // Number of rows in child to copy
  int n; // Number of cols in child to copy
  int ldd; // Leading dimension of dest (parent)
  int lds; // Leading dimension of src (child)
  double *dval; // Pointer to dest (parent)
  double *sval; // Pointer to src (child)
  long roffset; // Offset to rlist_direct
};

/* Copies delays from child to parent using one block per parent 
 * Note: src and dest pointers both contained in dinfo
 */
void __global__ add_delays(
      struct assemble_delay_type *dinfo, // information on each block
      const int *rlist_direct // children's rows indices in parents
      ) {
   dinfo += blockIdx.x;
   const int dskip = dinfo->dskip; // number of delays
   const int m = dinfo->m; // number of rows
   const int n = dinfo->n; // number of cols
   const int ldd = dinfo->ldd; // leading dimension of dest
   const int lds = dinfo->lds; // leading dimension of src

   double *const dest = dinfo->dval;
   const double *const src = dinfo->sval;
   rlist_direct += dinfo->roffset;

   for ( int y = threadIdx.y; y < n; y += blockDim.y ) {
      for ( int x = threadIdx.x; x < m; x += blockDim.x ) {
         if ( x < n ) {
            dest[x + y*ldd] = src[x + y*lds];
         }
         else {
            const int xt = dskip + rlist_direct[x - n] - 1;
            dest[xt + y*ldd] = src[x + y*lds];
         }
      }
   }
} 

} } // End namespace spral::ssids

/*******************************************************************************
 * Following routines are exported with C binding so can be called from Fortran
 ******************************************************************************/

using namespace spral::ssids;

extern "C" {

/* Invokes the add_delays<<<>>>() kernel */
void spral_ssids_add_delays( const cudaStream_t stream, const int ndblk,
      struct assemble_delay_type *const gpu_dinfo, int *const rlist_direct ) {
   if ( ndblk == 0 ) return; // Nothing to see here
   const dim3 threads(ADD_DELAYS_TX, ADD_DELAYS_TY);
   for ( int i = 0; i < ndblk; i += MAX_CUDA_BLOCKS ) {
      const int nb = min(MAX_CUDA_BLOCKS, ndblk - i);
      add_delays
         <<< nb, threads, 0, stream >>>
         ( gpu_dinfo + i, rlist_direct );
      CudaCheckError();
   }
}

/* Runs the kernel assemble<<<>>>() after setting up memory correctly. */
/* Requires gpu_next_sync[] to be of size >= (1+ncp)*sizeof(unsigned int) */
void spral_ssids_assemble(const cudaStream_t stream, const int nblk, const int blkoffset,
      struct assemble_blk_type *const blkdata, const int ncp,
      struct assemble_cp_type *const cpdata, double *const children,
      double *const parents, unsigned int *const gpu_next_sync) {
   /* Create and initialize synchronization objects using a single call:
      next_blk[1]
      sync[ncp]
    */
   CudaSafeCall(
         cudaMemsetAsync(gpu_next_sync,0,(1+ncp)*sizeof(unsigned int),stream)
         );
   /* Note, that we can only have at most 65535 blocks per dimn.
    * For some problems, nblk can exceed this, so we use more than one launch.
    * As the next block we look at is specified by next_blk this works fine.
    */
   const dim3 threads(HOGG_ASSEMBLE_NTX, HOGG_ASSEMBLE_NTY);
   for (int i=0; i<nblk; i+=MAX_CUDA_BLOCKS) {
      const int blocks = min(MAX_CUDA_BLOCKS, nblk-i);
      assemble
         <HOGG_ASSEMBLE_TX, HOGG_ASSEMBLE_TY,
          HOGG_ASSEMBLE_NTX, HOGG_ASSEMBLE_NTY>
         <<<blocks, threads, 0, stream>>>
         (&blkdata[blkoffset], cpdata, children, parents, &gpu_next_sync[0],
          &gpu_next_sync[1]);
      CudaCheckError();
   }
}

/* Same for solve phase */
extern "C"
void 
spral_ssids_assemble_solve_phase( 
      const cudaStream_t stream, 
      const int nrhs,
      const int nblk,
      struct assemble_blk_type *const blkdata, 
      const int ncp,
      struct assemble_cp_type *const cpdata,
      double *const pval,
      const double *const cval,
      unsigned int *const gpu_next_sync ) {
   cudaMemset(gpu_next_sync, 0, (1 + ncp)*sizeof(unsigned int));
   for ( int i = 0; i < nblk; i += MAX_CUDA_BLOCKS ) {
      const int blocks = min(MAX_CUDA_BLOCKS, nblk - i);
      if ( nrhs > 1 ) {
         const dim3 threads(HOGG_ASSEMBLE_NTX, HOGG_ASSEMBLE_NTY);
         cu_assemble4solve
           <HOGG_ASSEMBLE_TX, HOGG_ASSEMBLE_TY,
            HOGG_ASSEMBLE_NTX, HOGG_ASSEMBLE_NTY>
           <<< blocks, threads, 0, stream >>>
           ( nrhs, blkdata, cpdata, pval, cval, 
              &gpu_next_sync[0], &gpu_next_sync[1] );
         CudaCheckError();
      }
      else {
         const dim3 threads(HOGG_ASSEMBLE_NTX);
         cu_assemble4solve_onerhs
           < HOGG_ASSEMBLE_TX, HOGG_ASSEMBLE_NTX >
           <<< blocks, threads, 0, stream >>>
           ( blkdata, cpdata, pval, cval, 
              &gpu_next_sync[0], &gpu_next_sync[1] );
         CudaCheckError();
      }
   }
}

// Note: modified value lval is passed in via pointer in lndata, not as argument
void spral_ssids_load_nodes( const cudaStream_t stream, const int nblocks,
      const struct load_nodes_type *const lndata, const int *const list,
      const double *const mval ) {
  for ( int i = 0; i < nblocks; i += MAX_CUDA_BLOCKS ) {
    const int nb = min(MAX_CUDA_BLOCKS, nblocks - i);
    cu_load_nodes <<< nb, 128, 0, stream >>> ( lndata + i, list, mval );
    CudaCheckError();
  }
}

// Note: modified value lval is passed in via pointer in lndata, not as argument
void spral_ssids_load_nodes_sc( const cudaStream_t stream, const int nblocks,
      const struct load_nodes_type *const lndata, const int *const list, const int *const rlist,
      const double *const scale, const double *const mval ) {
  for ( int i = 0; i < nblocks; i += MAX_CUDA_BLOCKS ) {
    const int nb = min(MAX_CUDA_BLOCKS, nblocks - i);
    cu_load_nodes_sc <<< nb, 128, 0, stream >>> ( lndata + i, list, rlist, scale, mval );
    CudaCheckError();
  }
}

void spral_ssids_max_abs( const cudaStream_t stream, 
      const int nb, const long n, double *const u, double *const buff, double *const maxabs )
{
  cudaMemsetAsync(buff, 0, nb*sizeof(double), stream);
  cudaStreamSynchronize(stream);
  if ( n > 1024*nb )
    cu_max_abs< double, 256 ><<< nb, 256, 0, stream >>>( n, u, buff );
  else
    cu_max_abs< double, 32 ><<< nb, 32, 0, stream >>>( n, u, buff );
  CudaCheckError();
  cu_max_abs< double, 1024 ><<< 1, 1024, 0, stream >>>( nb, buff, maxabs );
  CudaCheckError();
}


} // end extern "C"
