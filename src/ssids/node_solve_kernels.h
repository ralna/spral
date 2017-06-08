#ifndef NODE_SOLVE_DATA
#define NODE_SOLVE_DATA

typedef size_t ptr_t;

template< typename ELEMENT_TYPE >
struct node_data {
  ELEMENT_TYPE* ptr_v;
  int ld;
  int nrows;
  int ncols;
};

template< typename ELEMENT_TYPE >
struct tile_presolve_data {
  ELEMENT_TYPE* ptr_diag;
  ELEMENT_TYPE* ptr_offd;
  int ldd;
  int ldo;
  int nrows;
  int ncols;
};

template< typename ELEMENT_TYPE >
struct node_solve_data {
  ELEMENT_TYPE* ptr_a;
  ELEMENT_TYPE* ptr_b;
  ELEMENT_TYPE* ptr_u;
  ELEMENT_TYPE* ptr_v;
  int lda;
  int ldb;
  int ldu;
  int ldv;
  int nrows;
  int ncols;
  int nrhs;
  int offb;
  long off_a;
  int off_b;
  int off_u;
  int off_v;
};

//////////////////////////////////////////////
namespace spral { namespace ssids {
//////////////////////////////////////////////

template< typename ELEMENT_TYPE >
struct l_inv_data {
  ELEMENT_TYPE* ptr_l;
  ELEMENT_TYPE* ptr_i;
  int ldl;
  int ldi;
  int n;
  int block;
  int nb;
};

void
cuda_init_presolve( 
   const cudaStream_t stream,
   const int nblocks,
   node_data< double > *const data
);

void
cuda_tile_presolve( 
   const cudaStream_t stream,
   const int nblocks,
   const int tile_size,
   tile_presolve_data< double > *const data,
   const int nud
);

void
cuda_multi_l_inv( 
   const cudaStream_t stream,
   const int nblocks,
   l_inv_data< double > *const data,
   const int tile_size,
   const int step
);

void
cuda_multi_l_inv_copy(
   const cudaStream_t stream, 
   const int nblocks,
   l_inv_data< double > *const data,
   const int tile_size
);

int
tile_presolve_ntiles(const int tile_size, const int nrows, const int ncols);

int
l_inv_nblocks(const int tile_size, const int n);

int
tile_presolve_setup(const int tile_size, 
                    const int nrows, const int ncols, 
                    double *const node, const int ldn,
                    double *const invd, const int ldi,
                    tile_presolve_data< double > *tile_data,
                    const int off);

//////////////////////////////////////////////
}} // namespace spral::ssids
//////////////////////////////////////////////

extern "C"
void
spral_ssids_multinode_dgemm_n( 
   const cudaStream_t stream, 
   const int nblocks,
   node_solve_data< double > *const data,
   const double a
);

extern "C"
void
spral_ssids_multinode_solve_n( 
  const cudaStream_t stream, 
  const int nblocks,
  const int nrhs,
  double *const a, double *const b, double *const u, double *const v,
  node_solve_data< double > *const data 
);

extern "C"
void
spral_ssids_multinode_solve_t(
  const cudaStream_t stream,
  const int nblocks,
  const int nrhs,
  double *const a, double *const b, double *const u, double *const v,
  node_solve_data< double > *const data
);

extern "C"
int
spral_ssids_multi_Ld_inv_init(
  const cudaStream_t stream,
  const int nnodes,
  node_data< double > *const data,
  const int tile_size,
  const int nud,
  double *const d_work 
);

extern "C"
int
spral_ssids_multi_Ld_inv( 
  const cudaStream_t stream, 
  const int nnodes,
  node_data< double > *const data,
  const int tile_size,
  double *const d_work
);

extern "C"
int
spral_ssids_node_solve_setup(
  const int nrows, const int ncols, const int nrhs,
  double *const a, const int lda, double *const b, const int ldb,
  double *const u, const int ldu, double *const v, const int ldv,
  node_solve_data< double > *data, const int off );

#endif

