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
   const cudaStream_t *stream, 
   int nblocks, 
   node_data< double >* data 
);

void
cuda_tile_presolve( 
   const cudaStream_t *stream, 
   int nblocks, 
   int tile_size, 
   tile_presolve_data< double >* data,
   int nud
);

void
cuda_multi_l_inv( 
   const cudaStream_t *stream, 
   int nblocks, 
   l_inv_data< double >* data, 
   int tile_size, 
   int step
);

void
cuda_multi_l_inv_copy( 
   const cudaStream_t *stream, 
   int nblocks, 
   l_inv_data< double >* data, 
   int tile_size
);

int
tile_presolve_ntiles(int tile_size, int nrows, int ncols);

int
l_inv_nblocks(int tile_size, int n);

int
tile_presolve_setup(int tile_size, 
                    int nrows, int ncols, 
                    double* node, int ldn,
                    double* invd, int ldi,
                    tile_presolve_data< double >* tile_data,
                    int off);

//////////////////////////////////////////////
}} // namespace spral::ssids
//////////////////////////////////////////////

extern "C"
void
spral_ssids_multinode_dgemm_n( 
   const cudaStream_t *stream, 
   int nblocks, 
   node_solve_data< double >* data, 
   double a 
);

extern "C"
void
spral_ssids_multinode_solve_n( 
  const cudaStream_t *stream, 
  int nblocks, 
  int nrhs, 
  double* a, double* b, double* u, double* v,
  node_solve_data< double >* data 
);

extern "C"
void
spral_ssids_multinode_solve_t( 
  const cudaStream_t *stream, 
  int nblocks, 
  int nrhs,
  double* a, double* b, double* u, double* v,
  node_solve_data< double >* data 
);

extern "C"
int
spral_ssids_multi_Ld_inv_init( 
  const cudaStream_t *stream, 
  int nnodes,
  node_data< double >* data,
  int tile_size,
  int nud,
  double* d_work 
);

extern "C"
int
spral_ssids_multi_Ld_inv( 
  const cudaStream_t *stream, 
  int nnodes,
  node_data< double >* data,
  int tile_size,
  double* d_work 
);

extern "C"
int
spral_ssids_node_solve_setup(
                  int nrows, int ncols, int nrhs,
                  double* a, int lda, double* b, int ldb,
                  double* u, int ldu, double* v, int ldv,
                  node_solve_data< double >* data,
                  int off );

#endif

