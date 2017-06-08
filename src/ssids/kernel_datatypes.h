#define MAX_CUDA_BLOCKS 65535

namespace spral { namespace ssids {

struct multinode_fact_type {
  int nrows;
  int ncols;
  double *lval;
  double *ldval;
  double *dval;
  int offp;
  int ib;
  int jb;
  int done;
  int rght;
  int lbuf;
};

struct cuda_stats {
  int num_two;
  int num_neg;
  int num_zero;
};

} } // end namespace spral::ssids
