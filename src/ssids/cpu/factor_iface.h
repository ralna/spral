#pragma once

namespace spral { namespace ssids { namespace cpu {

template <typename T>
struct cpu_node_data {
   /* Fixed data from analyse */
   const int nrow_expected; // Number of rows without delays
   const int ncol_expected; // Number of cols without delays
   struct cpu_node_data<T> *const first_child; // Pointer to our first child
   struct cpu_node_data<T> *const next_child; // Pointer to parent's next child
   const int *const rlist; // Pointer to row lists
   const bool even; // Indicates which stack we're using (odd or even distance
                    // from root)

   /* Data about A:
    * aval[i] goes into lcol[ amap[i] ] if there are no delays
    */
   int num_a; // Number of entries from A
   const int *const amap; // Map from A to node (length 2*num_a)

   /* Data that changes during factorize */
   int ndelay_in; // Number of delays arising from children
   int ndelay_out; // Number of delays arising to push into parent
   int nelim; // Number of columns succesfully eliminated
   T *lcol; // Pointer to start of factor data
   int *perm; // Pointer to permutation
   T *contrib; // Pointer to contribution block
};

struct cpu_factor_options {
   double small;
   double u;
   int print_level;
   int cpu_small_subtree_threshold;
   int cpu_task_block_size;
};

struct cpu_factor_stats {
   int flag;
   int num_delay;
   int num_neg;
   int num_two;
   int num_zero;
   int maxfront;
   int elim_at_pass[5];
   int elim_at_itr[5];
};

}}} /* namespaces spral::ssids::cpu */

extern "C" {

double *spral_ssids_smalloc_dbl(void *alloc, size_t sz);
int *spral_ssids_smalloc_int(void *alloc, size_t sz);

} /* end extern "C" */
