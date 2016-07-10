#pragma once

#include <cstddef>

namespace spral { namespace ssids { namespace cpu {

class SymbolicNode;

template <typename T>
struct cpu_node_data {
   /* Fixed data from analyse */
   struct cpu_node_data<T>* first_child; // Pointer to our first child
   struct cpu_node_data<T>* next_child; // Pointer to parent's next child
   SymbolicNode const* symb; // Symbolic node associated with this one

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
