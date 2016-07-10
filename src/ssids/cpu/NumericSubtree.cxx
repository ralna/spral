/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 */
#include "NumericSubtree.hxx"

#include <cstdio>

#include "NumericSubtree.hxx"
#include "SmallLeafSubtree.hxx"

using namespace spral::ssids::cpu;

/////////////////////////////////////////////////////////////////////////////
// anonymous namespace
namespace {

typedef double T;
const int BLOCK_SIZE = 16;
const int PAGE_SIZE = 16384;

} /* end of anon namespace */
//////////////////////////////////////////////////////////////////////////

extern "C"
void* spral_ssids_cpu_create_num_subtree_dbl(
      bool posdef,
      void const* symbolic_subtree_ptr,
      struct cpu_node_data<double>* nodes,
      const double *const aval, // Values of A
      const double *const scaling, // Scaling vector (NULL if none)
      void *const alloc,      // Pointer to Fortran allocator structure
      struct cpu_factor_options const* options, // Options in
      struct cpu_factor_stats* stats // Info out
      ) {
   auto const& symbolic_subtree = *static_cast<SymbolicSubtree const*>(symbolic_subtree_ptr);

   if(posdef) {
      try {
         auto* subtree = new NumericSubtree
            <true, BLOCK_SIZE, T, PAGE_SIZE>
            (symbolic_subtree, nodes, aval, scaling, alloc, options, stats);
         if(options->print_level > 9999) {
            printf("Final factors:\n");
            subtree->print();
         }
         return (void*) subtree;
      } catch(NotPosDefError npde) {
         stats->flag = SSIDS_ERROR_NOT_POS_DEF;
         return nullptr;
      }
   } else { /* indef */
      auto* subtree = new NumericSubtree
         <false, BLOCK_SIZE, T, PAGE_SIZE>
         (symbolic_subtree, nodes, aval, scaling, alloc, options, stats);
      if(options->print_level > 9999) {
         printf("Final factors:\n");
         subtree->print();
      }
      return (void*) subtree;
   }
}

extern "C"
void spral_ssids_cpu_destroy_num_subtree_dbl(bool posdef, void* target) {
   if(!target) return;

   if(posdef) {
      auto *subtree = static_cast<NumericSubtree<true, BLOCK_SIZE, T, PAGE_SIZE>*>(target);
      delete subtree;
   } else {
      auto *subtree = static_cast<NumericSubtree<false, BLOCK_SIZE, T, PAGE_SIZE>*>(target);
      delete subtree;
   }
}

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_cpu_subtree_solve_fwd_dbl(
      bool posdef,      // If true, performs A=LL^T, if false do pivoted A=LDL^T
      void* subtree_ptr,// pointer to relevant type of NumericSubtree
      int nrhs,         // number of right-hand sides
      double* x,        // ldx x nrhs array of right-hand sides
      int ldx           // leading dimension of x
      ) {

   // Call method
   if(posdef) { // Converting from runtime to compile time posdef value
      auto &subtree =
         *static_cast<NumericSubtree<true, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_fwd(nrhs, x, ldx);
   } else {
      auto &subtree =
         *static_cast<NumericSubtree<false, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_fwd(nrhs, x, ldx);
   }
}

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_cpu_subtree_solve_diag_dbl(
      bool posdef,      // If true, performs A=LL^T, if false do pivoted A=LDL^T
      void* subtree_ptr,// pointer to relevant type of NumericSubtree
      int nrhs,         // number of right-hand sides
      double* x,        // ldx x nrhs array of right-hand sides
      int ldx           // leading dimension of x
      ) {

   // Call method
   if(posdef) { // Converting from runtime to compile time posdef value
      auto &subtree =
         *static_cast<NumericSubtree<true, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_diag(nrhs, x, ldx);
   } else {
      auto &subtree =
         *static_cast<NumericSubtree<false, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_diag(nrhs, x, ldx);
   }
}

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_cpu_subtree_solve_diag_bwd_dbl(
      bool posdef,      // If true, performs A=LL^T, if false do pivoted A=LDL^T
      void* subtree_ptr,// pointer to relevant type of NumericSubtree
      int nrhs,         // number of right-hand sides
      double* x,        // ldx x nrhs array of right-hand sides
      int ldx           // leading dimension of x
      ) {

   // Call method
   if(posdef) { // Converting from runtime to compile time posdef value
      auto &subtree =
         *static_cast<NumericSubtree<true, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_diag_bwd(nrhs, x, ldx);
   } else {
      auto &subtree =
         *static_cast<NumericSubtree<false, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_diag_bwd(nrhs, x, ldx);
   }
}

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_cpu_subtree_solve_bwd_dbl(
      bool posdef,      // If true, performs A=LL^T, if false do pivoted A=LDL^T
      void* subtree_ptr,// pointer to relevant type of NumericSubtree
      int nrhs,         // number of right-hand sides
      double* x,        // ldx x nrhs array of right-hand sides
      int ldx           // leading dimension of x
      ) {

   // Call method
   if(posdef) { // Converting from runtime to compile time posdef value
      auto &subtree =
         *static_cast<NumericSubtree<true, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_bwd(nrhs, x, ldx);
   } else {
      auto &subtree =
         *static_cast<NumericSubtree<false, BLOCK_SIZE, T, PAGE_SIZE>*>(subtree_ptr);
      subtree.solve_bwd(nrhs, x, ldx);
   }
}

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_cpu_subtree_enquire_dbl(
      bool posdef,      // If true, performs A=LL^T, if false do pivoted A=LDL^T
      void const* subtree_ptr,// pointer to relevant type of NumericSubtree
      int* piv_order,   // pivot order, may be null, only used if indef
      double* d         // diagonal entries, may be null
      ) {

   // Call method
   if(posdef) { // Converting from runtime to compile time posdef value
      auto &subtree =
         *static_cast<NumericSubtree<true, BLOCK_SIZE, T, PAGE_SIZE> const*>(subtree_ptr);
      subtree.enquire(piv_order, d);
   } else {
      auto &subtree =
         *static_cast<NumericSubtree<false, BLOCK_SIZE, T, PAGE_SIZE> const*>(subtree_ptr);
      subtree.enquire(piv_order, d);
   }
}
