/* Copyright 2014-6 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 */
#include "factor.hxx"

#include <cstdio>

/////////////////////////////////////////////////////////////////////////////
// start of namespace spral::ssids::cpu
namespace spral { namespace ssids { namespace cpu {

/* FIXME: remove post debug */
template<typename T>
void print_node(bool posdef, int m, int n, int nelim, const int *perm, const int *rlist, const T *lcol, const T*d) {
   for(int i=0; i<m; i++) {
      if(i<n) printf("%d%s:", perm[i], (i<nelim)?"X":"D");
      else    printf("%d:", rlist[i-n]);
      for(int j=0; j<n; j++) printf(" %10.2e", lcol[j*m+i]);
      if(!posdef && i<nelim) printf("  d: %10.2e %10.2e\n", d[2*i+0], d[2*i+1]);
      else printf("\n");
   }
}

/* FIXME: remove post debug */
template<typename T>
void print_factors(
      bool posdef,
      int nnodes,
      struct cpu_node_data<T> *const nodes
      ) {
   for(int node=0; node<nnodes; node++) {
      printf("== Node %d ==\n", node);
      int m = nodes[node].nrow_expected + nodes[node].ndelay_in;
      int n = nodes[node].ncol_expected + nodes[node].ndelay_in;
      const int *rptr = &nodes[node].rlist[ nodes[node].ncol_expected ];
      print_node(posdef, m, n, nodes[node].nelim, nodes[node].perm,
            rptr, nodes[node].lcol, &nodes[node].lcol[m*n]);
   }
}

}}} /* end of namespace spral::ssids::cpu */
//////////////////////////////////////////////////////////////////////////

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_factor_cpu_dbl(
      bool posdef,     // If true, performs A=LL^T, if false do pivoted A=LDL^T
      int n,            // Maximum row index (+1)
      int nnodes,       // Number of nodes in assembly tree
      struct spral::ssids::cpu::cpu_node_data<double> *const nodes, // Data structure for node information
      const double *const aval, // Values of A
      const double *const scaling, // Scaling vector (NULL if none)
      void *const alloc,      // Pointer to Fortran allocator structure
      const struct spral::ssids::cpu::cpu_factor_options *const options, // Options in
      struct spral::ssids::cpu::cpu_factor_stats *const stats // Info out
      ) {

   // Initialize stats
   stats->flag = spral::ssids::cpu::SSIDS_SUCCESS;

   // Call relevant routine
   if(posdef) {
      try {
         spral::ssids::cpu::factor<true, double, 16, 16384, false>
            (n, nnodes, nodes, aval, scaling, alloc, options, stats);
      } catch(spral::ssids::cpu::NotPosDefError npde) {
         stats->flag = spral::ssids::cpu::SSIDS_ERROR_NOT_POS_DEF;
      }
   } else {
      spral::ssids::cpu::factor<false, double, 16, 16384, false>
         (n, nnodes, nodes, aval, scaling, alloc, options, stats);
   }

   // FIXME: Remove when done with debug
   if(options->print_level > 9999) {
      printf("Final factors:\n");
      spral::ssids::cpu::print_factors<double>(posdef, nnodes, nodes);
   }
}
