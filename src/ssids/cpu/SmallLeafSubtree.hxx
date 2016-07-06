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
#pragma once

/* SPRAL headers */
#include "kernels/assemble.hxx"
#include "factor.hxx"

namespace spral { namespace ssids { namespace cpu {

/** Handles the factorization of a small leaf subtree on a single core.
 *
 * This code uses supernodal working within the tree, and generates a
 * multifrontal-style contribution block above the tree. The analyse phase
 * generates internal data structures that guide the assembly process in an
 * efficient fashion, aiming to maximize vectorization.
 *
 * It is expected that the subtree will fit within L2 cache exclusively owned
 * by the executing thread.
 */
template <bool posdef,
          int BLOCK_SIZE,
          typename T,
          typename StackAllocator>
class SmallLeafSubtree {
public:
   /** Constructor performs analyse phase work */
   SmallLeafSubtree(int nnodes, struct cpu_node_data<T>* const nodes)
   : nnodes_(nnodes), nodes_(nodes)
   {}
   void factor(void* const alloc, StackAllocator &stalloc_odd, StackAllocator &stalloc_even, int* map, Workspace<T> &work, T const* aval, T const* scaling, struct cpu_factor_options const* options, struct cpu_factor_stats* stats) {
      /* Main loop: Iterate over nodes in order */
      for(int ni=0; ni<nnodes_; ++ni) {
         // Assembly
         assemble_node
            (posdef, ni, &nodes_[ni], alloc, &stalloc_odd, &stalloc_even, map,
             aval, scaling);
         // Update stats
         int nrow = nodes_[ni].nrow_expected + nodes_[ni].ndelay_in;
         stats->maxfront = std::max(stats->maxfront, nrow);
         // Factorization
         factor_node
            <posdef, T, BLOCK_SIZE>
            (ni, &nodes_[ni], options, stats);
         // Form update
         calculate_update<posdef>
            (&nodes_[ni], &stalloc_odd, &stalloc_even, &work);
      }
   }
   void solve_fwd() {
   }
   void solve_diag() {
   }
   void solve_bwd() {
   }
private:
   int nnodes_;
   struct cpu_node_data<T>* const nodes_;
};

}}} /* namespaces spral::ssids::cpu */
