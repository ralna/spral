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

#include "cpu_iface.hxx"
#include "factor.hxx"
#include "NumericNode.hxx"
#include "SmallLeafSymbolicSubtree.hxx"

/* SPRAL headers */

namespace spral { namespace ssids { namespace cpu {

template <typename T,
          int BLOCK_SIZE,
          typename FactorAllocator, // Allocator to use for factor storage
          typename ContribAllocator // Allocator to use for contribution blocks
          >
class SmallLeafNumericSubtree {
   SmallLeafNumericSubtree(SmallLeafSymbolicSubtree const& symb_, std::vector<NumericNode<T>>& old_nodes, T const* aval, T const* scaling, T* lcol, FactorAllocator& factor_alloc, ContribAllocator& contrib_alloc, Workspace& work, struct cpu_factor_options const* options, struct cpu_factor_stats &stats) 
      : old_nodes_(old_nodes)
   {
      for(int ni=symb_.sa_; ni<symb_.en_; ++ni) {
         // Assembly
         int* map = work.get_ptr<int>(symb_.symb_.n+1);
         assemble_node
            (true, ni, symb_.symb_[ni], &old_nodes_[ni], factor_alloc,
             contrib_alloc, map, aval, scaling);
         // Update stats
         int nrow = symb_.symb_[ni].nrow + old_nodes_[ni].ndelay_in;
         stats.maxfront = std::max(stats.maxfront, nrow);
         // Factorization
         factor_node
            <true, BLOCK_SIZE>
            (ni, symb_.symb_[ni], &old_nodes_[ni], options, stats);
         if(stats.flag<SSIDS_SUCCESS) return;
      }
   }

private:
   std::vector<NumericNode<T>>& old_nodes_;
   SmallLeafSymbolicSubtree const& symb_;
};

}}} /* namespaces spral::ssids::cpu */
