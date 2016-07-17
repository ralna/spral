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

#include <cstddef>
#include <vector>

#include "SmallLeafSymbolicSubtree.hxx"
#include "SymbolicNode.hxx"

namespace spral { namespace ssids { namespace cpu {

/** Symbolic factorization of a subtree to be factored on the CPU */
class SymbolicSubtree {
public:
   SymbolicSubtree(int n, int nnodes, int const* sptr, int const* sparent, long const* rptr, int const* rlist, int const* nptr, int const* nlist, struct cpu_factor_options const& options)
   : n(n), nnodes_(nnodes), nodes_(nnodes+1)
   {
      // FIXME: don't process nodes that are in small leaf subtrees
      /* Fill out basic details */
      for(int ni=0; ni<nnodes_; ++ni) {
         nodes_[ni].idx = ni;
         nodes_[ni].nrow = static_cast<int>(rptr[ni+1] - rptr[ni]);
         nodes_[ni].ncol = sptr[ni+1] - sptr[ni];
         nodes_[ni].first_child = nullptr;
         nodes_[ni].next_child = nullptr;
         nodes_[ni].rlist = &rlist[rptr[ni]-1]; // rptr is Fortran indexed
         nodes_[ni].num_a = nptr[ni+1] - nptr[ni];
         nodes_[ni].amap = &nlist[2*(nptr[ni]-1)]; // nptr is Fortran indexed
         nodes_[ni].parent = sparent[ni]-1; // sparent is Fortran indexed
         nodes_[ni].insmallleaf = false; // default to not in small leaf subtree
      }
      nodes_[nnodes_].first_child = nullptr; // List of roots
      /* Build child linked lists */
      for(int ni=0; ni<nnodes_; ++ni) {
         SymbolicNode *parent = &nodes_[ sparent[ni]-1 ];
         nodes_[ni].next_child = parent->first_child;
         parent->first_child = &nodes_[ni];
      }
      /* Count size of factors */
      nfactor_ = 0;
      for(int ni=0; ni<nnodes_; ++ni)
         nfactor_ += static_cast<size_t>(nodes_[ni].nrow)*nodes_[ni].ncol;
      /* Find small leaf subtrees */
      // Count flops below each node
      std::vector<long> flops(nnodes_+1, 0);
      for(int ni=0; ni<nnodes_; ++ni) {
         for(int k=0; k<nodes_[ni].ncol; ++k)
            flops[ni] += (nodes_[ni].nrow - k)*(nodes_[ni].nrow - k);
         flops[nodes_[ni].parent] += flops[ni];
      }
      // Start at least node and work way up using parents until too large
      for(int ni=0; ni<nnodes_; ) {
         if(nodes_[ni].first_child) { ++ni; continue; } // Not a leaf
         int last = ni;
         for(int current=ni; current<nnodes_; current=nodes_[current].parent) {
            if(flops[current] >= options.cpu_small_subtree_threshold) break;
            last = current;
         }
         if(last==ni) { ++ni; continue; } // No point for a single node
         // Nodes ni:last are in subtree
         small_leafs_.emplace_back(
               ni, last, sptr, sparent, rptr, rlist, nptr, nlist, *this
               );
         for(int i=ni; i<=last; ++i)
            nodes_[i].insmallleaf = true;
         ni = last+1; // Skip to next node not in this subtree
      }
   }

   SymbolicNode const& operator[](int idx) const {
      return nodes_[idx];
   }
   size_t get_factor_mem_est(double multiplier) const {
      size_t mem = n*sizeof(int) + (2*n+nfactor_)*sizeof(double);
      return std::max(mem, static_cast<size_t>(mem*multiplier));
   }
public:
   int const n; //< Maximum row index
private:
   int nnodes_;
   size_t nfactor_;
   std::vector<SymbolicNode> nodes_;
   std::vector<SmallLeafSymbolicSubtree> small_leafs_;

   template <bool posdef, size_t BLOCK_SIZE, typename T, size_t PAGE_SIZE, typename FactorAlloc, typename ContribAllocator>
   friend class NumericSubtree;
};

}}} /* end of namespace spral::ssids::cpu */
