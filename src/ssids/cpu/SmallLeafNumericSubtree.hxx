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

#include <memory>

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
   typedef typename std::allocator_traits<FactorAllocator>::template rebind_traits<double> FADoubleTraits;
   typedef typename std::allocator_traits<FactorAllocator>::template rebind_traits<int> FAIntTraits;
   typedef std::allocator_traits<ContribAllocator> CATraits;
public:
   SmallLeafNumericSubtree(SmallLeafSymbolicSubtree const& symb, std::vector<NumericNode<T>>& old_nodes, T const* aval, T const* scaling, FactorAllocator& factor_alloc, ContribAllocator& contrib_alloc, Workspace& work, struct cpu_factor_options const& options, struct cpu_factor_stats& stats) 
      : old_nodes_(old_nodes), symb_(symb), lcol_(FADoubleTraits::allocate(factor_alloc, symb.nfactor_))
   {
      /* Initialize nodes */
      for(int ni=symb_.sa_; ni<=symb_.en_; ++ni) {
         old_nodes_[ni].ndelay_in = 0;
         old_nodes_[ni].lcol = lcol_ + symb_[ni-symb_.sa_].lcol_offset;
      }
      memset(lcol_, 0, symb_.nfactor_*sizeof(T));

      /* Add aval entries */
      for(int ni=symb_.sa_; ni<=symb_.en_; ++ni)
         add_a(ni-symb_.sa_, symb_.symb_[ni], aval, scaling);

      /* Perform factorization */
      for(int ni=symb_.sa_; ni<=symb_.en_; ++ni) {
         // Assembly
         int* map = work.get_ptr<int>(symb_.symb_.n+1);
         assemble
            (ni-symb_.sa_, symb_.symb_[ni], &old_nodes_[ni], factor_alloc,
             contrib_alloc, map, aval, scaling);
         // Update stats
         int nrow = symb_.symb_[ni].nrow;
         stats.maxfront = std::max(stats.maxfront, nrow);
         // Factorization
         factor_node
            <true, BLOCK_SIZE>
            (ni, symb_.symb_[ni], &old_nodes_[ni], options, stats);
         if(stats.flag<SSIDS_SUCCESS) return;
      }
   }

private:
void add_a(
      int si,
      SymbolicNode const& snode,
      T const* aval,
      T const* scaling
      ) {
   double *lcol = lcol_ + symb_[si].lcol_offset;
   size_t ldl = align_lda<double>(snode.nrow);
   if(scaling) {
      /* Scaling to apply */
      for(int i=0; i<snode.num_a; i++) {
         long src  = snode.amap[2*i+0] - 1; // amap contains 1-based values
         long dest = snode.amap[2*i+1] - 1; // amap contains 1-based values
         int c = dest / snode.nrow;
         int r = dest % snode.nrow;
         T rscale = scaling[ snode.rlist[r]-1 ];
         T cscale = scaling[ snode.rlist[c]-1 ];
         size_t k = c*ldl + r;
         lcol[k] = rscale * aval[src] * cscale;
      }
   } else {
      /* No scaling to apply */
      for(int i=0; i<snode.num_a; i++) {
         long src  = snode.amap[2*i+0] - 1; // amap contains 1-based values
         long dest = snode.amap[2*i+1] - 1; // amap contains 1-based values
         int c = dest / snode.nrow;
         int r = dest % snode.nrow;
         size_t k = c*ldl + r;
         lcol[k] = aval[src];
      }
   }
}

void assemble(
      int si,
      SymbolicNode const& snode,
      NumericNode<T>* node,
      FactorAllocator& factor_alloc,
      ContribAllocator& contrib_alloc,
      int* map,
      T const* aval,
      T const* scaling
      ) {
   /* Rebind allocators */
   typename FAIntTraits::allocator_type factor_alloc_int(factor_alloc);

   /* Count incoming delays and determine size of node */
   int nrow = snode.nrow;
   int ncol = snode.ncol;

   /* Get space for contribution block + zero it */
   long contrib_dimn = snode.nrow - snode.ncol;
   node->contrib = (contrib_dimn > 0) ? CATraits::allocate(contrib_alloc, contrib_dimn*contrib_dimn) : nullptr;
   if(node->contrib)
      memset(node->contrib, 0, contrib_dimn*contrib_dimn*sizeof(T));

   /* Alloc + set perm */
   node->perm = FAIntTraits::allocate(factor_alloc_int, ncol); // ncol fully summed variables
   //node->perm = smalloc<int>(alloc, ncol); // ncol fully summed variables
   for(int i=0; i<snode.ncol; i++)
      node->perm[i] = snode.rlist[i];

   /* Add children */
   if(node->first_child != NULL) {
      /* Build lookup vector, allowing for insertion of delayed vars */
      /* Note that while rlist[] is 1-indexed this is fine so long as lookup
       * is also 1-indexed (which it is as it is another node's rlist[] */
      for(int i=0; i<snode.nrow; i++)
         map[ snode.rlist[i] ] = i;
      /* Loop over children adding contributions */
      for(auto* child=node->first_child; child!=NULL; child=child->next_child) {
         SymbolicNode const& csnode = *child->symb;
         /* Handle expected contributions (only if something there) */
         if(child->contrib) {
            int cm = csnode.nrow - csnode.ncol;
            for(int i=0; i<cm; i++) {
               int c = map[ csnode.rlist[csnode.ncol+i] ];
               T *src = &child->contrib[i*cm];
               if(c < snode.ncol) {
                  // Contribution added to lcol
                  int ldd = align_lda<double>(nrow);
                  T *dest = &node->lcol[c*ldd];
                  for(int j=i; j<cm; j++) {
                     int r = map[ csnode.rlist[csnode.ncol+j] ];
                     dest[r] += src[j];
                  }
               } else {
                  // Contribution added to contrib
                  // FIXME: Add after contribution block established?
                  int ldd = snode.nrow - snode.ncol;
                  T *dest = &node->contrib[(c-ncol)*ldd];
                  for(int j=i; j<cm; j++) {
                     int r = map[ csnode.rlist[csnode.ncol+j] ] - ncol;
                     dest[r] += src[j];
                  }
               }
            }
            /* Free memory from child contribution block */
            CATraits::deallocate(contrib_alloc, child->contrib, cm*cm);
         }
      }
   }
}

private:
   std::vector<NumericNode<T>>& old_nodes_;
   SmallLeafSymbolicSubtree const& symb_;
   T* lcol_;
};

}}} /* namespaces spral::ssids::cpu */
