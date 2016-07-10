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

namespace spral { namespace ssids { namespace cpu {

/** Symbolic representation of a node */
struct SymbolicNode {
   int idx; //< Index of node
   int nrow; //< Number of rows
   int ncol; //< Number of columns
   SymbolicNode* first_child; //< Pointer to first child in linked list
   SymbolicNode* next_child; //< Pointer to second child in linked list
   int const* rlist; //< Pointer to row lists
   int num_a; //< Number of entries mapped from A to L
   int const* amap; //< Pointer to map from A to L locations
   bool even; //< Indicates which stack we're using (odd or even distance
              //< from root)
};

/** Symbolic factorization of a subtree to be factored on the CPU */
class SymbolicSubtree {
public:
   SymbolicSubtree(int n, int nnodes, int const* sptr, int const* sparent, long const* rptr, int const* rlist, int const* nptr, int const* nlist)
   : n(n), nnodes_(nnodes), nodes_(nnodes+1)
   {
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
      }
      nodes_[nnodes_].first_child = nullptr; // List of roots
      /* Build child linked lists */
      for(int ni=0; ni<nnodes_; ++ni) {
         SymbolicNode *parent = &nodes_[ sparent[ni]-1 ];
         nodes_[ni].next_child = parent->first_child;
         parent->first_child = &nodes_[ni];
      }
      /* Setup odd/even distance from root */
      for(int ni=nnodes_-1; ni>=0; --ni) {
         if(sparent[ni] <= nnodes_) {
            SymbolicNode const& parent = nodes_[ sparent[ni]-1 ];
            nodes_[ni].even = !(parent.even);
         } else {
            nodes_[ni].even = true;
         }
      }
   }

   SymbolicNode const& operator[](int idx) const {
      return nodes_[idx];
   }
public:
   int const n; //< Maximum row index
private:
   int nnodes_;
   std::vector<SymbolicNode> nodes_;

   template <bool posdef, size_t BLOCK_SIZE, typename T, size_t PAGE_SIZE>
   friend class NumericSubtree;
};

}}} /* end of namespace spral::ssids::cpu */
