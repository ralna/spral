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

#include <vector>

namespace spral { namespace ssids { namespace cpu {

/** Symbolic representation of a node */
struct SymbolicNode {
   int nrow; //< Number of rows
   int ncol; //< Number of columns
   SymbolicNode* first_child; //< Pointer to first child in linked list
   SymbolicNode* next_child; //< Pointer to second child in linked list
   int const* rlist; //< Pointer to row lists
   bool even; //< Indicates which stack we're using (odd or even distance
              //< from root)
};

/** Symbolic factorization of a subtree to be factored on the CPU */
class SymbolicSubtree {
public:
   SymbolicSubtree(int nnodes, int const* sptr, long const* rptr, int const* rlist)
   : nnodes_(nnodes), nodes_(nnodes)
   {
      for(int ni=0; ni<nnodes_; ++ni) {
         nodes_[ni].nrow = static_cast<int>(rptr[ni+1] - rptr[ni]);
         nodes_[ni].ncol = sptr[ni+1] - sptr[ni];
         nodes_[ni].first_child = nullptr;
         nodes_[ni].next_child = nullptr;
         nodes_[ni].rlist = &rlist[rptr[ni]-1]; // NB: rptr is Fortran indexed
      }
   }

   SymbolicNode const& operator[](int idx) const {
      return nodes_[idx];
   }
private:
   int nnodes_;
   std::vector<SymbolicNode> nodes_;
};

}}} /* end of namespace spral::ssids::cpu */
