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
class SmallLeafSubtree {
private:
   class Node {
   public:
      int nrow;
      int ncol;
      int sparent;
      int* rlist;
   };

public:
   /** Constructor performs analyse phase work */
   SmallLeafSubtree(int sa, int en, int const* sptr, int const* sparent, long const* rptr, int const* rlist, int const* nptr, int const* nlist)
   : nnodes_(en-sa+1), nodes_(nnodes_), rlist_(nullptr), nptr_(nptr), nlist_(nlist)
   {
      /* Setup basic node information */
      nfactor_ = 0;
      for(int ni=sa; ni<en; ++ni) {
         nodes_[ni-sa].nrow = rptr[ni+1] - rptr[ni];
         nodes_[ni-sa].ncol = sptr[ni+1] - sptr[ni];
         nodes_[ni-sa].sparent = sparent[ni]-sa-1; // sparent is Fortran indexed
         nodes_[ni-sa].rlist = &rlist_[rptr[ni]-rptr[sa]];
         nfactor_ += nodes_[ni-sa].nrow * nodes_[ni-sa].ncol;
      }
      /* Construct rlist_ being offsets into parent node */
      rlist_ = new int[rptr[nnodes_] - rptr[0]];
      for(int ni=sa; ni<en; ++ni) {
         int const* ilist = &rlist[rptr[ni]-1]; // rptr is Fortran indexed
         int pnode = sparent[ni];
         int const* jlist = &rlist[rptr[pnode]-1]; // rptr is Fortran indexed
         int const* jstart = jlist;
         int *outlist = nodes_[ni-sa].rlist;
         for(int i=0; i<nodes_[ni-sa].nrow; ++i) {
            for(; *ilist != *jlist; ++jlist); // Finds match in jlist
            *(outlist++) = jlist - jstart;
         }
      }
   }
   ~SmallLeafSubtree() {
      delete[] rlist_;
   }
   template <typename T>
   T* factor(T const* aval, T const* scaling, T* lcol) {
      /* Zero factors */
      memset(lcol, 0, nfactor_*sizeof(T));
      /* Copy a values */
      for(int i=0; i<nptr_[nnodes_+1]; ++i) {
         // FIXME
      }
      /* Perform factorization */
      for(int ni=0; ni<nnodes_; ++ni) {
         // FIXME
      }
      // FIXME: Return contribution block
   }
private:
   int nnodes_;
   int nfactor_;
   std::vector<Node> nodes_;
   int* rlist_;
   int const* nptr_;
   int const* nlist_;
};

}}} /* namespaces spral::ssids::cpu */
