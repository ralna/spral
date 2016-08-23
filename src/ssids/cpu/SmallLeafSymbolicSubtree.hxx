/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

#include <memory>

#include "cpu_iface.hxx"
#include "SymbolicNode.hxx"

#include <cstdio> // FIXME debug only

namespace spral { namespace ssids { namespace cpu {

class SymbolicSubtree;

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
class SmallLeafSymbolicSubtree {
private:
   class Node {
   public:
      int nrow;
      int ncol;
      int sparent;
      int* rlist;
      int lcol_offset;
   };

public:
   /** Constructor performs analyse phase work */
   SmallLeafSymbolicSubtree(int sa, int en, int part_offset, int const* sptr, int const* sparent, long const* rptr, int const* rlist, int const* nptr, int const* nlist, SymbolicSubtree const& symb)
   : sa_(sa), en_(en), nnodes_(en-sa+1), parent_(sparent[part_offset+en]-1-part_offset),
     nodes_(nnodes_),
     rlist_(new int[rptr[part_offset+en+1]-rptr[part_offset+sa]], std::default_delete<int[]>()),
     nptr_(nptr), nlist_(nlist), symb_(symb)
   {
      /* Setup basic node information */
      nfactor_ = 0;
      int* newrlist = rlist_.get();
      for(int ni=sa; ni<=en; ++ni) {
         nodes_[ni-sa].nrow = rptr[part_offset+ni+1] - rptr[part_offset+ni];
         nodes_[ni-sa].ncol = sptr[part_offset+ni+1] - sptr[part_offset+ni];
         nodes_[ni-sa].sparent = sparent[part_offset+ni]-sa-1; // sparent is Fortran indexed
         // FIXME: subtract ncol off rlist for elim'd vars
         nodes_[ni-sa].rlist = &newrlist[rptr[part_offset+ni]-rptr[part_offset+sa]];
         nodes_[ni-sa].lcol_offset = nfactor_;
         size_t ldl = align_lda<double>(nodes_[ni-sa].nrow);
         nfactor_ += nodes_[ni-sa].ncol*ldl;
      }
      /* Construct rlist_ being offsets into parent node */
      for(int ni=sa; ni<=en; ++ni) {
         if(nodes_[ni-sa].ncol == nodes_[ni-sa].nrow) continue; // is root
         int const* ilist = &rlist[rptr[part_offset+ni]-1]; // rptr is Fortran indexed
         ilist += nodes_[ni-sa].ncol; // Skip eliminated vars
         int pnode = sparent[part_offset+ni]-1; //Fortran indexed
         int const* jlist = &rlist[rptr[pnode]-1]; // rptr is Fortran indexed
         int const* jstart = jlist;
         int *outlist = nodes_[ni-sa].rlist;
         for(int i=nodes_[ni-sa].ncol; i<nodes_[ni-sa].nrow; ++i) {
            for(; *ilist != *jlist; ++jlist); // Finds match in jlist
            *(outlist++) = jlist - jstart;
            ++ilist;
         }
      }
   }

   int get_parent() const { return parent_; }
   Node const& operator[](int idx) const { return nodes_[idx]; }
protected:
   int sa_;
   int en_;
   int nnodes_;
   int nfactor_;
   int parent_;
   std::vector<Node> nodes_;
   std::shared_ptr<int> rlist_;
   int const* nptr_;
   int const* nlist_;
   SymbolicSubtree const& symb_;
   
   template <bool posdef, typename T, typename FactorAllocator,
             typename PoolAllocator>
   friend class SmallLeafNumericSubtree;
};


}}} /* namespaces spral::ssids::cpu */
