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
#pragma once

#include <cstring>
#include "../smalloc.hxx"
#include "../SymbolicSubtree.hxx"

namespace spral { namespace ssids { namespace cpu {

template <typename T,
          typename StackAllocator>
void assemble_node(
      bool posdef,
      int ni, // FIXME: remove with debug
      SymbolicNode const& snode,
      struct cpu_node_data<T> *const node,
      void *const alloc,
      StackAllocator *stalloc_odd,
      StackAllocator *stalloc_even,
      int *const map,
      const T *const aval,
      const T *const scaling
      ) {
   /* Count incoming delays and determine size of node */
   node->ndelay_in = 0;
   for(struct cpu_node_data<T> *child=node->first_child; child!=NULL; child=child->next_child) {
      node->ndelay_in += child->ndelay_out;
   }
   int nrow = snode.nrow + node->ndelay_in;
   int ncol = snode.ncol + node->ndelay_in;

   /* Get space for node now we know it size using Fortran allocator + zero it*/
   // NB L is  nrow x ncol and D is 2 x ncol (but no D if posdef)
   size_t len = posdef ? ((size_t) nrow  ) * ncol  // posdef
                       : ((size_t) nrow+2) * ncol; // indef (includes D)
   if(!(node->lcol = smalloc<T>(alloc, len)))
      throw std::bad_alloc();
   memset(node->lcol, 0, len*sizeof(T));

   /* Get space for contribution block + zero it */
   long contrib_dimn = snode.nrow - snode.ncol;
   if(snode.even) {
      node->contrib = (contrib_dimn > 0) ? (T *) stalloc_even->alloc(contrib_dimn*contrib_dimn*sizeof(T)) : NULL;
   } else {
      node->contrib = (contrib_dimn > 0) ? (T *) stalloc_odd->alloc(contrib_dimn*contrib_dimn*sizeof(T)) : NULL;
   }
   if(node->contrib)
      memset(node->contrib, 0, contrib_dimn*contrib_dimn*sizeof(T));

   /* Alloc + set perm for expected eliminations at this node (delays are set
    * when they are imported from children) */
   node->perm = smalloc<int>(alloc, ncol); // ncol fully summed variables
   for(int i=0; i<snode.ncol; i++)
      node->perm[i] = snode.rlist[i];

   /* Add A */
   if(scaling) {
      /* Scaling to apply */
      for(int i=0; i<node->num_a; i++) {
         long src  = node->amap[2*i+0] - 1; // amap contains 1-based values
         long dest = node->amap[2*i+1] - 1; // amap contains 1-based values
         int c = dest / snode.nrow;
         int r = dest % snode.nrow;
         long k = c*nrow + r;
         if(r >= snode.ncol) k += node->ndelay_in;
         T rscale = scaling[ snode.rlist[r]-1 ];
         T cscale = scaling[ snode.rlist[c]-1 ];
         node->lcol[k] = rscale * aval[src] * cscale;
      }
   } else {
      /* No scaling to apply */
      for(int i=0; i<node->num_a; i++) {
         long src  = node->amap[2*i+0] - 1; // amap contains 1-based values
         long dest = node->amap[2*i+1] - 1; // amap contains 1-based values
         int c = dest / snode.nrow;
         int r = dest % snode.nrow;
         long k = c*nrow + r;
         if(r >= snode.ncol) k += node->ndelay_in;
         node->lcol[k] = aval[src];
      }
   }

   /* Add children */
   if(node->first_child != NULL) {
      /* Build lookup vector, allowing for insertion of delayed vars */
      /* Note that while rlist[] is 1-indexed this is fine so long as lookup
       * is also 1-indexed (which it is as it is another node's rlist[] */
      for(int i=0; i<snode.ncol; i++)
         map[ snode.rlist[i] ] = i;
      for(int i=snode.ncol; i<snode.nrow; i++)
         map[ snode.rlist[i] ] = i + node->ndelay_in;
      /* Loop over children adding contributions */
      int delay_col = snode.ncol;
      for(struct cpu_node_data<T> *child=node->first_child; child!=NULL; child=child->next_child) {
         SymbolicNode const& csnode = *child->symb;
         /* Handle delays - go to back of node
          * (i.e. become the last rows as in lower triangular format) */
         for(int i=0; i<child->ndelay_out; i++) {
            // Add delayed rows (from delayed cols)
            T *dest = &node->lcol[delay_col*(nrow+1)];
            int lds = csnode.nrow + child->ndelay_in;
            T *src = &child->lcol[(child->nelim+i)*(lds+1)];
            node->perm[delay_col] = child->perm[child->nelim+i];
            for(int j=0; j<child->ndelay_out-i; j++) {
               dest[j] = src[j];
            }
            // Add child's non-fully summed rows (from delayed cols)
            dest = node->lcol;
            src = &child->lcol[child->nelim*lds + child->ndelay_in +i*lds];
            for(int j=csnode.ncol; j<csnode.nrow; j++) {
               int r = map[ csnode.rlist[j] ];
               if(r < ncol) dest[r*nrow+delay_col] = src[j];
               else         dest[delay_col*nrow+r] = src[j];
            }
            delay_col++;
         }

         /* Handle expected contributions (only if something there) */
         if(child->contrib) {
            int cm = csnode.nrow - csnode.ncol;
            for(int i=0; i<cm; i++) {
               int c = map[ csnode.rlist[csnode.ncol+i] ];
               T *src = &child->contrib[i*cm];
               if(c < snode.ncol) {
                  // Contribution added to lcol
                  int ldd = nrow;
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
            if(csnode.even) {
               stalloc_even->free(child->contrib, cm*cm*sizeof(T));
            } else {
               stalloc_odd->free(child->contrib, cm*cm*sizeof(T));
            }
         }
      }
   }

   // FIXME: debug remove
   /*printf("Post asm node:\n");
   for(int i=0; i<nrow; i++) {
      for(int j=0; j<ncol; j++) printf(" %10.2e", node->lcol[j*nrow+i]);
      printf("\n");
   }*/
   /*printf("Post asm contrib:\n");
   int ldd = snode.nrow - snode.ncol;
   for(int i=0; i<ldd; i++) {
      for(int j=0; j<ldd; j++) printf(" %e", node->contrib[j*ldd+i]);
      printf("\n");
   }*/
}

}}} /* namespaces spral::ssids::cpu */
