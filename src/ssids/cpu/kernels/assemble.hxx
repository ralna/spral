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

#include<cstring>
#include<memory>

#include "../NumericNode.hxx"
#include "../SymbolicNode.hxx"

#include "../profile.hxx"

namespace spral { namespace ssids { namespace cpu {

template <typename T,
          typename FactorAlloc,
          typename ContribAlloc>
void assemble_pre(
      bool posdef,
      SymbolicNode const& snode,
      NumericNode<T>& node,
      FactorAlloc& factor_alloc,
      ContribAlloc& contrib_alloc,
      int* map,
      T const* aval,
      T const* scaling
      ) {
#ifdef PROFILE
   Profile::Task task_asm_pre("TA_ASM_PRE", omp_get_thread_num());
#endif
   /* Rebind allocators */
   typedef typename std::allocator_traits<FactorAlloc>::template rebind_traits<double> FADoubleTraits;
   typename FADoubleTraits::allocator_type factor_alloc_double(factor_alloc);
   typedef typename std::allocator_traits<FactorAlloc>::template rebind_traits<int> FAIntTraits;
   typename FAIntTraits::allocator_type factor_alloc_int(factor_alloc);
   typedef std::allocator_traits<ContribAlloc> CATraits;

   /* Count incoming delays and determine size of node */
   node.ndelay_in = 0;
   for(auto* child=node.first_child; child!=NULL; child=child->next_child) {
      node.ndelay_in += child->ndelay_out;
   }
   int nrow = snode.nrow + node.ndelay_in;
   int ncol = snode.ncol + node.ndelay_in;

   /* Get space for node now we know it size using Fortran allocator + zero it*/
   // NB L is  nrow x ncol and D is 2 x ncol (but no D if posdef)
   size_t ldl = align_lda<double>(nrow);
   size_t len = posdef ?  ldl    * ncol  // posdef
                       : (ldl+2) * ncol; // indef (includes D)
   node.lcol = FADoubleTraits::allocate(factor_alloc_double, len);
   //memset(node.lcol, 0, len*sizeof(T)); NOT REQUIRED as ContribAlloc is
   // required to ensure it is zero for us (i.e. uses calloc)

   /* Get space for contribution block + (explicitly do not zero it!) */
   long contrib_dimn = snode.nrow - snode.ncol;
   node.contrib = (contrib_dimn > 0) ? CATraits::allocate(contrib_alloc, contrib_dimn*contrib_dimn) : nullptr;

   /* Alloc + set perm for expected eliminations at this node (delays are set
    * when they are imported from children) */
   node.perm = FAIntTraits::allocate(factor_alloc_int, ncol); // ncol fully summed variables
   for(int i=0; i<snode.ncol; i++)
      node.perm[i] = snode.rlist[i];

   /* Add A */
   int const add_a_blk_sz = 256;
   for(int iblk=0; iblk<snode.num_a; iblk+=add_a_blk_sz) {
      #pragma omp task default(none) \
         firstprivate(iblk) \
         shared(snode, node, aval, scaling, ldl) \
         if(iblk+add_a_blk_sz < snode.num_a)
      if(scaling) {
         /* Scaling to apply */
         for(int i=iblk; i<std::min(iblk+add_a_blk_sz,snode.num_a); ++i) {
            long src  = snode.amap[2*i+0] - 1; // amap contains 1-based values
            long dest = snode.amap[2*i+1] - 1; // amap contains 1-based values
            int c = dest / snode.nrow;
            int r = dest % snode.nrow;
            long k = c*ldl + r;
            if(r >= snode.ncol) k += node.ndelay_in;
            T rscale = scaling[ snode.rlist[r]-1 ];
            T cscale = scaling[ snode.rlist[c]-1 ];
            node.lcol[k] = rscale * aval[src] * cscale;
         }
      } else {
         /* No scaling to apply */
         for(int i=iblk; i<std::min(iblk+add_a_blk_sz,snode.num_a); ++i) {
            long src  = snode.amap[2*i+0] - 1; // amap contains 1-based values
            long dest = snode.amap[2*i+1] - 1; // amap contains 1-based values
            int c = dest / snode.nrow;
            int r = dest % snode.nrow;
            long k = c*ldl + r;
            if(r >= snode.ncol) k += node.ndelay_in;
            node.lcol[k] = aval[src];
         }
      }
   }
   if(add_a_blk_sz < snode.num_a) {
      #pragma omp taskwait
   }
#ifdef PROFILE
   if(!node.first_child) task_asm_pre.done();
#endif

   /* Add children */
   if(node.first_child != NULL) {
      /* Build lookup vector, allowing for insertion of delayed vars */
      /* Note that while rlist[] is 1-indexed this is fine so long as lookup
       * is also 1-indexed (which it is as it is another node's rlist[] */
      for(int i=0; i<snode.ncol; i++)
         map[ snode.rlist[i] ] = i;
      for(int i=snode.ncol; i<snode.nrow; i++)
         map[ snode.rlist[i] ] = i + node.ndelay_in;
      /* Loop over children adding contributions */
      int delay_col = snode.ncol;
#ifdef PROFILE
      task_asm_pre.done();
#endif
      for(auto* child=node.first_child; child!=NULL; child=child->next_child) {
#ifdef PROFILE
         Profile::Task task_asm_pre("TA_ASM_PRE", omp_get_thread_num());
#endif
         SymbolicNode const& csnode = *child->symb;
         /* Handle delays - go to back of node
          * (i.e. become the last rows as in lower triangular format) */
         for(int i=0; i<child->ndelay_out; i++) {
            // Add delayed rows (from delayed cols)
            T *dest = &node.lcol[delay_col*(ldl+1)];
            int lds = align_lda<T>(csnode.nrow + child->ndelay_in);
            T *src = &child->lcol[(child->nelim+i)*(lds+1)];
            node.perm[delay_col] = child->perm[child->nelim+i];
            for(int j=0; j<child->ndelay_out-i; j++) {
               dest[j] = src[j];
            }
            // Add child's non-fully summed rows (from delayed cols)
            dest = node.lcol;
            src = &child->lcol[child->nelim*lds + child->ndelay_in +i*lds];
            for(int j=csnode.ncol; j<csnode.nrow; j++) {
               int r = map[ csnode.rlist[j] ];
               if(r < ncol) dest[r*ldl+delay_col] = src[j];
               else         dest[delay_col*ldl+r] = src[j];
            }
            delay_col++;
         }
#ifdef PROFILE
         task_asm_pre.done();
#endif

         /* Handle expected contributions (only if something there) */
         if(child->contrib) {
            int cm = csnode.nrow - csnode.ncol;
            int const block_size = 256; // FIXME: make configurable?
            for(int iblk=0; iblk<cm; iblk+=block_size) {
               #pragma omp task default(none) \
                  firstprivate(iblk) \
                  shared(map, child, snode, node, csnode, cm, nrow) \
                  if(iblk+block_size<cm)
               {
#ifdef PROFILE
                  Profile::Task task_asm_pre("TA_ASM_PRE", omp_get_thread_num());
#endif
                  for(int i=iblk; i<std::min(iblk+block_size,cm); i++) {
                     int c = map[ csnode.rlist[csnode.ncol+i] ];
                     T *src = &child->contrib[i*cm];
                     // NB: we handle contribution to contrib in assemble_post()
                     if(c < snode.ncol) {
                        // Contribution added to lcol
                        int ldd = align_lda<T>(nrow);
                        T *dest = &node.lcol[c*ldd];
                        for(int j=i; j<cm; j++) {
                           int r = map[ csnode.rlist[csnode.ncol+j] ];
                           dest[r] += src[j];
                        }
                     }
                  }
#ifdef PROFILE
                  task_asm_pre.done();
#endif
               } /* task */
            }
            if(cm > block_size) {
               // only wait if we've actually created tasks
               #pragma omp taskwait
            }
         }
      }
   }
}

template <typename T,
          typename ContribAlloc
          >
void assemble_post(
      SymbolicNode const& snode,
      NumericNode<T>& node,
      ContribAlloc& contrib_alloc,
      int* map
      ) {
   /* Rebind allocators */
   typedef std::allocator_traits<ContribAlloc> CATraits;

   /* Initialise variables */
   int ncol = snode.ncol + node.ndelay_in;

   /* Add children */
   if(node.first_child != NULL) {
      /* Build lookup vector, allowing for insertion of delayed vars */
      /* Note that while rlist[] is 1-indexed this is fine so long as lookup
       * is also 1-indexed (which it is as it is another node's rlist[] */
      for(int i=0; i<snode.ncol; i++)
         map[ snode.rlist[i] ] = i;
      for(int i=snode.ncol; i<snode.nrow; i++)
         map[ snode.rlist[i] ] = i + node.ndelay_in;
      /* Loop over children adding contributions */
      for(auto* child=node.first_child; child!=NULL; child=child->next_child) {
         SymbolicNode const& csnode = *child->symb;
         if(!child->contrib) continue;
         int cm = csnode.nrow - csnode.ncol;
         int const block_size = 256;
         for(int iblk=0; iblk<cm; iblk+=block_size) {
            #pragma omp task default(none) \
               firstprivate(iblk) \
               shared(map, child, snode, node, cm, csnode, ncol) \
               if(iblk+block_size<cm)
            {
#ifdef PROFILE
               Profile::Task task_asm("TA_ASM_POST", omp_get_thread_num());
#endif
               for(int i=iblk; i<std::min(iblk+block_size,cm); i++) {
                  int c = map[ csnode.rlist[csnode.ncol+i] ];
                  T *src = &child->contrib[i*cm];
                  // NB: only interested in contribution to generated element
                  if(c >= snode.ncol) {
                     // Contribution added to contrib
                     int ldd = snode.nrow - snode.ncol;
                     T *dest = &node.contrib[(c-ncol)*ldd];
                     for(int j=i; j<cm; j++) {
                        int r = map[ csnode.rlist[csnode.ncol+j] ] - ncol;
                        dest[r] += src[j];
                     }
                  }
               }
#ifdef PROFILE
               task_asm.done();
#endif
            } /* task */
         }
         if(cm > block_size) {
            // Only wait if we've actually lanuched tasks
            #pragma omp taskwait
         }
         /* Free memory from child contribution block */
         CATraits::deallocate(contrib_alloc, child->contrib, cm*cm);
      }
   }
}

}}} /* namespaces spral::ssids::cpu */
