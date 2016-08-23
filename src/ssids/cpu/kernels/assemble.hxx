/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

#include<cstring>
#include<memory>
#include<vector>

#include "ssids/contrib.h"
#include "../NumericNode.hxx"
#include "../SymbolicNode.hxx"
#include "../Workspace.hxx"

#include "../profile.hxx"

namespace spral { namespace ssids { namespace cpu {

/** Assemble a column.
 *
 * Performs the operation dest( idx(:) ) += src(:)
 */
template <typename T>
inline
void asm_col(int n, int const* idx, T const* src, T* dest) {
   int const nunroll = 4;
   int n2 = nunroll*(n/nunroll);
   for(int j=0; j<n2; j+=nunroll) {
      dest[ idx[j+0] ] += src[j+0];
      dest[ idx[j+1] ] += src[j+1];
      dest[ idx[j+2] ] += src[j+2];
      dest[ idx[j+3] ] += src[j+3];
   }
   for(int j=n2; j<n; j++)
      dest[ idx[j] ] += src[j];
}

template <typename T,
          typename FactorAlloc,
          typename PoolAlloc>
void assemble_pre(
      bool posdef,
      int n,
      SymbolicNode const& snode,
      void** child_contrib,
      NumericNode<T>& node,
      FactorAlloc& factor_alloc,
      PoolAlloc& pool_alloc,
      std::vector<Workspace>& work,
      T const* aval,
      T const* scaling
      ) {
#ifdef PROFILE
   Profile::Task task_asm_pre("TA_ASM_PRE");
#endif
   /* Rebind allocators */
   typedef typename std::allocator_traits<FactorAlloc>::template rebind_traits<double> FADoubleTraits;
   typename FADoubleTraits::allocator_type factor_alloc_double(factor_alloc);
   typedef typename std::allocator_traits<FactorAlloc>::template rebind_traits<int> FAIntTraits;
   typename FAIntTraits::allocator_type factor_alloc_int(factor_alloc);
   typedef std::allocator_traits<PoolAlloc> PATraits;
   typedef typename std::allocator_traits<PoolAlloc>::template rebind_traits<int> PAIntTraits;
   typename PAIntTraits::allocator_type pool_alloc_int(pool_alloc);

   /* Count incoming delays and determine size of node */
   node.ndelay_in = 0;
   for(auto* child=node.first_child; child!=NULL; child=child->next_child) {
      node.ndelay_in += child->ndelay_out;
   }
   for(int contrib_idx : snode.contrib) {
      int cn, ldcontrib, ndelay, lddelay;
      double const *cval, *delay_val;
      int const *crlist, *delay_perm;
      spral_ssids_contrib_get_data(
            child_contrib[contrib_idx], &cn, &cval, &ldcontrib, &crlist,
            &ndelay, &delay_perm, &delay_val, &lddelay
            );
      node.ndelay_in += ndelay;
   }
   int nrow = snode.nrow + node.ndelay_in;
   int ncol = snode.ncol + node.ndelay_in;

   /* Get space for node now we know it size using Fortran allocator + zero it*/
   // NB L is  nrow x ncol and D is 2 x ncol (but no D if posdef)
   size_t ldl = align_lda<double>(nrow);
   size_t len = posdef ?  ldl    * ncol  // posdef
                       : (ldl+2) * ncol; // indef (includes D)
   node.lcol = FADoubleTraits::allocate(factor_alloc_double, len);
   //memset(node.lcol, 0, len*sizeof(T)); NOT REQUIRED as PoolAlloc is
   // required to ensure it is zero for us (i.e. uses calloc)

   /* Get space for contribution block + (explicitly do not zero it!) */
   long contrib_dimn = snode.nrow - snode.ncol;
   node.contrib = (contrib_dimn > 0) ? PATraits::allocate(pool_alloc, contrib_dimn*contrib_dimn) : nullptr;

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
   int* map = nullptr;
   int delay_col = snode.ncol;
   if(node.first_child != NULL || snode.contrib.size() > 0) {
      /* Build lookup vector, allowing for insertion of delayed vars */
      /* Note that while rlist[] is 1-indexed this is fine so long as lookup
       * is also 1-indexed (which it is as it is another node's rlist[] */
      if(!map) map = PAIntTraits::allocate(pool_alloc_int, n+1);
      for(int i=0; i<snode.ncol; i++)
         map[ snode.rlist[i] ] = i;
      for(int i=snode.ncol; i<snode.nrow; i++)
         map[ snode.rlist[i] ] = i + node.ndelay_in;
      /* Loop over children adding contributions */
#ifdef PROFILE
      task_asm_pre.done();
#endif
      for(auto* child=node.first_child; child!=NULL; child=child->next_child) {
#ifdef PROFILE
         Profile::Task task_asm_pre("TA_ASM_PRE");
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
                  shared(map, child, snode, node, csnode, cm, nrow, work) \
                  if(iblk+block_size<cm)
               {
#ifdef PROFILE
                  Profile::Task task_asm_pre("TA_ASM_PRE");
#endif
                  int* cache = work[omp_get_thread_num()].get_ptr<int>(cm);
                  for(int j=iblk; j<cm; ++j)
                     cache[j] = map[ csnode.rlist[csnode.ncol+j] ];
                  for(int i=iblk; i<std::min(iblk+block_size,cm); i++) {
                     int c = cache[i];
                     T *src = &child->contrib[i*cm];
                     // NB: we handle contribution to contrib in assemble_post()
                     if(c < snode.ncol) {
                        // Contribution added to lcol
                        int ldd = align_lda<T>(nrow);
                        T *dest = &node.lcol[c*ldd];
                        asm_col(cm-i, &cache[i], &src[i], dest);
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
   /* Add any contribution block from other subtrees */
   for(int contrib_idx : snode.contrib) {
      int cn, ldcontrib, ndelay, lddelay;
      double const *cval, *delay_val;
      int const *crlist, *delay_perm;
      spral_ssids_contrib_get_data(
            child_contrib[contrib_idx], &cn, &cval, &ldcontrib, &crlist,
            &ndelay, &delay_perm, &delay_val, &lddelay
            );
      int* cache = work[omp_get_thread_num()].get_ptr<int>(cn);
      for(int j=0; j<cn; ++j)
         cache[j] = map[ crlist[j] ];
      /* Handle delays - go to back of node
       * (i.e. become the last rows as in lower triangular format) */
      for(int i=0; i<ndelay; i++) {
         // Add delayed rows (from delayed cols)
         T *dest = &node.lcol[delay_col*(ldl+1)];
         T const* src = &delay_val[i*(lddelay+1)];
         node.perm[delay_col] = delay_perm[i];
         for(int j=0; j<ndelay-i; j++) {
            dest[j] = src[j];
         }
         // Add child's non-fully summed rows (from delayed cols)
         dest = node.lcol;
         src = &delay_val[i*lddelay+ndelay];
         for(int j=0; j<cn; j++) {
            int r = cache[j];
            if(r < ncol) dest[r*ldl+delay_col] = src[j];
            else         dest[delay_col*ldl+r] = src[j];
         }
         delay_col++;
      }
      if(!cval) continue; // child was all delays, nothing more to do
      /* Handle expected contribution */
      for(int i=0; i<cn; ++i) {
         int c = cache[i];
         T const* src = &cval[i*ldcontrib];
         // NB: we handle contribution to contrib in assemble_post()
         if(c < snode.ncol) {
            // Contribution added to lcol
            int ldd = align_lda<T>(nrow);
            T *dest = &node.lcol[c*ldd];
            asm_col(cn-i, &cache[i], &src[i], dest);
         }
      }
   }
   if(map) PAIntTraits::deallocate(pool_alloc_int, map, n+1);
}

template <typename T,
          typename PoolAlloc
          >
void assemble_post(
      int n,
      SymbolicNode const& snode,
      void** child_contrib,
      NumericNode<T>& node,
      PoolAlloc& pool_alloc,
      std::vector<Workspace>& work
      ) {
   /* Rebind allocators */
   typedef std::allocator_traits<PoolAlloc> PATraits;
   typedef typename std::allocator_traits<PoolAlloc>::template rebind_traits<int> PAIntTraits;
   typename PAIntTraits::allocator_type pool_alloc_int(pool_alloc);

   /* Initialise variables */
   int ncol = snode.ncol + node.ndelay_in;

   /* Add children */
   int* map = nullptr;
   if(node.first_child != NULL || snode.contrib.size() > 0) {
      /* Build lookup vector, allowing for insertion of delayed vars */
      /* Note that while rlist[] is 1-indexed this is fine so long as lookup
       * is also 1-indexed (which it is as it is another node's rlist[] */
      if(!map) map = PAIntTraits::allocate(pool_alloc_int, n+1);
      // FIXME: probably don't need to worry about first ncol?
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
               shared(map, child, snode, node, cm, csnode, ncol, work) \
               if(iblk+block_size<cm)
            {
#ifdef PROFILE
               Profile::Task task_asm("TA_ASM_POST");
#endif
               int* cache = work[omp_get_thread_num()].get_ptr<int>(cm);
               for(int j=iblk; j<cm; ++j)
                  cache[j] = map[ csnode.rlist[csnode.ncol+j] ] - ncol;
               for(int i=iblk; i<std::min(iblk+block_size,cm); i++) {
                  int c = cache[i]+ncol;
                  T *src = &child->contrib[i*cm];
                  // NB: only interested in contribution to generated element
                  if(c >= snode.ncol) {
                     // Contribution added to contrib
                     int ldd = snode.nrow - snode.ncol;
                     T *dest = &node.contrib[(c-ncol)*ldd];
                     asm_col(cm-i, &cache[i], &src[i], dest);
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
         PATraits::deallocate(pool_alloc, child->contrib, cm*cm);
      }
   }
   /* Add any contribution block from other subtrees */
   for(int contrib_idx : snode.contrib) {
      int cn, ldcontrib, ndelay, lddelay;
      double const *cval, *delay_val;
      int const *crlist, *delay_perm;
      spral_ssids_contrib_get_data(
            child_contrib[contrib_idx], &cn, &cval, &ldcontrib, &crlist,
            &ndelay, &delay_perm, &delay_val, &lddelay
            );
      if(!cval) continue; // child was all delays, nothing to do
      int* cache = work[omp_get_thread_num()].get_ptr<int>(cn);
      for(int j=0; j<cn; ++j)
         cache[j] = map[ crlist[j] ] - ncol;
      for(int i=0; i<cn; ++i) {
         int c = cache[i]+ncol;
         T const* src = &cval[i*ldcontrib];
         // NB: only interested in contribution to generated element
         if(c >= snode.ncol) {
            // Contribution added to contrib
            int ldd = snode.nrow - snode.ncol;
            T *dest = &node.contrib[(c-ncol)*ldd];
            asm_col(cn-i, &cache[i], &src[i], dest);
         }
      }
      /* Free memory from child contribution block */
      spral_ssids_contrib_free_dbl(child_contrib[contrib_idx]);
   }
   if(map) PAIntTraits::deallocate(pool_alloc_int, map, n+1);
}

}}} /* namespaces spral::ssids::cpu */
