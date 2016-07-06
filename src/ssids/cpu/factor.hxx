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

/* Standard headers */
#include <cstddef>
#include <sstream>
#include <stdexcept>
/* SPRAL headers */
#include "cpu_iface.hxx"
#include "kernels/assemble.hxx"
#include "kernels/cholesky.hxx"
#include "kernels/CpuLDLT.cxx"

namespace spral { namespace ssids { namespace cpu {

const int SSIDS_SUCCESS = 0;
const int SSIDS_ERROR_NOT_POS_DEF = -6;

template <typename T>
class Workspace {
public:
   T *mem;
   size_t len;
   Workspace(size_t len)
      : mem(new T[len]), len(len)
      {}
   ~Workspace() {
      delete[] mem;
   }
   void ensure_length(size_t newlen) {
      if(len >= newlen) return; // Plenty big enough
      // Otherwise resize
      delete[] mem;
      mem = new T[newlen];
      len = newlen;
   }
};

/* Custom exceptions */
class NotPosDefError: public std::runtime_error {
public:
   int posn;

   NotPosDefError(int posn)
      : runtime_error("Matrix not positive definite"), posn(posn)
   {}

   virtual const char* what() const throw() {
      std::ostringstream cnvt;
      cnvt << std::runtime_error::what() << " (failed at column " << posn << ")";
      return cnvt.str().c_str();
   }
};

/* Factorize a node (indef) */
template <typename T, int BLOCK_SIZE>
void factor_node_indef(
      int ni, // FIXME: remove post debug
      struct cpu_node_data<T> *const node,
      const struct cpu_factor_options *const options,
      struct cpu_factor_stats *const stats
      ) {
   /* Extract useful information about node */
   int m = node->nrow_expected + node->ndelay_in;
   int n = node->ncol_expected + node->ndelay_in;
   T *lcol = node->lcol;
   T *d = &node->lcol[ ((long) m)*n ];
   int *perm = node->perm;

   /* Perform factorization */
   typedef CpuLDLT<T, BLOCK_SIZE> CpuLDLTSpec;
   //typedef CpuLDLT<T, BLOCK_SIZE, 5, true> CpuLDLTSpecDebug; // FIXME: debug remove
   struct CpuLDLTSpec::stat_type bubstats; // FIXME: not needed?
   node->nelim = CpuLDLTSpec(options->u, options->small).factor(m, n, perm, lcol, m, d, &bubstats);
   for(int i=0; i<5; i++) {
      stats->elim_at_pass[i] += bubstats.elim_at_pass[i];
   }
   int last_remain = n;
   for(int i=0; i<bubstats.nitr; i++) {
      stats->elim_at_itr[i] += last_remain - bubstats.remain[i];
      last_remain = bubstats.remain[i];
   }
   /*if(bubstats.nitr > 2) {
      printf("Node %d: %dx%d delay %d nitr %d\n", ni, m, n, n-node->nelim, bubstats.nitr);
      for(int i=0; i<bubstats.nitr; i++)
         printf("--> itr %d passes %d remain %d\n", i, bubstats.npass[i], bubstats.remain[i]);
   }*/

   /*for(int i=node->nelim; i<m; i++) {
      printf("%d:", i);
      for(int j=node->nelim; j<n; j++)
         printf(" %10.2e", lcol[j*m+i]);
      printf("\n");
   }*/

   /* Record information */
   node->ndelay_out = n - node->nelim;
   stats->num_delay += node->ndelay_out;
}
/* Factorize a node (posdef) */
template <typename T, int BLOCK_SIZE>
void factor_node_posdef(
      struct cpu_node_data<T> *const node,
      const struct cpu_factor_options *const options
      ) {
   /* Extract useful information about node */
   int m = node->nrow_expected;
   int n = node->ncol_expected;
   T *lcol = node->lcol;

   /* Perform factorization */
   int flag;
   #pragma omp parallel default(shared)
   {
      #pragma omp single
      cholesky_factor(m, n, lcol, m, options->cpu_task_block_size, &flag);
   } /* NB: implicit taskwait at end of parallel region */
   node->nelim = (flag!=-1) ? flag+1 : n;
   if(flag!=-1) throw NotPosDefError(flag);

   /* Record information */
   node->ndelay_out = 0;
}
/* Factorize a node (wrapper) */
template <bool posdef, typename T, int BLOCK_SIZE>
void factor_node(
      int ni,
      struct cpu_node_data<T> *const node,
      const struct cpu_factor_options *const options,
      struct cpu_factor_stats *const stats
      ) {
   if(posdef) factor_node_posdef<T, BLOCK_SIZE>(node, options);
   else       factor_node_indef <T, BLOCK_SIZE>(ni, node, options, stats);
}

/* Calculate update */
template <bool posdef, typename T, typename StackAllocator>
void calculate_update(
      struct cpu_node_data<T> *node,
      StackAllocator *stalloc_odd,
      StackAllocator *stalloc_even,
      Workspace<T> *work
      ) {
   // Check there is work to do
   int m = node->nrow_expected - node->ncol_expected;
   int n = node->nelim;
   if(n==0 && !node->first_child) {
      // If everything is delayed, and no contribution from children then
      // free contrib memory and mark as no contribution for parent's assembly
      // FIXME: actually loop over children and check one exists with contriub
      //        rather than current approach of just looking for children.
      if(node->even) {
         stalloc_even->free(node->contrib, m*m*sizeof(T));
      } else {
         stalloc_odd->free(node->contrib, m*m*sizeof(T));
      }
      node->contrib = NULL;
      return;
   }
   if(m==0 || n==0) return; // no-op

   if(posdef) {
      int ldl = node->nrow_expected;
      host_syrk<T>(FILL_MODE_LWR, OP_N, m, n,
            -1.0, &node->lcol[node->ncol_expected], ldl,
            1.0, node->contrib, m);
   } else {
      // Indefinte - need to recalculate LD before we can use it!

      // Calculate LD
      T *lcol = &node->lcol[node->ncol_expected+node->ndelay_in];
      int ldl = node->nrow_expected + node->ndelay_in;
      T *d = &node->lcol[ldl*(node->ncol_expected+node->ndelay_in)];
      work->ensure_length(m*n);
      T *ld = work->mem;
      for(int j=0; j<n;) {
         if(d[2*j+1] == 0.0) {
            // 1x1 pivot
            // (Actually stored as D^-1 so need to invert it again)
            if(d[2*j] == 0.0) {
               // Handle zero pivots with care
               for(int i=0; i<m; i++) {
                  ld[j*m+i] = 0.0;
               }
            } else {
               // Standard 1x1 pivot
               T d11 = 1/d[2*j];
               // And calulate ld
               for(int i=0; i<m; i++) {
                  ld[j*m+i] = d11*lcol[j*ldl+i];
               }
            }
            // Increment j
            j++;
         } else {
            // 2x2 pivot
            // (Actually stored as D^-1 so need to invert it again)
            T di11 = d[2*j]; T di21 = d[2*j+1]; T di22 = d[2*j+3];
            T det = di11*di22 - di21*di21;
            T d11 = di22 / det; T d21 = -di21 / det; T d22 = di11 / det;
            // And calulate ld
            for(int i=0; i<m; i++) {
               ld[j*m+i]     = d11*lcol[j*ldl+i] + d21*lcol[(j+1)*ldl+i];
               ld[(j+1)*m+i] = d21*lcol[j*ldl+i] + d22*lcol[(j+1)*ldl+i];
            }
            // Increment j
            j += 2;
         }
      }

      // Apply update to contrib block
      host_gemm<T>(OP_N, OP_T, m, m, n,
            -1.0, lcol, ldl, ld, m,
            1.0, node->contrib, m);
   }

   // FIXME: debug remove
   /*printf("Contrib = \n");
   for(int i=0; i<m; i++) {
      for(int j=0; j<m; j++) printf(" %e", node->contrib[j*m+i]);
      printf("\n");
   }*/
}

}}} /* end of namespace spral::ssids::cpu */
