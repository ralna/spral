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

#include "cpu_iface.hxx"
#include "factor.hxx"
#include "SymbolicSubtree.hxx"

namespace spral { namespace ssids { namespace cpu {

/** Class representing a subtree factorized on the CPU */
template <bool posdef, size_t BLOCK_SIZE, typename T, typename StackAllocator>
class NumericSubtree {
public:
   /** Constructor ties into symbolic_subtree */
   NumericSubtree(SymbolicSubtree const& symbolic_subtree, int nnodes, struct cpu_node_data<T>* nodes)
   : nnodes_(nnodes), nodes_(nodes), symb_(symbolic_subtree)
   {
      /* Associate symbolic nodes to numeric ones; copy tree structure */
      for(int ni=0; ni<nnodes_+1; ++ni) {
         nodes_[ni].symb = &symbolic_subtree[ni];
         auto* fc = symbolic_subtree[ni].first_child;
         nodes_[ni].first_child = fc ? &nodes_[fc->idx] : nullptr;
         auto* nc = symbolic_subtree[ni].next_child;
         nodes_[ni].next_child = nc ? &nodes_[nc->idx] :  nullptr;
      }
   }
   void factor(T const* aval, T const* scaling, void *const alloc, StackAllocator &stalloc_odd, StackAllocator &stalloc_even, Workspace<T> &work, int* map, struct cpu_factor_options const* options, struct cpu_factor_stats* stats) {

      /* Main loop: Iterate over nodes in order */
      for(int ni=0; ni<nnodes_; ++ni) {
         // Assembly
         assemble_node
            (posdef, ni, symb_[ni], &nodes_[ni], alloc, &stalloc_odd,
             &stalloc_even, map, aval, scaling);
         // Update stats
         int nrow = symb_[ni].nrow + nodes_[ni].ndelay_in;
         stats->maxfront = std::max(stats->maxfront, nrow);
         // Factorization
         factor_node
            <posdef, BLOCK_SIZE>
            (ni, symb_[ni], &nodes_[ni], options, stats);
         // Form update
         calculate_update<posdef>
            (symb_[ni], &nodes_[ni], &stalloc_odd, &stalloc_even, &work);
      }

      // Count stats
      // FIXME: gross hack for compat with bub (which needs to differentiate
      // between a natural zero and a 2x2 factor's second entry without counting)
      // SSIDS original data format [a11 a21 a22 xxx] seems more bizzare than
      // bub one [a11 a21 inf a22]
      if(posdef) {
         // no real changes to stats from zero initialization
      } else { // indefinite
         for(int ni=0; ni<nnodes_; ni++) {
            int m = symb_[ni].nrow + nodes_[ni].ndelay_in;
            int n = symb_[ni].ncol + nodes_[ni].ndelay_in;
            T *d = nodes_[ni].lcol + m*n;
            for(int i=0; i<nodes_[ni].nelim; i++)
               if(d[2*i] == std::numeric_limits<T>::infinity())
                  d[2*i] = d[2*i+1];
            for(int i=0; i<nodes_[ni].nelim; ) {
               T a11 = d[2*i];
               T a21 = d[2*i+1];
               if(a21 == 0.0) {
                  // 1x1 pivot (or zero)
                  if(a11 == 0.0) stats->num_zero++;
                  if(a11 < 0.0) stats->num_neg++;
                  i++;
               } else {
                  // 2x2 pivot
                  T a22 = d[2*(i+1)];
                  stats->num_two++;
                  T det = a11*a22 - a21*a21; // product of evals
                  T trace = a11 + a22; // sum of evals
                  if(det < 0) stats->num_neg++;
                  else if(trace < 0) stats->num_neg+=2;
                  i+=2;
               }
            }
         }
      }
   }

   void solve_fwd(int nrhs, double* x, int ldx) {
      /* Allocate memory */
      double* xlocal = new double[nrhs*symb_.n];
      int* map_alloc = (!posdef) ? new int[symb_.n] : nullptr; // only indef

      /* Main loop */
      for(int ni=0; ni<nnodes_; ++ni) {
         int m = symb_[ni].nrow;
         int n = symb_[ni].ncol;
         int nelim = (posdef) ? n
                              : nodes_[ni].nelim;
         int ndin = (posdef) ? 0
                             : nodes_[ni].ndelay_in;

         /* Build map (indef only) */
         int const *map;
         if(!posdef) {
            // indef need to allow for permutation and/or delays
            for(int i=0; i<n+ndin; ++i)
               map_alloc[i] = nodes_[ni].perm[i];
            for(int i=n; i<m; ++i)
               map_alloc[i+ndin] = symb_[ni].rlist[i];
            map = map_alloc;
         } else {
            // posdef there is no permutation
            map = symb_[ni].rlist;
         }

         /* Gather into dense vector xlocal */
         // FIXME: don't bother copying elements of x > m, just use beta=0
         //        in dgemm call and then add as we scatter
         for(int r=0; r<nrhs; ++r)
         for(int i=0; i<m+ndin; ++i)
            xlocal[r*symb_.n+i] = x[r*ldx + map[i]-1]; // Fortran indexed

         /* Perform dense solve */
         if(posdef) {
            cholesky_solve_fwd(m, n, nodes_[ni].lcol, m, nrhs, xlocal, symb_.n);
         } else { /* indef */
            ldlt_solve_fwd(m+ndin, nelim, nodes_[ni].lcol, m+ndin, nrhs, xlocal,
                  symb_.n);
         }

         /* Scatter result */
         for(int r=0; r<nrhs; ++r)
         for(int i=0; i<m+ndin; ++i)
            x[r*ldx + map[i]-1] = xlocal[r*symb_.n+i];
      }

      /* Cleanup memory */
      if(!posdef) delete[] map_alloc; // only used in indef case
      delete[] xlocal;
   }

   void solve_diag(int nrhs, double* x, int ldx) {
   }

   void solve_diag_bwd(int nrhs, double* x, int ldx) {
      /* Allocate memory */
      double* xlocal = new double[nrhs*symb_.n];
      int* map_alloc = (!posdef) ? new int[symb_.n] : nullptr; // only indef

      /* Perform solve */
      for(int ni=nnodes_-1; ni>=0; --ni) {
         int m = symb_[ni].nrow;
         int n = symb_[ni].ncol;
         int nelim = (posdef) ? n
                              : nodes_[ni].nelim;
         int ndin = (posdef) ? 0
                             : nodes_[ni].ndelay_in;

         /* Build map (indef only) */
         int const *map;
         if(!posdef) {
            // indef need to allow for permutation and/or delays
            for(int i=0; i<n+ndin; ++i)
               map_alloc[i] = nodes_[ni].perm[i];
            for(int i=n; i<m; ++i)
               map_alloc[i+ndin] = symb_[ni].rlist[i];
            map = map_alloc;
         } else {
            // posdef there is no permutation
            map = symb_[ni].rlist;
         }

         /* Gather into dense vector xlocal */
         for(int r=0; r<nrhs; ++r)
         for(int i=0; i<m+ndin; ++i)
            xlocal[r*symb_.n+i] = x[r*ldx + map[i]-1];

         /* Perform dense solve */
         if(posdef) {
            cholesky_solve_bwd(m, n, nodes_[ni].lcol, m, nrhs, xlocal, symb_.n);
         } else {
            ldlt_solve_diag(nelim, &nodes_[ni].lcol[(m+ndin)*(n+ndin)], xlocal);
            ldlt_solve_bwd(m+ndin, nelim, nodes_[ni].lcol, m+ndin, nrhs, xlocal,
                  symb_.n);
         }

         /* Scatter result (only first nelim entries have changed) */
         for(int r=0; r<nrhs; ++r)
         for(int i=0; i<nelim; ++i)
            x[r*ldx + map[i]-1] = xlocal[r*symb_.n+i];
      }

      /* Cleanup memory */
      if(!posdef) delete[] map_alloc; // only used in indef case
      delete[] xlocal;
   }

   void solve_bwd(int nrhs, double* x, int ldx) {
   }

   SymbolicSubtree const& get_symbolic_subtree() { return symb_; }
private:
   int nnodes_;
   struct cpu_node_data<T>* nodes_;
   SymbolicSubtree const& symb_;
};

}}} /* end of namespace spral::ssids::cpu */
