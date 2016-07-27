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
#include "NumericNode.hxx"
#include "SymbolicSubtree.hxx"
#include "SmallLeafNumericSubtree.hxx"

#include "profile.hxx"

namespace spral { namespace ssids { namespace cpu {

/** Class representing a subtree factorized on the CPU */
template <bool posdef,
          typename T,
          size_t PAGE_SIZE,
          typename FactorAllocator, // Allocator to use for factor storage
          typename ContribAllocator // Allocator to use for contribution blocks
          >
class NumericSubtree {
   typedef SmallLeafNumericSubtree<T, FactorAllocator, ContribAllocator> SLNS;
public:
   /* Delete copy constructors for safety re allocated memory */
   NumericSubtree(const NumericSubtree&) =delete;
   NumericSubtree& operator=(const NumericSubtree&) =delete;
   /** Performs factorization */
   NumericSubtree(
         SymbolicSubtree const& symbolic_subtree,
         T const* aval,
         T const* scaling,
         struct cpu_factor_options const& options,
         struct cpu_factor_stats& stats)
   : symb_(symbolic_subtree), nodes_(symbolic_subtree.nnodes_+1),
     small_leafs_(static_cast<SLNS*>(::operator new[](symb_.small_leafs_.size()*sizeof(SLNS)))),
     factor_alloc_(symbolic_subtree.get_factor_mem_est(options.multiplier))
   {
      Profile::init();

      /* Associate symbolic nodes to numeric ones; copy tree structure */
      for(int ni=0; ni<symb_.nnodes_+1; ++ni) {
         nodes_[ni].symb = &symbolic_subtree[ni];
         auto* fc = symbolic_subtree[ni].first_child;
         nodes_[ni].first_child = fc ? &nodes_[fc->idx] : nullptr;
         auto* nc = symbolic_subtree[ni].next_child;
         nodes_[ni].next_child = nc ? &nodes_[nc->idx] :  nullptr;
      }

      /* Allocate workspaces */
      int num_threads = omp_get_max_threads();
      std::vector<struct cpu_factor_stats> thread_stats(num_threads);
      std::vector<Workspace> work;
      work.reserve(num_threads);
      for(int i=0; i<num_threads; ++i)
         work.emplace_back(PAGE_SIZE);

      #pragma omp parallel default(shared)
      {
         /* Initalise stats */
         int this_thread = omp_get_thread_num();
         thread_stats[this_thread].flag = SSIDS_SUCCESS;
         thread_stats[this_thread].num_delay = 0;
         thread_stats[this_thread].num_neg = 0;
         thread_stats[this_thread].num_two = 0;
         thread_stats[this_thread].num_zero = 0;
         thread_stats[this_thread].maxfront = 0;
         thread_stats[this_thread].not_first_pass = 0;
         thread_stats[this_thread].not_second_pass = 0;

         #pragma omp taskgroup
         {
         /* Loop over small leaf subtrees */
         #pragma omp single
         if(posdef)
         for(unsigned int si=0; si<symb_.small_leafs_.size(); ++si) {
            auto* parent_lcol = &nodes_[symb_.small_leafs_[si].get_parent()];
            #pragma omp task default(none) \
               firstprivate(si) \
               shared(aval, options, scaling, thread_stats, work) \
               depend(inout: parent_lcol[0:1])
            {
               int this_thread = omp_get_thread_num();
#ifdef PROFILE
               Profile::Task task_subtree("TA_SUBTREE", this_thread);
#endif
               auto const& leaf = symb_.small_leafs_[si];
               new (&small_leafs_[si]) SLNS(leaf, nodes_, aval, scaling,
                     factor_alloc_, contrib_alloc_, work[this_thread],
                     options, thread_stats[this_thread]);
               #pragma omp cancel taskgroup \
                  if(thread_stats[this_thread].flag<SSIDS_SUCCESS)
#ifdef profile
               task_subtree.done();
#endif
            }
         }
         #pragma omp taskwait // FIXME: remove

         /* Loop over singleton nodes in order */
         #pragma omp single
         for(int ni=0; ni<symb_.nnodes_; ++ni) {
            if(posdef && symb_[ni].insmallleaf) continue; // already handled
            // FIXME: depending inout on parent is not a good way to represent
            //        the dependency.
            auto* this_lcol = &nodes_[ni]; // for depend
            auto* parent_lcol = &nodes_[symb_[ni].parent]; // for depend
            #pragma omp task default(none) \
               firstprivate(ni) \
               shared(aval, options, scaling, thread_stats, work) \
               depend(inout: this_lcol[0:1]) \
               depend(inout: parent_lcol[0:1])
            {
               /*printf("%d: Node %d parent %d (of %d) size %d x %d\n",
                     omp_get_thread_num(), ni, symb_[ni].parent, symb_.nnodes_,
                     symb_[ni].nrow, symb_[ni].ncol);*/
               int this_thread = omp_get_thread_num();
               // Assembly
#ifdef PROFILE
               Profile::Task task_assemble("TA_ASSEMBLE", this_thread);
#endif
               int* map = work[this_thread].get_ptr<int>(symb_.n+1);
               assemble_node
                  (posdef, ni, symb_[ni], &nodes_[ni], factor_alloc_,
                   contrib_alloc_, map, aval, scaling);
#ifdef PROFILE
               task_assemble.done();
               Profile::Task task_factor("TA_FACTOR", this_thread);
#endif
               // Update stats
               int nrow = symb_[ni].nrow + nodes_[ni].ndelay_in;
               thread_stats[this_thread].maxfront = std::max(thread_stats[this_thread].maxfront, nrow);
               // Factorization
               factor_node
                  <posdef>
                  (ni, symb_[ni], &nodes_[ni], options,
                   thread_stats[this_thread], contrib_alloc_);
               #pragma omp cancel taskgroup \
                  if(thread_stats[this_thread].flag<SSIDS_SUCCESS)
            }
         }
         } // taskgroup

         // Reduce thread_stats
         #pragma omp single 
         {
            stats.flag = SSIDS_SUCCESS;
            stats.num_delay = 0;
            stats.num_neg = 0;
            stats.num_two = 0;
            stats.num_zero = 0;
            stats.maxfront = 0;
            stats.not_first_pass = 0;
            stats.not_second_pass = 0;
            for(auto tstats : thread_stats) {
               stats.flag =
                  (stats.flag == SSIDS_SUCCESS) ? tstats.flag
                                                : std::min(stats.flag, tstats.flag); 
               stats.num_delay += tstats.num_delay;
               stats.num_neg += tstats.num_neg;
               stats.num_two += tstats.num_two;
               stats.num_zero += tstats.num_zero;
               stats.maxfront = std::max(stats.maxfront, tstats.maxfront);
               stats.not_first_pass += tstats.not_first_pass;
               stats.not_second_pass += tstats.not_second_pass;
            }
         }
      }

      // Count stats
      // FIXME: Do this as we go along...
      if(posdef) {
         // no real changes to stats from zero initialization
      } else { // indefinite
         for(int ni=0; ni<symb_.nnodes_; ni++) {
            int m = symb_[ni].nrow + nodes_[ni].ndelay_in;
            int n = symb_[ni].ncol + nodes_[ni].ndelay_in;
            int ldl = align_lda<T>(m);
            T *d = nodes_[ni].lcol + n*ldl;
            for(int i=0; i<nodes_[ni].nelim; ) {
               T a11 = d[2*i];
               T a21 = d[2*i+1];
               if(i+1==nodes_[ni].nelim || std::isfinite(d[2*i+2])) {
                  // 1x1 pivot (or zero)
                  if(a11 == 0.0) stats.num_zero++;
                  if(a11 < 0.0) stats.num_neg++;
                  i++;
               } else {
                  // 2x2 pivot
                  T a22 = d[2*(i+1)];
                  stats.num_two++;
                  T det = a11*a22 - a21*a21; // product of evals
                  T trace = a11 + a22; // sum of evals
                  if(det < 0) stats.num_neg++;
                  else if(trace < 0) stats.num_neg+=2;
                  i+=2;
               }
            }
         }
      }
      Profile::end();
   }
   ~NumericSubtree() {
      delete[] small_leafs_;
   }

   void solve_fwd(int nrhs, double* x, int ldx) const {
      /* Allocate memory */
      double* xlocal = new double[nrhs*symb_.n];
      int* map_alloc = (!posdef) ? new int[symb_.n] : nullptr; // only indef

      /* Main loop */
      for(int ni=0; ni<symb_.nnodes_; ++ni) {
         int m = symb_[ni].nrow;
         int n = symb_[ni].ncol;
         int nelim = (posdef) ? n
                              : nodes_[ni].nelim;
         int ndin = (posdef) ? 0
                             : nodes_[ni].ndelay_in;
         int ldl = align_lda<T>(m+ndin);

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
            cholesky_solve_fwd(m, n, nodes_[ni].lcol, ldl, nrhs, xlocal, symb_.n);
         } else { /* indef */
            ldlt_app_solve_fwd(m+ndin, nelim, nodes_[ni].lcol, ldl, nrhs,
                  xlocal, symb_.n);
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

   template <bool do_diag, bool do_bwd>
   void solve_diag_bwd_inner(int nrhs, double* x, int ldx) const {
      if(posdef && !do_bwd) return; // diagonal solve is a no-op for posdef

      /* Allocate memory - map only needed for indef bwd/diag_bwd solve */
      double* xlocal = new double[nrhs*symb_.n];
      int* map_alloc = (!posdef && do_bwd) ? new int[symb_.n]
                                           : nullptr;

      /* Perform solve */
      for(int ni=symb_.nnodes_-1; ni>=0; --ni) {
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
            if(do_bwd) {
               for(int i=0; i<n+ndin; ++i)
                  map_alloc[i] = nodes_[ni].perm[i];
               for(int i=n; i<m; ++i)
                  map_alloc[i+ndin] = symb_[ni].rlist[i];
               map = map_alloc;
            } else { // if only doing diagonal, only need first nelim<=n+ndin
               map = nodes_[ni].perm;
            }
         } else {
            // posdef there is no permutation
            map = symb_[ni].rlist;
         }

         /* Gather into dense vector xlocal */
         int blkm = (do_bwd) ? m+ndin
                             : nelim;
         int ldl = align_lda<T>(blkm);
         for(int r=0; r<nrhs; ++r)
         for(int i=0; i<blkm; ++i)
            xlocal[r*symb_.n+i] = x[r*ldx + map[i]-1];

         /* Perform dense solve */
         if(posdef) {
            cholesky_solve_bwd(m, n, nodes_[ni].lcol, ldl, nrhs, xlocal, symb_.n);
         } else {
            if(do_diag) ldlt_app_solve_diag(
                  nelim, &nodes_[ni].lcol[(n+ndin)*ldl], xlocal
                  );
            if(do_bwd) ldlt_app_solve_bwd(
                  m+ndin, nelim, nodes_[ni].lcol, ldl, nrhs, xlocal, symb_.n
                  );
         }

         /* Scatter result (only first nelim entries have changed) */
         for(int r=0; r<nrhs; ++r)
         for(int i=0; i<nelim; ++i)
            x[r*ldx + map[i]-1] = xlocal[r*symb_.n+i];
      }

      /* Cleanup memory */
      if(!posdef && do_bwd) delete[] map_alloc; // only used in indef case
      delete[] xlocal;
   }

   void solve_diag(int nrhs, double* x, int ldx) const {
      solve_diag_bwd_inner<true, false>(nrhs, x, ldx);
   }

   void solve_diag_bwd(int nrhs, double* x, int ldx) const {
      solve_diag_bwd_inner<true, true>(nrhs, x, ldx);
   }

   void solve_bwd(int nrhs, double* x, int ldx) const {
      solve_diag_bwd_inner<false, true>(nrhs, x, ldx);
   }

   /** Returns information on diagonal entries and/or pivot order.
    * Note that piv_order is only set in indefinite case.
    * One of piv_order or d may be null in indefinite case.
    */
   void enquire(int *piv_order, double* d) const {
      if(posdef) {
         for(int ni=0; ni<symb_.nnodes_; ++ni) {
            int blkm = symb_[ni].nrow;
            int nelim = symb_[ni].ncol;
            int ldl = align_lda<T>(blkm);
            for(int i=0; i<nelim; ++i)
               *(d++) = nodes_[ni].lcol[i*(ldl+1)];
         }
      } else { /*indef*/
         for(int ni=0, piv=0; ni<symb_.nnodes_; ++ni) {
            int blkm = symb_[ni].nrow + nodes_[ni].ndelay_in;
            int blkn = symb_[ni].ncol + nodes_[ni].ndelay_in;
            int ldl = align_lda<T>(blkm);
            int nelim = nodes_[ni].nelim;
            double const* dptr = &nodes_[ni].lcol[blkn*ldl];
            for(int i=0; i<nelim; ) {
               if(i+1==nelim || std::isfinite(dptr[2*i+2])) {
                  /* 1x1 pivot */
                  if(piv_order) {
                     piv_order[nodes_[ni].perm[i]-1] = (piv++);
                  }
                  if(d) {
                     *(d++) = dptr[2*i+0];
                     *(d++) = 0.0;
                  }
                  i+=1;
               } else {
                  /* 2x2 pivot */
                  if(piv_order) {
                     piv_order[nodes_[ni].perm[i]-1] = -(piv++);
                     piv_order[nodes_[ni].perm[i+1]-1] = -(piv++);
                  }
                  if(d) {
                     *(d++) = dptr[2*i+0];
                     *(d++) = dptr[2*i+1];
                     *(d++) = dptr[2*i+2];
                     *(d++) = 0.0;
                  }
                  i+=2;
               }
            }
         }
      }
   }

   /** Allows user to alter D values, indef case only. */
   void alter(double const* d) {
      for(int ni=0; ni<symb_.nnodes_; ++ni) {
         int blkm = symb_[ni].nrow + nodes_[ni].ndelay_in;
         int blkn = symb_[ni].ncol + nodes_[ni].ndelay_in;
         int ldl = align_lda<T>(blkm);
         int nelim = nodes_[ni].nelim;
         double* dptr = &nodes_[ni].lcol[blkn*ldl];
         for(int i=0; i<nelim; ++i) {
            dptr[2*i+0] = *(d++);
            dptr[2*i+1] = *(d++);
         }
      }
   }

	void print() const {
		for(int node=0; node<symb_.nnodes_; node++) {
			printf("== Node %d ==\n", node);
			int m = symb_[node].nrow + nodes_[node].ndelay_in;
			int n = symb_[node].ncol + nodes_[node].ndelay_in;
         int ldl = align_lda<T>(m);
         int nelim = nodes_[node].nelim;
			int const* rlist = &symb_[node].rlist[ symb_[node].ncol ];
			for(int i=0; i<m; ++i) {
				if(i<n) printf("%d%s:", nodes_[node].perm[i], (i<nelim)?"X":"D");
				else    printf("%d:", rlist[i-n]);
				for(int j=0; j<n; j++) printf(" %10.2e", nodes_[node].lcol[j*ldl+i]);
            T const* d = &nodes_[node].lcol[n*ldl];
				if(!posdef && i<nelim)
               printf("  d: %10.2e %10.2e", d[2*i+0], d[2*i+1]);
		      printf("\n");
			}
		}
	}

   SymbolicSubtree const& get_symbolic_subtree() { return symb_; }

private:
   SymbolicSubtree const& symb_;
   std::vector<NumericNode<T>> nodes_;
   SLNS *small_leafs_; // Apparently emplace_back isn't threadsafe, so
      // std::vector is out. So we use placement new instead.
   FactorAllocator factor_alloc_;
   ContribAllocator contrib_alloc_;
};

}}} /* end of namespace spral::ssids::cpu */
