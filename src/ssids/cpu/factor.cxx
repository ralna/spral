/* Copyright 2014 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 */

/* Standard headers */
#include <atomic>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <limits> // FIXME: remove when done with debug if unneeded
#include <new>
#include <queue>
#include <sstream>
#include <stdexcept>
/* SPRAL headers */
#include "factor_iface.h"
#include "AlignedAllocator.hxx"
#include "kernels/CpuLDLT.cxx"
#include "kernels/cholesky.hxx"

/////////////////////////////////////////////////////////////////////////////
// start of namespace spral::ssids::cpu
namespace spral { namespace ssids { namespace cpu {

const int SSIDS_SUCCESS = 0;
const int SSIDS_ERROR_NOT_POS_DEF = -6;

/** Class for stack-based allocation.
 *
 * Designed to make quick allocation and deallocation of fronts in a stack
 * based fashion more efficient
 */
template<size_t PAGE_SIZE>
class StackAllocator {
private:
   class Page {
   private:
      char mem[PAGE_SIZE];
   public:
      Page *prev, *next;
      size_t top;
      Page(Page *prev)
         : prev(prev), next(NULL), top(0)
         {}
      void *alloc(size_t size) {
         if(top+size > PAGE_SIZE) return NULL;
         void *ptr = mem + top;
         top += size;
         return ptr;
      }
      int free(size_t size) {
         top -= size;
         return top;
      }
   };
   Page *current, *last;
public:
   StackAllocator(void)
      : current(new Page(NULL)), last(current)
      {}
   ~StackAllocator() {
      // Work backwards freeing Pages
      while(last) {
         Page *prev = last->prev;
         delete last;
         last = prev;
      }
   }
   
   void *alloc(size_t size) {
      if(size > PAGE_SIZE/2) {
         // Requesting a large block, use standard allocator
         return new char[size];
      } else {
         // Requesting a small block
         void *ptr = current->alloc(size);
         if(!ptr) {
            // Ran out of space on page, find a new one
            if(current->next) {
               // Reuse existing page
               current = current->next;
            } else {
               // Allocate a new one
               last = new Page(current);
               current->next = last;
               current = last;
            }
            // Allocation will defintely work this time
            ptr = current->alloc(size);
         }
         // Return allocation
         return ptr;
      }
   }

   void free(void *ptr, size_t size) {
      if(size > PAGE_SIZE/2) {
         // Was directly allocated
         delete[] (char *) ptr;
      } else {
         // Call page's free function
         int remain = current->free(size);
         if(remain==0 && current->prev) {
            // Page now clear (and not first page), iterate backwards
            current = current->prev;
         }
      }
   }
};

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

/* Generic wrapper around Fortran-defined smalloc calls */
template<typename T>
T *smalloc(void *alloc, size_t len);
template<>
double *smalloc(void *alloc, size_t len) {
   return spral_ssids_smalloc_dbl(alloc, len);
}
template<>
int *smalloc(void *alloc, size_t len) {
   return spral_ssids_smalloc_int(alloc, len);
}

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

template <typename T>
struct cpu_node_data {
   /* Fixed data from analyse */
   const int nrow_expected; // Number of rows without delays
   const int ncol_expected; // Number of cols without delays
   struct cpu_node_data<T> *const first_child; // Pointer to our first child
   struct cpu_node_data<T> *const next_child; // Pointer to parent's next child
   const int *const rlist; // Pointer to row lists
   const bool even; // Indicates which stack we're using (odd or even distance
                    // from root)

   /* Data about A:
    * aval[i] goes into lcol[ amap[i] ] if there are no delays
    */
   int num_a; // Number of entries from A
   const int *const amap; // Map from A to node (length 2*num_a)

   /* Data that changes during factorize */
   int ndelay_in; // Number of delays arising from children
   int ndelay_out; // Number of delays arising to push into parent
   int nelim; // Number of columns succesfully eliminated
   T *lcol; // Pointer to start of factor data
   int *perm; // Pointer to permutation
   T *contrib; // Pointer to contribution block
};

struct cpu_factor_options {
   double small;
   double u;
   int print_level;
   int cpu_small_subtree_threshold;
   int cpu_task_block_size;
};

struct cpu_factor_stats {
   int flag;
   int num_delay;
   int num_neg;
   int num_two;
   int num_zero;
   int maxfront;
   int elim_at_pass[5];
   int elim_at_itr[5];
};

template <typename T,
          size_t PAGE_SIZE>
void assemble_node(
      bool posdef,
      int ni, // FIXME: remove with debug
      struct cpu_node_data<T> *const node,
      void *const alloc,
      StackAllocator<PAGE_SIZE> *stalloc_odd,
      StackAllocator<PAGE_SIZE> *stalloc_even,
      int *const map,
      const T *const aval,
      const T *const scaling
      ) {
   /* Count incoming delays and determine size of node */
   node->ndelay_in = 0;
   for(struct cpu_node_data<T> *child=node->first_child; child!=NULL; child=child->next_child) {
      node->ndelay_in += child->ndelay_out;
   }
   int nrow = node->nrow_expected + node->ndelay_in;
   int ncol = node->ncol_expected + node->ndelay_in;

   /* Get space for node now we know it size using Fortran allocator + zero it*/
   // NB L is  nrow x ncol and D is 2 x ncol (but no D if posdef)
   size_t len = posdef ? ((size_t) nrow  ) * ncol  // posdef
                       : ((size_t) nrow+2) * ncol; // indef (includes D)
   if(!(node->lcol = smalloc<T>(alloc, len)))
      throw std::bad_alloc();
   memset(node->lcol, 0, len*sizeof(T));

   /* Get space for contribution block + zero it */
   long contrib_dimn = node->nrow_expected - node->ncol_expected;
   if(node->even) {
      node->contrib = (contrib_dimn > 0) ? (T *) stalloc_even->alloc(contrib_dimn*contrib_dimn*sizeof(T)) : NULL;
   } else {
      node->contrib = (contrib_dimn > 0) ? (T *) stalloc_odd->alloc(contrib_dimn*contrib_dimn*sizeof(T)) : NULL;
   }
   if(node->contrib)
      memset(node->contrib, 0, contrib_dimn*contrib_dimn*sizeof(T));

   /* Alloc + set perm for expected eliminations at this node (delays are set
    * when they are imported from children) */
   node->perm = smalloc<int>(alloc, ncol); // ncol fully summed variables
   for(int i=0; i<node->ncol_expected; i++)
      node->perm[i] = node->rlist[i];

   /* Add A */
   if(scaling) {
      /* Scaling to apply */
      for(int i=0; i<node->num_a; i++) {
         long src  = node->amap[2*i+0] - 1; // amap contains 1-based values
         long dest = node->amap[2*i+1] - 1; // amap contains 1-based values
         int c = dest / node->nrow_expected;
         int r = dest % node->nrow_expected;
         long k = c*nrow + r;
         if(r >= node->ncol_expected) k += node->ndelay_in;
         T rscale = scaling[ node->rlist[r]-1 ];
         T cscale = scaling[ node->rlist[c]-1 ];
         node->lcol[k] = rscale * aval[src] * cscale;
      }
   } else {
      /* No scaling to apply */
      for(int i=0; i<node->num_a; i++) {
         long src  = node->amap[2*i+0] - 1; // amap contains 1-based values
         long dest = node->amap[2*i+1] - 1; // amap contains 1-based values
         int c = dest / node->nrow_expected;
         int r = dest % node->nrow_expected;
         long k = c*nrow + r;
         if(r >= node->ncol_expected) k += node->ndelay_in;
         node->lcol[k] = aval[src];
      }
   }

   /* Add children */
   if(node->first_child != NULL) {
      /* Build lookup vector, allowing for insertion of delayed vars */
      /* Note that while rlist[] is 1-indexed this is fine so long as lookup
       * is also 1-indexed (which it is as it is another node's rlist[] */
      for(int i=0; i<node->ncol_expected; i++)
         map[ node->rlist[i] ] = i;
      for(int i=node->ncol_expected; i<node->nrow_expected; i++)
         map[ node->rlist[i] ] = i + node->ndelay_in;
      /* Loop over children adding contributions */
      int delay_col = node->ncol_expected;
      for(struct cpu_node_data<T> *child=node->first_child; child!=NULL; child=child->next_child) {
         /* Handle delays - go to back of node
          * (i.e. become the last rows as in lower triangular format) */
         for(int i=0; i<child->ndelay_out; i++) {
            // Add delayed rows (from delayed cols)
            T *dest = &node->lcol[delay_col*(nrow+1)];
            int lds = child->nrow_expected + child->ndelay_in;
            T *src = &child->lcol[(child->nelim+i)*(lds+1)];
            node->perm[delay_col] = child->perm[child->nelim+i];
            for(int j=0; j<child->ndelay_out-i; j++) {
               dest[j] = src[j];
            }
            // Add child's non-fully summed rows (from delayed cols)
            dest = node->lcol;
            src = &child->lcol[child->nelim*lds + child->ndelay_in +i*lds];
            for(int j=child->ncol_expected; j<child->nrow_expected; j++) {
               int r = map[ child->rlist[j] ];
               if(r < ncol) dest[r*nrow+delay_col] = src[j];
               else         dest[delay_col*nrow+r] = src[j];
            }
            delay_col++;
         }

         /* Handle expected contributions (only if something there) */
         if(child->contrib) {
            int cm = child->nrow_expected - child->ncol_expected;
            for(int i=0; i<cm; i++) {
               int c = map[ child->rlist[child->ncol_expected+i] ];
               T *src = &child->contrib[i*cm];
               if(c < node->ncol_expected) {
                  // Contribution added to lcol
                  int ldd = nrow;
                  T *dest = &node->lcol[c*ldd];
                  for(int j=i; j<cm; j++) {
                     int r = map[ child->rlist[child->ncol_expected+j] ];
                     dest[r] += src[j];
                  }
               } else {
                  // Contribution added to contrib
                  // FIXME: Add after contribution block established?
                  int ldd = node->nrow_expected - node->ncol_expected;
                  T *dest = &node->contrib[(c-ncol)*ldd];
                  for(int j=i; j<cm; j++) {
                     int r = map[ child->rlist[child->ncol_expected+j] ] - ncol;
                     dest[r] += src[j];
                  }
               }
            }
            /* Free memory from child contribution block */
            if(child->even) {
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
   int ldd = node->nrow_expected - node->ncol_expected;
   for(int i=0; i<ldd; i++) {
      for(int j=0; j<ldd; j++) printf(" %e", node->contrib[j*ldd+i]);
      printf("\n");
   }*/
}
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

/* FIXME: remove post debug */
template<typename T>
void print_node(bool posdef, int m, int n, int nelim, const int *perm, const int *rlist, const T *lcol, const T*d) {
   for(int i=0; i<m; i++) {
      if(i<n) printf("%d%s:", perm[i], (i<nelim)?"X":"D");
      else    printf("%d:", rlist[i-n]);
      for(int j=0; j<n; j++) printf(" %10.2e", lcol[j*m+i]);
      if(!posdef && i<nelim) printf("  d: %10.2e %10.2e\n", d[2*i+0], d[2*i+1]);
      else printf("\n");
   }
}
/* FIXME: remove post debug */
template<typename T>
void print_factors(
      bool posdef,
      int nnodes,
      struct cpu_node_data<T> *const nodes
      ) {
   for(int node=0; node<nnodes; node++) {
      printf("== Node %d ==\n", node);
      int m = nodes[node].nrow_expected + nodes[node].ndelay_in;
      int n = nodes[node].ncol_expected + nodes[node].ndelay_in;
      const int *rptr = &nodes[node].rlist[ nodes[node].ncol_expected ];
      print_node(posdef, m, n, nodes[node].nelim, nodes[node].perm,
            rptr, nodes[node].lcol, &nodes[node].lcol[m*n]);
   }
}

/* Calculate update */
template <bool posdef, typename T, size_t PAGE_SIZE>
void calculate_update(
      struct cpu_node_data<T> *node,
      StackAllocator<PAGE_SIZE> *stalloc_odd,
      StackAllocator<PAGE_SIZE> *stalloc_even,
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

/* Simplistic multifrontal factorization */
template <bool posdef,
          typename T,
          size_t BLOCK_SIZE,
          size_t PAGE_SIZE,
          bool timing
          >
void factor(
      int n,            // Maximum row index (+1)
      int nnodes,       // Number of nodes in assembly tree
      struct cpu_node_data<T> *const nodes, // Data structure for node information
      const T *const aval, // Values of A
      const T *const scaling, // Scaling vector (NULL if no scaling)
      void *const alloc,      // Fortran allocator pointer
      const struct cpu_factor_options *const options, // Options in
      struct cpu_factor_stats *const stats // Info out
      ) {

   // Allocate workspaces
   StackAllocator<PAGE_SIZE> stalloc_odd;
   StackAllocator<PAGE_SIZE> stalloc_even;
   Workspace<T> work(PAGE_SIZE);
   int *map = new int[n+1]; // +1 to allow for indexing with 1-indexed array

   // Initialize statistics
   stats->num_delay = 0;
   stats->num_neg = 0;
   stats->num_two = 0;
   stats->num_zero = 0;
   stats->maxfront = 0;
   for(int i=0; i<5; i++) stats->elim_at_pass[i] = 0;
   for(int i=0; i<5; i++) stats->elim_at_itr[i] = 0;

   /* Main loop: Iterate over nodes in order */
   for(int ni=0; ni<nnodes; ++ni) {
      // Assembly
      assemble_node
         <T, PAGE_SIZE>
         (posdef, ni, &nodes[ni], alloc, &stalloc_odd, &stalloc_even, map, aval,
          scaling);
      // Update stats
      int nrow = nodes[ni].nrow_expected + nodes[ni].ndelay_in;
      stats->maxfront = std::max(stats->maxfront, nrow);
      // Factorization
      factor_node
         <posdef, T, BLOCK_SIZE>
         (ni, &nodes[ni], options, stats);
      // Form update
      calculate_update
         <posdef, T, PAGE_SIZE>
         (&nodes[ni], &stalloc_odd, &stalloc_even, &work);
   }

   // Count stats
   // FIXME: gross hack for compat with bub (which needs to differentiate
   // between a natural zero and a 2x2 factor's second entry without counting)
   // SSIDS original data format [a11 a21 a22 xxx] seems more bizzare than
   // bub one [a11 a21 inf a22]
   if(posdef) {
      // no real changes to stats from zero initialization
   } else { // indefinite
      for(int ni=0; ni<nnodes; ni++) {
         int m = nodes[ni].nrow_expected + nodes[ni].ndelay_in;
         int n = nodes[ni].ncol_expected + nodes[ni].ndelay_in;
         T *d = nodes[ni].lcol + m*n;
         for(int i=0; i<nodes[ni].nelim; i++)
            if(d[2*i] == std::numeric_limits<T>::infinity())
               d[2*i] = d[2*i+1];
         for(int i=0; i<nodes[ni].nelim; ) {
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

   /* Free memory */
   delete[] map;

}


}}} /* end of namespace spral::ssids::cpu */
//////////////////////////////////////////////////////////////////////////

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_factor_cpu_dbl(
      bool posdef,     // If true, performs A=LL^T, if false do pivoted A=LDL^T
      int n,            // Maximum row index (+1)
      int nnodes,       // Number of nodes in assembly tree
      struct spral::ssids::cpu::cpu_node_data<double> *const nodes, // Data structure for node information
      const double *const aval, // Values of A
      const double *const scaling, // Scaling vector (NULL if none)
      void *const alloc,      // Pointer to Fortran allocator structure
      const struct spral::ssids::cpu::cpu_factor_options *const options, // Options in
      struct spral::ssids::cpu::cpu_factor_stats *const stats // Info out
      ) {

   // Initialize stats
   stats->flag = spral::ssids::cpu::SSIDS_SUCCESS;

   // Call relevant routine
   if(posdef) {
      try {
         spral::ssids::cpu::factor<true, double, 16, 16384, false>
            (n, nnodes, nodes, aval, scaling, alloc, options, stats);
      } catch(spral::ssids::cpu::NotPosDefError npde) {
         stats->flag = spral::ssids::cpu::SSIDS_ERROR_NOT_POS_DEF;
      }
   } else {
      spral::ssids::cpu::factor<false, double, 16, 16384, false>
         (n, nnodes, nodes, aval, scaling, alloc, options, stats);
   }

   // FIXME: Remove when done with debug
   if(options->print_level > 9999) {
      printf("Final factors:\n");
      spral::ssids::cpu::print_factors<double>(posdef, nnodes, nodes);
   }
}
