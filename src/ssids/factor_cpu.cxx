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
/* External library headers */
#include "bub/bub.hxx"
#include <hwloc.h>
#include <papi.h> // FIXME: remove or make dependency in autoconf
#include <pthread.h>
/* SPRAL headers */
#include "factor_cpu_iface.h"

/////////////////////////////////////////////////////////////////////////////
// start of namespace spral::ssids::internal
namespace spral { namespace ssids { namespace internal {

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
   typedef bub::CpuLDLT<T, BLOCK_SIZE> CpuLDLTSpec;
   typedef bub::CpuLDLT<T, BLOCK_SIZE, 5, true> CpuLDLTSpecDebug; // FIXME: debug remove
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
   typedef bub::CpuLLT<T, BLOCK_SIZE> CpuLLTSpec;
   //typedef bub::CpuLLT<T, 4, true> CpuLLTSpecDebug; //FIXME: remove
   int flag = CpuLLTSpec().factor(m, n, lcol, m);
   node->nelim = (flag) ? flag : n;
   if(flag) throw NotPosDefError(flag);

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
      host_syrk<T>(bub::FILL_MODE_LWR, bub::OP_N, m, n,
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
      host_gemm<T>(bub::OP_N, bub::OP_T, m, m, n,
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

/** Class representing a unit of work (i.e. a (potentially splitable) task) */
class WorkUnit {
public:
   /** Returns estimated execution time in ms. A value of -1 indicates no idea */
   virtual long estTime(void)=0;
   /** Signal work unit that it should try and split itself to produce more work */
   virtual void split(void)=0;
   /** Start performing work */
   virtual void exec(void)=0;
};

/** Class representing a memory resource. Acts as an address, data only
 * accesible once actualise() method has been called. */
class MemResource {
   /** Returns true if memory is available locally (i.e. actualise() called) */
   virtual bool isActualised(void)=0;
   /** Called to bring current copy of data into local pool */
   virtual void actualise(void)=0;
};

template <bool posdef,
          typename T,
          size_t BLOCK_SIZE,
          int PAGE_SIZE,
          bool timing
          >
class NodeWork : public WorkUnit {
public: // FIXME: do we care enough to make this private?
   /** Node index */
   int index;
   /** Pointer to Fortran data structure */
   struct cpu_node_data<T> *data;
   /** Pointer to options structure */
   const struct cpu_factor_options *const options;
   /** Pointer to Fortran allocator */
   void *const alloc;
   /** Pointer to values of A */
   const T *const aval;
   /** Pointer to scaling arrays */
   const T *const scaling;

   /* Pointers for executiong unit specific resources */
   int *map;
   StackAllocator<PAGE_SIZE> *stalloc_odd;
   StackAllocator<PAGE_SIZE> *stalloc_even;
   Workspace<T> *work;
   struct cpu_factor_stats *stats;

   /* Timing results */
   long_long atime, ftime, ctime;
public:
   NodeWork(int index, struct cpu_node_data<T> *data,
         const struct cpu_factor_options *options, void *alloc,
         const T *aval, const T* scaling)
      : index(index), data(data), options(options), alloc(alloc), aval(aval),
        scaling(scaling)
      {}
   void check(int m, int n, int ndelay) {
      int nrow = data->nrow_expected + data->ndelay_in;
      int ncol = data->ncol_expected + data->ndelay_in;
      int del = data->ndelay_out;
      //if(nrow != m || ncol != n || del != ndelay)
      //if(nrow == m && ncol == n && del != ndelay)
      if(ndelay > 1.5*ndelay)
         printf("Node %d differs from ma97\n   SSIDS %dx%d delay %d\n   ma97 %dx%d delay %d\n", index, nrow, ncol, del, m, n, ndelay);
   }
   void setResource(
         int *const map,
         StackAllocator<PAGE_SIZE> *stalloc_odd,
         StackAllocator<PAGE_SIZE> *stalloc_even,
         Workspace<T> *work,
         struct cpu_factor_stats *stats
         ) {
      this->map = map;
      this->stalloc_odd = stalloc_odd;
      this->stalloc_even = stalloc_even;
      this->work = work;
      this->stats = stats;
   }
   long estTime(void) { return -1; }
   void split(void) { /* Noop */ }
   void exec(void) { 
      // Assembly
      if(timing) atime = PAPI_get_real_usec();
      assemble_node
         <T, PAGE_SIZE>
         (posdef, index, data, alloc, stalloc_odd, stalloc_even, map, aval,
          scaling);
      if(timing) atime = PAPI_get_real_usec() - atime;
      // Update stats
      int nrow = data->nrow_expected + data->ndelay_in;
      if(nrow > stats->maxfront) stats->maxfront = nrow;
      // Factorization
      if(timing) ftime = PAPI_get_real_usec();
      factor_node
         <posdef, T, BLOCK_SIZE>
         (index, data, options, stats);
      if(timing) ftime = PAPI_get_real_usec() - ftime;
      // Form update
      if(timing) ctime = PAPI_get_real_usec();
      calculate_update
         <posdef, T, PAGE_SIZE>
         (data, stalloc_odd, stalloc_even, work);
      if(timing) ctime = PAPI_get_real_usec() - ctime;
   }
};

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

   std::ifstream nodefile("ma97_nodes.dat");

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

   if(timing) {
      if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
         exit(1);
   }

   /* Set up work units */
   typedef NodeWork<posdef, T, BLOCK_SIZE, PAGE_SIZE, timing> NodeWorkSpec;
   std::vector<WorkUnit*> wUnits;
   wUnits.reserve(nnodes);
   for(int ni=0; ni<nnodes; ni++) {
      wUnits.push_back(
            new NodeWorkSpec(ni, &nodes[ni], options, alloc, aval, scaling)
            );
   }

   /* Main loop: Iterate over nodes in order */
   for(auto itr=wUnits.begin(); itr!=wUnits.end(); itr++) {
      NodeWorkSpec *nwu = dynamic_cast<NodeWorkSpec*> (*itr);
      nwu->setResource(map, &stalloc_odd, &stalloc_even, &work, stats);
      int idx, m, n, ndelay;
      nodefile >> idx >> m >> n >> ndelay;
      nwu->exec();
      nwu->check(m, n, ndelay);
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


   /* Extract and display timing information */
   if(timing) {
      const int TRES = 250;
      float *atime, *ftime, *ctime;
      long ldt = (stats->maxfront-1)/TRES + 1;
      atime = new float[ldt*ldt];
      for(long i=0; i<ldt*ldt; i++) atime[i] = 0.0;
      ftime = new float[ldt*ldt];
      for(long i=0; i<ldt*ldt; i++) ftime[i] = 0.0;
      ctime = new float[ldt*ldt];
      for(long i=0; i<ldt*ldt; i++) ctime[i] = 0.0;

      for(auto itr=wUnits.begin(); itr!=wUnits.end(); itr++) {
         NodeWorkSpec *nwu = dynamic_cast<NodeWorkSpec*> (*itr);
         const struct cpu_node_data<T> *data = nwu->data;
         int nrow = data->nrow_expected + data->ndelay_in;
         int ncol = data->ncol_expected + data->ndelay_in;
         atime[ncol/TRES*ldt + nrow/TRES] += nwu->atime*1e-6;
         ftime[ncol/TRES*ldt + nrow/TRES] += nwu->ftime*1e-6;
         ctime[ncol/TRES*ldt + nrow/TRES] += nwu->ctime*1e-6;
         if(nrow>3000)
            printf("%dx%d fact took %e\n", nrow, ncol, nwu->ftime*1e-6);
      }
      // Find maximum nonzero entries
      int maxr=0, maxc=0;
      float atotal = 0.0;
      for(int j=0; j<ldt; j++)
      for(int i=0; i<ldt; i++) {
         if(atime[j*ldt+i] > 0) {
            if(i>=maxr) maxr=i+1;
            if(j>=maxc) maxc=j+1;
            atotal += atime[j*ldt+i];
         }
      }
      
      // Print times
      printf("atime: %e\n", atotal);
      for(int i=0; i<maxr; i++) {
         printf("%6d:", i*TRES);
         atotal = 0.0;
         for(int j=0; j<maxc; j++) {
            if(atime[j*ldt+i] > 0) printf(" %8.2e", atime[j*ldt+i]);
            else                   printf(" %8s", "");
            atotal += atime[j*ldt+i];
         }
         printf(" =%8.2e\n", atotal);
      }
      float ftotal = 0.0;
      for(int j=0; j<maxc; j++)
      for(int i=0; i<maxr; i++)
         ftotal += ftime[j*ldt+i];
      printf("ftime: %e\n", ftotal);
      for(int i=0; i<maxr; i++) {
         printf("%6d:", i*TRES);
         ftotal = 0.0;
         for(int j=0; j<maxc; j++) {
            if(ftime[j*ldt+i] > 0) printf(" %8.2e", ftime[j*ldt+i]);
            else                   printf(" %8s", "");
            ftotal += ftime[j*ldt+i];
         }
         printf(" =%8.2e\n", ftotal);
      }
      float ctotal = 0.0;
      for(int j=0; j<maxc; j++)
      for(int i=0; i<maxr; i++)
         ctotal += ctime[j*ldt+i];
      printf("ctime: %e\n", ctotal);
      for(int i=0; i<maxr; i++) {
         printf("%6d:", i*TRES);
         ctotal = 0.0;
         for(int j=0; j<maxc; j++) {
            if(ctime[j*ldt+i] > 0) printf(" %8.2e", ctime[j*ldt+i]);
            else                   printf(" %8s", "");
            ctotal += ctime[j*ldt+i];
         }
         printf(" =%8.2e\n", ctotal);
      }

      delete[] atime;
      delete[] ftime;
      delete[] ctime;
   }
}

/** Represents an actual execution unit (i.e. context for a thread) */
class ExecUnit {
private:
   /** hwloc info for associated PU */
   hwloc_topology_t topology;
   hwloc_obj_t pu;
   /** pthread object associated with this ExecUnit */
   pthread_t thread;
   /** flag when to stop */
   enum stop_flags {
      KEEP_GOING,
      ABORT_IMMEDIATLY,
      STOP_WHEN_DONE
   };
   std::atomic<enum stop_flags> stopflag;
   /** Pipeline of work */
   std::queue<WorkUnit *> work;
   /** static function we can target with pthread_create to call exec()
    *  NB: This is where binding is actually performed. */
   static
   void *exec_this(void *eu_ptr) {
      // FIXME: Error checking
      ExecUnit *eu = (ExecUnit *) eu_ptr;
      hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
      hwloc_bitmap_set(cpuset, eu->pu->os_index);
      hwloc_set_cpubind(eu->topology, cpuset, HWLOC_CPUBIND_THREAD);
      hwloc_bitmap_free(cpuset);
      eu->exec();
      return NULL;
   }
public:
   ExecUnit (
         hwloc_topology_t topology,
         hwloc_obj_t pu
         )
      : topology(topology), pu(pu), stopflag(KEEP_GOING)
      {}
   /** Function to add work to be performed */
   void enqueue(WorkUnit *item) {
      // FIXME: thread protect
      work.push(item);
   }
   /** Function to signal an abort to execution immediately */
   void abort(void) {
      stopflag = ABORT_IMMEDIATLY;
   }
   /** Function to signal that exec() should return when it runs out of work */
   void end(void) {
      stopflag = STOP_WHEN_DONE;
   }
   /** Set the current thread into a spin loop running this exec unit's work */
   void exec(void) {
      printf("Hello from (hopefully) %d\n", pu->logical_index);
      //while(stopflag==KEEP_GOING) {}
   }
   /** Launch a thread, bind its affinity and have it call exec() */
   void launch_thread(void) {
      // NB binding actually performed by exec_this() function
      // FIXME: Error check
      pthread_create(&thread, NULL, exec_this, this);
   }
};

class ExecContext {
private:
   std::vector<ExecContext *> eContexts;
   std::vector<ExecUnit *> eUnits;
public:
   /** Construct ourselves from a given topology */
   ExecContext(
         hwloc_topology_t topology,
         hwloc_obj_t root) {
      // Strip extraneous layers
      while(root->arity == 1 && root->type != HWLOC_OBJ_PU)
         root = root->children[0];
      // Extract children
      for (unsigned int i = 0; i < root->arity; i++) {
         hwloc_obj_t obj = root->children[i];
         while(obj->arity == 1 && obj->type != HWLOC_OBJ_PU)
            obj = obj->children[0];
         if(obj->type == HWLOC_OBJ_PU)
            eUnits.push_back(new ExecUnit(topology, obj));
         else
            eContexts.push_back(new ExecContext(topology, obj));
      }
      
   }
   void abort(void) {
      for(auto itr=eContexts.begin(); itr!=eContexts.end(); itr++)
         (*itr)->abort();
      for(auto itr=eUnits.begin(); itr!=eUnits.end(); itr++)
         (*itr)->abort();
   }
   void end(void) {
      for(auto itr=eContexts.begin(); itr!=eContexts.end(); itr++)
         (*itr)->abort();
      for(auto itr=eUnits.begin(); itr!=eUnits.end(); itr++)
         (*itr)->abort();
   }
   void launch_threads(void) {
      for(auto itr=eContexts.begin(); itr!=eContexts.end(); itr++)
         (*itr)->launch_threads();
      for(auto itr=eUnits.begin(); itr!=eUnits.end(); itr++)
         (*itr)->launch_thread();
   }
};

static void print_children(hwloc_topology_t topology, hwloc_obj_t obj, int depth) {
   char string[128];
   hwloc_obj_snprintf(string, sizeof(string), topology, obj, "#", 0);
   printf("%*s%s\n", 2*depth, "", string);
   for (unsigned int i = 0; i < obj->arity; i++) {
      print_children(topology, obj->children[i], depth + 1);
   }
}

}}} /* end of namespace spral::ssids::internal */
//////////////////////////////////////////////////////////////////////////

/* Double precision wrapper around templated routines */
extern "C"
void spral_ssids_factor_cpu_dbl(
      bool posdef,     // If true, performs A=LL^T, if false do pivoted A=LDL^T
      int n,            // Maximum row index (+1)
      int nnodes,       // Number of nodes in assembly tree
      struct spral::ssids::internal::cpu_node_data<double> *const nodes, // Data structure for node information
      const double *const aval, // Values of A
      const double *const scaling, // Scaling vector (NULL if none)
      void *const alloc,      // Pointer to Fortran allocator structure
      const struct spral::ssids::internal::cpu_factor_options *const options, // Options in
      struct spral::ssids::internal::cpu_factor_stats *const stats // Info out
      ) {

   // Initialize execution topology
   hwloc_topology_t topology;
   hwloc_topology_init(&topology);
   hwloc_topology_load(topology);
   //print_children(topology, hwloc_get_root_obj(topology), 0);
   spral::ssids::internal::ExecContext eContext(topology, hwloc_get_root_obj(topology));
   //eContext.launch_threads();

   // Initialize stats
   stats->flag = spral::ssids::internal::SSIDS_SUCCESS;

   // Call relevant routine
   if(posdef) {
      try {
         spral::ssids::internal::factor<true, double, 16, 16384, false>
            (n, nnodes, nodes, aval, scaling, alloc, options, stats);
      } catch(spral::ssids::internal::NotPosDefError npde) {
         stats->flag = spral::ssids::internal::SSIDS_ERROR_NOT_POS_DEF;
      }
   } else {
      spral::ssids::internal::factor<false, double, 16, 16384, false>
         (n, nnodes, nodes, aval, scaling, alloc, options, stats);
   }

   // FIXME: Remove when done with debug
   if(options->print_level > 9999) {
      printf("Final factors:\n");
      spral::ssids::internal::print_factors<double>(posdef, nnodes, nodes);
   }
}
