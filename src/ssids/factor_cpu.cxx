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
#include <cstring>
#include <cstdio>
#include <limits> // FIXME: remove when done with debug if unneeded
#include <sstream>
#include <stdexcept>
/* External library headers */
#include "bub/bub.hxx"
/* SPRAL headers */
#include "factor_cpu_iface.h"

/////////////////////////////////////////////////////////////////////////////
// start of namespace spral::ssids::internal
namespace spral { namespace ssids { namespace internal {

const int SSIDS_SUCCESS = 0;
const int SSIDS_ERROR_NOT_POS_DEF = -6;

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
};

struct cpu_factor_stats {
   int flag;
};

template <typename T>
void assemble_node(
      struct cpu_node_data<T> *const node,
      void *const alloc,
      int *const map,
      const T *const aval,
      const T *const scaling
      ) {
   /* Count incoming delays and determine size of node */
   node->ndelay_in = 0;
   for(struct cpu_node_data<T> *child=node->first_child; child!=NULL; child=child->next_child)
      node->ndelay_in += child->ndelay_out;
   int nrow = node->nrow_expected + node->ndelay_in;
   int ncol = node->ncol_expected + node->ndelay_in;

   /* Get space for node now we know it size using Fortran allocator */
   size_t len = ((size_t) nrow+2) * ncol; // L is nrow x ncol and D is 2 x ncol
   node->lcol = smalloc<T>(alloc, len);
   node->perm = smalloc<int>(alloc, ncol); // ncol fully summed variables

   /* Get space for contribution block */
   long contrib_dimn = node->nrow_expected - node->ncol_expected;
   node->contrib = (contrib_dimn > 0) ? new T[contrib_dimn*contrib_dimn] : NULL;

   /* Zero node L and D */
   memset(node->lcol, 0, len*sizeof(T));

   /* Set perm for expected eliminations at this node */
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
         long k = node->ndelay_in*(nrow+1) + c*nrow + r;
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
         long k = node->ndelay_in*(nrow+1) + c*nrow + r;
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
      /* Loop over children adding contrubutions */
      for(struct cpu_node_data<T> *child=node->first_child; child!=NULL; child=child->next_child) {
         /* Handle delays - go to back of node */
         for(int i=0; i<child->ndelay_out; i++) {
            // Add delayed rows (in delayed cols)
            int delay_col = node->ncol_expected + i;
            T *dest = &node->lcol[delay_col*(nrow+1)];
            int lds = child->nrow_expected + child->ndelay_in;
            T *src = &child->lcol[child->nelim*(lds+1)];
            node->perm[delay_col] = child->perm[child->nelim+i];
            for(int j=i; j<child->ndelay_out; j++)
               dest[j] = src[i];
            // Add child's non-fully summed rows (in delayed cols)
            dest = &node->lcol[delay_col*nrow];
            src = &child->lcol[child->nelim*lds + child->ndelay_in];
            for(int j=child->ncol_expected; j<child->nrow_expected; j++) {
               int r = map[ child->rlist[j] ];
               dest[r] = src[j];
            }
         }
         /* Handle expected contributions */
         for(int i=child->ncol_expected; i<child->nrow_expected; i++) {
            int c = map[ child->rlist[i] ];
            int lds = child->nrow_expected - child->ncol_expected;
            T *src = &child->contrib[(i-child->ncol_expected)*(lds+1)];
            T *dest; int ldd;
            if(c < node->ncol_expected) {
               // Contribution added to lcol
               ldd = nrow;
               dest = &node->lcol[c*ldd];
               for(int j=child->ncol_expected; j<child->nrow_expected; j++) {
                  int r = map[ child->rlist[j] ];
                  dest[r] += src[r];
               }
            } else {
               // Contribution added to contrib
               // FIXME: Add after contribution block established?
               ldd = node->nrow_expected - node->ncol_expected;
               dest = &node->contrib[(c-node->ncol_expected)*lds];
               for(int j=child->ncol_expected; j<child->nrow_expected; j++) {
                  int r = map[ child->rlist[j] ] - child->ncol_expected;
                  dest[r] = src[r];
               }
            }
         }
         
         /* Free memory from child contribution block */
         delete[] child->contrib;
      }
   }
}
/* Factorize a node (indef) */
template <typename T, int BLOCK_SIZE>
void factor_node_indef(
      struct cpu_node_data<T> *const node,
      const struct cpu_factor_options *const options
      ) {
   /* Extract useful information about node */
   int m = node->nrow_expected + node->ndelay_in;
   int n = node->ncol_expected + node->ndelay_in;
   T *lcol = node->lcol;
   T *d = &node->lcol[ ((long) m)*n ];
   int *perm = node->perm;

   /* Perform factorization */
   typedef bub::CpuLDLT<T, BLOCK_SIZE> CpuLDLTSpec;
   //typedef bub::CpuLDLT<T, 5, true> CpuLDLTSpec; // FIXME: remove
   node->nelim = CpuLDLTSpec(options->u, options->small).factor(m, n, perm, lcol, m, d);

   /* Record information */
   node->ndelay_out = n - node->nelim;
}
/* Factorize a node (psdef) */
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
   //typedef bub::CpuLLT<T, 4, true> CpuLLTSpec; //FIXME: remove
   int flag = CpuLLTSpec().factor(m, n, lcol, m);
   node->nelim = (flag) ? flag : n;
   if(flag) throw NotPosDefError(flag);
}
/* Factorize a node (wrapper) */
template <bool posdef, typename T, int BLOCK_SIZE>
void factor_node(
      struct cpu_node_data<T> *const node,
      const struct cpu_factor_options *const options
      ) {
   if(posdef) factor_node_posdef<T, BLOCK_SIZE>(node, options);
   else       factor_node_indef <T, BLOCK_SIZE>(node, options);
}

/* FIXME: remove post debug */
template<typename T>
void print_node(int m, int n, const int *perm, const T *lcol, const T*d) {
   for(int i=0; i<m; i++) {
      printf("%d:", perm[i]);
      for(int j=0; j<n; j++) printf(" %10.2e", lcol[j*m+i]);
      printf("  d: %10.2e %10.2e\n", d[2*i+0], d[2*i+1]);
   }
}
/* FIXME: remove post debug */
template<typename T>
void print_factors(
      int nnodes,
      struct cpu_node_data<T> *const nodes
      ) {
   for(int node=0; node<nnodes; node++) {
      printf("== Node %d ==\n", node);
      int m = nodes[node].nrow_expected + nodes[node].ndelay_in;
      int n = nodes[node].ncol_expected + nodes[node].ndelay_in;
      print_node(m, n, nodes[node].perm, nodes[node].lcol, &nodes[node].lcol[m*n]);
   }
}

/* Calculate update */
template <bool posdef, typename T>
void calculate_update(
      struct cpu_node_data<T> *node
      ) {
   if(posdef) {
      int m = node->nrow_expected - node->ncol_expected;
      if(m==0) return; // no-op
      int k = node->nelim;
      int ldl = node->nrow_expected;
      host_syrk<T>(bub::FILL_MODE_LWR, bub::OP_N, m, k,
            -1.0, &node->lcol[node->ncol_expected], ldl,
            1.0, node->contrib, m);
   } else {
   }
}

/* Simplistic multifrontal factorization */
template <bool posdef, typename T>
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

   int *map = new int[n+1]; // +1 to allow for indexing with 1-indexed array

   /* Main loop: Iterate over nodes in order */
   for(int ni=0; ni<nnodes; ni++) {
      // Assembly
      assemble_node<T>(&nodes[ni], alloc, map, aval, scaling);
      // Factorization
      factor_node<posdef, T, 16>(&nodes[ni], options);
      // Form update
      calculate_update<posdef>(&nodes[ni]);
   }

   // FIXME: gross hack for compat with bub (which needs to differentiate
   // between a natural zero and a 2x2 factor's second entry without counting)
   // SSIDS original data format [a11 a21 a22 xxx] seems more bizzare than
   // bub one [a11 a21 inf a22]
   for(int ni=0; ni<nnodes; ni++) {
      int m = nodes[ni].nrow_expected + nodes[ni].ndelay_in;
      int n = nodes[ni].ncol_expected + nodes[ni].ndelay_in;
      T *d = nodes[ni].lcol + m*n;
      for(int i=0; i<2*nodes[ni].nelim; i++)
         if(d[i] == std::numeric_limits<T>::infinity())
            d[i] = d[i+1];
   }

   /* Free memory */
   delete[] map;
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

   const bool debug = false;

   // Initialize stats
   stats->flag = spral::ssids::internal::SSIDS_SUCCESS;

   // Call relevant routine
   if(posdef) {
      try {
         spral::ssids::internal::factor<true, double>
            (n, nnodes, nodes, aval, scaling, alloc, options, stats);
      } catch(spral::ssids::internal::NotPosDefError npde) {
         stats->flag = spral::ssids::internal::SSIDS_ERROR_NOT_POS_DEF;
      }
   } else {
      spral::ssids::internal::factor<false, double>
         (n, nnodes, nodes, aval, scaling, alloc, options, stats);
   }

   // FIXME: Remove when done with debug
   if(debug) {
      printf("Final factors:\n");
      spral::ssids::internal::print_factors<double>(nnodes, nodes);
   }
}
