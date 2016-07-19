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

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <utility>
#include <vector>

#include <omp.h>

#include "../AlignedAllocator.hxx"
#include "block_ldlt.hxx"
#include "common.hxx"
#include "wrappers.hxx"

namespace spral { namespace ssids { namespace cpu {

template<typename T,
         int BLOCK_SIZE,
         bool debug=false,
         typename Alloc=AlignedAllocator<T>
         >
class CpuLDLT {
private:
   /// u specifies the pivot tolerance
   const T u;
   /// small specifies the threshold for entries to be considered zero
   const T small;
   /// Allocator
   mutable Alloc alloc;

   /** Workspace allocated on a per-thread basis */
   struct ThreadWork {
      T *ld;

      ThreadWork() {
         ld = new T[BLOCK_SIZE*BLOCK_SIZE];
      }
      ~ThreadWork() {
         delete[] ld;
      }
   };
   /** Workspace allocated on a global basis, along with functions to aquire
    *  memory from it.
    *  FIXME: This can almost certainly be implemented better. */
   class GlobalWork {
   public:
      GlobalWork(int num_blocks)
      : num_blocks(num_blocks), block_pool(nullptr), block_free(nullptr),
        nextBlock(0)
      {
         block_pool = new T[num_blocks*BLOCK_SIZE*BLOCK_SIZE];
         block_free = new bool[num_blocks];
         for(int i=0; i<num_blocks; i++) block_free[i] = true;
         omp_init_lock(&lock);
      }
      ~GlobalWork() {
         delete[] block_pool;
         delete[] block_free;
         omp_destroy_lock(&lock);
      }

      /** Get next free block in a thread-safe fashion.
       *  Return nullptr if it can't find a free block.
       */
      T *acquire_block() {
         T *ptr = nullptr; // Return value
         omp_set_lock(&lock);
         
         /* We literally just loop over the blocks using nextBlock until
          * we find a free one or conclude it doesn't exist */
         for(int end_idx = nextBlock + num_blocks; nextBlock < end_idx; ) {
            int idx = (nextBlock++) % num_blocks;
            if(block_free[idx]) {
               block_free[idx] = false;
               ptr = &block_pool[idx*BLOCK_SIZE*BLOCK_SIZE];
               break;
            }
         }
         omp_unset_lock(&lock);
         return ptr;
      }
      /** Get next free block in a thread-safe fashion.
       *  Keep trying until it suceeds, use taskyield after each try.
       *  NB: This may deadlock if there are no other tasks releasing blocks.
       */
      T *acquire_block_wait() {
         while(true) {
            T *ptr = acquire_block();
            if(ptr) {
               return ptr;
            }
            #pragma omp taskyield
         }
      }
      /** Marks a block obtained by acquire_block() as availabel for reuse */
      void release_block(const T *const ptr) {
         int idx = (ptr-block_pool) / (BLOCK_SIZE*BLOCK_SIZE);
         // NB: I'm not sure we need to acquire a lock here, but best be safe
         omp_set_lock(&lock);
         block_free[idx] = true;
         omp_unset_lock(&lock);
      }
   private:
      int num_blocks; //< Number of blocks available in block_pool
      T *block_pool; //< Pointer to block pool
      bool *block_free; //< block_free[i] is true if block i is free
      int nextBlock; //< next block to check (use with modulus to get index)
      omp_lock_t lock; //< lock used for safe access to nextBlock, block_free
   };


   // FIXME: Horrendously inefficient, rework code to cope without
   void permute_a_d_to_elim_order(int n, const int *order, T *a, int lda, T *d, int *perm) {
      // Allocate memory
      T *acopy = new T[n*n];
      int *permcopy = new int[n];

      // Copy a (lwr -> lwr+upr)
      for(int j=0; j<n; j++) {
         for(int i=0; i<j; i++)
            acopy[j*n+i] = a[i*lda+j];
         for(int i=j; i<n; i++)
            acopy[j*n+i] = a[j*lda+i];
      }

      // And put a back in permuted form
      for(int jj=0; jj<n; jj++) {
         int j = order[jj];
         for(int ii=jj; ii<n; ii++) {
            int i = order[ii];
            a[jj*lda+ii] = acopy[j*n+i];
         }
      }

      // Copy d and perm
      for(int j=0; j<2*n; j++)
         acopy[j] = d[j];
      for(int i=0; i<n; i++)
         permcopy[i] = perm[i];

      // And put d and perm back in permuted form
      for(int jj=0; jj<n; jj++) {
         int j = order[jj];
         d[jj*2+0] = acopy[2*j+0];
         d[jj*2+1] = acopy[2*j+1];
         perm[jj] = permcopy[j];
      }

      // Free memory
      delete[] acopy;
      delete[] permcopy;
   }

   struct col_data {
      int oldelim;
      int nelim;
      int npass;
      omp_lock_t lock;
      int perm[BLOCK_SIZE]; // FIXME: Probably not needed, but really good for debug
      int elim_order[BLOCK_SIZE];
      T d[2*BLOCK_SIZE];

      col_data() : 
         nelim(0)
      {
         omp_init_lock(&lock);
         // FIXME: remove
         for(int i=0; i<BLOCK_SIZE; i++) elim_order[i] = -1;
      }
      ~col_data() {
         omp_destroy_lock(&lock);
      }

      void permuted_store(T *user_d, int *user_perm) {
         for(int j=0; j<BLOCK_SIZE; j++) {
            int col = elim_order[j];
            if(col < 0 ) continue; // Padding column
            user_d[2*col+0] = d[2*j+0];
            user_d[2*col+1] = d[2*j+1];
            user_perm[col] = perm[j];
         }
      }

   };

   class BlockData {
      public:
      T *aglobal;
      /// Latest accepted value of A or L
      T __attribute__((aligned(32))) aval[BLOCK_SIZE*BLOCK_SIZE];
      /// Trial value of L
      T *lwork;

      BlockData()
      : aglobal(nullptr), lwork(nullptr)
      {}

      void init(int nrow, int ncol, T *a, int lda) {
         aglobal = a;
         load(nrow, ncol, aglobal, lda);
      }

      /// Loads block of a into this->aval.
      //  If too few rows/cols, pad with NaNs.
      //  Treat too few cols specially: it will be marked as already eliminated
      //  and pushed to the left to avoid trying to eliminate it.
      // FIXME: Is NaN init actually required if count as already elim?
      void load(int nrow, int ncol, T *a, int lda) {
         int roffset = BLOCK_SIZE-nrow;
         int coffset = BLOCK_SIZE-ncol;
         // Pad missing columns with NaNs
         T *av = aval;
         for(int j=0; j<coffset; j++)
            for(int i=0; i<BLOCK_SIZE; i++)
               aval[j*BLOCK_SIZE+i] = std::numeric_limits<T>::quiet_NaN();
         av += coffset*BLOCK_SIZE;
         // Pad missing rows with NaNs
         for(int j=0; j<ncol; j++)
            for(int i=0; i<roffset; i++)
               av[j*BLOCK_SIZE+i] = std::numeric_limits<T>::quiet_NaN();
         av += roffset;
         // Load actual data
         for(int j=0; j<ncol; j++)
            for(int i=0; i<nrow; i++)
               av[j*BLOCK_SIZE+i] = a[j*lda+i];
      }

      void create_restore_point() {
         for(int j=0; j<BLOCK_SIZE; j++)
         for(int i=0; i<BLOCK_SIZE; i++)
            lwork[j*BLOCK_SIZE+i] = aval[j*BLOCK_SIZE+i];
      }

      void create_restore_point_with_row_perm(const int *lperm) {
         for(int j=0; j<BLOCK_SIZE; j++)
         for(int i=0; i<BLOCK_SIZE; i++) {
            int r = lperm[i];
            lwork[j*BLOCK_SIZE+i] = aval[j*BLOCK_SIZE+r];
         }
         for(int j=0; j<BLOCK_SIZE; j++)
         for(int i=0; i<BLOCK_SIZE; i++)
            aval[j*BLOCK_SIZE+i] = lwork[j*BLOCK_SIZE+i];
      }

      void create_restore_point_with_col_perm(const int *lperm) {
         for(int j=0; j<BLOCK_SIZE; j++) {
            int c = lperm[j];
            for(int i=0; i<BLOCK_SIZE; i++)
               lwork[j*BLOCK_SIZE+i] = aval[c*BLOCK_SIZE+i];
         }
         for(int j=0; j<BLOCK_SIZE; j++)
         for(int i=0; i<BLOCK_SIZE; i++)
            aval[j*BLOCK_SIZE+i] = lwork[j*BLOCK_SIZE+i];
      }

      /** Restores any columns that have failed back to their previous
       *  values stored in lwork[] */
      void restore_part(int rfrom, int cfrom) {
         for(int j=cfrom; j<BLOCK_SIZE; j++)
         for(int i=rfrom; i<BLOCK_SIZE; i++)
            aval[j*BLOCK_SIZE+i] = lwork[j*BLOCK_SIZE+i];
      }

      /** Restores any columns that have failed back to their previous
       *  values stored in lwork[]. Applies a symmetric permutation while
       *  doing so. */
      void restore_part_with_sym_perm(int from, const int *lperm) {
         for(int j=from; j<BLOCK_SIZE; j++) {
            int c = lperm[j];
            for(int i=from; i<BLOCK_SIZE; i++) {
               int r = lperm[i];
               aval[j*BLOCK_SIZE+i] = (r>c) ? lwork[c*BLOCK_SIZE+r]
                                            : lwork[r*BLOCK_SIZE+c];
            }
         }
      }

      void permuted_store(T *user_a, int lda, const struct col_data *idata, const struct col_data *jdata) {
         for(int j=0; j<BLOCK_SIZE; j++) {
            int col = jdata->elim_order[j];
            if(col < 0 ) continue; // Padding column
            for(int i=0; i<BLOCK_SIZE; i++) {
               int row = idata->elim_order[i];
               if(row < 0) continue; // Padding row
               if(row > col)
                  user_a[col*lda+row] = aval[j*BLOCK_SIZE+i];
               else
                  user_a[row*lda+col] = aval[j*BLOCK_SIZE+i];
            }
         }
      }

      void permuted_store_diag(T *user_a, int lda, const struct col_data *cdata) {
         for(int j=0; j<BLOCK_SIZE; j++) {
            int col = cdata->elim_order[j];
            if(col < 0 ) continue; // Padding column
            for(int i=j; i<BLOCK_SIZE; i++) {
               int row = cdata->elim_order[i];
               if(row < 0) continue; // Padding row
               if(row > col)
                  user_a[col*lda+row] = aval[j*BLOCK_SIZE+i];
               else
                  user_a[row*lda+col] = aval[j*BLOCK_SIZE+i];
            }
         }
      }

      void col_permuted_store(int nrow, T *user_a, int lda, const struct col_data *jdata) {
         T *av = aval + (BLOCK_SIZE-nrow);
         for(int j=0; j<BLOCK_SIZE; j++) {
            int col = jdata->elim_order[j];
            if(col < 0 ) continue; // Padding column
            for(int i=0; i<nrow; i++)
               user_a[col*lda+i] = av[j*BLOCK_SIZE+i];
         }
      }

      /** Check if a block satisifies pivot threshold (colwise version) */
      int check_threshold(int rfrom, int cfrom, T u) {
         // Perform thrshold test for each uneliminated row/column
         for(int j=cfrom; j<BLOCK_SIZE; j++)
         for(int i=rfrom; i<BLOCK_SIZE; i++)
            if(fabs(aval[j*BLOCK_SIZE+i]) > 1.0/u) {
               if(debug) printf("Failed %d,%d:%e\n", i, j, fabs(aval[j*BLOCK_SIZE+i]));
               return j;
            }
         // If we get this far, everything is good
         return BLOCK_SIZE;
      }

      /** Check if a block satisifies pivot threshold (rowwise version) */
      int check_threshold_trans(int rfrom, int cfrom, T u) {
         // Perform thrshold test for each uneliminated row/column
         for(int j=rfrom; j<BLOCK_SIZE; j++)
         for(int i=cfrom; i<BLOCK_SIZE; i++)
            if(fabs(aval[i*BLOCK_SIZE+j]) > 1.0/u) {
               if(debug) printf("Failed %d,%d:%e\n", i, j, fabs(aval[i*BLOCK_SIZE+j]));
               return j;
            }
         // If we get this far, everything is good
         return BLOCK_SIZE;
      }

      /** Performs solve with diagonal block \f$L_{21} = A_{21} L_{11}^{-T} D_1^{-1}\f$. Designed for below diagonal. */
      void apply_pivot(int rfrom, int cfrom, const T *diag, int ldd, const T *d, const T small) {
         if(rfrom >= BLOCK_SIZE || cfrom >= BLOCK_SIZE) return; // no-op
         // Perform solve L_11^-T
         host_trsm<T>(SIDE_RIGHT, FILL_MODE_LWR, OP_T, DIAG_UNIT, BLOCK_SIZE-rfrom, BLOCK_SIZE-cfrom, 1.0, &diag[cfrom*BLOCK_SIZE+cfrom], ldd, &aval[cfrom*BLOCK_SIZE+rfrom], BLOCK_SIZE);
         // Perform solve L_21 D^-1
         for(int i=cfrom; i<BLOCK_SIZE; ) {
            if(d[2*i+1]==0.0) {
               // 1x1 pivot
               T d11 = d[2*i];
               if(d11 == 0.0) {
                  // Handle zero pivots carefully
                  for(int j=rfrom; j<BLOCK_SIZE; j++) {
                     T v = aval[i*BLOCK_SIZE+j];
                     aval[i*BLOCK_SIZE+j] = 
                        (fabs(v)<small) ? 0.0
                                        : std::numeric_limits<T>::infinity()*v;
                     // NB: *v above handles NaNs correctly
                  }
               } else {
                  // Non-zero pivot, apply in normal fashion
                  for(int j=rfrom; j<BLOCK_SIZE; j++)
                     aval[i*BLOCK_SIZE+j] *= d11;
               }
               i++;
            } else {
               // 2x2 pivot
               T d11 = d[2*i];
               T d21 = d[2*i+1];
               T d22 = d[2*i+3];
               for(int j=rfrom; j<BLOCK_SIZE; j++) {
                  T a1 = aval[i*BLOCK_SIZE+j];
                  T a2 = aval[(i+1)*BLOCK_SIZE+j];
                  aval[i*BLOCK_SIZE+j]     = d11*a1 + d21*a2;
                  aval[(i+1)*BLOCK_SIZE+j] = d21*a1 + d22*a2;
               }
               i += 2;
            }
         }
      }

      /** Performs solve with diagonal block \f$L_{21}^T = D_1^{-T} L_{11}^{-1} A_{21}^T \f$. Designed for left of diagonal.  */
      void apply_pivot_trans(int rfrom, int cfrom, const T *diag, int ldd, const T *d, const T small) {
         if(rfrom >= BLOCK_SIZE || cfrom >= BLOCK_SIZE) return; // no uneliminated columns
         // Perform solve L_11^-1
         host_trsm<T>(SIDE_LEFT, FILL_MODE_LWR, OP_N, DIAG_UNIT, BLOCK_SIZE-rfrom, BLOCK_SIZE-cfrom, 1.0, &diag[rfrom*ldd+rfrom], ldd, &aval[cfrom*BLOCK_SIZE+rfrom], BLOCK_SIZE);
         // Perform solve D^-T L_21^T
         for(int i=rfrom; i<BLOCK_SIZE; ) {
            if(d[2*i+1]==0.0) {
               // 1x1 pivot
               T d11 = d[2*i];
               if(d11 == 0.0) {
                  // Handle zero pivots carefully
                  for(int j=cfrom; j<BLOCK_SIZE; j++) {
                     T v = aval[j*BLOCK_SIZE+i];
                     aval[j*BLOCK_SIZE+i] = 
                        (fabs(v)<small) ? 0.0 // *v handles NaNs
                                        : std::numeric_limits<T>::infinity()*v;
                     // NB: *v above handles NaNs correctly
                  }
               } else {
                  // Non-zero pivot, apply in normal fashion
                  for(int j=cfrom; j<BLOCK_SIZE; j++) {
                     aval[j*BLOCK_SIZE+i] *= d11;
                  }
               }
               i++;
            } else {
               // 2x2 pivot
               T d11 = d[2*i];
               T d21 = d[2*i+1];
               T d22 = d[2*i+3];
               for(int j=cfrom; j<BLOCK_SIZE; j++) {
                  T a1 = aval[j*BLOCK_SIZE+i];
                  T a2 = aval[j*BLOCK_SIZE+(i+1)];
                  aval[j*BLOCK_SIZE+i]     = d11*a1 + d21*a2;
                  aval[j*BLOCK_SIZE+(i+1)] = d21*a1 + d22*a2;
               }
               i += 2;
            }
         }
      }

      /** Apply successful pivot update to all uneliminated columns 
       *  (this.aval in non-transpose) */
      void update_n(int oldelim, int nelim, const T *l, int ldl, const T *ld, int ldld, int rfrom, int cfrom) {
         host_gemm(OP_N, OP_T,
               BLOCK_SIZE-rfrom, BLOCK_SIZE-cfrom, nelim-oldelim,
               -1.0, &ld[oldelim*ldld+rfrom], ldld, &l[oldelim*ldl+cfrom], ldl,
               1.0, &aval[cfrom*BLOCK_SIZE+rfrom], BLOCK_SIZE);
         /*for(int j=cfrom; j<BLOCK_SIZE; j++)
            for(int i=rfrom; i<BLOCK_SIZE; i++)
               for(int k=oldelim; k<nelim; k++)
                  aval[j*BLOCK_SIZE+i] -= l[k*ldl+j] * ld[k*ldld+i];*/
      }

      /** Apply successful pivot update to all uneliminated columns
       * (this.aval in transpose) */
      void update_t(int oldelim, int nelim, const T *l, int ldl, const T *ld, int ldld, int rfrom, int cfrom) {
         host_gemm(OP_N, OP_N,
               BLOCK_SIZE-rfrom, BLOCK_SIZE-cfrom, nelim-oldelim,
               -1.0, &ld[oldelim*ldld+rfrom], ldld, &l[cfrom*ldl+oldelim], ldl,
               1.0, &aval[cfrom*BLOCK_SIZE+rfrom], BLOCK_SIZE);
         /*for(int j=cfrom; j<BLOCK_SIZE; j++)
            for(int i=rfrom; i<BLOCK_SIZE; i++)
               for(int k=oldelim; k<nelim; k++)
                  aval[j*BLOCK_SIZE+i] -= l[j*ldl+k] * ld[k*ldld+i];*/
      }

      void rcmax_unelim(int rfrom, int cfrom, int *ridx, T *rval, int *cidx, T *cval) {
         for(int j=cfrom; j<BLOCK_SIZE; j++) {
            for(int i=rfrom; i<BLOCK_SIZE; i++) {
               T v = fabs(aval[j*BLOCK_SIZE+i]);
               if(v > cval[j]) {
                  cidx[j] = i;
                  cval[j] = v;
               }
               if(v > rval[i]) {
                  ridx[i] = j;
                  rval[i] = v;
               }
            }
         }
      }

      void rcmax_unelim_diag(int from, int *idx, T *val) {
         for(int j=from; j<BLOCK_SIZE; j++) {
            for(int i=j; i<BLOCK_SIZE; i++) {
               T v = fabs(aval[j*BLOCK_SIZE+i]);
               if(v > val[j]) {
                  idx[j] = i;
                  val[j] = v;
               }
               if(v > val[i]) {
                  idx[i] = j;
                  val[i] = v;
               }
            }
         }
      }

      /** Returns column idx in vector xch */
      void get_col(int idx, T xch[BLOCK_SIZE]) const {
         for(int i=0; i<BLOCK_SIZE; i++)
            xch[i] = aval[idx*BLOCK_SIZE+i];
      }
      /** Sets column idx to vector xch */
      void set_col(int idx, const T xch[BLOCK_SIZE]) {
         for(int i=0; i<BLOCK_SIZE; i++)
            aval[idx*BLOCK_SIZE+i] = xch[i];
      }
      /** Returns row idx in vector xch */
      void get_row(int idx, T xch[BLOCK_SIZE]) const {
         for(int i=0; i<BLOCK_SIZE; i++)
            xch[i] = aval[i*BLOCK_SIZE+idx];
      }
      /** Sets row idx to vector xch */
      void set_row(int idx, const T xch[BLOCK_SIZE]) {
         for(int i=0; i<BLOCK_SIZE; i++)
            aval[i*BLOCK_SIZE+idx] = xch[i];
      }
      /** Returns row/column idx in vector xch (for diagonal blocks) */
      void get_rc_diag(int idx, T xch[BLOCK_SIZE]) const {
         for(int i=0; i<idx; i++)
            xch[i] = aval[i*BLOCK_SIZE+idx];
         for(int i=idx; i<BLOCK_SIZE; i++)
            xch[i] = aval[idx*BLOCK_SIZE+i];
      }
      /** Sets row/column idx to vector xch (for diagonal blocks) */
      void set_rc_diag(int idx, const T xch[BLOCK_SIZE]) {
         for(int i=0; i<idx; i++)
            aval[i*BLOCK_SIZE+idx] = xch[i];
         for(int i=idx; i<BLOCK_SIZE; i++)
            aval[idx*BLOCK_SIZE+i] = xch[i];
      }
   };

   /** Calculates LD from L and D */
   template <enum operation op>
   static
   void calcLD(int m, int n, const T *l, int ldl, const T *d, T *ld, int ldld) {
      for(int col=0; col<n; ) {
         if(d[2*col+1]==0.0) {
            // 1x1 pivot
            T d11 = d[2*col];
            if(d11 != 0.0) d11 = 1/d11; // Zero pivots just cause zeroes
            for(int row=0; row<m; row++)
               ld[col*ldld+row] = d11 * ((op==OP_N) ? l[col*ldl+row]
                                                    : l[row*ldl+col]);
            col++;
         } else {
            // 2x2 pivot
            T d11 = d[2*col];
            T d21 = d[2*col+1];
            T d22 = d[2*col+3];
            T det = d11*d22 - d21*d21;
            d11 = d11/det;
            d21 = d21/det;
            d22 = d22/det;
            for(int row=0; row<m; row++) {
               T a1 = (op==OP_N) ? l[col*ldl+row]     : l[row*ldl+col];
               T a2 = (op==OP_N) ? l[(col+1)*ldl+row] : l[row*ldl+(col+1)];
               ld[col*ldld+row]     =  d22*a1 - d21*a2;
               ld[(col+1)*ldld+row] = -d21*a1 + d11*a2;
            }
            col += 2;
         }
      }
   }


   bool run_elim(int &next_elim, const int mblk, const int nblk, struct col_data *cdata, BlockData *blkdata, GlobalWork &global_work, ThreadWork all_thread_work[]) {
      bool changed = false;

      // FIXME: is global_lperm really the best way?
      int *global_lperm = new int[nblk*BLOCK_SIZE];

      // Record how many were eliminated prior to this pass of blk col
      for(int blk=0; blk<nblk; blk++)
         cdata[blk].oldelim = cdata[blk].nelim;

      /* Inner loop - iterate over block columns */
      for(int blk=0; blk<nblk; blk++) {
         // Don't bother adding tasks if we eliminated everything already
         if(cdata[blk].oldelim>=BLOCK_SIZE) continue;

         if(debug) {
            printf("Bcol %d:\n", blk);
            print_mat(mblk, nblk, blkdata, cdata);
         }

         // Factor diagonal: depend on cdata[blk] as we do some init here
         #pragma omp task default(none) \
            firstprivate(blk) \
            shared(blkdata, cdata, global_lperm, global_work, all_thread_work)\
            depend(inout: blkdata[blk*mblk+blk:1]) \
            depend(inout: cdata[blk:1])
         {
            int thread_num = omp_get_thread_num();
            ThreadWork &thread_work = all_thread_work[thread_num];
            blkdata[blk*mblk+blk].lwork = global_work.acquire_block_wait();
            BlockData &dblk = blkdata[blk*mblk+blk];
            int *lperm = &global_lperm[blk*BLOCK_SIZE];
            for(int i=0; i<BLOCK_SIZE; i++)
               lperm[i] = i;
            dblk.create_restore_point();
            block_ldlt<T, BLOCK_SIZE>(cdata[blk].oldelim, cdata[blk].perm, dblk.aval, cdata[blk].d, thread_work.ld, u, small, lperm);
            // Initialize threshold check (no lock required becuase task depend)
            cdata[blk].npass = BLOCK_SIZE;
         }
         
         // Loop over off-diagonal blocks applying pivot
         for(int jblk=0; jblk<blk; jblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, jblk) \
               shared(blkdata, cdata, global_lperm, global_work) \
               depend(in: blkdata[blk*mblk+blk:1]) \
               depend(inout: blkdata[jblk*mblk+blk:1]) \
               depend(in: cdata[blk:1])
            {
               typedef void (BlockData::*BlockDataApplyTPtr)(int, int, const T*, int, const T*, const T);
               BlockDataApplyTPtr bdatp = &BlockData::apply_pivot_trans;
               BlockData &dblk = blkdata[blk*mblk+blk];
               BlockData &cblk = blkdata[jblk*mblk+blk];
               const int *lperm = &global_lperm[blk*BLOCK_SIZE];
               // Perform necessary operations
               cblk.lwork = global_work.acquire_block_wait();
               cblk.create_restore_point_with_row_perm(lperm);
               cblk.apply_pivot_trans(cdata[blk].oldelim, cdata[jblk].nelim, dblk.aval, BLOCK_SIZE, cdata[blk].d, small);
               // Update threshold check
               int blkpass = cblk.check_threshold_trans(cdata[blk].oldelim, cdata[jblk].nelim, u);
               omp_set_lock(&cdata[blk].lock);
               if(blkpass < cdata[blk].npass)
                  cdata[blk].npass = blkpass;
               omp_unset_lock(&cdata[blk].lock);
            }
         }
         for(int iblk=blk+1; iblk<mblk; iblk++) {
            #pragma omp task default(none) \
               firstprivate(blk, iblk) \
               shared(blkdata, cdata, global_lperm, global_work) \
               depend(in: blkdata[blk*mblk+blk:1]) \
               depend(inout: blkdata[blk*mblk+iblk:1]) \
               depend(in: cdata[blk:1])
            {
               BlockData &dblk = blkdata[blk*mblk+blk];
               BlockData &rblk = blkdata[blk*mblk+iblk];
               const int *lperm = &global_lperm[blk*BLOCK_SIZE];
               // Perform necessary operations
               rblk.lwork = global_work.acquire_block_wait();
               rblk.create_restore_point_with_col_perm(lperm);
               int rfrom = (iblk < nblk) ? cdata[iblk].nelim : 0;
               rblk.apply_pivot(rfrom, cdata[blk].oldelim, dblk.aval, BLOCK_SIZE, cdata[blk].d, small);
               // Update threshold check
               int blkpass = rblk.check_threshold(rfrom, cdata[blk].oldelim, u);
               omp_set_lock(&cdata[blk].lock);
               if(blkpass < cdata[blk].npass)
                  cdata[blk].npass = blkpass;
               omp_unset_lock(&cdata[blk].lock);
            }
         }

         // Adjust column once all applys have finished and we know final
         // number of passed columns.
         #pragma omp task default(none) \
            firstprivate(blk) \
            shared(blkdata, cdata, changed, next_elim, global_work) \
            depend(inout: cdata[blk:1])
         {
            // Adjust to avoid splitting 2x2 pivots
            if(cdata[blk].npass>cdata[blk].oldelim) {
               T d11 = cdata[blk].d[2*(cdata[blk].npass-1) + 0];
               T d21 = cdata[blk].d[2*(cdata[blk].npass-1) + 1];
               if(d21!=0.0 && d11!=std::numeric_limits<T>::infinity()) {
                  // last passed entry was first part of 2x2
                  cdata[blk].npass--; 
               }
            }
            if(debug) printf("Adjusted to %d\n", cdata[blk].npass);
            // Count threshold
            for(int i=cdata[blk].oldelim; i<cdata[blk].npass; i++) {
               cdata[blk].elim_order[i] = next_elim++;
               cdata[blk].nelim++;
            }

            // Record if we eliminated anything
            changed = changed || (cdata[blk].oldelim != cdata[blk].nelim);
         }

         // Update uneliminated columns
         for(int jblk=0; jblk<blk; jblk++) {
            for(int iblk=jblk; iblk<mblk; iblk++) {
               // Calculate block index we depend on for i
               // (we only work with lower half of matrix)
               int iblk_idx = (blk < iblk) ? blk*mblk+iblk
                                           : iblk*mblk+blk;
               #pragma omp task default(none) \
                  firstprivate(blk, iblk, jblk) \
                  shared(cdata, blkdata, all_thread_work, global_work) \
                  depend(inout: blkdata[jblk*mblk+iblk:1]) \
                  depend(in: cdata[blk:1]) \
                  depend(in: blkdata[jblk*mblk+blk:1]) \
                  depend(in: blkdata[iblk_idx:1])
               {
                  // If we're on the block row we've just eliminated, restore
                  // any failed rows and release resources storing checkpoint
                  if(iblk==blk) {
                     if(cdata[blk].nelim < BLOCK_SIZE)
                        blkdata[jblk*mblk+iblk]
                           .restore_part(cdata[blk].nelim, cdata[jblk].nelim);
                     global_work.release_block(blkdata[jblk*mblk+blk].lwork);
                  }
                  // Perform actual update (if required)
                  if(cdata[blk].oldelim != cdata[blk].nelim) {
                     int thread_num = omp_get_thread_num();
                     ThreadWork &thread_work = all_thread_work[thread_num];
                     int oldelim = cdata[blk].oldelim;
                     int nelim = cdata[blk].nelim;
                     int rfrom = (iblk < nblk) ? cdata[iblk].nelim : 0;
                     if(blk <= iblk) {
                        calcLD<OP_N>(BLOCK_SIZE-rfrom, nelim-oldelim,
                           &blkdata[blk*mblk+iblk].aval[oldelim*BLOCK_SIZE+rfrom], BLOCK_SIZE,
                           &cdata[blk].d[2*oldelim],
                           &thread_work.ld[oldelim*BLOCK_SIZE+rfrom], BLOCK_SIZE);
                     } else {
                        calcLD<OP_T>(BLOCK_SIZE-rfrom, nelim-oldelim,
                           &blkdata[iblk*mblk+blk].aval[rfrom*BLOCK_SIZE+oldelim], BLOCK_SIZE,
                           &cdata[blk].d[2*oldelim],
                           &thread_work.ld[oldelim*BLOCK_SIZE+rfrom], BLOCK_SIZE);
                     }
                     blkdata[jblk*mblk+iblk].update_t(oldelim, nelim,
                           blkdata[jblk*mblk+blk].aval, BLOCK_SIZE,
                           thread_work.ld, BLOCK_SIZE,
                           rfrom, cdata[jblk].nelim);
                  }
               }
            }
         }
         for(int jblk=blk; jblk<nblk; jblk++) {
            for(int iblk=jblk; iblk<mblk; iblk++) {
               #pragma omp task default(none) \
                  firstprivate(blk, iblk, jblk) \
                  shared(cdata, blkdata, all_thread_work, global_lperm, \
                         global_work) \
                  depend(inout: blkdata[jblk*mblk+iblk:1]) \
                  depend(in: cdata[blk:1]) \
                  depend(in: blkdata[blk*mblk+iblk:1]) \
                  depend(in: blkdata[blk*mblk+jblk:1])
               {
                  // If we're on the block col we've just eliminated, restore
                  // any failed cols and release checkpoint resources
                  if(jblk==blk) {
                     if(cdata[blk].nelim < BLOCK_SIZE) {
                        if(iblk==blk) {
                           // Diagonal block needs to apply a permutation
                           const int *lperm = &global_lperm[blk*BLOCK_SIZE];
                           blkdata[jblk*mblk+iblk].restore_part_with_sym_perm(
                                 cdata[blk].nelim, lperm
                                 );
                        } else {
                           int rfrom = (iblk < nblk) ? cdata[iblk].nelim : 0;
                           blkdata[jblk*mblk+iblk]
                              .restore_part(rfrom, cdata[blk].nelim);
                        }
                     }
                     global_work.release_block(blkdata[jblk*mblk+iblk].lwork);
                  }
                  // Perform actual update (if required)
                  if(cdata[blk].oldelim != cdata[blk].nelim) {
                     int thread_num = omp_get_thread_num();
                     ThreadWork &thread_work = all_thread_work[thread_num];
                     int oldelim = cdata[blk].oldelim;
                     int nelim = cdata[blk].nelim;
                     int rfrom = (iblk < nblk) ? cdata[iblk].nelim : 0;
                     calcLD<OP_N>(BLOCK_SIZE-rfrom, nelim-oldelim,
                           &blkdata[blk*mblk+iblk].aval[oldelim*BLOCK_SIZE+rfrom],
                           BLOCK_SIZE, &cdata[blk].d[2*oldelim],
                           &thread_work.ld[oldelim*BLOCK_SIZE+rfrom], BLOCK_SIZE);
                     blkdata[jblk*mblk+iblk].update_n(oldelim, nelim,
                           blkdata[blk*mblk+jblk].aval, BLOCK_SIZE,
                           thread_work.ld, BLOCK_SIZE,
                           rfrom, cdata[jblk].nelim);
                  }
               }
            }
         }
      }
      #pragma omp taskwait

      delete[] global_lperm;

      return changed;
   }

   static
   void print_mat(int m, int n, const int *perm, const bool *eliminated, const T *a, int lda) {
      for(int row=0; row<m; row++) {
         if(row < n)
            printf("%d%s:", perm[row], eliminated[row]?"X":" ");
         else
            printf("%d%s:", row, "U");
         for(int col=0; col<n; col++)
            printf(" %10.4f", a[col*lda+row]);
         printf("\n");
      }
   }

   static
   void print_mat(int mblk, int nblk, const BlockData *blkdata, const struct col_data *cdata) {
      for(int rblk=0; rblk<mblk; rblk++) {
         for(int row=0; row<BLOCK_SIZE; row++) {
            int r = rblk*BLOCK_SIZE+row;
            if(r < nblk*BLOCK_SIZE)
               printf("%d%s:", cdata[rblk].perm[row], (row<cdata[rblk].nelim)?"X":" ");
            else
               printf("%d%s:", r, "U");
            for(int cblk=0; cblk<nblk; cblk++) {
               const BlockData &blk = blkdata[cblk*mblk+rblk];
               for(int col=0; col<BLOCK_SIZE; col++)
                  printf(" %10.4f", blk.aval[col*BLOCK_SIZE+row]);
            }
            printf("\n");
         }
      }
   }

public:
   CpuLDLT(T u, T small)
   : u(u), small(small)
   {}

   /** Factorize an entire matrix */
   int factor(int m, int n, int *perm, T *a, int lda, T *d) {
      /* Sanity check arguments */
      if(m < n) return -1;
      if(lda < n) return -4;

      /* Initialize useful quantities:
       * If we have m > n, then need to separate diag block and rect part to
       * make handling easier - hence the funny calculation for mblk. */
      int nblk = (n-1) / BLOCK_SIZE + 1;
      int mblk = (m>n) ? nblk + (m-n-1) / BLOCK_SIZE + 1 : nblk;
      int next_elim = 0;

      /* Load data block-wise */
      typedef typename std::allocator_traits<Alloc>::template rebind_alloc<BlockData> BlockDataAlloc;
      BlockDataAlloc bdalloc(alloc);
      BlockData *blkdata = std::allocator_traits<BlockDataAlloc>::allocate(
            bdalloc, mblk*nblk
            );
      for(int i=0; i<mblk*nblk; i++)
         std::allocator_traits<BlockDataAlloc>::construct(
               bdalloc, &blkdata[i]
               );
      for(int jblk=0; jblk<nblk; jblk++) {
         // Diagonal block part
         for(int iblk=0; iblk<nblk; iblk++)
            blkdata[jblk*mblk+iblk].init(
                  ((iblk+1)*BLOCK_SIZE<n) ? BLOCK_SIZE : (n - iblk*BLOCK_SIZE),
                  ((jblk+1)*BLOCK_SIZE<n) ? BLOCK_SIZE : (n - jblk*BLOCK_SIZE),
                  &a[jblk*BLOCK_SIZE*lda + iblk*BLOCK_SIZE], lda
                  );
         // Rectangular block below it
         int m2 = m-n;
         T *arect = &a[n];
         for(int iblk=0; iblk<mblk-nblk; iblk++)
            blkdata[jblk*mblk+nblk+iblk].init(
                  ((iblk+1)*BLOCK_SIZE<m2) ? BLOCK_SIZE : (m2 - iblk*BLOCK_SIZE),
                  ((jblk+1)*BLOCK_SIZE<n) ? BLOCK_SIZE : (n - jblk*BLOCK_SIZE),
                  &arect[jblk*BLOCK_SIZE*lda + iblk*BLOCK_SIZE], lda
                  );
      }

      /* Temporary workspaces */
      struct col_data *cdata = new struct col_data[nblk];

      /* Load column data */
      for(int blk=0; blk<nblk-1; blk++)
         for(int i=0; i<BLOCK_SIZE; i++)
            cdata[blk].perm[i] = perm[blk*BLOCK_SIZE+i];
      {
         // Handle last block specially to allow for undersize
         int coffset = nblk*BLOCK_SIZE - n;
         for(int i=0; i<BLOCK_SIZE-coffset; i++)
            cdata[nblk-1].perm[coffset+i] = perm[(nblk-1)*BLOCK_SIZE+i];
      }
      if(n < nblk*BLOCK_SIZE) {
         // Account for extra cols as "already eliminated"
         cdata[nblk-1].nelim = nblk*BLOCK_SIZE - n;
         // perm is too small for extra "NaN" rows, adjust to acommodate
         // Fill out "extra" entries with negative numbers for debugging sanity
         int coffset = nblk*BLOCK_SIZE - n;
         for(int i=0,j=-1; i<coffset; i++, j--) {
            cdata[nblk-1].perm[i] = j;
            cdata[nblk-1].elim_order[i] = -1;
         }
      }

      /* Main loop
       *    - Each pass leaves any failed pivots in place and keeps everything
       *      up-to-date.
       *    - If no pivots selected across matrix, perform swaps to get large
       *      entries into diagonal blocks
       */
      int num_threads = omp_get_max_threads();
      ThreadWork all_thread_work[num_threads];
      // FIXME: Following line is a maximum! Make smaller?
      GlobalWork global_work((nblk*(nblk+1))/2+mblk*nblk);
      run_elim(next_elim, mblk, nblk, cdata, blkdata, global_work, all_thread_work);

      // Calculate number of successful eliminations (removing any dummy cols)
      int num_elim = next_elim;

      // Complete elimination order for anything that wasn't eliminated
      for(int blk=0; blk<nblk; blk++)
         for(int i=cdata[blk].nelim; i<BLOCK_SIZE; i++)
            cdata[blk].elim_order[i] = next_elim++;

      // Store data back in correct permutation
      for(int jblk=0; jblk<nblk; jblk++) {
         // Diagonal block part
         cdata[jblk].permuted_store(d, perm);
         blkdata[jblk*mblk+jblk].permuted_store_diag(a, lda, &cdata[jblk]);
         for(int iblk=jblk+1; iblk<nblk; iblk++)
            blkdata[jblk*mblk+iblk].permuted_store(a, lda, &cdata[iblk], &cdata[jblk]);
         // Rectangular part
         int m2 = m-n;
         T *arect = &a[n];
         for(int iblk=0; iblk<mblk-nblk; iblk++)
            blkdata[jblk*mblk+nblk+iblk].col_permuted_store(
               ((iblk+1)*BLOCK_SIZE<m2) ? BLOCK_SIZE : (m2 - iblk*BLOCK_SIZE),
               &arect[iblk*BLOCK_SIZE], lda, &cdata[jblk]);
      }

      // FIXME: debug? calculate elimianted array
      bool *eliminated = new bool[nblk*BLOCK_SIZE];
      for(int i=0; i<n; i++) eliminated[i] = false;
      for(int blk=0; blk<nblk; blk++) {
         for(int j=0; j<cdata[blk].nelim; j++)
            eliminated[blk*BLOCK_SIZE+j] = true;
      }

      // Construct elim_order
      // FIXME: just use cdata in store
      int *elim_order = new int[nblk*BLOCK_SIZE];
      int j=0;
      for(int blk=0; blk<nblk; blk++)
         for(int i=0; i<BLOCK_SIZE; i++)
            if(cdata[blk].elim_order[i] >= 0)
               elim_order[cdata[blk].elim_order[i]] = j++;

      if(debug) {
         printf("FINAL:\n");
         print_mat(m, n, perm, eliminated, a, lda);
         printf("elim order:");
         for(int i=0; i<n; i++) printf(" %d", elim_order[i]);
         printf("\n");
      }

      // Permute to match elim_order [FIXME inefficient!]
      //permute_a_d_to_elim_order(n, elim_order, a, lda, d, perm);
      
      // Free memory
      delete[] eliminated;
      delete[] elim_order;
      delete[] cdata;
      for(int i=0; i<mblk*nblk; i++)
         std::allocator_traits<BlockDataAlloc>::destroy(
               bdalloc, &blkdata[i]
               );
      std::allocator_traits<BlockDataAlloc>::deallocate(
            bdalloc, blkdata, mblk*nblk
            );

      return num_elim;
   }
};

template <typename T>
void ldlt_solve_fwd(int m, int n, T const* l, int ldl, int nrhs, T* x, int ldx) {
   if(nrhs==1) {
      host_trsv(FILL_MODE_LWR, OP_N, DIAG_UNIT, n, l, ldl, x, 1);
      if(m > n)
         gemv(OP_N, m-n, n, -1.0, &l[n], ldl, x, 1, 1.0, &x[n], 1);
   } else {
      host_trsm(SIDE_LEFT, FILL_MODE_LWR, OP_N, DIAG_UNIT, n, nrhs, 1.0, l, ldl, x, ldx);
      if(m > n)
         host_gemm(OP_N, OP_N, m-n, nrhs, n, -1.0, &l[n], ldl, x, ldx, 1.0, &x[n], ldx);
   }
}

template <typename T>
void ldlt_solve_diag(int n, T const* d, T* x) {
   for(int i=0; i<n; ) {
      if(d[2*i+1]==0.0) {
         // 1x1 pivot
         T d11 = d[2*i];
         x[i] *= d11;
         i++;
      } else {
         // 2x2 pivot
         T d11 = d[2*i];
         T d21 = d[2*i+1];
         T d22 = d[2*i+3];
         T x1 = x[i];
         T x2 = x[i+1];
         x[i]   = d11*x1 + d21*x2;
         x[i+1] = d21*x1 + d22*x2;
         i += 2;
      }
   }
}

template <typename T>
void ldlt_solve_bwd(int m, int n, T const* l, int ldl, int nrhs, T* x, int ldx) {
   if(nrhs==1) {
      if(m > n)
         gemv(OP_T, m-n, n, -1.0, &l[n], ldl, &x[n], 1, 1.0, x, 1);
      host_trsv(FILL_MODE_LWR, OP_T, DIAG_UNIT, n, l, ldl, x, 1);
   } else {
      if(m > n)
         host_gemm(OP_T, OP_N, n, nrhs, m-n, -1.0, &l[n], ldl, &x[n], ldx, 1.0, x, ldx);
      host_trsm(SIDE_LEFT, FILL_MODE_LWR, OP_T, DIAG_UNIT, n, nrhs, 1.0, l, ldl, x, ldx);
   }
}

}}} /* namespaces spral::ssids::cpu */
