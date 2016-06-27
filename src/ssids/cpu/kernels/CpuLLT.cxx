/* COPYRIGHT (c) 2014-5 The Science and Technology Facilities Council (STFC)
 *
 * This code is NOT under the BSD licence. (We are currently evaluating the
 * release path for this code).
 * To licence this code please contact hsl@stfc.ac.uk.
 */

#pragma once

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <fstream>
#include <ostream>
#include <sstream>

#include <immintrin.h>
//#include "/usr/local/src/iaca-lin64/include/iacaMarks.h" // FIXME: remove debug

#include "common.hxx"
#include "CpuLog.hxx"
#include "wrappers.hxx"

namespace spral { namespace ssids { namespace cpu {

template<typename T,
         enum cpu_arch arch=CPU_BEST_ARCH,
         int BLOCK_SIZE=16,
         bool debug=false,
         bool LOG=false
         >
class CpuLLT {
public:
   CpuLLT()
   : log(LOG?10000:0)
   {}

   /** Factorize an entire matrix */
   int factor(int m, int n, T *a, int lda) {
      /* Sanity check arguments */
      if(m < n) return -1;
      if(lda < n) return -4;

      /* Initialize useful quantities */
      int mblk = (m-1) / BLOCK_SIZE + 1;
      int nblk = (n-1) / BLOCK_SIZE + 1;

      /* Load data block-wise */
      void *aligned_ptr;
      if(posix_memalign(&aligned_ptr, 32, mblk*nblk*sizeof(struct block_data)))
         throw new std::bad_alloc();
      struct block_data *blkdata = new(aligned_ptr) struct block_data[mblk*nblk];
      for(int jblk=0; jblk<nblk; jblk++) {
         for(int iblk=0; iblk<mblk; iblk++)
            blkdata[jblk*mblk+iblk].init(
                  ((iblk+1)*BLOCK_SIZE<m) ? BLOCK_SIZE : (m - iblk*BLOCK_SIZE),
                  ((jblk+1)*BLOCK_SIZE<n) ? BLOCK_SIZE : (n - jblk*BLOCK_SIZE),
                  &a[jblk*BLOCK_SIZE*lda + iblk*BLOCK_SIZE], lda
                  );
      }

      /* Main call */
      int flag;
      #pragma omp parallel default(shared)
      #pragma omp single
      {
         flag = run_elim(m, n, blkdata);
      }

      if(LOG) {
         std::ostringstream filename;
         filename << "dpotrf" << nblk << ".prof.fig";
         std::ofstream proffile;
         proffile.open(filename.str().c_str());
         log.writeFig(proffile, getTaskName);
         proffile.close();
      }

      if(debug) {
         printf("on return from run_elim():\n");
         print_mat(mblk, nblk, flag, blkdata);
      }

      // Store data back in correct permutation
      for(int jblk=0; jblk<nblk; jblk++) {
         blkdata[jblk*mblk+jblk].store_diag(
               ((jblk+1)*BLOCK_SIZE<m) ? BLOCK_SIZE : (m - jblk*BLOCK_SIZE),
               ((jblk+1)*BLOCK_SIZE<n) ? BLOCK_SIZE : (n - jblk*BLOCK_SIZE),
               &a[jblk*BLOCK_SIZE*(lda+1)], lda);
         for(int iblk=jblk+1; iblk<mblk; iblk++) {
            blkdata[jblk*mblk+iblk].store(
                  ((iblk+1)*BLOCK_SIZE<m) ? BLOCK_SIZE : (m - iblk*BLOCK_SIZE),
                  ((jblk+1)*BLOCK_SIZE<n) ? BLOCK_SIZE : (n - jblk*BLOCK_SIZE),
                  &a[jblk*BLOCK_SIZE*lda+iblk*BLOCK_SIZE], lda);
         }
      }

      if(debug) {
         printf("FINAL:\n");
         print_mat(m, n, flag, a, lda);
      }
      
      // Free memory
      //delete[] blkdata;
      for(int i=0; i<mblk*nblk; i++) blkdata[i].~block_data();
      free(aligned_ptr);

      return flag;
   }

   /** Factorize a square block without restricting pivots */
   int factor_block(int m, int n, T *__restrict__ a) {
      /* Main loop */
      for(int p=0; p<n; p++) {
         /* Handle diagonal pivot */
         T d11 = a[p*BLOCK_SIZE+p];
         if(d11 <= 0.0) return p+1; // NB use fortran indexing of failed
                                    // pivot to distinguish first pivot
                                    // failure from success!
         d11 = sqrt(d11);
         a[p*BLOCK_SIZE+p] = d11;

         /* Apply pivot */
         d11 = 1.0 / d11;
         T *acol = &a[p*BLOCK_SIZE];
         for(int r=p+1; r<m; r++)
            acol[r] *= d11;

         /* Perform update */
         for(int c=p+1; c<n; c++)
            for(int r=c; r<m; r++)
               a[c*BLOCK_SIZE+r] -= acol[c]*acol[r];
      }

      return 0;
   }

private:
   CpuLog log;

   struct block_data {
      /// Latest accepted value of A or L
      T aval[BLOCK_SIZE*BLOCK_SIZE] __attribute__ ((aligned(32)));
      /// Trial value of L
      T *lwork;

      block_data() {}

      void init(int nrow, int ncol, T *a, int lda) {
         load(nrow, ncol, a, lda);
      }

      /// Loads block of a into this->aval.
      //  If too few rows/cols, pad with NaNs.
      //  Treat too few cols specially: it will be marked as already eliminated
      //  and pushed to the left to avoid trying to eliminate it.
      // FIXME: Is NaN init actually required?
      void load(int nrow, int ncol, T *a, int lda) {
         // Pad missing columns with NaNs
         for(int j=ncol; j<BLOCK_SIZE; j++)
            for(int i=0; i<BLOCK_SIZE; i++)
               aval[j*BLOCK_SIZE+i] = std::numeric_limits<T>::quiet_NaN();
         // Pad missing rows with NaNs
         for(int j=0; j<ncol; j++)
            for(int i=nrow; i<BLOCK_SIZE; i++)
               aval[j*BLOCK_SIZE+i] = std::numeric_limits<T>::quiet_NaN();
         // Load actual data
         for(int j=0; j<ncol; j++)
            for(int i=0; i<nrow; i++)
               aval[j*BLOCK_SIZE+i] = a[j*lda+i];
      }

      void store(int nrow, int ncol, T *__restrict__ a, int lda) __restrict__ const {
         // Store actual data
         for(int j=0; j<ncol; j++)
            for(int i=0; i<nrow; i++)
                a[j*lda+i] = aval[j*BLOCK_SIZE+i];
      }

      void store_diag(int nrow, int ncol, T *__restrict__ a, int lda) __restrict__ const {
         // Store actual data
         for(int j=0; j<ncol; j++)
            for(int i=j; i<nrow; i++)
                a[j*lda+i] = aval[j*BLOCK_SIZE+i];
      }

      /** Performs solve with diagonal block \f$L_{21} = A_{21} L_{11}^{-T}\f$. Designed for below diagonal. */
      void apply_pivot(const T *diag) {
         // Perform solve L_11^-T
         host_trsm<T>(SIDE_RIGHT, FILL_MODE_LWR, OP_T, DIAG_NON_UNIT,
               BLOCK_SIZE, BLOCK_SIZE, 1.0,
               diag, BLOCK_SIZE,
               aval, BLOCK_SIZE
               );
      }

      /** Apply successful pivot update to all uneliminated columns 
       *  (this.aval in non-transpose) */
      void update(const T *__restrict__ l1, const T *__restrict__ l2) {
         host_gemm<T>(OP_N, OP_T, BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE, -1.0, l1, BLOCK_SIZE, l2, BLOCK_SIZE, 1.0, aval, BLOCK_SIZE);
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

   static
   std::string getTaskName(int id) {
      switch(id) {
         case 0: return "Factor";
         case 1: return "Apply";
         case 2: return "Update";
         default: return "Unknown";
      }
   }

   int run_elim(const int m, const int n, struct block_data *blkdata) {
      const int mblk = (m-1) / BLOCK_SIZE + 1;
      const int nblk = (n-1) / BLOCK_SIZE + 1;
      // Calculate rows/cols in final block row/col
      const int mfinal = (m%BLOCK_SIZE==0) ? BLOCK_SIZE : m % BLOCK_SIZE;
      const int nfinal = (n%BLOCK_SIZE==0) ? BLOCK_SIZE : n %BLOCK_SIZE;
      int npass = n;
      for(int blk=0; blk<nblk; blk++) {

         if(debug) {
            printf("Bcol %d:\n", blk);
            print_mat(mblk, nblk, blk*BLOCK_SIZE, blkdata);
         }

         // Factor diagonal 
         struct block_data *dblk = &blkdata[blk*mblk+blk];
         #pragma omp task default(none) \
            firstprivate(dblk, blk) \
            shared(npass, blkdata) \
            depend(inout: blkdata[blk*mblk+blk:1])
         if (npass>=n) {
            CpuLog::LogTask *t;
            if(LOG) t = &log.tstart(0, blk, blk);
            int flag = factor_block(
                  (blk==mblk-1) ? mfinal : BLOCK_SIZE,
                  (blk==nblk-1) ? nfinal : BLOCK_SIZE,
                  dblk->aval
                  );
            // NB: No need for atomic on next line as we will detect failure
            //     before start of each task
            if(flag) npass = blk*BLOCK_SIZE + flag;
            if(LOG) t->end();
         }

         // Loop over off-diagonal blocks applying pivot
         for(int iblk=blk+1; iblk<mblk; iblk++) {
            #pragma omp task default(none) \
               firstprivate(dblk, blk, iblk) \
               shared(blkdata, npass) \
               depend(in: blkdata[blk*mblk+blk:1]) \
               depend(inout: blkdata[blk*mblk+iblk:1])
            if (npass>=n) {
               CpuLog::LogTask *t;
               if(LOG) t = &log.tstart(1, iblk, blk);
               struct block_data *rblk = &blkdata[blk*mblk+iblk];
               // Perform necessary operations
               rblk->apply_pivot(dblk->aval);
               if(LOG) t->end();
            }
         }
         
         // Update uneliminated columns
         for(int jblk=blk+1; jblk<nblk; jblk++) {
            for(int iblk=jblk; iblk<mblk; iblk++) {
               #pragma omp task default(none) \
                  firstprivate(blk, iblk, jblk) \
                  shared(blkdata, npass) \
                  depend(in: blkdata[blk*mblk+iblk:1]) \
                  depend(in: blkdata[blk*mblk+jblk:1]) \
                  depend(inout: blkdata[jblk*mblk+iblk:1])
               if (npass>=n) {
                  CpuLog::LogTask *t;
                  if(LOG) t = &log.tstart(2, iblk, jblk, blk);
                  blkdata[jblk*mblk+iblk].update(
                        blkdata[blk*mblk+iblk].aval,
                        blkdata[blk*mblk+jblk].aval
                        );
                  if(LOG) t->end();
               }
            }
         }
      }
      #pragma omp taskwait

      return npass; // Everything passed
   }

   static
   void print_mat(int m, int n, int nelim, const T *a, int lda) {
      for(int row=0; row<m; row++) {
         if(row < n)
            printf("%d%s:", row, (row<nelim)?"X":" ");
         else
            printf("%d%s:", row, "U");
         for(int col=0; col<n; col++)
            printf(" %10.4f", a[col*lda+row]);
         printf("\n");
      }
   }

   static
   void print_mat(int mblk, int nblk, int nelim, const struct block_data *blkdata) {
      for(int rblk=0; rblk<mblk; rblk++) {
         for(int row=0; row<BLOCK_SIZE; row++) {
            int r = rblk*BLOCK_SIZE+row;
            if(r < nblk*BLOCK_SIZE)
               printf("%d%s:", r, (r<nelim)?"X":" ");
            else
               printf("%d%s:", r, "U");
            for(int cblk=0; cblk<nblk; cblk++) {
               const struct block_data *blk = &blkdata[cblk*mblk+rblk];
               for(int col=0; col<BLOCK_SIZE; col++)
                  printf(" %10.4f", blk->aval[col*BLOCK_SIZE+row]);
            }
            printf("\n");
         }
      }
   }


};

}}} /* namespaces spral::ssids::cpu */
