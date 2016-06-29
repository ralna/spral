#include "cholesky.hxx"

#include "wrappers.hxx"

namespace spral { namespace ssids { namespace cpu {

/** Perform Cholesky factorization of lower triangular matrix a[] in place.
 *
 * Note that calling this routine just enqueues tasks, these will not
 * necessarily be complete until an omp taskwait construct is encountered.
 *
 * \param m the number of rows
 * \param n the number of columns
 * \param a the matrix to be factorized, only lower triangle is used, however
 *    upper triangle may get overwritten with rubbish
 * \param lda the leading dimension of a
 * \param blksz the block size to use for parallelization. Blocks are aimed to
 *    contain at most blksz**2 entries.
 * \param info is initialized to -1, and will be changed to the index of any
 *    column where a non-zero column is encountered.
 */
void cholesky_factor(int m, int n, double* a, int lda, int blksz, int &info) {
   if(n < blksz) {
      // Adjust so blocks have blksz**2 entries
      blksz = int((long(blksz)*blksz) / n);
   }

   info = -1;

   /* FIXME: Would this be better row-wise to ensure critical path, rather than
    * its current col-wise implementation ensuring maximum work available??? */
   for(int j=0; j<n; j+=blksz) {
      int blkn = std::min(blksz, n-j);
      /* Diagonal Block Factorization Task */
      #pragma omp task default(none) \
         firstprivate(j, blkn) \
         shared(m, a, lda, blksz, info) \
         depend(inout: a[j*(lda+1):1])
      if(info==-1) {
         int blkm = std::min(blksz, m-j);
         int flag = lapack_potrf(FILL_MODE_LWR, blkn, &a[j*(lda+1)], lda);
         if(flag > 0) {
            // Matrix was not positive definite
            info = flag-1; // flag uses Fortran indexing
         } else if(blkm>blkn) {
            // Diagonal block factored OK, handle some rectangular part of block
            host_trsm(SIDE_RIGHT, FILL_MODE_LWR, OP_T, DIAG_NON_UNIT,
                  blkm-blkn, blkn, 1.0, &a[j*(lda+1)], lda,
                  &a[j*(lda+1)+blkn], lda);
         }
      }
      /* Column Solve Tasks */
      for(int i=j+blksz; i<blksz; i+=blksz) {
         int blkm = std::min(blksz, m-i);
         #pragma omp task default(none) \
            firstprivate(i, j, blkn, blkm) \
            shared(a, lda, info) \
            depend(in: a[j*(lda+1):1]) \
            depend(inout: a[j*lda + i:1])
         if(info==-1) {
            host_trsm(SIDE_RIGHT, FILL_MODE_LWR, OP_T, DIAG_NON_UNIT,
                  blkm-blkn, blkn, 1.0, &a[j*(lda+1)], lda, &a[j*lda+i], lda);
         }
      }
      /* Schur Update Tasks */
      for(int k=j+blksz; k<n; k+=blksz) {
         int blkk = std::min(blksz, n-k);
         for(int i=k; i<m; i+=blksz) {
            #pragma omp task default(none) \
               firstprivate(i, j, k, blkn, blkk) \
               shared(n, a, lda, blksz, info) \
               depend(in: a[j*lda+k:1]) \
               depend(in: a[j*lda+i:1]) \
               depend(inout: a[k*lda+i:1])
            if(info==-1) {
               int blkm = std::min(blksz, n-i);
               host_gemm(OP_N, OP_T, blkm, blkk, blkn, -1.0, &a[j*lda+k], lda,
                     &a[j*lda+i], lda, 1.0, &a[k*lda+i], lda);
            }
         }
      }
   }
}

}}} /* namespaces spral::ssids::cpu */
