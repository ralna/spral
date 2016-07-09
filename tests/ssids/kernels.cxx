#include "kernels/framework.hxx"
#include "ssids/cpu/kernels/cholesky.hxx"
#include "ssids/cpu/kernels/ldlt_nopiv.hxx"
#include "ssids/cpu/kernels/wrappers.hxx"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

using namespace spral::ssids::cpu;

int test_cholesky(int m, int n, int blksz, bool debug=false) {
   /* Generate random dense posdef matrix of size m */
   int lda = m;
   double *a = new double[m*lda];
   gen_posdef(m, a, lda);
   /* Take a copy */
   double *l = new double[m*lda];
   memcpy(l, a, m*lda*sizeof(double));
   
   /* Factor first m x n part with our code */
   if(debug) { printf("PRE:\n"); print_mat(" %e", m, l, lda); }
   int info;
   #pragma omp parallel default(shared)
   {
      #pragma omp single
      {
         cholesky_factor(m, n, l, lda, blksz, &info);
      }
   } /* implicit task wait on exit from parallel region */
   if(debug) { printf("POST:\n"); print_mat(" %e", m, l, lda); }
   if(m>n) {
      /* Schur complement update remainder block */
      host_gemm<double>(OP_N, OP_T, m-n, m-n, n, -1.0, &l[n], lda, &l[n], m, 1.0, &l[n*lda+n], lda);
      if(debug) { printf("post schur:\n"); print_mat(" %e", m, l, lda); }
      /* Factor remaining part using LAPACK Cholesky factorization */
      lapack_potrf<double>(FILL_MODE_LWR, m-n, &l[n*lda+n], lda);
      if(debug) { printf("post potrf:\n"); print_mat(" %e", m, l, lda); }
   }

   /* Generate a rhs corresponding to x=1.0 */
   double *rhs = new double[m];
   memset(rhs, 0, m*sizeof(double));
   for(int j=0; j<m; ++j) {
      rhs[j] += a[j*lda+j] * 1.0;
      for(int i=j+1; i<m; ++i) {
         rhs[j] += a[j*lda+i] * 1.0;
         rhs[i] += a[j*lda+i] * 1.0;
      }
   }

   /* Perform solves with 1 and 3 rhs */
   int nrhs = 3;
   int ldsoln = m+3;
   double *soln = new double[(nrhs+1)*ldsoln];
   for(int r=0; r<nrhs+1; ++r)
      memcpy(&soln[r*ldsoln], rhs, m*sizeof(double));
   if(debug) { printf("rhs ="); print_vec(" %e", m, soln); }
   cholesky_solve_fwd(m, n, l, lda, 1, soln, ldsoln);
   cholesky_solve_fwd(m, n, l, lda, nrhs, &soln[ldsoln], ldsoln);
   if(debug) { printf("post fwd ="); print_vec(" %e", m, soln); }
   host_trsm<double>(SIDE_LEFT, FILL_MODE_LWR, OP_N, DIAG_NON_UNIT, m-n, nrhs+1, 1.0, &l[n*lda+n], lda, &soln[n], ldsoln);
   host_trsm<double>(SIDE_LEFT, FILL_MODE_LWR, OP_T, DIAG_NON_UNIT, m-n, nrhs+1, 1.0, &l[n*lda+n], lda, &soln[n], ldsoln);
   cholesky_solve_bwd(m, n, l, lda, 1, soln, ldsoln);
   cholesky_solve_bwd(m, n, l, lda, nrhs, &soln[ldsoln], ldsoln);
   if(debug) { printf("post bwd ="); print_vec(" %e", m, soln); }

   double fwderr = forward_error(m, nrhs, soln, ldsoln);
   double bwderr = backward_error(m, a, lda, rhs, nrhs, soln, ldsoln);

   if(debug) printf("fwderr = %e\nbwderr = %e\n", fwderr, bwderr);

   /* Cleanup memory */
   delete[] a;
   delete[] l;
   delete[] rhs;
   delete[] soln;

   if(bwderr >= 1e-14 || std::isnan(bwderr)) return -1; // Failed accuracy test

   return 0; // Test passed
}

int test_ldlt(int m, int n, bool debug=false) {
   /* Generate random dense posdef matrix of size m */
   int lda = m;
   double *a = new double[m*lda];
   gen_posdef(m, a, lda);
   /* Take a copy */
   double *l = new double[m*lda];
   memcpy(l, a, m*lda*sizeof(double));
   
   /* Factor first m x n part with our code */
   double *work = new double[2*m];
   if(debug) { printf("PRE:\n"); print_mat(" %e", m, l, lda); }
   ldlt_nopiv_factor(m, n, l, lda, work);
   if(debug) { printf("POST:\n"); print_mat(" %e", m, l, lda); }
   if(m>n) {
      /* Schur complement update remainder block */
      for(int j=0; j<n-1; j+=2) {
         /* Calculate D (we store D^-1) */
         double a11 = l[    j*lda+  j];
         double a21 = l[    j*lda+j+1];
         double a22 = l[(j+1)*lda+j+1];
         double det = a11*a22 - a21*a21;
         double d11 = a22/det;
         double d21 = -a21/det;
         double d22 = a11/det;
         /* Find LD for current 2x2 column*/
         for(int i=n; i<m; ++i) {
            work[  i] = d11*l[j*lda+i] + d21*l[(j+1)*lda+i];
            work[m+i] = d21*l[j*lda+i] + d22*l[(j+1)*lda+i];
         }
         /* Apply schur complement from this 2x2 block */
         host_gemm<double>(OP_N, OP_T, m-n, m-n, 2, -1.0, &l[j*lda+n], lda, &work[n], m, 1.0, &l[n*lda+n], lda);
      }
      if(n%2!=0) {
         /* n is odd, last column is a 1x1 pivot */
         int j = n-1;
         double d11 = 1/l[j*lda+j];
         for(int i=n; i<m; ++i)
            work[i] = d11*l[j*lda+i];
         /* Apply schur complement from this 1x1 block */
         host_gemm<double>(OP_N, OP_T, m-n, m-n, 1, -1.0, &l[j*lda+n], lda, &work[n], m, 1.0, &l[n*lda+n], lda);
      }
      if(debug) { printf("post schur:\n"); print_mat(" %e", m, l, lda); }
      /* Factor remaining part using LAPACK Cholesky factorization */
      lapack_potrf<double>(FILL_MODE_LWR, m-n, &l[n*lda+n], lda);
      if(debug) { printf("post potrf:\n"); print_mat(" %e", m, l, lda); }
   }

   /* Generate a rhs corresponding to x=1.0 */
   double *rhs = new double[m];
   memset(rhs, 0, m*sizeof(double));
   for(int j=0; j<m; ++j) {
      rhs[j] += a[j*lda+j] * 1.0;
      for(int i=j+1; i<m; ++i) {
         rhs[j] += a[j*lda+i] * 1.0;
         rhs[i] += a[j*lda+i] * 1.0;
      }
   }

   /* Perform a solve */
   double *soln = new double[m];
   memcpy(soln, rhs, m*sizeof(double));
   if(debug) { printf("rhs ="); print_vec(" %e", m, soln); }
   ldlt_nopiv_solve_fwd(m, n, l, lda, soln);
   if(debug) { printf("post fwd ="); print_vec(" %e", m, soln); }
   host_trsv<double>(FILL_MODE_LWR, OP_N, DIAG_NON_UNIT, m-n, &l[n*lda+n], lda, &soln[n], 1);
   host_trsv<double>(FILL_MODE_LWR, OP_T, DIAG_NON_UNIT, m-n, &l[n*lda+n], lda, &soln[n], 1);
   ldlt_nopiv_solve_diag(m, n, l, lda, soln);
   if(debug) { printf("post diag ="); print_vec(" %e", m, soln); }
   ldlt_nopiv_solve_bwd(m, n, l, lda, soln);
   if(debug) { printf("post bwd ="); print_vec(" %e", m, soln); }

   double fwderr = forward_error(m, 1, soln, m);
   double bwderr = backward_error(m, a, lda, rhs, 1, soln, m);

   if(debug) printf("fwderr = %e\nbwderr = %e\n", fwderr, bwderr);

   /* Cleanup memory */
   delete[] a;
   delete[] l;
   delete[] work;
   delete[] rhs;
   delete[] soln;

   if(bwderr >= 1e-14 || std::isnan(bwderr)) return -1; // Failed accuracy test

   return 0; // Test passed
}

int main(void) {
   int nerr = 0;

   /* Cholesky tests (m, n, blksz) */
#if 1
   TEST(test_cholesky(1, 1, 1));
   TEST(test_cholesky(2, 2, 2));
   TEST(test_cholesky(2, 2, 4));
   TEST(test_cholesky(2, 2, 1));
   TEST(test_cholesky(3, 3, 1));
   TEST(test_cholesky(4, 4, 2));
   TEST(test_cholesky(5, 5, 2));
   TEST(test_cholesky(2, 1, 1));
   TEST(test_cholesky(3, 2, 1));
   TEST(test_cholesky(6, 4, 2));
   TEST(test_cholesky(500, 234, 32));
   TEST(test_cholesky(733, 231, 32));
   TEST(test_cholesky(733, 231, 19));
#endif

   /* LDL^T no pivoting tests (m, n) */
#if 1
   TEST(test_ldlt(1, 1));
   TEST(test_ldlt(2, 2));
   TEST(test_ldlt(3, 3));
   TEST(test_ldlt(4, 4));
   TEST(test_ldlt(5, 5));
   TEST(test_ldlt(4, 2));
   TEST(test_ldlt(5, 3));
   TEST(test_ldlt(6, 4));
   TEST(test_ldlt(500, 234));
   TEST(test_ldlt(733, 231));
#endif

   return nerr;
}
