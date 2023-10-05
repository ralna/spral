/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * Licence: BSD licence, see LICENCE file for details
 *
 */
#include "ldlt_nopiv.hxx"

#include <cmath>
#include <cstring>

#include "framework.hxx"
#include "ssids/cpu/kernels/ldlt_nopiv.hxx"
#include "ssids/cpu/kernels/wrappers.hxx"

using namespace spral::ssids::cpu;

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

   if(bwderr >= 1e-14 || std::isnan(bwderr)) return 1; // Failed accuracy test

   return 0; // Test passed
}

int run_ldlt_nopiv_tests() {
   int nerr = 0;

   /* LDL^T no pivoting tests (m, n) */
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

   return nerr;
}
