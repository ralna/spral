/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * Licence: BSD licence, see LICENCE file for details
 *
 */
#include "cholesky.hxx"

#include <cmath>
#include <cstring>

#include "framework.hxx"
#include "ssids/cpu/kernels/cholesky.hxx"
#include "ssids/cpu/kernels/wrappers.hxx"

using namespace spral::ssids::cpu;

int test_cholesky(int m, int n, int blksz, bool debug=false) {
   /* Generate random dense posdef matrix of size m */
   int lda = m;
   double *a = new double[m*lda];
   gen_posdef(m, a, lda);
   /* Take a copy */
   double *l = new double[m*lda];
   memcpy(l, a, m*lda*sizeof(double));

   /* Factor first m x n part with our code, generate contrib block */
   if(debug) { printf("PRE:\n"); print_mat(" %e", m, l, lda); }
   int info;
   #pragma omp parallel default(shared)
   {
      #pragma omp single
      {
         cholesky_factor(m, n, l, lda, 1.0, &l[n*lda+n], lda, blksz, &info);
      }
   } /* implicit task wait on exit from parallel region */
   if(debug) { printf("POST:\n"); print_mat(" %e", m, l, lda); }
   if(m>n) {
      /* Schur complement update remainder block */
      //host_gemm<double>(OP_N, OP_T, m-n, m-n, n, -1.0, &l[n], lda, &l[n], m, 1.0, &l[n*lda+n], lda);
      if(debug) { printf("post schur:\n"); print_mat(" %e", m, l, lda); }
      /* Factor remaining part using LAPACK Cholesky factorization */
      lapack_potrf<double>(FILL_MODE_LWR, m-n, &l[n*lda+n], lda);
      if(debug) { printf("post potrf:\n"); print_mat(" %e", m, l, lda); }
   }

   /* Generate a rhs corresponding to x=1.0 */
   double *rhs = new double[m];
   gen_rhs(m, a, lda, rhs);

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

   if(bwderr >= 1e-14 || std::isnan(bwderr)) return 1; // Failed accuracy test

   return 0; // Test passed
}

int run_cholesky_tests() {
   int nerr=0;

   /* Cholesky tests (m, n, blksz) */
   TEST(test_cholesky(1, 1, 1));
   TEST(test_cholesky(2, 2, 2));
   TEST(test_cholesky(2, 2, 4));
   TEST(test_cholesky(2, 2, 1));
   TEST(test_cholesky(3, 3, 1));
   TEST(test_cholesky(4, 4, 2));
   TEST(test_cholesky(5, 5, 2));
   TEST(test_cholesky(2, 1, 1));
   TEST(test_cholesky(3, 2, 1));
   TEST(test_cholesky(5, 1, 2));
   TEST(test_cholesky(5, 3, 2));
   TEST(test_cholesky(6, 4, 2));
   TEST(test_cholesky(500, 234, 32));
   TEST(test_cholesky(733, 231, 32));
   TEST(test_cholesky(733, 231, 19));
   TEST(test_cholesky(1668, 204, 256));

   return nerr;
}
