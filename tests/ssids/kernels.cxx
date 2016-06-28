#include "ssids/cpu/kernels/ldlt_nopiv.hxx"
#include "ssids/cpu/kernels/wrappers.hxx"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>

using namespace spral::ssids::cpu;

void gen_posdef(int n, double* a, int lda) {
   /* Fill matrix with random numbers from Unif [-1.0,1.0] */
   for(int j=0; j<n; ++j)
   for(int i=j; i<n; ++i)
      a[j*lda+i] = 1.0 - (2.0*rand()) / RAND_MAX ;
   /* Make diagonally dominant */
   for(int i=0; i<n; ++i) a[i*lda+i] = fabs(a[i*lda+i]);
   for(int j=0; j<n; ++j)
   for(int i=j; i<n; ++i) {
      a[j*lda+j] += fabs(a[j*lda+i]);
      a[i*lda+i] += fabs(a[j*lda+i]);
   }
   /* Fill upper triangle with NaN */
   for(int j=0; j<n; ++j)
   for(int i=0; i<j; ++i)
      a[j*lda+i] = std::numeric_limits<double>::signaling_NaN();
}

int test_ldlt(int m, int n) {
   /* Generate random dense posdef matrix of size m */
   int lda = m;
   double *a = new double[m*lda];
   gen_posdef(m, a, lda);
   /* Take a copy */
   double *l = new double[m*lda];
   memcpy(l, a, m*lda*sizeof(double));
   
   /* Factor first m x n part with our code */
   double *work = new double[2*m];
   ldlt_nopiv_factor(m, n, l, lda, work);
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
      host_gemm<double>(OP_N, OP_T, m-n, m-n, 2, -1.0, &l[j*lda+n], lda, work, m, 1.0, &l[n*lda+n], lda);
   }
   if(n%2!=0) {
      /* n is odd, last column is a 1x1 pivot */
      int j = n-1;
      double d11 = 1/l[j*lda+j];
      for(int i=n; n<m; ++i)
         work[i] = d11*l[j*lda+i];
      /* Apply schur complement from this 1x1 block */
      host_gemm<double>(OP_N, OP_T, m-n, m-n, 1, -1.0, &l[j*lda+n], lda, work, m, 1.0, &l[n*lda+n], lda);
   }
   /* Factor remaining part using LAPACK Cholesky factorization */
   lapack_potrf<double>(FILL_MODE_LWR, m-n, &l[n*lda+n], lda);

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
   ldlt_nopiv_solve_fwd(m, n, l, lda, soln);
   host_trsv<double>(FILL_MODE_LWR, OP_N, DIAG_NON_UNIT, m-n, &l[n*lda+n], lda, &soln[n], 1);
   host_trsv<double>(FILL_MODE_LWR, OP_T, DIAG_NON_UNIT, m-n, &l[n*lda+n], lda, &soln[n], 1);
   ldlt_nopiv_solve_diag(m, n, l, lda, soln);
   ldlt_nopiv_solve_bwd(m, n, l, lda, soln);

   /* Calculate residual vector and anorm*/
   double *resid = new double[m];
   double *rowsum = new double[m];
   memcpy(resid, rhs, m*sizeof(double));
   memset(rowsum, 0, m*sizeof(double));
   for(int j=0; j<m; ++j) {
      resid[j] -= a[j*lda+j] * soln[j];
      rowsum[j] += fabs(a[j*lda+j]);
      for(int i=j+1; i<m; ++i) {
         resid[j] -= a[j*lda+i] * soln[i];
         resid[i] -= a[j*lda+i] * soln[j];
         rowsum[j] += fabs(a[j*lda+i]);
         rowsum[i] += fabs(a[j*lda+i]);
      }
   }
   double anorm = 0.0;
   for(int i=0; i<m; ++i)
      anorm = std::max(anorm, rowsum[i]);

   /* Check scaled backwards error */
   double rhsnorm=0.0, residnorm=0.0, solnnorm=0.0, fwderr=0.0;
   for(int i=0; i<m; ++i) {
      rhsnorm = std::max(rhsnorm, fabs(rhs[i]));
      residnorm = std::max(residnorm, fabs(resid[i]));
      solnnorm = std::max(solnnorm, fabs(soln[i]));
      fwderr = std::max(fwderr, fabs(soln[i] - 1.0));
   }
   //printf("%e / %e %e %e\n", residnorm, anorm, solnnorm, rhsnorm);
   double bwderr = residnorm / (anorm*solnnorm + rhsnorm);

   printf("fwderr = %e\nbwderr = %e\n", fwderr, bwderr);

   if(bwderr >= 1e-14 || isnan(bwderr)) return -1; // Failed accuracy test

   /* Cleanup memory */
   delete[] a;
   delete[] l;
   delete[] work;
   delete[] rhs;
   delete[] soln;
   delete[] resid;
   delete[] rowsum;

   return 0; // Test passed
}

#define ANSI_COLOR_RED     "\x1b[31;1m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_BLUE    "\x1b[34;1m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define TEST(func) \
   if(!func) { \
      printf(ANSI_COLOR_BLUE "%-20s" ANSI_COLOR_RESET " [ " ANSI_COLOR_GREEN "ok" ANSI_COLOR_RESET " ]\n", #func); \
   } else { \
      printf(ANSI_COLOR_BLUE "%-20s" ANSI_COLOR_RESET " [" ANSI_COLOR_RED "fail" ANSI_COLOR_RESET "]\n", #func); \
      nerr++; \
   }

int main(void) {
   int nerr = 0;

   /* LDL^T no pivoting tests */
   TEST(test_ldlt(3, 3));
#if 0
   TEST(test_ldlt(1, 1));
   TEST(test_ldlt(2, 2));
   TEST(test_ldlt(4, 4));
   TEST(test_ldlt(5, 5));
   TEST(test_ldlt(4, 2));
#endif

   return nerr;
}
