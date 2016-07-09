#include "framework.hxx"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

/** Generates a random dense positive definte matrix. Off diagonal entries are
 * Unif[-1,1]. Each diagonal entry a_ii = Unif[0.1,1.1] + sum_{i!=j} |a_ij|. */
void gen_posdef(int n, double* a, int lda) {
   /* Fill matrix with random numbers from Unif [-1.0,1.0] */
   for(int j=0; j<n; ++j)
   for(int i=j; i<n; ++i)
      a[j*lda+i] = 1.0 - (2.0*rand()) / RAND_MAX ;
   /* Make diagonally dominant */
   for(int i=0; i<n; ++i) a[i*lda+i] = fabs(a[i*lda+i]) + 0.1;
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

void print_vec(char const* format, int n, double const* vec) {
   for(int i=0; i<n; ++i)
      printf(format, vec[i]);
   printf("\n");
}

void print_mat(char const* format, int n, double const* a, int lda) {
   for(int i=0; i<n; ++i) {
      printf("%d:", i);
      for(int j=0; j<=i; ++j)
         printf(format, a[j*lda+i]);
      printf("\n");
   }
}

/** Calculate scaled backward error ||Ax-b|| / ( ||A|| ||x|| + ||b|| ).
 * All norms are infinity norms. */
double backward_error(int n, double const* a, int lda, double const* rhs, int nrhs, double const* soln, int ldsoln) {
   /* Allocate memory */
   double *resid = new double[n];
   double *rowsum = new double[n];

   /* Calculate residual vector and anorm*/
   double worstbwderr = 0.0;
   for(int r=0; r<nrhs; ++r) {
      memcpy(resid, rhs, n*sizeof(double));
      memset(rowsum, 0, n*sizeof(double));
      for(int j=0; j<n; ++j) {
         resid[j] -= a[j*lda+j] * soln[r*ldsoln+j];
         rowsum[j] += fabs(a[j*lda+j]);
         for(int i=j+1; i<n; ++i) {
            resid[j] -= a[j*lda+i] * soln[r*ldsoln+i];
            resid[i] -= a[j*lda+i] * soln[r*ldsoln+j];
            rowsum[j] += fabs(a[j*lda+i]);
            rowsum[i] += fabs(a[j*lda+i]);
         }
      }
      double anorm = 0.0;
      for(int i=0; i<n; ++i)
         anorm = std::max(anorm, rowsum[i]);

      /* Check scaled backwards error */
      double rhsnorm=0.0, residnorm=0.0, solnnorm=0.0;
      for(int i=0; i<n; ++i) {
         rhsnorm = std::max(rhsnorm, fabs(rhs[i]));
         residnorm = std::max(residnorm, fabs(resid[i]));
         solnnorm = std::max(solnnorm, fabs(r*ldsoln+soln[i]));
      }

      worstbwderr = std::max(worstbwderr, residnorm/(anorm*solnnorm + rhsnorm));
   }

   /* Cleanup */
   delete[] resid;
   delete[] rowsum;

   /* Return result */
   //printf("%e / %e %e %e\n", residnorm, anorm, solnnorm, rhsnorm);
   return worstbwderr;
}

/** Calculates forward error ||soln-x||_inf assuming x=1.0 */
double forward_error(int n, int nrhs, double const* soln, int ldx) {
   /* Check scaled backwards error */
   double fwderr=0.0;
   for(int r=0; r<nrhs; ++r)
   for(int i=0; i<n; ++i)
      fwderr = std::max(fwderr, fabs(soln[r*ldx+i] - 1.0));
   return fwderr;
}
