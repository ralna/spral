/* examples/C/ssmfe/hermitian.c - Example code for SPRAL_SSMFE package */
/* Hermitian operator example */
#include "spral.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>

/* central differences for i d/dx */
void apply_idx(int n, int m, const double complex *x_ptr, double complex *y_ptr) {
   /* Use "variable-modified types" to simplify matrix indexing */
   const double complex (*x)[m][n] = (const double complex (*)[m][n]) x_ptr;
   double complex (*y)[m][n] = (double complex (*)[m][n]) y_ptr;
   for(int j=0; j<m; j++) {
      for(int i=0; i<n; i++) {
         int il = (i==0)   ? n-1 : i-1;
         int ir = (i==n-1) ? 0   : i+1;
         (*y)[j][i] = _Complex_I * ((*x)[j][ir] - (*x)[j][il]);
      }
   }
}

/* main routine */
int main(void) {
   const int n   = 80;                 /* problem size */
   const int nep = 5;                  /* eigenpairs wanted */

   double *lambda = malloc(nep * sizeof(*lambda));   /* eigenvalues */
   double complex (*X)[nep][n] = malloc(sizeof(*X)); /* eigenvectors */
   struct spral_ssmfe_rciz rci;        /* reverse communication data */
   struct spral_ssmfe_options options; /* options */
   void *keep;                         /* private data */
   struct spral_ssmfe_inform inform;   /* information */

   /* Initialize options to default values */
   spral_ssmfe_default_options(&options);

   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      spral_ssmfe_standard_double_complex(&rci, nep, nep, lambda, n,
         &(*X)[0][0], n, &keep, &options, &inform);
      switch ( rci.job ) {
      case 1:
         apply_idx(n, rci.nx, rci.x, rci.y);
         break;
      case 2:
         // No preconditioning
         break;
      default:
         goto finished;
      }
   }
finished:
   printf("%d eigenpairs converged in %d iterations\n", inform.left, inform.iteration);
   for(int i=0; i<inform.left; i++)
      printf(" lambda[%1d] = %13.7e\n", i, lambda[i]);
   spral_ssmfe_free_double_complex(&keep, &inform);
   free(lambda);
   free(X);

   /* Success */
   return 0;
}
