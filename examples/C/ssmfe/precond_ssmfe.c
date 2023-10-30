/* examples/C/ssmfe/precond_ssmfe.c */
/* Laplacian on a square grid (using SPRAL_SSMFE routines) */
#include "spral.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>

/* Header that implements Laplacian and preconditioners */
#include "laplace2d.h"

int main(void) {
   const int m   = 20;     /* grid points along each side */
   const int n   = m*m;    /* problem size */
   const int nep = 5;      /* eigenpairs wanted */

   double *lambda = malloc(2*nep * sizeof(*lambda)); /* eigenvalues */
   double (*X)[2*nep][n] = malloc(sizeof(*X));       /* eigenvectors */
   struct spral_ssmfe_rcid rci;           /* reverse communication data */
   struct spral_ssmfe_options options;    /* options */
   void *keep;                            /* private data */
   struct spral_ssmfe_inform inform;      /* information */

   /* Initialize options to default values */
   spral_ssmfe_default_options(&options);
   /* gap between the last converged eigenvalue and the rest of the spectrum
    * must be at least 0.1 times average gap between computed eigenvalues */
   options.left_gap = -0.1;

   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      spral_ssmfe_standard_double(&rci, nep, 2*nep, lambda, n, &(*X)[0][0], n,
         &keep, &options, &inform);
      switch ( rci.job ) {
      case 1:
         apply_laplacian(m, m, rci.nx, rci.x, rci.y);
         break;
      case 2:
         apply_gauss_seidel_step(m, m, rci.nx, rci.x, rci.y);
         break;
      default:
         goto finished;
      }
   }
finished:
   printf("%d eigenpairs converged in %d iterations\n", inform.left, inform.iteration);
   for(int i=0; i<inform.left; i++)
      printf(" lambda[%1d] = %13.7e\n", i, lambda[i]);
   spral_ssmfe_free_double(&keep, &inform);
   free(lambda);
   free(X);

   /* Success */
   return 0;
}
