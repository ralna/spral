/* examples/C/ssmfe/shift_invert.c */
/* Laplacian on a rectangular grid by shift-invert via LDLT factorization */
#include "spral.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>

/* Headers that implements Laplacian and preconditioners and LDLT support */
#include "laplace2d.h"
#include "ldltf.h"

int main(void) {
   const int nx = 8;          /* grid points along x */
   const int ny = 8;          /* grid points along y */
   const int n = nx*ny;       /* problem size */
   const double sigma = 1.0;  /* shift */

   int *ipiv = malloc(n * sizeof(*ipiv));        /* LDLT pivot index */
   double *lambda = malloc(n * sizeof(*lambda)); /* eigenvalues */
   double (*X)[n][n] = malloc(sizeof(*X));       /* eigenvectors */
   double (*A)[n][n] = malloc(sizeof(*A));       /* matrix */
   double (*LDLT)[n][n] = malloc(sizeof(*LDLT)); /* factors */
   double (*work)[n][n] = malloc(sizeof(*work)); /* work array for dsytrf */
   struct spral_ssmfe_options options;    /* eigensolver options */
   struct spral_ssmfe_inform inform;      /* information */
   struct spral_ssmfe_rcid rci;           /* reverse communication data */
   void *keep;                            /* private data */

   /* Initialize options to default values */
   spral_ssmfe_default_options(&options);

   /* Set up then perform LDLT factorization of the shifted matrix */
   set_laplacian_matrix(nx, ny, n, A);
   for(int j=0; j<n; j++)
   for(int i=0; i<n; i++)
      (*LDLT)[j][i] = (i==j) ? (*A)[j][i] - sigma
                             : (*A)[j][i];
   cwrap_dsytrf('L', n, &(*LDLT)[0][0], n, ipiv, &(*work)[0][0], n*n);

   /* Main loop */
   int left = num_neg_D(n, n, LDLT, ipiv);   /* all evalues to left of sigma */
   int right = 5;                            /* 5 evalues to right of sigma */
   rci.job = 0; keep = NULL;
   while(true) {
      spral_ssmfe_standard_shift_double(&rci, sigma, left, right, n, lambda,
            n, &(*X)[0][0], n, &keep, &options, &inform);
      switch( rci.job ) {
      case 1:
         cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, n,
            1.0, &(*A)[0][0], n, rci.x, n, 0.0, rci.y, n);
         break;
      case 2:
         // No preconditioning
         break;
      case 9:
         cblas_dcopy(n*rci.nx, rci.x, 1, rci.y, 1);
         cwrap_dsytrs('L', n, rci.nx, &(*LDLT)[0][0], n, ipiv, rci.y, n);
         break;
      default:
         goto finished;
      }
   }
finished:
   printf("Eigenvalues near %e (took %d iterations)\n", sigma, inform.iteration);
   for(int i=0; i<inform.left+inform.right; i++)
      printf(" lambda[%1d] = %13.7e\n", i, lambda[i]);
   spral_ssmfe_free_double(&keep, &inform);
   free(ipiv);
   free(lambda);
   free(X);
   free(A);
   free(LDLT);
   free(work);

   /* Success */
   return 0;
}
