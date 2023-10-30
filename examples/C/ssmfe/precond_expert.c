/* examples/C/ssmfe/precond_expert.c */
/* Laplacian on a square grid (using SPRAL_SSMFE_EXPERT routines) */
#include "spral.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>

/* Header that implements Laplacian and preconditioners */
#include "laplace2d.h"

int main(void) {
   const int ngrid = 20;          /* grid points along each side */
   const int n     = ngrid*ngrid; /* problem size */
   const int nep   = 5;           /* eigenpairs wanted */
   const int m     = 3;           /* dimension of the iterated subspace */

   int state = SPRAL_RANDOM_INITIAL_SEED; /* PRNG state */

   int *ind = malloc(m * sizeof(*ind));          /* permutation index */
   double *lambda = malloc(n * sizeof(*lambda)); /* eigenvalues */
   double (*X)[n][n] = malloc(sizeof(*X));       /* eigenvectors */
   /* Work arrays */
   double (*rr)[3][2*m][2*m] = malloc(sizeof(*rr));
   double (*W)[6][m][n] = malloc(sizeof(*W));
   double (*U)[m][n] = malloc(sizeof(*U));

   /* Derived types */
   struct spral_ssmfe_rcid rci;            /* reverse communication data */
   struct spral_ssmfe_options options;     /* options */
   void *keep;                             /* private data */
   struct spral_ssmfe_inform inform;       /* information */

   /* Initialize options to default values */
   spral_ssmfe_default_options(&options);
   /* gap between the last converged eigenvalue and the rest of the spectrum
    * must be at least 0.1 times average gap between computed eigenvalues */
   options.left_gap = -0.1;
   /* block size is small, convergence may be slow, allow more iterations */
   options.max_iterations = 200;

   /* Initialize W to lin indep vectors by randomizing */
   for(int i=0; i<n; i++)
   for(int j=0; j<m; j++)
      (*W)[0][j][i] = spral_random_real(&state, true);

   int ncon = 0;   /* number of converged eigenpairs */
   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      spral_ssmfe_expert_standard_double(&rci, nep, n, lambda, m, &(*rr)[0][0][0],
         ind, &keep, &options, &inform);
      switch ( rci.job ) {
      case 1:
         apply_laplacian(
            ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
            );
         break;
      case 2:
         apply_gauss_seidel_step (
            ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
            );
         break;
      case 5:
         if ( rci.i < 0 ) continue;
         for(int k=0; k<rci.nx; k++) {
           int j = ncon + k;
           cblas_dcopy( n, &(*W)[0][rci.jx+k][0], 1, &(*X)[j][0], 1 );
         }
         ncon += rci.nx;
         break;
      case 11:
         if ( rci.i == 0 ) {
            if ( rci.kx != rci.ky || rci.jx > rci.jy ) {
               cblas_dcopy(n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1);
            } else if ( rci.jx < rci.jy ) {
               for(int j=rci.nx-1; j>=0; j--)
                  cblas_dcopy(n, &(*W)[rci.kx][rci.jx+j][0], 1, &(*W)[rci.ky][rci.jy+j][0], 1);
            }
         } else {
            for(int i=0; i<n; i++) {
               for(int j=0; j<rci.nx; j++)
                  (*U)[j][i] = (*W)[rci.kx][ind[j]][i];
               for(int j=0; j<rci.nx; j++)
                  (*W)[rci.kx][j][i] = (*U)[j][i];
               if(rci.ky != rci.kx) {
                  for(int j=0; j<rci.nx; j++)
                    (*U)[j][i] = (*W)[rci.ky][ind[j]][i];
                  for(int j=0; j<rci.nx; j++)
                    (*W)[rci.ky][j][i] = (*U)[j][i];
               }
            }
         }
         break;
      case 12:
         for(int i=0; i<rci.nx; i++)
           (*rr)[rci.k][rci.j+i][rci.i+i] =
             cblas_ddot(n, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1);
         break;
      case 13:
         for(int i=0; i<rci.nx; i++) {
            if( rci.kx == rci.ky ) {
               double s = cblas_dnrm2(n, &(*W)[rci.kx][rci.jx+i][0], 1);
               if( s > 0 )
                  cblas_dscal(n, 1/s, &(*W)[rci.kx][rci.jx+i][0], 1);
            } else {
               double s = sqrt(fabs(cblas_ddot(
                  n, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1)
                  ));
               if ( s > 0 ) {
                  cblas_dscal(n, 1/s, &(*W)[rci.kx][rci.jx+i][0], 1);
                  cblas_dscal(n, 1/s, &(*W)[rci.ky][rci.jy+i][0], 1);
               } else {
                  for(int j=0; j<n; j++)
                     (*W)[rci.ky][rci.jy+i][j] = 0.0;
               }
            }
         }
         break;
      case 14:
         for(int i=0; i<rci.nx; i++) {
           double s = -(*rr)[rci.k][rci.j+i][rci.i+i];
           cblas_daxpy(n, s, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1);
         }
         break;
      case 15:
         if ( rci.nx > 0 && rci.ny > 0 )
            cblas_dgemm(
               CblasColMajor, CblasTrans, CblasNoTrans, rci.nx, rci.ny, n,
               rci.alpha, &(*W)[rci.kx][rci.jx][0], n, &(*W)[rci.ky][rci.jy][0], n,
               rci.beta, &(*rr)[rci.k][rci.j][rci.i], 2*m
               );
         break;
      case 16: // Fall through to 17
      case 17:
         if( rci.ny < 1 ) continue;
         if( rci.nx < 1 ) {
            if( rci.job == 17 ) continue;
            if( rci.beta == 1.0 ) continue;
            for(int j=rci.jy; j<rci.jy+rci.ny; j++)
               cblas_dscal(n, rci.beta, &(*W)[rci.ky][j][0], 1);
            continue;
         }
         if( rci.job == 17 ) {
            cblas_dgemm(
               CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.ny, rci.nx,
               1.0, &(*W)[rci.kx][rci.jx][0], n, &(*rr)[rci.k][rci.j][rci.i], 2*m,
               0.0, &(*W)[rci.ky][rci.jy][0], n
               );
            cblas_dcopy(n*rci.ny, &(*W)[rci.ky][rci.jy][0], 1, &(*W)[rci.kx][rci.jx][0], 1);
         } else {
            cblas_dgemm(
               CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.ny, rci.nx,
               rci.alpha, &(*W)[rci.kx][rci.jx][0], n, &(*rr)[rci.k][rci.j][rci.i],
               2*m, rci.beta, &(*W)[rci.ky][rci.jy][0], n
               );
         }
         break;
      case 21: // Fall through to 22
      case 22:
         if( ncon > 0 ) {
            cblas_dgemm(
               CblasColMajor, CblasTrans, CblasNoTrans, ncon, rci.nx, n,
               1.0, &(*X)[0][0], n, &(*W)[rci.ky][rci.jy][0], n, 0.0, &(*U)[0][0], n
               );
            cblas_dgemm(
               CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, ncon,
               -1.0, &(*X)[0][0], n, &(*U)[0][0], n, 1.0, &(*W)[rci.kx][rci.jx][0], n
               );
         }
         break;
      case 999:
         if( rci.k == 0 ) {
            if( rci.jx > 1 ) {
               for(int j=0; j<rci.jx; j++)
               for(int i=0; i<n; i++)
                  (*W)[0][j][i] = spral_random_real(&state, true);
            }
         }
         break;
      default:
         goto finished;
      }
   }
finished:
   if(inform.flag != 0) printf("inform.flag = %d\n", inform.flag);
   printf("%3d eigenpairs converged in %d iterations\n", ncon, inform.iteration);
   for(int i=0; i<ncon; i++)
      printf(" lambda[%1d] = %13.7e\n", i, lambda[i]);
   spral_ssmfe_expert_free(&keep, &inform);
   free(ind);
   free(lambda);
   free(X);
   free(rr);
   free(W);
   free(U);

   /* Success */
   return 0;
}
