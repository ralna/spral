/* examples/C/ssmfe/precond_expert.c */
/* Laplacian on a square grid (using SPRAL_SSMFE_EXPERT routines) */
#include "spral.h"
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>
#include <math.h>

/* Header that implements Laplacian and preconditioners */
#include "laplace2d.h"
#include "ldltf.h"

int test_core(void);
int test_core_d(int, int);
int test_core_z(int, int);
int test_expert(void);
int test_expert_d(int);
int test_expert_z(int);
int test_ssmfe(void);
int test_ssmfe_d(int);
int test_ssmfe_z(int);

void zdcopy( int n, double complex* x, double* y );
void dzcopy( int n, double* y, double complex* x );
void vector_operations_d
  ( struct spral_ssmfe_rcid rci,
    int n, int m, int kw, int *ind,
    double (*W)[kw][m][n],
    double (*rr)[3][2*m][2*m],
    double *U);

void vector_operations_z
  ( struct spral_ssmfe_rciz rci,
    int n, int m, int kw, int *ind,
    double complex (*W)[kw][m][n],
    double complex (*rr)[3][2*m][2*m],
    double complex *U);

int main(void) {

  int errors = 0;
  int err;

  fprintf(stdout, "testing ssmfe_core...\n");
  err = test_core();
  errors += err;
  fprintf(stdout, "%d errors\n", err);

  fprintf(stdout, "testing ssmfe_expert...\n");
  err = test_expert();
  errors += err;
  fprintf(stdout, "%d errors\n", err);

  fprintf(stdout, "testing ssmfe...\n");
  err = test_ssmfe();
  errors += err;
  fprintf(stdout, "%d errors\n", err);

  fprintf(stdout, "=============================\n");
  fprintf(stdout, "Total number of errors = %d\n", errors);

  return errors;
}

int test_core(void) {

  int errors = 0;
  int err;

  fprintf(stdout, "testing standard_double...");
  err = test_core_d(0, 0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_double...");
  err = test_core_d(1, 0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing largest_double...");
  err = test_core_d(0, 1);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing standard_double_complex...");
  err = test_core_z(0, 0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_double_complex...");
  err = test_core_z(1, 0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing largest_double_complex...");
  err = test_core_z(0, 1);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  return errors;

}

int test_core_d(int problem, int largest) {
   const int ngrid = 20;          /* grid points along each side */
   const int n     = ngrid*ngrid; /* problem size */
   const int nep   = 5;           /* eigenpairs wanted */
   const int m     = 3;           /* dimension of the iterated subspace */

   const int ncon_x = 5;
   const double lambda_x[] = { /* expected eigenvalues */
      4.4676695e-02,
      1.1119274e-01,
      1.1119274e-01,
      1.7770878e-01,
      2.2040061e-01
   };
   const double lambda_lx[] = { /* expected eigenvalues in 'largest' mode */
      -7.9553233,
      -7.8888073,
      -7.8888073,
      -7.8222912,
      -7.7795994
   };
   int errors;

   int state = SPRAL_RANDOM_INITIAL_SEED; /* PRNG state */

   int *ind = malloc(m * sizeof(*ind));          /* permutation index */
   double *lambda = malloc(n * sizeof(*lambda)); /* eigenvalues storage */
   double (*X)[n][n] = malloc(sizeof(*X));       /* eigenvectors storage */
   /* Work arrays */
   double *lmd = malloc(m * sizeof(*lmd));
   double (*rr)[3][2*m][2*m] = malloc(sizeof(*rr));
   double (*W)[8][m][n] = malloc(sizeof(*W));
   double (*U)[m][n] = malloc(sizeof(*U));
   double *V = malloc(m * (sizeof(*V)));

   /* Derived types */
   struct spral_ssmfe_rcid rci;            /* reverse communication data */
   struct spral_ssmfe_core_options options;/* options */
   void *keep;                             /* private data */
   struct spral_ssmfe_inform inform;       /* information */

   /* Initialize options to default values */
   spral_ssmfe_core_default_options(&options);

   /* Initialize W to lin indep vectors by randomizing */
   for(int i=0; i<n; i++)
   for(int j=0; j<m; j++)
      (*W)[0][j][i] = spral_random_real(&state, true);

   int ncon = 0;   /* number of converged eigenpairs */
   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      if ( largest != 0 )
        spral_ssmfe_largest_double(&rci, problem, nep, m, lmd,
           &(*rr)[0][0][0], ind, &keep, &options, &inform);
      else
        spral_ssmfe_double(&rci, problem, nep, 0, m, lmd,
           &(*rr)[0][0][0], ind, &keep, &options, &inform);
      switch ( rci.job ) {
      case 1:
         apply_laplacian(
            ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
            );
         if ( largest != 0 )
            cblas_dscal( n*rci.nx, -1.0, &(*W)[rci.ky][rci.jy][0], 1 );
         break;
      case 2:
         if ( largest != 0 )
            cblas_dcopy(
               n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1
               );
         else
            apply_gauss_seidel_step (
               ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
               );
         break;
      case 3:
         cblas_dcopy(
            n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1
            );
         break;
      case 4:
         for ( int i = 0; i < m; i++ ) {
            if ( inform.converged[i] != 0 )
               continue;
            if ( inform.err_X[i] > 0 && inform.err_X[i] < 1e-6 )
               inform.converged[i] = 1;
         }
         break;
      case 5:
         if ( rci.i < 0 ) continue;
         for(int k=0; k<rci.nx; k++) {
           int j = ncon + k;
           lambda[j] = lmd[rci.jx + k];
           cblas_dcopy( n, &(*W)[0][rci.jx + k][0], 1, &(*X)[j][0], 1 );
         }
         ncon += rci.nx;
         if ( ncon >= nep || inform.iteration > 300 )
            goto finished;
         break;
      case 11:
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
      case 17:
        vector_operations_d( rci, n, m, 8, ind, W, rr, V );
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
            if ( problem != 0 )
              cblas_dgemm(
                 CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, ncon,
                 -1.0, &(*X)[0][0], n, &(*U)[0][0], n, 1.0, &(*W)[rci.ky][rci.jy][0], n
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
   errors = 0;
   if ( inform.flag != 0 ) errors++;
   if ( ncon != ncon_x ) errors++;
   if ( largest != 0 ) {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_lx[i]) > 1e-6 ) errors++;
   }
   else {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_x[i]) > 1e-6 ) errors++;
   }
   spral_ssmfe_core_free(&keep, &inform);
   free(ind);
   free(lambda);
   free(X);
   free(lmd);
   free(rr);
   free(W);
   free(U);
   free(V);

   return errors;
}

int test_core_z(int problem, int largest) {
   const int ngrid = 20;          /* grid points along each side */
   const int n     = ngrid*ngrid; /* problem size */
   const int nep   = 5;           /* eigenpairs wanted */
   const int m     = 3;           /* dimension of the iterated subspace */

   const int ncon_x = 5;       /* expected number of converged eigenpairs */
   const double lambda_x[] = { /* expected eigenvalues */
      4.4676695e-02,
      1.1119274e-01,
      1.1119274e-01,
      1.7770878e-01,
      2.2040061e-01
   };
   const double lambda_lx[] = { /* expected eigenvalues in 'largest' mode */
      -7.9553233,
      -7.8888073,
      -7.8888073,
      -7.8222912,
      -7.7795994
   };

   const double complex ZERO = 0.0;
   const double complex ONE = 1.0;
   const double complex MINUS_ONE = -1.0;

   int errors;

   int state = SPRAL_RANDOM_INITIAL_SEED; /* PRNG state */

   int *ind = malloc(m * sizeof(*ind));            /* permutation index */
   double *lambda = malloc(n * sizeof(*lambda));   /* eigenvalues storage */
   double complex (*X)[n][n] = malloc(sizeof(*X)); /* eigenvectors storage */
   /* Work arrays */
   double *lmd = malloc(m * sizeof(*lmd));
   double complex (*rr)[3][2*m][2*m] = malloc(sizeof(*rr));
   double complex (*W)[8][m][n] = malloc(sizeof(*W));
   double complex (*U)[m][n] = malloc(sizeof(*U));
   double complex *V = malloc(m * sizeof(*V));

   /* Derived types */
   struct spral_ssmfe_rciz rci;            /* reverse communication data */
   struct spral_ssmfe_core_options options;/* options */
   void *keep;                             /* private data */
   struct spral_ssmfe_inform inform;       /* information */

   /* Initialize options to default values */
   spral_ssmfe_core_default_options(&options);

   /* Initialize W to lin indep vectors by randomizing */
   for(int i=0; i<n; i++)
   for(int j=0; j<m; j++)
      (*W)[0][j][i] = spral_random_real(&state, true);

   int ncon = 0;   /* number of converged eigenpairs */
   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      if ( largest != 0 )
        spral_ssmfe_largest_double_complex(&rci, problem, nep, m, lmd,
           &(*rr)[0][0][0], ind, &keep, &options, &inform);
      else
        spral_ssmfe_double_complex(&rci, problem, nep, 0, m, lmd,
           &(*rr)[0][0][0], ind, &keep, &options, &inform);
      switch ( rci.job ) {
      case 1:
         apply_laplacian_z(
            ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
            );
         if ( largest != 0 )
            cblas_zscal( n*rci.nx, &MINUS_ONE, &(*W)[rci.ky][rci.jy][0], 1 );
         break;
      case 2:
         if ( largest != 0 )
            cblas_zcopy(
               n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1
               );
         else
            apply_gauss_seidel_step_z (
               ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
               );
         break;
      case 3:
         cblas_zcopy(
            n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1
            );
         break;
      case 4:
         for ( int i = 0; i < m; i++ ) {
            if ( inform.converged[i] != 0 )
               continue;
            if ( inform.err_X[i] > 0 && inform.err_X[i] < 1e-6 )
               inform.converged[i] = 1;
         }
         break;
      case 5:
         if ( rci.i < 0 ) continue;
         for(int k=0; k<rci.nx; k++) {
           int j = ncon + k;
           lambda[j] = lmd[rci.jx + k];
           cblas_zcopy( n, &(*W)[0][rci.jx + k][0], 1, &(*X)[j][0], 1 );
         }
         ncon += rci.nx;
         if ( ncon >= nep || inform.iteration > 300 )
            goto finished;
         break;
      case 11:
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
      case 17:
        vector_operations_z( rci, n, m, 8, ind, W, rr, V );
         break;
      case 21: // Fall through to 22
      case 22:
         if( ncon > 0 ) {
            cblas_zgemm(
               CblasColMajor, CblasTrans, CblasNoTrans, ncon, rci.nx, n,
               &ONE, &(*X)[0][0], n, &(*W)[rci.ky][rci.jy][0], n, &ZERO, &(*U)[0][0], n
               );
            cblas_zgemm(
               CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, ncon,
               &MINUS_ONE, &(*X)[0][0], n, &(*U)[0][0], n, &ONE,
               &(*W)[rci.kx][rci.jx][0], n
               );
            if ( problem != 0 )
              cblas_zgemm(
                 CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, ncon,
                 &MINUS_ONE, &(*X)[0][0], n, &(*U)[0][0], n, &ONE,
                 &(*W)[rci.ky][rci.jy][0], n
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
   errors = 0;
   if ( inform.flag != 0 ) errors++;
   if ( ncon != ncon_x ) errors++;
   if ( largest != 0 ) {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_lx[i]) > 1e-6 ) errors++;
   }
   else {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_x[i]) > 1e-6 ) errors++;
   }
   spral_ssmfe_core_free(&keep, &inform);
   free(ind);
   free(lambda);
   free(X);
   free(lmd);
   free(rr);
   free(W);
   free(U);
   free(V);

   return errors;
}

int test_expert(void) {

  int errors = 0;
  int err;

  fprintf(stdout, "testing standard_double...");
  err = test_expert_d(0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_double...");
  err = test_expert_d(1);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing standard_shift_double...");
  err = test_expert_d(2);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_shift_double...");
  err = test_expert_d(3);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing buckling_double...");
  err = test_expert_d(4);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing standard_double_complex...");
  err = test_expert_z(0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_double_complex...");
  err = test_expert_z(1);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing standard_shift_double_complex...");
  err = test_expert_z(2);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_shift_double_complex...");
  err = test_expert_z(3);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing buckling_double_complex...");
  err = test_expert_z(4);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  return errors;
}

int test_expert_d(int problem) {

   const int ngrid = 20;          /* grid points along each side */
   const int n     = ngrid*ngrid; /* problem size */
   const int m     = 3;           /* dimension of the iterated subspace */

   const int ncon_x = 6;
   const double lambda_x[] = {
      4.4676695e-02,
      1.1119274e-01,
      1.1119274e-01,
      1.7770878e-01,
      2.2040061e-01,
      2.2040061e-01
   };
   const double lambda_six[] = {
      9.5108266e-01,
      9.5108266e-01,
      8.8141871e-01,
      8.8141871e-01,
      8.4187478e-01,
      8.4187478e-01
   };

   int errors;
   int nep = 5;

   int state = SPRAL_RANDOM_INITIAL_SEED; /* PRNG state */

   int *ind = malloc(m * sizeof(*ind));          /* permutation index */
   int *ipiv = malloc(n * sizeof(*ipiv));        /* LDLT pivot index */
   double sigma;                                 /* shift */
   double (*A)[n][n] = malloc(sizeof(*A));       /* matrix */
   double (*LDLT)[n][n] = malloc(sizeof(*LDLT)); /* factors */
   double *lambda = malloc(n * sizeof(*lambda)); /* eigenvalues */
   double (*X)[n][n] = malloc(sizeof(*X));       /* eigenvectors */
   /* Work arrays */
   double (*work)[n][n] = malloc(sizeof(*work)); /* work array for dsytrf */
   double (*rr)[3][2*m][2*m] = malloc(sizeof(*rr));
   double (*W)[8][m][n] = malloc(sizeof(*W));
   double (*U)[m][n] = malloc(sizeof(*U));
   double *V = malloc(m * sizeof(*V));

   /* Derived types */
   struct spral_ssmfe_rcid rci;            /* reverse communication data */
   struct spral_ssmfe_options options;     /* options */
   void *keep;                             /* private data */
   struct spral_ssmfe_inform inform;       /* information */

   if ( problem > 1 ) {
     sigma = 1.0;
     set_laplacian_matrix(ngrid, ngrid, n, A);
     for(int j=0; j<n; j++)
     for(int i=0; i<n; i++)
        (*LDLT)[j][i] = (i==j) ? (*A)[j][i] - sigma
                               : (*A)[j][i];
     cwrap_dsytrf('L', n, &(*LDLT)[0][0], n, ipiv, &(*work)[0][0], n*n);
     int i = num_neg_D(n, n, LDLT, ipiv);   /* all evalues to left of sigma */
     if ( i < nep ) nep = i;
   }

   /* Initialize options to default values */
   spral_ssmfe_default_options(&options);
   /* gap between the last converged eigenvalue and the rest of the spectrum
    * must be at least 0.1 times average gap between computed eigenvalues */
   options.left_gap = -0.1;
   /* block size is small, convergence may be slow, allow more iterations */
   options.max_iterations = 300;

   /* Initialize W to lin indep vectors by randomizing */
   for(int i=0; i<n; i++)
   for(int j=0; j<m; j++)
      (*W)[0][j][i] = spral_random_real(&state, true);

   int ncon = 0;   /* number of converged eigenpairs */
   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      switch ( problem ) {
      case 0:
        spral_ssmfe_expert_standard_double(&rci, nep, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 1:
        spral_ssmfe_expert_generalized_double(&rci, nep, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 2:
        spral_ssmfe_expert_standard_shift_double(&rci, sigma, nep, 0, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 3:
        spral_ssmfe_expert_generalized_shift_double(
           &rci, sigma, nep, 0, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 4:
        spral_ssmfe_expert_buckling_double(
           &rci, sigma, nep, 0, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      }
      switch ( rci.job ) {
      case 1:
         if ( problem < 4 )
           apply_laplacian(
              ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
              );
         else
           cblas_dcopy
              ( n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1 );
         break;
      case 2:
         if ( problem < 2 )
           apply_gauss_seidel_step (
              ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
              );
         else
            cblas_dcopy
               ( n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1 );
         break;
      case 3:
         if ( problem < 4 )
           cblas_dcopy
              ( n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1 );
         else
           apply_laplacian(
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
      case 9:
         cblas_dcopy(
            n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1
            );
         cwrap_dsytrs(
            'L', n, rci.nx, &(*LDLT)[0][0], n, ipiv, &(*W)[rci.ky][rci.jy][0], n
            );
         break;
      case 11:
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
      case 17:
        vector_operations_d( rci, n, m, 8, ind, W, rr, V );
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
            if ( problem == 1 || problem == 3 )
              cblas_dgemm(
                 CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, ncon,
                 -1.0, &(*X)[0][0], n, &(*U)[0][0], n, 1.0, &(*W)[rci.ky][rci.jy][0], n
                 );
            else if ( problem == 4 )
               apply_laplacian(
                  ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0],
                  &(*W)[rci.ky][rci.jy][0]
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
   errors = 0;
   if ( inform.flag != 0 ) errors++;
   if ( ncon != ncon_x ) errors++;
   if ( problem < 2 ) {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_x[i]) > 1e-6 ) errors++;
   }
   else {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_six[i]) > 1e-6 ) errors++;
   }
   spral_ssmfe_expert_free(&keep, &inform);
   free(ind);
   free(ipiv);
   free(A);
   free(LDLT);
   free(lambda);
   free(X);
   free(work);
   free(rr);
   free(W);
   free(U);
   free(V);

   return errors;
}

int test_expert_z(int problem) {

   const int ngrid = 20;          /* grid points along each side */
   const int n     = ngrid*ngrid; /* problem size */
   const int m     = 3;           /* dimension of the iterated subspace */

   const int ncon_x = 6;
   const double lambda_x[] = {
      4.4676695e-02,
      1.1119274e-01,
      1.1119274e-01,
      1.7770878e-01,
      2.2040061e-01,
      2.2040061e-01
   };
   const double lambda_six[] = {
      9.5108266e-01,
      9.5108266e-01,
      8.8141871e-01,
      8.8141871e-01,
      8.4187478e-01,
      8.4187478e-01
   };

   const double complex ZERO = 0.0;
   const double complex ONE = 1.0;
   const double complex NONE = -1.0;

   int errors;
   int nep = 5;

   int state = SPRAL_RANDOM_INITIAL_SEED; /* PRNG state */

   int *ind = malloc(m * sizeof(*ind));                    /* permutation index */
   int *ipiv = malloc(n * sizeof(*ipiv));                   /* LDLT pivot index */
   double sigma;                  /* shift */
   double *lambda = malloc(n * sizeof(*lambda));              /* eigenvalues */
   double complex (*X)[n][n] = malloc(sizeof(*X));        /* eigenvectors */
   /* Work arrays */
   double (*B)[n][n] = malloc(sizeof(*B));
   double (*LDLT)[n][n] = malloc(sizeof(*LDLT));
   double (*work)[n][n] = malloc(sizeof(*work));
   double complex (*rr)[3][2*m][2*m] = malloc(sizeof(*rr));
   double complex (*W)[8][m][n] = malloc(sizeof(*W));
   double complex (*U)[m][n] = malloc(sizeof(*U));
   double complex *V = malloc(m * sizeof(*V));

   /* Derived types */
   struct spral_ssmfe_rciz rci;            /* reverse communication data */
   struct spral_ssmfe_options options;     /* options */
   void *keep;                             /* private data */
   struct spral_ssmfe_inform inform;       /* information */

   if ( problem > 1 ) {
     sigma = 1.0;
//     set_laplacian_matrix_z(ngrid, ngrid, n, A);
     set_laplacian_matrix(ngrid, ngrid, n, B);
     for(int j=0; j<n; j++)
     for(int i=0; i<n; i++)
        (*LDLT)[j][i] = (i==j) ? (*B)[j][i] - sigma
                               : (*B)[j][i];
//        (*LDLT)[j][i] = (i==j) ? (*A)[j][i] - sigma
//                               : (*A)[j][i];
     cwrap_dsytrf('L', n, &(*LDLT)[0][0], n, ipiv, &(*work)[0][0], n*n);
   }

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
      switch ( problem ) {
      case 0:
        spral_ssmfe_expert_standard_double_complex(&rci, nep, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 1:
        spral_ssmfe_expert_generalized_double_complex(&rci, nep, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 2:
        spral_ssmfe_expert_standard_shift_double_complex(
           &rci, sigma, nep, 0, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 3:
        spral_ssmfe_expert_generalized_shift_double_complex(
           &rci, sigma, nep, 0, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      case 4:
        spral_ssmfe_expert_buckling_double_complex(
           &rci, sigma, nep, 0, n, lambda,
           m, &(*rr)[0][0][0], ind, &keep, &options, &inform);
        break;
      }
      switch ( rci.job ) {
      case 1:
         if ( problem < 4 )
           apply_laplacian_z(
              ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
              );
         else
            cblas_zcopy
               ( n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1 );
         break;
      case 2:
         if ( problem < 2 )
            apply_gauss_seidel_step_z (
               ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
               );
         else
            cblas_zcopy
               ( n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1 );
         break;
      case 3:
         if ( problem < 4 )
           cblas_zcopy
              ( n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1 );
         else
           apply_laplacian_z(
              ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0], &(*W)[rci.ky][rci.jy][0]
              );
         break;
      case 5:
         if ( rci.i < 0 ) continue;
         for(int k=0; k<rci.nx; k++) {
           int j = ncon + k;
           cblas_zcopy( n, &(*W)[0][rci.jx+k][0], 1, &(*X)[j][0], 1 );
         }
         ncon += rci.nx;
         break;
      case 9:
         for(int i=0; i<n; i++)
           for(int j=0; j<rci.nx; j++)
              (*work)[j][i] = (*W)[rci.kx][rci.jx + j][i];
         cwrap_dsytrs(
            'L', n, rci.nx, &(*LDLT)[0][0], n, ipiv, &(*work)[0][0], n
            );
         for(int i=0; i<n; i++)
           for(int j=0; j<rci.nx; j++)
              (*W)[rci.ky][rci.jy + j][i] = (*work)[j][i];
         break;
      case 11:
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
      case 17:
        vector_operations_z( rci, n, m, 8, ind, &(*W), &(*rr), V );
         break;
      case 21: // Fall through to 22
      case 22:
         if( ncon > 0 ) {
            cblas_zgemm(
               CblasColMajor, CblasTrans, CblasNoTrans, ncon, rci.nx, n,
               &ONE, &(*X)[0][0], n, &(*W)[rci.ky][rci.jy][0], n, &ZERO, &(*U)[0][0], n
               );
            cblas_zgemm(
               CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, ncon,
               &NONE, &(*X)[0][0], n, &(*U)[0][0], n, &ONE, &(*W)[rci.kx][rci.jx][0], n
               );
            if ( problem == 1 || problem == 3 )
              cblas_zgemm(
                 CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.nx, ncon,
                 &NONE, &(*X)[0][0], n, &(*U)[0][0], n, &ONE, &(*W)[rci.ky][rci.jy][0], n
                 );
            else if ( problem == 4 )
               apply_laplacian_z(
                  ngrid, ngrid, rci.nx, &(*W)[rci.kx][rci.jx][0],
                  &(*W)[rci.ky][rci.jy][0]
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
   errors = 0;
   if ( inform.flag != 0 ) errors++;
   if ( ncon != ncon_x ) errors++;
   if ( problem < 2 ) {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_x[i]) > 1e-6 ) errors++;
   }
   else {
     for ( int i = 0; i < ncon && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_six[i]) > 1e-6 ) errors++;
   }
   spral_ssmfe_expert_free(&keep, &inform);
   free(ind);
   free(ipiv);
   free(lambda);
   free(X);
   free(B);
   free(LDLT);
   free(work);
   free(rr);
   free(W);
   free(U);
   free(V);

   return errors;
}

int test_ssmfe(void) {

  int errors = 0;
  int err;

  fprintf(stdout, "testing standard_double...");
  err = test_ssmfe_d(0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_double...");
  err = test_ssmfe_d(1);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing standard_shift_double...");
  err = test_ssmfe_d(2);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_shift_double...");
  err = test_ssmfe_d(3);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing buckling_double...");
  err = test_ssmfe_d(4);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing standard_double_complex...");
  err = test_ssmfe_z(0);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_double_complex...");
  err = test_ssmfe_z(1);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing standard_shift_double_complex...");
  err = test_ssmfe_z(2);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing generalized_shift_double_complex...");
  err = test_ssmfe_z(3);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  fprintf(stdout, "testing buckling_double_complex...");
  err = test_ssmfe_z(4);
  if ( err != 0 )
    fprintf(stdout, "%d errors\n", err);
  else
    fprintf(stdout, "ok\n");
  errors += err;

  return errors;
}

int test_ssmfe_d(int problem) {

   const int m   = 20;     /* grid points along each side */
   const int n   = m*m;    /* problem size */

   const int ncon_x = 6;
   const double sigma = 1.0;  /* shift */
   const double lambda_x[] = {
      4.4676695e-02,
      1.1119274e-01,
      1.1119274e-01,
      1.7770878e-01,
      2.2040061e-01,
      2.2040061e-01
   };
   const double lambda_six[] = {
      8.4187478e-01,
      8.4187478e-01,
      8.8141871e-01,
      8.8141871e-01,
      9.5108266e-01,
      9.5108266e-01
   };

   int errors;
   int nep = 5;

   int *ipiv = malloc(n * sizeof(*ipiv));        /* LDLT pivot index */

   double *lambda = malloc(n * sizeof(*lambda)); /* eigenvalues */
   double (*X)[n][n] = malloc(sizeof(*X));       /* eigenvectors */
   double (*A)[n][n] = malloc(sizeof(*A));       /* matrix */
   double (*LDLT)[n][n] = malloc(sizeof(*LDLT)); /* factors */
   double (*work)[n][n] = malloc(sizeof(*work)); /* work array for dsytrf */

   struct spral_ssmfe_rcid rci;           /* reverse communication data */
   struct spral_ssmfe_options options;    /* options */
   void *keep;                            /* private data */
   struct spral_ssmfe_inform inform;      /* information */

   if ( problem > 1 ) {
     set_laplacian_matrix( m, m, n, A );
     for(int j=0; j<n; j++)
     for(int i=0; i<n; i++)
        (*LDLT)[j][i] = (i==j) ? (*A)[j][i] - sigma
                               : (*A)[j][i];
     cwrap_dsytrf( 'L', n, &(*LDLT)[0][0], n, ipiv, &(*work)[0][0], n*n );
   }

   /* Initialize options to default values */
   spral_ssmfe_default_options(&options);
   /* gap between the last converged eigenvalue and the rest of the spectrum
    * must be at least 0.1 times average gap between computed eigenvalues */
   options.left_gap = -0.1;

   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      switch ( problem ) {
      case 0:
         spral_ssmfe_standard_double( &rci, nep, 2*nep, lambda, n, &(*X)[0][0], n,
            &keep, &options, &inform );
         break;
      case 1:
         spral_ssmfe_standard_double( &rci, nep, 2*nep, lambda, n, &(*X)[0][0], n,
            &keep, &options, &inform );
         break;
      case 2:
         spral_ssmfe_standard_shift_double( &rci, sigma, nep, 0, n, lambda,
            n, &(*X)[0][0], n, &keep, &options, &inform );
         break;
      case 3:
         spral_ssmfe_generalized_shift_double( &rci, sigma, nep, 0, n, lambda,
            n, &(*X)[0][0], n, &keep, &options, &inform );
         break;
      case 4:
         spral_ssmfe_buckling_double( &rci, sigma, nep, 0, n, lambda,
            n, &(*X)[0][0], n, &keep, &options, &inform );
         break;
      }
      switch ( rci.job ) {
      case 1:
         if ( problem < 4 )
            apply_laplacian(m, m, rci.nx, rci.x, rci.y);
         else
            cblas_dcopy( n*rci.nx, rci.x, 1, rci.y, 1 );
         break;
      case 2:
         if ( problem < 2 )
            apply_gauss_seidel_step(m, m, rci.nx, rci.x, rci.y);
         break;
      case 3:
         if ( problem < 4 )
            cblas_dcopy( n*rci.nx, rci.x, 1, rci.y, 1 );
         else
            apply_laplacian(m, m, rci.nx, rci.x, rci.y);
         break;
      case 9:
         cblas_dcopy( n*rci.nx, rci.x, 1, rci.y, 1 );
         cwrap_dsytrs( 'L', n, rci.nx, &(*LDLT)[0][0], n, ipiv, rci.y, n );
         break;
      default:
         goto finished;
      }
   }
finished:
   errors = 0;
   if ( inform.flag != 0 ) errors++;
   if ( inform.left != ncon_x ) errors++;
   if ( problem < 2 ) {
     for ( int i = 0; i < inform.left && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_x[i]) > 1e-6 ) errors++;
   }
   else {
     for ( int i = 0; i < inform.left && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_six[i]) > 1e-6 ) errors++;
   }
   spral_ssmfe_free_double(&keep, &inform);
   free(ipiv);
   free(lambda);
   free(X);
   free(A);
   free(LDLT);
   free(work);
   return errors;
}

int test_ssmfe_z(int problem) {

   const int m   = 20;     /* grid points along each side */
   const int n   = m*m;    /* problem size */

   const int ncon_x = 6;
   const double sigma = 1.0;  /* shift */
   const double lambda_x[] = {
      4.4676695e-02,
      1.1119274e-01,
      1.1119274e-01,
      1.7770878e-01,
      2.2040061e-01,
      2.2040061e-01
   };
   const double lambda_six[] = {
      8.4187478e-01,
      8.4187478e-01,
      8.8141871e-01,
      8.8141871e-01,
      9.5108266e-01,
      9.5108266e-01
   };

   int errors;
   int nep = 5;

   int *ipiv = malloc(n * sizeof(*ipiv));          /* LDLT pivot index */

   double *lambda = malloc(n * sizeof(*lambda));   /* eigenvalues */
   double complex (*X)[n][n] = malloc(sizeof(*X)); /* eigenvectors */
   double (*A)[n][n] = malloc(sizeof(*A));         /* matrix */
   double (*LDLT)[n][n] = malloc(sizeof(*LDLT));   /* factors */
   double (*work)[n][n] = malloc(sizeof(*work));   /* work array for dsytrf */

   struct spral_ssmfe_rciz rci;           /* reverse communication data */
   struct spral_ssmfe_options options;    /* options */
   void *keep;                            /* private data */
   struct spral_ssmfe_inform inform;      /* information */

   if ( problem > 1 ) {
     set_laplacian_matrix( m, m, n, A );
     for(int j=0; j<n; j++)
     for(int i=0; i<n; i++)
        (*LDLT)[j][i] = (i==j) ? (*A)[j][i] - sigma
                               : (*A)[j][i];
     cwrap_dsytrf( 'L', n, &(*LDLT)[0][0], n, ipiv, &(*work)[0][0], n*n );
   }

   /* Initialize options to default values */
   spral_ssmfe_default_options(&options);
   /* gap between the last converged eigenvalue and the rest of the spectrum
    * must be at least 0.1 times average gap between computed eigenvalues */
   options.left_gap = -0.1;

   rci.job = 0; keep = NULL;
   while(true) { /* reverse communication loop */
      switch ( problem ) {
      case 0:
         spral_ssmfe_standard_double_complex(
            &rci, nep, 2*nep, lambda, n, &(*X)[0][0], n,
            &keep, &options, &inform );
         break;
      case 1:
         spral_ssmfe_standard_double_complex(
            &rci, nep, 2*nep, lambda, n, &(*X)[0][0], n,
            &keep, &options, &inform );
         break;
      case 2:
         spral_ssmfe_standard_shift_double_complex(
            &rci, sigma, nep, 0, n, lambda,
            n, &(*X)[0][0], n, &keep, &options, &inform );
         break;
      case 3:
         spral_ssmfe_generalized_shift_double_complex(
            &rci, sigma, nep, 0, n, lambda,
            n, &(*X)[0][0], n, &keep, &options, &inform );
         break;
      case 4:
         spral_ssmfe_buckling_double_complex(
            &rci, sigma, nep, 0, n, lambda,
            n, &(*X)[0][0], n, &keep, &options, &inform );
         break;
      }
      switch ( rci.job ) {
      case 1:
         if ( problem < 4 )
            apply_laplacian_z( m, m, rci.nx, rci.x, rci.y );
         else
            cblas_zcopy( n*rci.nx, rci.x, 1, rci.y, 1 );
         break;
      case 2:
         if ( problem < 2 )
            apply_gauss_seidel_step_z( m, m, rci.nx, rci.x, rci.y );
         break;
      case 3:
         if ( problem < 4 )
            cblas_zcopy( n*rci.nx, rci.x, 1, rci.y, 1 );
         else
            apply_laplacian_z( m, m, rci.nx, rci.x, rci.y );
         break;
      case 9:
         zdcopy( n*rci.nx, rci.x, &(*work)[0][0] );
         cwrap_dsytrs( 'L', n, rci.nx, &(*LDLT)[0][0], n, ipiv, &(*work)[0][0], n );
         dzcopy( n*rci.nx, &(*work)[0][0], rci.y );
         break;
      default:
         goto finished;
      }
   }
finished:
   errors = 0;
   if ( inform.flag != 0 ) errors++;
   if ( inform.left != ncon_x ) errors++;
   if ( problem < 2 ) {
     for ( int i = 0; i < inform.left && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_x[i]) > 1e-6 ) errors++;
   }
   else {
     for ( int i = 0; i < inform.left && i < ncon_x; i++ )
        if ( fabs(lambda[i] - lambda_six[i]) > 1e-6 ) errors++;
   }
   if ( problem == -2 ) {
   printf("%d eigenpairs converged in %d iterations\n", inform.left, inform.iteration);
   for(int i=0; i<inform.left; i++)
      printf(" lambda[%1d] = %13.7e\n", i, lambda[i]);
   }
   spral_ssmfe_free_double_complex(&keep, &inform);
   free(ipiv);
   free(lambda);
   free(X);
   free(A);
   free(LDLT);
   free(work);
   return errors;
}

void zdcopy( int n, double complex* x, double* y ) {
   for ( int i = 0; i < n; i++ )
      y[i] = x[i];
}

void dzcopy( int n, double* x, double complex* y ) {
   for ( int i = 0; i < n; i++ )
      y[i] = x[i];
}

void vector_operations_d
  ( struct spral_ssmfe_rcid rci,
    int n, int m, int kw, int *ind,
    double (*W)[kw][m][n],
    double (*rr)[3][2*m][2*m],
    double *U) {

   switch ( rci.job ) {
      case 11:
         if ( rci.i == 0 ) {
            if ( rci.kx != rci.ky || rci.jx > rci.jy ) {
               cblas_dcopy(
                  n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1
                  );
            } else if ( rci.jx < rci.jy ) {
               for(int j=rci.nx-1; j>=0; j--)
                  cblas_dcopy(
                     n, &(*W)[rci.kx][rci.jx+j][0], 1, &(*W)[rci.ky][rci.jy+j][0], 1
                     );
            }
         } else {
            for(int i=0; i<n; i++) {
               for(int j=0; j<rci.nx; j++)
                  U[j] = (*W)[rci.kx][ind[j]][i];
               for(int j=0; j<rci.nx; j++)
                  (*W)[rci.kx][j][i] = U[j];
               if(rci.ky != rci.kx) {
                  for(int j=0; j<rci.nx; j++)
                    U[j] = (*W)[rci.ky][ind[j]][i];
                  for(int j=0; j<rci.nx; j++)
                    (*W)[rci.ky][j][i] = U[j];
               }
            }
         }
         break;
      case 12:
         for(int i=0; i<rci.nx; i++)
           (*rr)[rci.k][rci.j+i][rci.i+i] =
             cblas_ddot(
                n, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1
                );
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
           cblas_daxpy(
              n, s, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1
              );
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
         if( rci.job == 17 ) {
            cblas_dgemm(
               CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.ny, rci.nx,
               1.0, &(*W)[rci.kx][rci.jx][0], n, &(*rr)[rci.k][rci.j][rci.i], 2*m,
               0.0, &(*W)[rci.ky][rci.jy][0], n
               );
            cblas_dcopy(
               n*rci.ny, &(*W)[rci.ky][rci.jy][0], 1, &(*W)[rci.kx][rci.jx][0], 1
               );
         } else {
            cblas_dgemm(
               CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.ny, rci.nx,
               rci.alpha, &(*W)[rci.kx][rci.jx][0], n, &(*rr)[rci.k][rci.j][rci.i],
               2*m, rci.beta, &(*W)[rci.ky][rci.jy][0], n
               );
         }
         break;
   }
}

void vector_operations_z
  ( struct spral_ssmfe_rciz rci,
    int n, int m, int kw, int *ind,
    double complex (*W)[kw][m][n],
    double complex (*rr)[3][2*m][2*m],
    double complex *U) {

   const double complex ZERO = 0.0;
   const double complex ONE = 1.0;

   double complex z;

    switch ( rci.job ) {
    case 11:
       if ( rci.i == 0 ) {
          if ( rci.kx != rci.ky || rci.jx > rci.jy ) {
             cblas_zcopy(
                n*rci.nx, &(*W)[rci.kx][rci.jx][0], 1, &(*W)[rci.ky][rci.jy][0], 1
                );
          } else if ( rci.jx < rci.jy ) {
             for(int j=rci.nx-1; j>=0; j--)
                cblas_zcopy(
                   n, &(*W)[rci.kx][rci.jx+j][0], 1, &(*W)[rci.ky][rci.jy+j][0], 1
                   );
          }
       } else {
          for(int i=0; i<n; i++) {
             for(int j=0; j<rci.nx; j++)
                U[j] = (*W)[rci.kx][ind[j]][i];
             for(int j=0; j<rci.nx; j++)
                (*W)[rci.kx][j][i] = U[j];
             if(rci.ky != rci.kx) {
                for(int j=0; j<rci.nx; j++)
                  U[j] = (*W)[rci.ky][ind[j]][i];
                for(int j=0; j<rci.nx; j++)
                  (*W)[rci.ky][j][i] = U[j];
             }
          }
       }
       break;
    case 12:
       for(int i=0; i<rci.nx; i++) {
          cblas_zdotc_sub(
              n, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1, &z
              );
          (*rr)[rci.k][rci.j+i][rci.i+i] = z;
       }
       break;
    case 13:
       for(int i=0; i<rci.nx; i++) {
          if( rci.kx == rci.ky ) {
             double s = cblas_dznrm2(n, &(*W)[rci.kx][rci.jx+i][0], 1);
             if( s > 0 ) {
                z = ONE/s;
                cblas_zscal(n, &z, &(*W)[rci.kx][rci.jx+i][0], 1);
             }
          } else {
             cblas_zdotc_sub(
                n, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1, &z
                );
             double s = sqrt(fabs(z));
             if ( s > 0 ) {
                z = ONE/s;
                cblas_zscal(n, &z, &(*W)[rci.kx][rci.jx+i][0], 1);
                cblas_zscal(n, &z, &(*W)[rci.ky][rci.jy+i][0], 1);
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
         z = s;
         cblas_zaxpy(
            n, &z, &(*W)[rci.kx][rci.jx+i][0], 1, &(*W)[rci.ky][rci.jy+i][0], 1
            );
       }
       break;
    case 15:
       if ( rci.nx > 0 && rci.ny > 0 )
          cblas_zgemm(
             CblasColMajor, CblasTrans, CblasNoTrans, rci.nx, rci.ny, n,
             &rci.alpha, &(*W)[rci.kx][rci.jx][0], n, &(*W)[rci.ky][rci.jy][0], n,
             &rci.beta, &(*rr)[rci.k][rci.j][rci.i], 2*m
             );
       break;
    case 16: // Fall through to 17
    case 17:
       if( rci.job == 17 ) {
          cblas_zgemm(
             CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.ny, rci.nx,
             &ONE, &(*W)[rci.kx][rci.jx][0], n, &(*rr)[rci.k][rci.j][rci.i], 2*m,
             &ZERO, &(*W)[rci.ky][rci.jy][0], n
             );
          cblas_zcopy(
             n*rci.ny, &(*W)[rci.ky][rci.jy][0], 1, &(*W)[rci.kx][rci.jx][0], 1
             );
       } else {
          cblas_zgemm(
             CblasColMajor, CblasNoTrans, CblasNoTrans, n, rci.ny, rci.nx,
             &rci.alpha, &(*W)[rci.kx][rci.jx][0], n, &(*rr)[rci.k][rci.j][rci.i],
             2*m, &rci.beta, &(*W)[rci.ky][rci.jy][0], n
             );
       }
       break;
     }
}
