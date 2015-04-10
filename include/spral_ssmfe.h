#ifndef SPRAL_SSMFE_H
#define SPRAL_SSMFE_H

#include <stdbool.h>
#include <complex.h>

/************************************
 * Derived types
 ************************************/

struct spral_ssmfe_rcid {
   int job;
   int nx;
   int jx;
   int kx;
   int ny;
   int jy;
   int ky;
   int i;
   int j;
   int k;
   double alpha;
   double beta;
   double *x;
   double *y;
   char unused[80]; // Allow for future expansion
};

struct spral_ssmfe_rciz {
   int job;
   int nx;
   int jx;
   int kx;
   int ny;
   int jy;
   int ky;
   int i;
   int j;
   int k;
   double complex alpha;
   double complex beta;
   double complex *x;
   double complex *y;
   char unused[80]; // Allow for future expansion
};

struct spral_ssmfe_core_options {
   int array_base; // Not in Fortran type
   double cf_max;
   int err_est;
   int extra_left;
   int extra_right;
   double min_gap;
   bool minAprod;
   bool minBprod;
   char unused[80]; // Allow for future expansion
};

struct spral_ssmfe_inform {
   int flag;
   int stat;
   int non_converged;
   int iteration;
   int left;
   int right;
   int *converged;
   double next_left;
   double next_right;
   double *residual_norms;
   double *err_lambda;
   double *err_X;
   char unused[80]; // Allow for future expansion
};

/************************************
 * SSMFE subroutines 
 ************************************/

/************************************
 * SSMFE_EXPERT (additional) subroutines 
 ************************************/

/************************************
 * SSMFE_CORE (additional) subroutines 
 ************************************/

/* Initialize options to defaults */
void spral_ssmfe_core_default_options(struct spral_ssmfe_core_options *options);
/* Core driver routine for (real) standard and generalized eigenproblems */
void spral_ssmfe_ssmfe_double(struct spral_ssmfe_rcid *rci, int problem,
      int left, int right, int m, double *lambda, double *rr, int *ind,
      void **keep, const struct spral_ssmfe_core_options *options,
      struct spral_ssmfe_inform *inform);
void spral_ssmfe_ssmfe_double_complex(struct spral_ssmfe_rciz *rci, int problem,
      int left, int right, int m, double *lambda, double complex *rr, int *ind,
      void **keep, const struct spral_ssmfe_core_options *options,
      struct spral_ssmfe_inform *inform);
/* Free memory */
int spral_ssmfe_core_free(void **keep, struct spral_ssmfe_inform *inform);

#endif // SPRAL_SSMFE_H
