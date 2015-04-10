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

struct spral_ssmfe_options {
   int array_base;
   int print_level;
   int unit_error;
   int unit_warning;
   int unit_diagnostic;
   int max_iterations;
   int user_x;
   int err_est;
   double abs_tol_lambda;
   double rel_tol_lambda;
   double abs_tol_residual;
   double rel_tol_residual;
   double tol_x;
   double left_gap;
   double right_gap;
   int extra_left;
   int extra_right;
   int max_left;
   int max_right;
   bool minAprod;
   bool minBprod;
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

/* Initialize options to defaults */
void spral_ssmfe_default_options(struct spral_ssmfe_options *options);

/************************************
 * SSMFE_EXPERT subroutines (additional to those shared with SSMFE)
 ************************************/

/* Find leftmost eigenpairs of std eigenvalue problem, optional precond */
void spral_ssmfe_standard_double(struct spral_ssmfe_rcid *rci, int left,
      int mep, double *lambda, int m, double *rr, int *ind,
      void **keep, const struct spral_ssmfe_options *optoins,
      struct spral_ssmfe_inform *inform);
/* Free memory */
void spral_ssmfe_expert_free(void **keep, struct spral_ssmfe_inform *inform);

/************************************
 * SSMFE_CORE subroutines 
 ************************************/

/* Initialize options to defaults */
void spral_ssmfe_core_default_options(struct spral_ssmfe_core_options *options);
/* Core driver routine for left and right values of (real) eigenproblems */
void spral_ssmfe_double(struct spral_ssmfe_rcid *rci, int problem,
      int left, int right, int m, double *lambda, double *rr, int *ind,
      void **keep, const struct spral_ssmfe_core_options *options,
      struct spral_ssmfe_inform *inform);
/* Core driver routine for left and right values of (complex) eigenproblems */
void spral_ssmfe_double_complex(struct spral_ssmfe_rciz *rci, int problem,
      int left, int right, int m, double *lambda, double complex *rr, int *ind,
      void **keep, const struct spral_ssmfe_core_options *options,
      struct spral_ssmfe_inform *inform);
/* Core driver routine for largest values of (real) eigenproblems */
void spral_ssmfe_largest_double(struct spral_ssmfe_rcid *rci, int problem,
      int nep, int m, double *lambda, double *rr, int *ind,
      void **keep, const struct spral_ssmfe_core_options *options,
      struct spral_ssmfe_inform *inform);
/* Core driver routine for largest values of (complex) eigenproblems */
void spral_ssmfe_largest_double_complex(struct spral_ssmfe_rciz *rci,
      int problem, int nep, int m, double *lambda, double complex *rr, int *ind,
      void **keep, const struct spral_ssmfe_core_options *options,
      struct spral_ssmfe_inform *inform);
/* Free memory */
void spral_ssmfe_core_free(void **keep, struct spral_ssmfe_inform *inform);

#endif // SPRAL_SSMFE_H
