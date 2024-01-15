/* examples/C/lsmr.c - Example code for SPRAL_LSMR package */
#include "spral.h"

#include <stdlib.h>
#include <stdio.h>

void matrix_mult(int m, int n, int const *ptr, int const *row,
      double const *val, double const *v, double *u);
void matrix_mult_trans(int m, int n, int const *ptr, int const *row,
      double const *val, double const *u, double *v);

int main(void) {
   /* Derived types */
   void *keep;
   struct spral_lsmr_options options;
   struct spral_lsmr_inform inform;

   /* Initialize derived types */
   keep = NULL;
   spral_lsmr_default_options(&options);
   options.unit_diagnostics = -1; /* Switch off diagnostic printing */

   /* Data for matrix:
    * ( 1.0     -1.0 )
    * (     2.0      )
    * ( 2.0      2.0 )
    * ( 5.0 3.0 -2.0 )
    * (          6.0 ) */
   const int m = 5, n = 3;
   int ptr[]    = {   0,             3,         5,                 9 };
   int row[]    = {   0,   2,   3,   1,   3,    0,   2,    3,   4 };
   double val[] = { 1.0, 2.0, 5.0, 2.0, 3.0, -1.0, 2.0, -2.0, 6.0 };
   /* Data for rhs b */
   double b[] = { 1.0, 1.0, 1.0, 1.0, 1.0 };

   /* prepare for LSMR calls (using no preconditioning) */
   int action = 0;
   double u[m], v[n], x[n];
   for(int i=0; i<m; ++i) u[i] = b[i];

   bool done = false;
   while(!done) {
      spral_lsmr_solve(&action, m, n, u, v, x, &keep, &options, &inform, NULL);

      switch(action) {
      case 0: /* we are done */
         printf("Exit LSMR with inform.flag = %d and inform.itn = %d\n",
               inform.flag, inform.itn);
         printf("LS solution is:\n");
         for(int i=0; i<n; ++i) printf(" %10.2f", x[i]);
         printf("\n");
         done = true;
         break;

      case 1: /* Compute v = v + A'*u without altering u */
         matrix_mult_trans(m, n, ptr, row, val, u, v);
         break;

      case 2: /* Compute u = u + A*v  without altering v */
         matrix_mult(m, n, ptr, row, val, v, u);
         break;
      }
   }

   spral_lsmr_free(&keep);
   return 0; // Success
}

/* Takes b and computes u = u + A * v (A in CSC format) */
void matrix_mult(int m, int n, int const *ptr, int const *row,
      double const *val, double const *v, double *u) {
   for(int j=0; j<n; ++j) {
      for(int k=ptr[j]; k<ptr[j+1]; ++k) {
         int i = row[k];
         u[i] += val[k]*v[j];
      }
   }
}

/* Takes b and computes v = v + A^T * u (A in CSC format) */
void matrix_mult_trans(int m, int n, int const *ptr, int const *row,
      double const *val, double const *u, double *v) {
   for(int j=0; j<n; ++j) {
      for(int k=ptr[j]; k<ptr[j+1]; ++k) {
         int i = row[k];
         v[j] += val[k]*u[i];
      }
   }
}
