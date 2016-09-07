/* examples/C/scaling/hungarian_sym.c - Example code for SPRAL_SCALING */
#include <stdlib.h>
#include <stdio.h>
#include "spral.h"

int main(void) {
   /* Derived types */
   struct spral_scaling_hungarian_options options;
   struct spral_scaling_hungarian_inform inform;

   /* Other variables */
   int match[5];
   double scaling[5];

   /* Data for symmetric matrix:
    * ( 2  1         )
    * ( 1  4  1    8 )
    * (    1  3  2   )
    * (       2      )
    * (    8       2 ) */
   int n = 5;
   int ptr[]    = { 0,        2,             5,      7,7,   8 };
   int row[]    = { 0,   1,   1,   2,   4,   2,   3,   4   };
   double val[] = { 2.0, 1.0, 4.0, 1.0, 8.0, 3.0, 2.0, 2.0 };
   printf("Initial matrix:\n");
   spral_print_matrix(-1, SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, val, 0);

   /* Perform symmetric scaling */
   spral_scaling_hungarian_default_options(&options);
   spral_scaling_hungarian_sym(n, ptr, row, val, scaling, match, &options, &inform);
   if(inform.flag<0) {
      printf("spral_scaling_hungarian_sym() returned with error %5d", inform.flag);
      exit(1);
   }

   /* Print scaling and matching */
   printf("Matching:");
   for(int i=0; i<n; i++) printf(" %10d", match[i]);
   printf("\nScaling: ");
   for(int i=0; i<n; i++) printf(" %10.2le", scaling[i]);
   printf("\n");

   /* Calculate scaled matrix and print it */
   for(int i=0; i<n; i++) {
      for(int j=ptr[i]; j<ptr[i+1]; j++)
         val[j] = scaling[row[j]] * val[j] * scaling[i];
   }
   printf("Scaled matrix:\n");
   spral_print_matrix(-1, SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, val, 0);

   /* Success */
   return 0;
}
