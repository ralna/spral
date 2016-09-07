/* examples/C/scaling/auction_unsym.c - Example code for SPRAL_SCALING */
#include <stdlib.h>
#include <stdio.h>
#include "spral.h"

int main(void) {
   /* Derived types */
   struct spral_scaling_auction_options options;
   struct spral_scaling_auction_inform inform;

   /* Other variables */
   int match[5];
   double rscaling[5], cscaling[5];

   /* Data for unsymmetric matrix:
    * ( 2  5         )
    * ( 1  4       7 )
    * (    1     2   )
    * (       3      )
    * (    8       2 ) */
   int m = 5, n = 5;
   int ptr[]    = { 0,        2,                  6,   7,   8,       10 };
   int row[]    = { 0,   1,   0,   1,   2,   4,   3,   2,   1,   4   };
   double val[] = { 2.0, 1.0, 5.0, 4.0, 1.0, 8.0, 3.0, 2.0, 7.0, 2.0 };
   printf("Initial matrix:\n");
   spral_print_matrix(-1, SPRAL_MATRIX_REAL_UNSYM, m, n, ptr, row, val, 0);

   /* Perform symmetric scaling */
   spral_scaling_auction_default_options(&options);
   spral_scaling_auction_unsym(m, n, ptr, row, val, rscaling, cscaling, match,
         &options, &inform);
   if(inform.flag<0) {
      printf("spral_scaling_auction_unsym() returned with error %5d", inform.flag);
      exit(1);
   }

   /* Print scaling and matching */
   printf("Matching:");
   for(int i=0; i<n; i++) printf(" %10d", match[i]);
   printf("\nRow Scaling: ");
   for(int i=0; i<m; i++) printf(" %10.2le", rscaling[i]);
   printf("\nCol Scaling: ");
   for(int i=0; i<n; i++) printf(" %10.2le", cscaling[i]);
   printf("\n");

   /* Calculate scaled matrix and print it */
   for(int i=0; i<n; i++) {
      for(int j=ptr[i]; j<ptr[i+1]; j++)
         val[j] = rscaling[row[j]] * val[j] * cscaling[i];
   }
   printf("Scaled matrix:\n");
   spral_print_matrix(-1, SPRAL_MATRIX_REAL_UNSYM, m, n, ptr, row, val, 0);

   /* Success */
   return 0;
}
