/* examples/C/random_matrix.c - Example code for SPRAL_RANDOM_MATRIX package */
#include "spral.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
   int state = SPRAL_RANDOM_INITIAL_SEED;

   int m=4, n=5, nnz=8;
   int ptr[n+1], row[nnz];
   double val[nnz];

   /* Generate matrix */
   printf("Generating a %d x %d non-singular matrix with %d non-zeroes\n",
         m, n, nnz);
   if(spral_random_matrix_generate(&state, 0, m, n, nnz, ptr, row, val,
            SPRAL_RANDOM_MATRIX_NONSINGULAR)) {
      printf("Error return from spral_random_matrix_generate()\n");
      exit(1);
   }

   /* Print matrix using utility routine from SPRAL_MATRIX_UTILS package */
   printf("Generated matrix:\n");
   spral_print_matrix(-1, 0, m, n, ptr, row, val, 0);

   /* Success */
   return 0;
}
