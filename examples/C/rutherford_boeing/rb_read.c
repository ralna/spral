/* examples/C/rutherford_boeing/rb_read.c
 * Example code for SPRAL_RUTHERFORD_BOEING */
#include <stdint.h>
#include <stdio.h>
#include "spral.h"

int main(void) {
   /* Matrix data */
   void *handle;
   char title[73]; // Maximum length 72 character + '\0'
   int matrix_type, m, n;
   int64_t *ptr;
   int *row;
   double *val;

   /* Initalise options */
   struct spral_rb_read_options options;
   spral_rb_default_read_options(&options);

   /* Read matrix */
   spral_rb_read("matrix.rb", &handle, &matrix_type, &m, &n, &ptr, &row, &val,
         &options, title, NULL, NULL);

   /* Print matrix */
   printf("Matrix '%s'\n", title);
   spral_print_matrix_i64d(-1, matrix_type, m, n, ptr, row, val, 0);

   /* Free handle */
   spral_rb_free_handle(&handle);
}
