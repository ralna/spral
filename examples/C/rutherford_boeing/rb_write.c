/* examples/C/rutherford_boeing/rb_write.c
 * Example code for SPRAL_RUTHERFORD_BOEING */
#include <stdint.h>
#include <stdio.h>
#include "spral.h"

int main(void) {
   /* Initialize options */
   struct spral_rb_write_options options;
   spral_rb_default_write_options(&options);

   /* Data for symmetric matrix
    * ( 2  1         )
    * ( 1  4  1    8 )
    * (    1  3  2   )
    * (       2      )
    * (    8       2 ) */
   int n = 5;
   int64_t ptr[]   = { 0,        2,             5,      7,7,   8 };
   int row[]    = { 0,   1,   1,   2,   4,   2,   3,   4   };
   double val[] = { 2.0, 1.0, 4.0, 1.0, 8.0, 3.0, 2.0, 2.0 };

   /* Write matrix */
   spral_rb_write("matrix.rb", SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row,
      val, &options, "SPRAL_RUTHERFORD_BOEING test matrix", NULL);
}
