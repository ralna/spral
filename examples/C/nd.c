/* examples/Fortran/nd.f90 - Example code for SPRAL_ND package */
#include "spral.h"
#include <stdio.h>

int main(void) {
   /* Derived types */
   struct spral_nd_options options;
   struct spral_nd_inform inform;

   /* Initialize options */
   spral_nd_default_options(&options);

   /* Data for matrix:
    *     0  1  2  3  4  5  6  7
    * 0 (    X     X  X  X       )
    * 1 ( X        X     X  X    )
    * 2 (             X  X       )
    * 3 ( X  X              X    )
    * 4 ( X     X        X     X )
    * 5 ( X  X  X     X     X  X )
    * 6 (    X     X     X       )
    * 7 (             X  X       ) */
   const int n = 8;
   int ptr[] = { 0,          4,       7,    9, 10,  12,   14, 14, 14 };
   int row[] = { 1, 3, 4, 5, 3, 5, 6, 4, 5, 6, 5, 7, 6, 7 };

   /* Find and print nested dissection ordering */
   int perm[n];
   spral_nd_order(0, n, ptr, row, perm, &options, &inform);
   printf("Permutation:");
   for(int i=0; i<n; i++) printf(" %5d", perm[i]);
   printf("\n");

   return 0;
}
