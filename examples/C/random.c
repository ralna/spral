/* examples/C/random.c - Example code for SPRAL_RANDOM package */
#include "spral.h"
#include <limits.h>
#include <stdint.h>
#include <stdio.h>

#define bool2str(VAL) ( (VAL) ? "true" : "false" )

int main(void) {
   int state = SPRAL_RANDOM_INITIAL_SEED;

   // Store initial random seed so we can reuse it later
   int initial_seed = state;

   // Generate some random values
   printf("\nSome random values\n");
   printf("Sample Unif(-1,1)               = %16.12f\n",
         spral_random_real(&state, false));
   printf("Sample Unif(0,1)                = %16.12f\n",
         spral_random_real(&state, true));
   printf("Sample Unif(1, ..., 20)         = %16d\n",
         spral_random_integer(&state, 20));
   printf("Sample Unif(1, ..., 20*INT_MAX) = %16ld\n",
         spral_random_long(&state,((int64_t)20)*INT_MAX));
   printf("Sample B(1,0.5)                 = %16s\n",
         bool2str(spral_random_logical(&state)));

   // Restore initial seed
   state = initial_seed;

   // Generate the same random values
   printf("\nThe same random values again\n");
   printf("Sample Unif(-1,1)               = %16.12f\n",
         spral_random_real(&state, false));
   printf("Sample Unif(0,1)                = %16.12f\n",
         spral_random_real(&state, true));
   printf("Sample Unif(1, ..., 20)         = %16d\n",
         spral_random_integer(&state, 20));
   printf("Sample Unif(1, ..., 20*INT_MAX) = %16ld\n",
         spral_random_long(&state, ((int64_t)20)*INT_MAX));
   printf("Sample B(1,0.5)                 = %16s\n",
         bool2str(spral_random_logical(&state)));

   /* Success */
   return 0;
}
