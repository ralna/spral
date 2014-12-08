#ifndef SPRAL_SCALING_H
#define SPRAL_SCALING_H

#include <stdbool.h>

struct spral_scaling_hungarian_options {
   int array_base;
   bool scale_if_singular;
};

struct spral_scaling_hungarian_inform {
   int flag;
   int stat;
};

/* Set default values for hungarian_options */
void spral_scaling_hungarian_default_options(
      struct spral_scaling_hungarian_options *options);
/* Scale a symmetric matrix using Hungarian algorithm */
void spral_scaling_hungarian_sym(int n, const int *ptr, const int *row,
      const double *val, int *match, double *scaling,
      const struct spral_scaling_hungarian_options *options,
      struct spral_scaling_hungarian_inform *inform);

#endif // SPRAL_SCALING_H
