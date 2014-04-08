#ifndef SPRAL_SSIDS_H
#define SPRAL_SSIDS_H

#include <stdbool.h>

struct spral_ssids_options {
};

struct spral_ssids_inform {
};

/* Perform analysis phase */
void spral_ssids_analyse(bool check, int n, int *order, const int *ptr,
      const int *row, const double *val, void **akeep,
      const struct spral_ssids_options *options,
      struct spral_ssids_inform *inform);

#endif // SPRAL_SSIDS_H
