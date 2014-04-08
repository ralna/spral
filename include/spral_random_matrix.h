#ifndef SPRAL_RANDOM_MATRIX_H
#define SPRAL_RANDOM_MATRIX_H

#include <stdbool.h>

#define SPRAL_RANDOM_MATRIX_FINDEX        1
#define SPRAL_RANDOM_MATRIX_NONSINGULAR   2
#define SPRAL_RANDOM_MATRIX_SORT          4

/* Generate an m x n random matrix with nnz non-zero entries */
int spral_random_matrix_generate(int *state, int matrix_type, int m, int n,
      int nnz, int *ptr, int *row, double *val, int flags);

#endif // SPRAL_RANDOM_MATRIX_H
