#ifndef SPRAL_MATRIX_UTIL_H
#define SPRAL_MATRIX_UTIL_H

/* Note: At present, interface is only partially defined! */

enum spral_matrix_type {
   SPRAL_MATRIX_UNDEFINED=0,
   SPRAL_MATRIX_REAL_RECT=1,        SPRAL_MATRIX_CPLX_RECT=-1,
   SPRAL_MATRIX_REAL_UNSYM=2,       SPRAL_MATRIX_CPLX_UNSYM=-2,
   SPRAL_MATRIX_REAL_SYM_PSDEF=3,   SPRAL_MATRIX_CPLX_HERM_PSDEF=-3,
   SPRAL_MATRIX_REAL_SYM_INDEF=4,   SPRAL_MATRIX_CPLX_HERM_INDEF=-4,
                                    SPRAL_MATRIX_CPLX_SYM=-5,
   SPRAL_MATRIX_REAL_SKEW=6,        SPRAL_MATRIX_CPLX_SKEW=-6
};

void spral_print_matrix(int lines, enum spral_matrix_type matrix_type, int m,
      int n, int *ptr, int *row, double *val, int base);

#endif // SPRAL_MATRIX_UTIL_H
