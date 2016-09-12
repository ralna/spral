#ifndef SPRAL_RUTHERFORD_BOEING_H
#define SPRAL_RUTHERFORD_BOEING_H

#ifdef __cplusplus
extern "C" {
#endif

struct spral_rb_read_options {
   int array_base;
   bool add_diagonal;
   float extra_space;
   int format;
   int lwr_upr_full;
   int values;
};

void spral_rb_default_read_options(struct spral_rb_read_options *options);
int spral_rb_peek(const char *filename, int *m, int *n, long *nelt, long *nvar,
      long *nval, enum spral_matrix_type *matrix_type, char *type_code,
      char *title, char *identifier);
int spral_rb_read_i32d(const char *filename, void **handle, int *m, int *n,
      int **ptr, int **row, double **val,
      const struct spral_rb_read_options *options, char *type_code, char *title,
      char *identifier);
int spral_rb_read_i64d(const char *filename, void **handle, int *m, int *n,
      long **ptr, int **row, double **val,
      const struct spral_rb_read_options *options, char *type_code, char *title,
      char *identifier);
void spral_rb_free_handle(void **handle);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SPRAL_RUTHERFORD_BOEING_H
