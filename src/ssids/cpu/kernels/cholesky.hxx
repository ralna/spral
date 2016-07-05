namespace spral { namespace ssids { namespace cpu {

void cholesky_factor(int m, int n, double* a, int lda, int blksz, int *info);
void cholesky_solve_fwd(int m, int n, double const* a, int lda, double* x);
void cholesky_solve_bwd(int m, int n, double const* a, int lda, double* x);

}}} /* namespaces spral::ssids::cpu */
