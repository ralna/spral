/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 */
namespace spral { namespace ssids { namespace cpu {

void cholesky_factor(int m, int n, double* a, int lda, int blksz, int *info);
void cholesky_solve_fwd(int m, int n, double const* a, int lda, int nrhs, double* x, int ldx);
void cholesky_solve_bwd(int m, int n, double const* a, int lda, int nrhs, double* x, int ldx);

}}} /* namespaces spral::ssids::cpu */
