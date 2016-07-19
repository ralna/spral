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
#pragma once

namespace spral { namespace ssids { namespace cpu {

int ldlt_tpp_factor(int m, int n, int* perm, double* a, int lda, double* d, double* ld, int ldld, double u, double small, int nleft=0, double *aleft=nullptr, int ldleft=0);
void ldlt_tpp_solve_fwd(int m, int n, double const* l, int ldl, int nrhs, double* x, int ldx);
void ldlt_tpp_solve_diag(int n, double const* d, double* x);
void ldlt_tpp_solve_bwd(int m, int n, double const* l, int ldl, int nrhs, double* x, int ldx);

}}} /* end of namespace spral::ssids::cpu */
