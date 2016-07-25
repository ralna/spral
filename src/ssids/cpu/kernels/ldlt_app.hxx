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

template<typename T>
int ldlt_app_factor(int m, int n, int *perm, T *a, int lda, T *d, double u, double small);

template <typename T>
void ldlt_app_solve_fwd(int m, int n, T const* l, int ldl, int nrhs, T* x, int ldx);

template <typename T>
void ldlt_app_solve_diag(int n, T const* d, T* x);

template <typename T>
void ldlt_app_solve_bwd(int m, int n, T const* l, int ldl, int nrhs, T* x, int ldx);

}}} /* namespaces spral::ssids::cpu */
