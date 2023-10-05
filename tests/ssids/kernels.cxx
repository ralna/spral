/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * Licence: BSD licence, see LICENCE file for details
 *
 */

#include <cstdio>

#include <fenv.h>

#include "kernels/framework.hxx"

#include "kernels/block_ldlt.hxx"
#include "kernels/cholesky.hxx"
#include "kernels/ldlt_app.hxx"
#include "kernels/ldlt_nopiv.hxx"
#include "kernels/ldlt_tpp.hxx"

int main(void) {
   int nerr = 0;

   // Enable trapping of bad numerics (NB: can give false positives
   // eg in upper triangle if it's deliberately allowed to contain rubbish)
#if 0
   feenableexcept(FE_INVALID | // NaNs
                  FE_OVERFLOW | // Infs
                  FE_DIVBYZERO); // divide by zero
#endif

   nerr += run_cholesky_tests();
   nerr += run_ldlt_nopiv_tests();
   nerr += run_ldlt_tpp_tests();
   nerr += run_block_ldlt_tests();
   nerr += run_ldlt_app_tests();

   if(nerr==0) {
      printf(ANSI_COLOR_BLUE "\n====================================\n"
             ANSI_COLOR_GREEN  "   All tests passed sucessfully\n"
             ANSI_COLOR_BLUE   "====================================\n"
             ANSI_COLOR_RESET);
   } else {
      printf(ANSI_COLOR_BLUE "\n====================================\n"
             ANSI_COLOR_RED    "   %d tests FAILED!\n"
             ANSI_COLOR_BLUE  "====================================\n"
             ANSI_COLOR_RESET, nerr);
   }

   return nerr;
}
