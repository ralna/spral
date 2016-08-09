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
#include "omp.hxx"

#include <cstdio>

/* This file wraps the C interface for OpenMP in C++ for style/safety */
namespace spral { namespace omp {

bool cancel_support() {
   return omp_get_cancellation();
}

static int cancel_warning_issued = 0;
void warn_if_no_cancel() {
   if(!cancel_support()) {
      int warn;
      #pragma omp atomic capture
      warn = cancel_warning_issued++;
      if(!warn)
         printf("\nCancellation support is not enabled.\n"
                "For best performance, enable cancellation support by setting\n"
                "the environment variable OMP_CANCELLATION=true.\n");
   }
}

}} /* namepsace spral::omp */
