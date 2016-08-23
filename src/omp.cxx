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

int get_global_thread_num() {
   int nbelow = 1;
   int thread_num = 0;
   for(int level=omp_get_level(); level>0; --level) {
      thread_num += nbelow * omp_get_ancestor_thread_num(level);
      nbelow *= omp_get_team_size(level);
   }
   return thread_num;
}

bool cancel_support() {
   return omp_get_cancellation();
}
bool nested_support() {
   return omp_get_nested();
}
bool proc_bind_support() {
   return (omp_get_proc_bind() != omp_proc_bind_false);
}

static int cancel_warning_issued = 0;
void warn_if_no_cancel() {
   if(!cancel_support()) {
      int warn;
      #pragma omp atomic capture
      warn = cancel_warning_issued++;
      if(!warn)
         printf("\nOpenMP cancellation support is not enabled.\n"
                "For best performance, enable cancellation support by setting\n"
                "the environment variable OMP_CANCELLATION=true.\n");
   }
}

static int nested_warning_issued = 0;
void warn_if_no_nested() {
   if(!nested_support()) {
      int warn;
      #pragma omp atomic capture
      warn = nested_warning_issued++;
      if(!warn)
         printf("\nOpenMP nested parallelism support is not enabled.\n"
                "For best performance, enable nested support by setting\n"
                "the environment variable OMP_NESTED=true.\n");
   }
}

static int proc_bind_warning_issued = 0;
void warn_if_no_proc_bind() {
   if(!proc_bind_support()) {
      int warn;
      #pragma omp atomic capture
      warn = proc_bind_warning_issued++;
      if(!warn)
         printf("\nOpenMP processor binding support is not enabled.\n"
                "For best performance, enable support by setting\n"
                "the environment variable OMP_PROC_BIND=true.\n");
   }
}

}} /* namepsace spral::omp */
