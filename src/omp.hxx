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

#include <omp.h>

/* This file wraps the C interface for OpenMP in C++ for style/safety */
namespace spral { namespace omp {

/** Light-weight wrapper around omp_lock_t. Disabled move semantics. */
class Lock {
public:
   Lock(Lock const&) =delete;
   Lock& operator=(Lock const&) =delete;
   Lock() {
      omp_init_lock(&lock_);
   }
   ~Lock() {
      omp_destroy_lock(&lock_);
   }
   inline
   void set() {
      omp_set_lock(&lock_);
   }
   inline
   void unset() {
      omp_unset_lock(&lock_);
   }
   inline
   bool test() {
      return omp_test_lock(&lock_);
   }
private:
   omp_lock_t lock_;
};

/// Return global thread number (=thread number if not nested)
int get_global_thread_num();

/// Returns true if omp cancel is supported
bool cancel_support();
/// Returns true if omp nesting is supported
bool nested_support();
/// Returns true if omp processor binding is supported
bool proc_bind_support();

/// Prints a warning message if cancel is not supported
void warn_if_no_cancel();
/// Prints a warning message if nesting is not supported
void warn_if_no_nested();
/// Prints a warning message if processor binding is not supported
void warn_if_no_proc_bind();

}} /* end of namespace spral::omp */
