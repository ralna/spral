/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 *
 *  \brief
 *  Additional support functions and wrappers for OpenMP.
 */
#pragma once

#include <omp.h>

/* This file wraps the C interface for OpenMP in C++ for style/safety */
namespace spral { namespace omp {

/**
 * \brief Safe wrapper around omp_lock_t ensuring init/cleanup.
 *        See AcquiredLock for locking functionality.
 *
 * This acts as an underlying resource that may be aquired by instantiating
 * an AcquiredLock with this as an argument.
 *
 * \sa AcquiredLock
 */
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
private:
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

   omp_lock_t lock_;

   friend class AcquiredLock;
};

/**
 * \brief RAII lock. Acquires lock on construction, releases on destruction.
 */
class AcquiredLock {
public:
   AcquiredLock(Lock& lock)
   : lock_(lock)
   {
      lock_.set();
   }
   ~AcquiredLock() {
      lock_.unset();
   }
private:
   Lock& lock_; ///< Underlying lock.
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
