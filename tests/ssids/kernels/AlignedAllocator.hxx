/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

#include <cstddef>
#include <cstdlib>
#ifdef _WIN32
# include <malloc.h>  // _aligned_malloc, _aligned_free
#endif
#include <new>

namespace spral { namespace test {

template <class T>
class AlignedAllocator {
public:
  // Number of bytes boundary we align to
#if defined(__AVX512F__)
  const int alignment = 64;
#elif defined(__AVX__)
  const int alignment = 32;
#else
  const int alignment = 16;
#endif

   typedef T value_type;

   AlignedAllocator() = default;
   template <class U>
   AlignedAllocator(const AlignedAllocator<U>&) {}

   T* allocate(size_t n) {
      // NB: size allocated must be multiple of alignment
      size_t sz = alignment*( ( n*sizeof(T)-1 )/alignment + 1);
#ifdef _WIN32
      auto ptr = _aligned_malloc(sz, alignment);
#elif defined(_ISOC11_SOURCE) /* FIXME: is this ever true in C++? */
      auto ptr = aligned_alloc(alignment, sz);
#else /* _POSIX_C_SOURCE > 200112L || _XOPEN_SOURCE >= 600 */
      void *ptr;
      int test = posix_memalign(&ptr, alignment, sz);
      if(test) throw std::bad_alloc();
#endif
      if(!ptr) throw std::bad_alloc();
      return static_cast<T*>(ptr);
   }

   void deallocate(T *p, size_t n) {
#ifdef _WIN32
      _aligned_free(p);
#else
      free(p);
#endif
   }
};

template<typename T, typename U>
bool operator== (const AlignedAllocator<T> &, const AlignedAllocator<U> &) {
   return true;
}

template<typename T, typename U>
bool operator!= (const AlignedAllocator<T> &lhs, const AlignedAllocator<U> &rhs) {
   return !(lhs==rhs);
}

}} /* namespaces spral::test */
