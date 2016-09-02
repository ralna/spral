/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

#include <cstddef>
#include <cstdlib>
#include <new>

namespace spral { namespace test {

template <class T>
class AlignedAllocator {
public:
   const int alignment = 32; // Number of bytes boundary we align to

   typedef T value_type;

   AlignedAllocator() = default;
   template <class U>
   AlignedAllocator(const AlignedAllocator<U>&) {}

   T* allocate(size_t n) {
      // NB: size allocated must be multiple of alignment
      size_t sz = alignment*( ( n*sizeof(T)-1 )/alignment + 1);
#ifdef _ISOC11_SOURCE /* FIXME: is this ever true in C++? */
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
      free(p);
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
