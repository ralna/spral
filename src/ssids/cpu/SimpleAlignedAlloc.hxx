/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

#include <memory>

#include "omp.hxx"

namespace spral { namespace ssids { namespace cpu {

/**
 * \brief Simple allocator of aligned memory.
 *
 * Real pointer address to free is located at ptr - 1.
 */
template <typename T>
class SimpleAlignedAllocator {
#if defined(__AVX512F__)
  int const align = 64;
#elif defined(__AVX__)
  int const align = 32;
#else
  int const align = 16;
#endif
public:
   typedef T value_type;

   SimpleAlignedAllocator(size_t sz) {}

   template<typename U>
   SimpleAlignedAllocator(SimpleAlignedAllocator<U> const& other)
   {}

   T* allocate(std::size_t n)
   {
      // alloc sufficient space to allow for alignment and storage of base addr
      size_t size = n*sizeof(T) + align;
      void *mem = malloc( size+sizeof(void*) );
      // align it, allowing for base addr beforehand
      void *ptr = reinterpret_cast<void*>(
            reinterpret_cast<char*>(mem) + sizeof(void*)
            );
      if(!std::align(align, n*sizeof(T), ptr, size)) throw std::bad_alloc();
      // store base address to left of returned pointer
      void **vptr = reinterpret_cast<void**>(ptr) - 1;
      *vptr = mem;
      // finally return pointer
      return static_cast<T*>(ptr);
   }

   void deallocate(T* ptr, std::size_t n)
   {
      // extract base address
      void **vptr = reinterpret_cast<void**>(ptr) - 1;
      // free it
      free(*vptr);
   }
};

}}} /* namespaces spral::ssids::cpu */
