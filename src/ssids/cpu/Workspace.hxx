/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

#include <memory>

#include "compat.hxx" // in case std::align not defined

namespace spral { namespace ssids { namespace cpu {

/** A Workspace is a chunk of memory that can be reused. The get_ptr<T>(len)
 * function provides a pointer to it after ensuring it is of at least the
 * given size. */
class Workspace {
   static int const align = 32;
public:
   Workspace(size_t sz)
   {
      alloc_and_align(sz);
   }
   ~Workspace() {
      ::operator delete(mem_);
   }
   void alloc_and_align(size_t sz) {
      sz_ = sz+align;
      mem_ = ::operator new(sz_);
      mem_aligned_ = mem_;
      std::align(align, sz, mem_aligned_, sz_);
   }
   template <typename T>
   T* get_ptr(size_t len) {
      if(sz_ < len*sizeof(T)) {
         // Need to resize
         ::operator delete(mem_);
         alloc_and_align(len*sizeof(T));
      }
      return static_cast<T*>(mem_aligned_);
   }
private:
   void* mem_;
   void* mem_aligned_;
   size_t sz_;
};

}}} /* end of namespace spral::ssids::cpu */
