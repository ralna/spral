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

#include <memory>

#ifdef __INTEL_COMPILER
// Define our own std::align as intel (at least <=16.0.3) don't
namespace std {
void* align(std::size_t alignment, std::size_t size, void*& ptr, std::size_t& space) {
   auto cptr = reinterpret_cast<uintptr_t>(ptr);
   auto pad = cptr % alignment;
   if(pad == 0) return ptr;
   pad = alignment - pad;
   cptr += pad;
   space -= pad;
   ptr = reinterpret_cast<void*>(cptr);
   return ptr;
}
} /* namespace std */
#endif

namespace spral { namespace ssids { namespace cpu {

namespace pool_alloc_internal {

/** A single fixed size page of memory with allocate function.
 * We are required to guaruntee it is zero'd, so use calloc rather than anything
 * else for the allocation.
 * Deallocation is not supported.
 */
class Page {
   static const int align = 32; // 32 byte alignment
public:
   Page(size_t sz, Page* next=nullptr)
   : next(next), mem_(calloc(sz+align, 1)), ptr_(mem_), space_(sz+align)
   {}
   ~Page() {
      free(mem_);
   }
   void* allocate(size_t sz) {
      if(!std::align(32, sz, ptr_, space_)) return nullptr;
      void* ret = ptr_;
      ptr_ = (char*)ptr_ + sz;
      space_ -= sz;
      return ret;
   }
public:
   Page* const next;
private:
   void *const mem_; // Pointer to memory so we can free it
   void *ptr_; // Next address to return
   size_t space_; // Amount of free memory
};

/** A memory allocation pool consisting of one or more pages.
 * Deallocation is not supported.
 */
class Pool {
   const size_t PAGE_SIZE = 8*1024*1024; // 8MB
public:
   Pool(size_t initial_size)
   : top_page_(new Page(std::max(PAGE_SIZE, initial_size)))
   {}
   Pool(const Pool&) =delete; // Not copyable
   Pool& operator=(const Pool&) =delete; // Not copyable
   ~Pool() {
      /* Iterate over linked list deleting pages */
      for(Page* page=top_page_; page; ) {
         Page* next = page->next;
         delete page;
         page = next;
      }
   }
   void* allocate(size_t sz) {
      void* ptr;
      #pragma omp critical
      {
         ptr = top_page_->allocate(sz);
         if(!ptr) { // Insufficient space on current top page, make a new one
            top_page_ = new Page(std::max(PAGE_SIZE, sz), top_page_);
            ptr = top_page_->allocate(sz);
         }
      }
      return ptr;
   }
private:
   Page* top_page_;
};

} /* namespace spral::ssids::cpu::pool_alloc_internal */

/** An allocator built on top of a pool of pages.
 * Deallocation is not supported.
 */
template <typename T>
class PoolAlloc {
public :
   typedef T               value_type;

   PoolAlloc(size_t initial_size)
   : pool_(new pool_alloc_internal::Pool(initial_size))
   {}

   /** Rebind a type T to a type U PoolAlloc */
   template <typename U>
   PoolAlloc(PoolAlloc<U> &other)
   : pool_(other.pool_)
   {}

   T* allocate(std::size_t n) {
      return static_cast<T*>(pool_->allocate(n*sizeof(T)));
   }
   void deallocate(T* p, std::size_t n) {
      throw std::runtime_error("Deallocation not supported on PoolAlloc");
   }
   template<class U>
   bool operator==(PoolAlloc<U> const& rhs) {
      return true;
   }
   template<class U>
   bool operator!=(PoolAlloc<U> const& rhs) {
      return !(*this==rhs);
   }
protected:
   std::shared_ptr<pool_alloc_internal::Pool> pool_;
   template <typename U> friend class PoolAlloc;
};

}}} /* namepsace spral::ssids::cpu */
