/* Copyright 2014-6 The Science and Technology Facilities Council (STFC)
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

namespace spral { namespace ssids { namespace cpu {

/** Class for stack-based allocation.
 *
 * Designed to make quick allocation and deallocation of fronts in a stack
 * based fashion more efficient
 */
template<size_t PAGE_SIZE>
class StackAllocator {
private:
   class Page {
   private:
      char mem[PAGE_SIZE];
   public:
      Page *prev, *next;
      size_t top;
      Page(Page *prev)
         : prev(prev), next(NULL), top(0)
         {}
      void *alloc(size_t size) {
         if(top+size > PAGE_SIZE) return NULL;
         void *ptr = mem + top;
         top += size;
         return ptr;
      }
      int free(size_t size) {
         top -= size;
         return top;
      }
   };
   Page *current, *last;
public:
   StackAllocator(void)
      : current(new Page(NULL)), last(current)
      {}
   ~StackAllocator() {
      // Work backwards freeing Pages
      while(last) {
         Page *prev = last->prev;
         delete last;
         last = prev;
      }
   }
   
   void *alloc(size_t size) {
      if(size > PAGE_SIZE/2) {
         // Requesting a large block, use standard allocator
         return new char[size];
      } else {
         // Requesting a small block
         void *ptr = current->alloc(size);
         if(!ptr) {
            // Ran out of space on page, find a new one
            if(current->next) {
               // Reuse existing page
               current = current->next;
            } else {
               // Allocate a new one
               last = new Page(current);
               current->next = last;
               current = last;
            }
            // Allocation will defintely work this time
            ptr = current->alloc(size);
         }
         // Return allocation
         return ptr;
      }
   }

   void free(void *ptr, size_t size) {
      if(size > PAGE_SIZE/2) {
         // Was directly allocated
         delete[] (char *) ptr;
      } else {
         // Call page's free function
         int remain = current->free(size);
         if(remain==0 && current->prev) {
            // Page now clear (and not first page), iterate backwards
            current = current->prev;
         }
      }
   }
};

}}} /* namespaces spral::ssids::cpu */
