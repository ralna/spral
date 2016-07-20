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

#include <memory>
#include <type_traits>
#include <vector>

namespace spral { namespace ssids { namespace cpu {

/** Pool of blocks that can be aquired/released, providing a way to cap
 *  memory usage whilst avoiding fragmentation.
 *  Further, guaruntees blocks are aligned to 32-byte boundaries so are
 *  suitable for AVX usage. */
template <typename T, std::size_t BLOCK_SIZE>
class BlockPool {
   static const std::size_t align_ = 32; //< Alignement for AVX is 32 bytes
   /// add sufficient padding to each block for alignment through typedef
   typedef typename std::aligned_storage<BLOCK_SIZE*BLOCK_SIZE*sizeof(T), align_>::type block_type;
public:
   /* Not copyable */
   BlockPool(BlockPool const&) =delete;
   BlockPool& operator=(BlockPool const&) =delete;
   /** Constructor allocates memory pool */
   BlockPool(std::size_t num_blocks)
   : mem_(new block_type[num_blocks])
   {
      // Set up stack of free blocks such that we issue them in order
      pool_.reserve(num_blocks);
      for(int i=num_blocks-1; i>=0; --i) {
         pool_.push_back(reinterpret_cast<T*>(mem_ + i));
      }
      // Initialize lock
      omp_init_lock(&lock);
   }
   ~BlockPool() {
      // FIXME: Throw an exception if we've not had all memory returned? 
      omp_destroy_lock(&lock);
      delete[] mem_;
   }

   /** Get next free block in a thread-safe fashion.
    *  Return nullptr if it can't find a free block.
    */
   T *get_nowait() {
      T *ptr = nullptr;
      omp_set_lock(&lock);
         if(pool_.size() > 0) {
            ptr = pool_.back();
            pool_.pop_back();
         }
      omp_unset_lock(&lock);
      return ptr;
   }
   /** Get next free block in a thread-safe fashion.
    *  Keep trying until it suceeds, use taskyield after each try.
    *  NB: This may deadlock if there are no other tasks releasing blocks.
    */
   T *get_wait() {
      while(true) {
         T *ptr = get_nowait();
         if(ptr) {
            return ptr;
         }
         #pragma omp taskyield
      }
   }
   /** Release block acquired using get_*() function for reuse */
   void release(T *const ptr) {
      omp_set_lock(&lock);
         pool_.push_back(ptr);
      omp_unset_lock(&lock);
   }
private:
   block_type* const mem_; //< pointer to allocated memory
   std::vector<T*> pool_; //< stack of free blocks
   omp_lock_t lock; //< lock used for safe access to pool_
};

}}} /* namespace spral::ssids::cpu */
