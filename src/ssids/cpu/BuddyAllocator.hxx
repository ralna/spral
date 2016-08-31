/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

//#define MEM_STATS

#include <memory>

#include "omp.hxx"

namespace spral { namespace ssids { namespace cpu {

namespace buddy_alloc_internal {

template <typename CharAllocator=std::allocator<char>>
class Page {
   typedef typename std::allocator_traits<CharAllocator>::template rebind_traits<int> IntAllocTraits;
   static int const nlevel=16;
   static int const align=32;
   static int const ISSUED_FLAG = -2; // flag value as issued
public:
   Page(Page const&) =delete; // not copyable
   Page& operator=(Page const&) =delete; // not copyable
   Page(size_t size, CharAllocator const &alloc=CharAllocator())
   : alloc_(alloc)
   {
      min_size_ = std::max(size_t(1), (size-1) / (1<<(nlevel-1)) + 1);
      min_size_ = align * ((min_size_-1)/align + 1); // make muliple of align
      size_ = min_size_<<(nlevel-1);
      /* Allocate memory of sufficient size and align it */
      mem_ = std::allocator_traits<CharAllocator>::allocate(alloc_, size_+align);
      size_t space = size_+align; 
      void* to_align = mem_;
      std::align(align, size, to_align, space);
      base_ = static_cast<char*>(to_align);
      typename IntAllocTraits::allocator_type intAlloc(alloc_);
      next_ = IntAllocTraits::allocate(intAlloc, 1<<(nlevel-1));
      /* Initialize data structures */
      head_[nlevel-1] = 0; next_[0] = -1; // a single free block at top level
      for(int i=0; i<nlevel-1; ++i)
         head_[i] = -1; // ... and no free blocks at other levels
#ifdef MEM_STATS
      printf("BuddyAllocator: Creating new page %p size %ld\n", mem_, size_);
#endif /* MEM_STATS */
   }
   // Move constructor
   Page(Page&& other) noexcept
   : alloc_(other.alloc_), min_size_(other.min_size_), size_(other.size_),
     mem_(other.mem_), base_(other.base_), next_(other.next_)
   {
      other.mem_ = nullptr;
      other.base_ = nullptr;
      other.next_ = nullptr;
      for(int i=0; i<nlevel; ++i)
         head_[i] = other.head_[i];
#ifdef MEM_STATS
      used_ = other.used_;
      max_used_ = other.max_used_;
#endif /* MEM_STATS */
   }
   ~Page() noexcept(false) {
      if(next_ && head_[nlevel-1] != 0)
         throw std::runtime_error("outstanding allocations on cleanup\n");
      if(next_) {
         typename IntAllocTraits::allocator_type intAlloc(alloc_);
         IntAllocTraits::deallocate(intAlloc, next_, 1<<(nlevel-1));
#ifdef MEM_STATS
         printf("BuddyAllocator: Allocated %16ld (%.2e GB)\n", 
               size_, 1e-9*size_);
         printf("BuddyAllocator: Max Used  %16ld (%.2e GB)\n",
               max_used_, 1e-9*max_used_);
#endif /* MEM_STATS */
      }
      if(mem_)
         std::allocator_traits<CharAllocator>::deallocate(
               alloc_, mem_, size_+align
               );
   }
   void* allocate(std::size_t sz) {
      if(sz > size_) return nullptr; // too big: don't even try 
      // Determine which level of block we're trying to find
      int level = sz_to_level(sz);
      void* ptr = addr_to_ptr(get_next_ptr(level));
#ifdef MEM_STATS
      if(ptr) {
         used_ += sz;
         max_used_ = std::max(max_used_, used_);
      }
#endif /* MEM_STATS */
      return ptr;
   }
   void deallocate(void* ptr, std::size_t sz) {
      int idx = ptr_to_addr(ptr);
      int level = sz_to_level(sz);
      mark_free(idx, level);
#ifdef MEM_STATS
      used_ -= sz;
#endif /* MEM_STATS */
   }
   /** Return true if this Page owners given pointer */
   bool is_owner(void* ptr) {
      int idx = ptr_to_addr(ptr);
      return (idx>=0 && idx<(1<<(nlevel-1)));
   }
   size_t count_free() const {
      size_t free=0;
      for(int i=0; i<nlevel; ++i) {
         for(int p=head_[i]; p!=-1; p=next_[p])
            free += (1<<i) * min_size_;
      }
      return free;
   }
   void print() const {
      printf("Page %p size %ld free %ld\n", mem_, size_, count_free());
   }
private:
   /** Returns next ptr at given level, creating one if required.
    *  If we cannot create one, return -1 */
   int get_next_ptr(int level) {
      if(level<0 || level>=nlevel) return -1; // invalid level
      if(head_[level] == -1) {
         // Need to split next level up to get one
         int above = get_next_ptr(level+1);
         if(above==-1) return -1; // couldn't find one
         split_block(level+1, above);
      }
      int p = head_[level];
      head_[level] = next_[p];
      next_[p] = ISSUED_FLAG;
      return p;
   }

   /** Marks given block as free, tries to merge with partner if possible */
   void mark_free(int idx, int level) {
      if(level < nlevel-1) {
         // There exists a partner, see if we can merge with it
         int partner = get_partner(idx, level);
         if(next_[partner] != ISSUED_FLAG) {
            // Partner is free in *some* list, not necessarily this level
            if(remove_from_free_list(partner, level)) {
               // It was this level - we can merge
               mark_free(std::min(idx, partner), level+1);
               return;
            }
         }
      }
      // Otherwise, can't merge, add to free list
      next_[idx] = head_[level];
      head_[level] = idx;
   }

   /** Finds the given address in free list for level and removes it.
    *  Returns false if it cannot be found, true otherwise.
    */
   bool remove_from_free_list(int idx, int level) {
      int prev = -1;
      int current = head_[level];
      while(current!=-1 && current != idx) {
         prev = current;
         current = next_[current];
      }
      if(current != idx) return false; // can't find it
      if(prev==-1) {
         // at the head
         head_[level] = next_[idx];
      } else {
         // in the middle
         next_[prev] = next_[idx];
      }
      return true; // success
   }

   /** Splits the given block */
   void split_block(int level, int block) {
      int left = block;
      int right = get_partner(block, level-1);
      next_[right] = head_[level-1];
      next_[left] = right;
      head_[level-1] = left;
   }

   /** Given address location, return pointer */
   void* addr_to_ptr(int idx) {
      return (idx==-1) ? nullptr : base_ + idx*min_size_;
   }

   /** Given pointer, return address */
   int ptr_to_addr(void* ptr) {
      return
         static_cast<uintptr_t>(static_cast<char*>(ptr)-base_) / min_size_;
   }

   /** Given a size, find the relevant level */
   int sz_to_level(std::size_t sz) {
      int val = sz / min_size_;
      // Find next power of 2 higher than val
      int level = 0;
      while((val>>level) > 0) ++level;
      return level;
   }

   /** Given an index find its partner at given level */
   int get_partner(int idx, int level) {
      return idx ^ (1<<level);
   }

   CharAllocator alloc_;
   size_t min_size_; ///< size of smallest block we can allocate
   size_t size_; ///< maximum size
   char* mem_; ///< pointer to memory allocation
   char* base_; ///< aligned memory base
   int head_[nlevel]; ///< first free block at each level's size
   int *next_; ///< next free block at given level
#ifdef MEM_STATS
   size_t used_ = 0; ///< Total amount used
   size_t max_used_ = 0; ///< High water mark of used_
#endif /* MEM_STATS */
};

template <typename CharAllocator>
class Table {
   typedef Page<CharAllocator> PageSpec;
   typedef typename std::allocator_traits<CharAllocator>::template rebind_alloc<PageSpec> PageAlloc;
public:
   Table(const Table&) =delete;
   Table& operator=(const Table&) =delete;
   Table(std::size_t sz, CharAllocator const& alloc=CharAllocator())
   : alloc_(alloc), max_sz_(sz), pages_(PageAlloc(alloc))
   {
      pages_.emplace_back(max_sz_, alloc_);
   }

   void* allocate(std::size_t sz) {
      // Try allocating in existing pages
      spral::omp::AcquiredLock scopeLock(lock_);
      void* ptr;
      for(auto& page: pages_) {
         ptr = page.allocate(sz);
         if(ptr) break; // allocation suceeded
      }
      if(!ptr) {
         // Failed to alloc on existing page: make a bigger page and use it
#ifdef MEM_STATS
         printf("Failed to allocate %ld on existing page...\n", sz);
         for(auto& page: pages_)
            page.print();
#endif /* MEM_STATS */
         size_t old_max_sz = max_sz_;
         try {
            max_sz_ = std::max(2*max_sz_, sz);
            pages_.emplace_back(max_sz_, alloc_);
         } catch(std::bad_alloc const&) {
            // Failed to alloc block twice as big, try one the same size
            max_sz_ = old_max_sz;
            try {
               max_sz_ = std::max(max_sz_, sz);
               pages_.emplace_back(max_sz_, alloc_);
            } catch(std::bad_alloc const&) {
               // That didn't work either, try one of just big enough for sz
               pages_.emplace_back(sz, alloc_);
               // If this fails, we just give up and propogate std::bad_alloc 
            }
         }
         ptr = pages_.back().allocate(sz);
      }
      return ptr;
   }

   void deallocate(void* ptr, std::size_t sz) {
      spral::omp::AcquiredLock scopeLock(lock_);
      for(auto& page: pages_) {
         if(page.is_owner(ptr)) {
            page.deallocate(ptr, sz);
            break;
         }
      }
   }

private:
   CharAllocator alloc_;
   std::size_t max_sz_;
   std::vector<PageSpec, PageAlloc> pages_;
   spral::omp::Lock lock_;
};

} /* namespace buddy_alloc_internal */

/** Simple buddy-system allocator */
template <typename T, typename BaseAllocator>
class BuddyAllocator {
   typedef typename std::allocator_traits<BaseAllocator>::template rebind_alloc<char> CharAllocator;
public:
   typedef T value_type;

   BuddyAllocator(size_t size, BaseAllocator const& base=BaseAllocator())
   : table_(new buddy_alloc_internal::Table<CharAllocator>(size*sizeof(T), base))
   {}
   template<typename U, typename UBaseAllocator>
   BuddyAllocator(BuddyAllocator<U, UBaseAllocator> const& other)
   : table_(other.table_)
   {}

   T* allocate(std::size_t n)
   {
      return static_cast<T*>(table_.get()->allocate(n*sizeof(T)));
   }

   void deallocate(T* ptr, std::size_t n)
   {
      table_.get()->deallocate(ptr, n*sizeof(T));
   }
private:
   std::shared_ptr<buddy_alloc_internal::Table<CharAllocator>> table_;
   template<typename U, typename UAlloc>
   friend class BuddyAllocator;
};

}}} /* namespaces spral::ssids::cpu */
