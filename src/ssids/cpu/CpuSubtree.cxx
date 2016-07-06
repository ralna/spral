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
#include "CpuSubtree.hxx"

#include "StackAllocator.hxx"

using namespace spral::ssids::cpu;

extern "C"
void* spral_ssids_create_cpu_subtree_dbl(bool posdef, void const* symbolic_subtree_ptr, int nnodes, struct cpu_node_data<double>* nodes) {
   const int BLOCK_SIZE = 16;
   const int PAGE_SIZE = 16384;
   typedef double T;

   auto const& symbolic_subtree = *static_cast<SymbolicSubtree const*>(symbolic_subtree_ptr);

   return (posdef) ? (void*) new CpuSubtree
                     <true, BLOCK_SIZE, T, StackAllocator<PAGE_SIZE>>
                     (symbolic_subtree, nnodes, nodes)
                   : (void*) new CpuSubtree
                     <false, BLOCK_SIZE, T, StackAllocator<PAGE_SIZE>>
                     (symbolic_subtree, nnodes, nodes);
}

extern "C"
void spral_ssids_destroy_cpu_subtree_dbl(bool posdef, void* subtree_ptr) {
   const int BLOCK_SIZE = 16;
   const int PAGE_SIZE = 16384;
   typedef double T;

   if(posdef) {
      auto *subtree = static_cast<CpuSubtree<true, BLOCK_SIZE, T, StackAllocator<PAGE_SIZE>>*>(subtree_ptr);
      delete subtree;
   } else {
      auto *subtree = static_cast<CpuSubtree<false, BLOCK_SIZE, T, StackAllocator<PAGE_SIZE>>*>(subtree_ptr);
      delete subtree;
   }
}
