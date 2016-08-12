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
#include "SymbolicSubtree.hxx"

using namespace spral::ssids::cpu;

extern "C"
void* spral_ssids_cpu_create_symbolic_subtree(
      int n, int sa, int en, int const* sptr, int const* sparent,
      long const* rptr, int const* rlist, int const* nptr, int const* nlist,
      int ncontrib, int const* contrib_idx,
      struct cpu_factor_options const* options) {
   return (void*) new SymbolicSubtree(
         n, sa, en, sptr, sparent, rptr, rlist, nptr, nlist, ncontrib,
         contrib_idx, *options
         );
}

extern "C"
void spral_ssids_cpu_destroy_symbolic_subtree(void* target) {
   if(!target) return;

   auto *subtree = static_cast<SymbolicSubtree*>(target);
   delete subtree;
}
