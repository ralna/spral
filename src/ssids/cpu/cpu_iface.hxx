/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

#include <cstddef>

namespace spral { namespace ssids { namespace cpu {

enum struct PivotMethod : int {
   app_aggressive = 1,
   app_block      = 2,
   tpp            = 3
};

struct cpu_factor_options {
   int print_level;
   bool action;
   double small;
   double u;
   double multiplier;
   long cpu_small_subtree_threshold;
   int cpu_task_block_size;
   PivotMethod pivot_method;
};

/** Return nearest value greater than supplied lda that is multiple of alignment */
template<typename T>
size_t align_lda(size_t lda) {
   int const align = 32;
   static_assert(align % sizeof(T) == 0, "Can only align if T divides align");
   int const Talign = align / sizeof(T);
   return Talign*((lda-1)/Talign + 1);
}

}}} /* namespaces spral::ssids::cpu */
