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

#include <cstddef>

namespace spral { namespace ssids { namespace cpu {

enum struct PivotMethod : int {
   app_aggressive = 0,
   app_block      = 1,
   tpp            = 2
};

struct cpu_factor_options {
   double multiplier;
   double small;
   double u;
   int print_level;
   int cpu_small_subtree_threshold;
   int cpu_task_block_size;
   PivotMethod pivot_method;
};

struct cpu_factor_stats {
   int flag;
   int num_delay;
   int num_neg;
   int num_two;
   int num_zero;
   int maxfront;
   int not_first_pass;
   int not_second_pass;
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
