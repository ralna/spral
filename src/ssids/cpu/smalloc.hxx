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

/* Generic wrapper around Fortran-defined smalloc calls */
template<typename T>
T *smalloc(void *alloc, size_t len);

}}} /* end of namespace spral::ssids::cpu */
