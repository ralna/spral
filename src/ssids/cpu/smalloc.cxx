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
#include "smalloc.hxx"

#include "factor_iface.hxx"

namespace spral { namespace ssids { namespace cpu {

template<>
double *smalloc<double>(void *alloc, size_t len) {
   return spral_ssids_smalloc_dbl(alloc, len);
}
template<>
int *smalloc<int>(void *alloc, size_t len) {
   return spral_ssids_smalloc_int(alloc, len);
}

}}} /* end of namespace spral::ssids::cpu */
