/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#include "ssids/profile.hxx"

#ifdef PROFILE
struct timespec spral::ssids::Profile::tstart;
#endif

using namespace spral::ssids;

extern "C"
void spral_ssids_profile_begin() {
   Profile::init();
}

extern "C"
void spral_ssids_profile_end() {
   Profile::end();
}
