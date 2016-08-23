/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#include "ssids/cpu/profile.hxx"

#ifdef PROFILE
struct timespec spral::ssids::cpu::Profile::tstart;
#endif

using namespace spral::ssids::cpu;

extern "C"
void spral_ssids_cpu_profile_begin() {
   Profile::init();
}

extern "C"
void spral_ssids_cpu_profile_end() {
   Profile::end();
}
