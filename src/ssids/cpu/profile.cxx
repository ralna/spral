#include "profile.hxx"

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
