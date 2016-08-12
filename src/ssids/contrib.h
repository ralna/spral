#ifndef SPRAL_SSIDS_CONTRIB_H
#define SPRAL_SSIDS_CONTRIB_H

#ifdef __cplusplus
extern "C" {
#endif

void spral_ssids_contrib_get_data(void const* contrib, int* n,
      double const** val, int const** rlist, int* ndelay,
      int const** delay_perm, double const** delay_val, int* lddelay);

void spral_ssids_contrib_free_dbl(void* contrib);

#ifdef __cplusplus
}
#endif

#endif /* SPRAL_SSIDS_CONTRIB_H */
