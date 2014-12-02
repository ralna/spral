! GPU version
! Specifies exact implementation of user facing types
module spral_ssids_type_select
   use spral_ssids_akeep_gpu, only : &
      ssids_akeep => ssids_akeep_gpu
   use spral_ssids_fkeep_gpu, only : &
      ssids_fkeep => ssids_fkeep_gpu
   use spral_ssids_inform_gpu, only : &
      ssids_inform => ssids_inform_gpu
   implicit none
end module spral_ssids_type_select
