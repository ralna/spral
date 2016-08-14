! GPU version
! Specifies exact implementation of user facing types
module spral_ssids_type_select
   use spral_ssids_akeep, only : &
      ssids_akeep => ssids_akeep_base
   use spral_ssids_fkeep, only : &
      ssids_fkeep => ssids_fkeep_base
   use spral_ssids_gpu_inform, only : &
      ssids_inform => ssids_inform_gpu
   implicit none
contains
   logical function detect_gpu()
      ! FIXME: actually implement some gpu detection...
      detect_gpu = .true.
   end function detect_gpu
end module spral_ssids_type_select
