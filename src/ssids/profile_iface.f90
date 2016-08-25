!> \file
!> \copyright 2016 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Jonathan Hogg
module spral_ssids_profile
   implicit none

   private
   public :: profile_begin, &
             profile_end

   interface
      subroutine profile_begin() &
            bind(C, name="spral_ssids_profile_begin")
      end subroutine profile_begin
      subroutine profile_end() &
            bind(C, name="spral_ssids_profile_end")
      end subroutine profile_end
   end interface

end module spral_ssids_profile
