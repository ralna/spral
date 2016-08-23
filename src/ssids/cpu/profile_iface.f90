module spral_ssids_cpu_profile
   implicit none

   private
   public :: cpu_profile_begin, &
             cpu_profile_end

   interface
      subroutine cpu_profile_begin() &
            bind(C, name="spral_ssids_cpu_profile_begin")
      end subroutine cpu_profile_begin
      subroutine cpu_profile_end() &
            bind(C, name="spral_ssids_cpu_profile_end")
      end subroutine cpu_profile_end
   end interface

end module spral_ssids_cpu_profile
