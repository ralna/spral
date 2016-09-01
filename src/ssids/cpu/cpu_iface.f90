!> \file
!> \copyright 2016 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Jonathan Hogg
module spral_ssids_cpu_iface
   use, intrinsic :: iso_c_binding
   use spral_ssids_datatypes, only : ssids_options
   use spral_ssids_inform, only : ssids_inform
   implicit none

   private
   public :: cpu_factor_options, cpu_factor_stats
   public :: cpu_copy_options_in, cpu_copy_stats_out

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> @brief Interoperable subset of ssids_options
   !> @details Interoperates with cpu_factor_options C++ type
   !> @sa spral_ssids_datatypes::ssids_options
   !> @sa spral::ssids::cpu::cpu_factor_options
   type, bind(C) :: cpu_factor_options
      real(C_DOUBLE) :: multiplier
      real(C_DOUBLE) :: small
      real(C_DOUBLE) :: u
      integer(C_INT) :: print_level
      integer(C_LONG) :: cpu_small_subtree_threshold
      integer(C_INT) :: cpu_task_block_size
      integer(C_INT) :: pivot_method
   end type cpu_factor_options

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> @brief Interoperable subset of ssids_inform
   !> @details Interoperates with ThreadStats C++ type
   !> @sa spral_ssids_inform::ssids_inform
   !> @sa spral::ssids::cpu::ThreadStats
   type, bind(C) :: cpu_factor_stats
      integer(C_INT) :: flag
      integer(C_INT) :: num_delay
      integer(C_INT) :: num_neg
      integer(C_INT) :: num_two
      integer(C_INT) :: num_zero
      integer(C_INT) :: maxfront
      integer(C_INT) :: not_first_pass
      integer(C_INT) :: not_second_pass
   end type cpu_factor_stats

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Copy subset of ssids_options to interoperable type
subroutine cpu_copy_options_in(foptions, coptions)
   type(ssids_options), intent(in) :: foptions
   type(cpu_factor_options), intent(out) :: coptions

   coptions%multiplier     = foptions%multiplier
   coptions%small          = foptions%small
   coptions%u              = foptions%u
   coptions%print_level    = foptions%print_level
   coptions%cpu_small_subtree_threshold = foptions%cpu_small_subtree_threshold
   coptions%cpu_task_block_size         = foptions%cpu_task_block_size
   coptions%pivot_method   = foptions%pivot_method
end subroutine cpu_copy_options_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Copy subset of ssids_inform from interoperable type
subroutine cpu_copy_stats_out(n, cstats, finform)
   integer, intent(in) :: n
   type(cpu_factor_stats), intent(in) :: cstats
   type(ssids_inform), intent(inout) :: finform

   ! Copy stats
   finform%flag         = cstats%flag
   finform%num_delay    = cstats%num_delay
   finform%num_neg      = cstats%num_neg
   finform%num_two      = cstats%num_two
   finform%matrix_rank  = n - cstats%num_zero
   finform%maxfront     = cstats%maxfront
   finform%not_first_pass = cstats%not_first_pass
   finform%not_second_pass = cstats%not_second_pass
end subroutine cpu_copy_stats_out

end module spral_ssids_cpu_iface
