!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! General utility functions needed by all nd modules
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module spral_nd_util
   use spral_nd_types
   implicit none

   ! NB: Everything is public

contains

!
! Prints out error associated with flag (if options is set to print errors)
!
subroutine nd_print_error(flag, options, context)
   integer, intent(in) :: flag ! Error flag to print message for
   type(nd_options), intent(in) :: options
   character(len=*), intent(in) :: context ! context to print with message

   if(options%print_level.lt.0) return ! Print level too low, don't print
   if(options%unit_error.le.0) return ! Invalid unit, don't print

   ! Otherwise print the error message
   call nd_print_message(flag, options%unit_error, context)
end subroutine nd_print_error

subroutine nd_print_diagnostic(diag_level, options, text)
   integer, intent(in) :: diag_level ! level of diagnostic this is
   type(nd_options), intent(in) :: options
   character(len=*), intent(in) :: text ! text to print

   if(options%print_level.lt.diag_level) return ! Above print level, don't print
   if(options%unit_diagnostics.le.0) return ! Invalid unit, don't print

   write (options%unit_diagnostics,'(a)') text
end subroutine nd_print_diagnostic

!
! Prints out errors and warnings according to value of flag
!
subroutine nd_print_message(flag,unit,context)
   integer, intent(in) :: flag ! Error flag to print message for
   integer, intent(in) :: unit ! Fortran unit to print message to
   character(len=*), intent(in) :: context ! context to print with message

   if (unit.lt.0) return ! No output

   if (flag<0) write (unit,fmt='('' ERROR: '')')
   write (unit,advance='no',fmt='('' '', a,'':'')') &
      context(1:len_trim(context))

   select case (flag)
   case (0)
      write (unit,'(A)') ' successful completion'
   case (ND_ERR_MEMORY_ALLOC)
      write (unit,'(A)') ' memory allocation failure'
   case (ND_ERR_MEMORY_DEALLOC)
      write (unit,'(A)') ' memory deallocation failure'
   case (ND_ERR_N)
      write (unit,'(A)') ' n<1'
   end select
end subroutine nd_print_message

!
! Returns ptr(idx) if idx.le.n, or ne+1 otherwise
!
integer function nd_get_ptr(idx, n, ne, ptr)
   integer, intent(in) :: idx, n, ne
   integer, dimension(n), intent(in) :: ptr

   if(idx.le.n) then
      nd_get_ptr = ptr(idx)
   else
      nd_get_ptr = ne+1
   endif
end function nd_get_ptr

!
! Calculates the cost function used to determine whether a particular
! partition is better than another or not.
!
subroutine cost_function(a_weight_1, a_weight_2, a_weight_sep, sumweight, &
      balance_tol, imbal, costf, tau)
   integer, intent(in) :: a_weight_1, a_weight_2, a_weight_sep ! Weighted
      ! size of partitions and separator
   integer, intent(in) :: sumweight
   real(wp), intent(in) :: balance_tol
   logical, intent(in) :: imbal ! Use penalty function?
   integer, intent(in) :: costf
   real(wp), intent(out) :: tau

   real(wp) :: a_wgt1, a_wgt2, balance 

   real(wp), parameter :: beta = 0.5

   a_wgt1 = max(1,a_weight_1)
   a_wgt2 = max(1,a_weight_2)
   balance = real(max(a_wgt1,a_wgt2)) / min(a_wgt1,a_wgt2)

   if (costf.le.1) then
      tau = a_weight_sep / real(a_wgt1 * a_wgt2)
      if (imbal .and. balance.ge.balance_tol) &
         tau = tau + sumweight-2
   else
      tau = a_weight_sep * (1.0_wp + (beta*abs(a_wgt1-a_wgt2))/sumweight)
      if (imbal .and. balance.ge.balance_tol) &
         tau = tau + sumweight * (1.0+beta)
   end if

end subroutine cost_function

end module spral_nd_util
