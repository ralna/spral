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

! ---------------------------------------------------
! nd_convert_partition_flags
! ---------------------------------------------------
! Given a partition array, convert the partition into a flag array
subroutine nd_convert_partition_flags(a_n,a_n1,a_n2,partition,flag_1, &
    flag_2,flag_sep,flags)
  integer, intent(in) :: a_n ! order of matrix
  integer, intent(in) :: a_n1 ! Size of partition 1
  integer, intent(in) :: a_n2 ! Size of partition 2
  integer, intent(in) :: partition(a_n) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end.
  integer, intent(in) :: flag_1 ! flag for rows in partition 1
  integer, intent(in) :: flag_2 ! flag for rows in partition 2
  integer, intent(in) :: flag_sep ! flag for rows in separator
  integer, intent(out) :: flags(a_n) ! flags(i) contains flag for row i
  ! and indicates which partition it is in

  integer :: j, k

  do j = 1, a_n1
    k = partition(j)
    flags(k) = flag_1
  end do
  do j = a_n1 + 1, a_n1 + a_n2
    k = partition(j)
    flags(k) = flag_2
  end do
  do j = a_n1 + a_n2 + 1, a_n
    k = partition(j)
    flags(k) = flag_sep
  end do

end subroutine nd_convert_partition_flags

! ---------------------------------------------------
! nd_flags_partition
! ---------------------------------------------------
! Given a partition array, convert the partition into a flag array
subroutine nd_convert_flags_partition(a_n,a_n1,a_n2,flags,flag_1, &
    flag_2,partition)
  integer, intent(in) :: a_n ! order of matrix
  integer, intent(in) :: a_n1 ! Size of partition 1
  integer, intent(in) :: a_n2 ! Size of partition 2
  integer, intent(in) :: flags(a_n) ! flags(i) contains flag for row i
  ! and indicates which partition it is in
  integer, intent(in) :: flag_1 ! flag for rows in partition 1
  integer, intent(in) :: flag_2 ! flag for rows in partition 2
  integer, intent(out) :: partition(a_n) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end.

  integer :: i, j, k, l

  i = 1
  j = a_n1 + 1
  k = a_n1 + a_n2 + 1
  do l = 1, a_n
    if (flags(l).eq.flag_1) then
      partition(i) = l
      i = i + 1
    else
      if (flags(l).eq.flag_2) then
        partition(j) = l
        j = j + 1
      else
        partition(k) = l
        k = k + 1
      end if
    end if
  end do

end subroutine nd_convert_flags_partition

end module spral_nd_util
