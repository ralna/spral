! examples/Fortran/ssids.f90 - Example code for SPRAL_SSIDS package
program ssids_example
   use spral_ssids
   implicit none

   integer, parameter :: long = selected_int_kind(18)

   ! Derived types
   type (ssids_akeep)   :: akeep
   type (ssids_fkeep)   :: fkeep
   type (ssids_options) :: options
   type (ssids_inform)  :: inform

   ! Parameters
   integer, parameter :: wp = kind(0.0d0)

   ! Matrix data
   logical :: posdef
   integer :: n, row(9)
   integer(long) :: ptr(6)
   real(wp) :: val(9)

   ! Other variables
   integer :: piv_order(5), cuda_error
   real(wp) :: x(5)
   logical :: check

   ! Data for matrix:
   ! ( 2  1         )
   ! ( 1  4  1    1 )
   ! (    1  3  2   )
   ! (       2 -1   )
   ! (    1       2 )
   posdef = .false.
   n = 5
   ptr(1:n+1)        = (/ 1,        3,             6,        8,    9,  10 /)
   row(1:ptr(n+1)-1) = (/ 1,   2,   2,   3,   5,   3,   4,   4,    5   /)
   val(1:ptr(n+1)-1) = (/ 2.0, 1.0, 4.0, 1.0, 1.0, 3.0, 2.0, -1.0, 2.0 /)

   ! The right-hand side with solution (1.0, 2.0, 3.0, 4.0, 5.0)
   x(1:n) = (/ 4.0, 17.0, 19.0, 2.0, 12.0 /)

   ! Perform analyse and factorise with data checking
   check = .true.
   call ssids_analyse(check, n, ptr, row, akeep, options, inform)
   if(inform%flag<0) go to 100
   call ssids_factor(posdef, val, akeep, fkeep, options, inform)
   if(inform%flag<0) go to 100

   ! Solve
   call ssids_solve(x,akeep,fkeep,options,inform)
   if(inform%flag<0) go to 100
   write(*,'(a,/,(3es18.10))') ' The computed solution is:', x(1:n)

   ! Determine and print the pivot order
   call ssids_enquire_indef(akeep, fkeep, options, inform, piv_order=piv_order)
   write(*,"(a,10i5)") ' Pivot order:', piv_order(1:n)

   100 continue
   call ssids_free(akeep, fkeep, cuda_error)

end program ssids_example
