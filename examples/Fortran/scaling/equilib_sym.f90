! examples/Fortran/scaling/equilib_sym.f90 - Example code for SPRAL_SCALING
program equilib_scale_sym_example
   use spral_scaling
   use spral_matrix_util, only : print_matrix, &
                                 SPRAL_MATRIX_REAL_SYM_INDEF
   implicit none

   ! Derived types
   type (equilib_options)  :: options
   type (equilib_inform)   :: inform

   ! Parameters
   integer, parameter :: wp = kind(0.0d0)

   ! Matrix data
   integer :: n, ptr(6), row(8)
   real(wp) :: val(8)

   ! Other variables
   integer :: i, j
   real(wp) :: scaling(5)

   ! Data for symmetric matrix:
   ! ( 2  1         )
   ! ( 1  4  1    8 )
   ! (    1  3  2   )
   ! (       2      )
   ! (    8       2 )
   n = 5
   ptr(1:n+1)        = (/ 1,        3,             6,      8,8,   9 /)
   row(1:ptr(n+1)-1) = (/ 1,   2,   2,   3,   5,   3,   4,   5   /)
   val(1:ptr(n+1)-1) = (/ 2.0, 1.0, 4.0, 1.0, 8.0, 3.0, 2.0, 2.0 /)
   write(*, "(a)") "Initial matrix:"
   call print_matrix(6, -1, SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, val)

   ! Perform symmetric scaling
   call equilib_scale_sym(n, ptr, row, val, scaling, options, inform)
   if(inform%flag<0) then
      write(*, "(a, i5)") "equilib_scale_sym() returned with error ", &
         inform%flag
      stop
   endif

   ! Print scaling and matching
   write(*,"(a,10es10.2)") 'Scaling: ', scaling(1:n)

   ! Calculate scaled matrix and print it
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         val(j) = scaling(i) * val(j) * scaling(row(j))
      end do
   end do
   write(*, "(a)") "Scaled matrix:"
   call print_matrix(6, -1, SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, val)

end program equilib_scale_sym_example
