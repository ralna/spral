! examples/Fortran/scaling/equilib_unsym.f90 - Example code for SPRAL_SCALING
program equilib_scale_unsym_example
   use spral_scaling
   use spral_matrix_util, only : print_matrix, &
                                 SPRAL_MATRIX_REAL_UNSYM
   implicit none

   ! Derived types
   type (equilib_options)  :: options
   type (equilib_inform)   :: inform

   ! Parameters
   integer, parameter :: wp = kind(0.0d0)

   ! Matrix data
   integer :: m, n, ptr(6), row(10)
   real(wp) :: val(10)

   ! Other variables
   integer :: i, j
   real(wp) :: rscaling(5), cscaling(5)

   ! Data for unsymmetric matrix:
   ! ( 2  5         )
   ! ( 1  4       7 )
   ! (    1     2   )
   ! (       3      )
   ! (    8       2 )
   m = 5; n = 5
   ptr(1:n+1)        = (/ 1,        3,                  7,   8,   9,       11 /)
   row(1:ptr(n+1)-1) = (/ 1,   2,   1,   2,   3,   5,   4,   3,   2,   5   /)
   val(1:ptr(n+1)-1) = (/ 2.0, 1.0, 5.0, 4.0, 1.0, 8.0, 3.0, 2.0, 7.0, 2.0 /)
   write(*, "(a)") "Initial matrix:"
   call print_matrix(6, -1, SPRAL_MATRIX_REAL_UNSYM, m, n, ptr, row, val)

   ! Perform unsymmetric scaling
   call equilib_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, options, &
      inform)
   if(inform%flag<0) then
      write(*, "(a, i5)") "equilib_scale_unsym() returned with error ", &
         inform%flag
      stop
   endif

   ! Print scaling and matching
   write(*,"(a,10es10.2)") 'Row Scaling: ', rscaling(1:m)
   write(*,"(a,10es10.2)") 'Col Scaling: ', cscaling(1:n)

   ! Calculate scaled matrix and print it
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         val(j) = rscaling(row(j)) * val(j) * cscaling(i)
      end do
   end do
   write(*, "(a)") "Scaled matrix:"
   call print_matrix(6, -1, SPRAL_MATRIX_REAL_UNSYM, m, n, ptr, row, val)

end program equilib_scale_unsym_example
