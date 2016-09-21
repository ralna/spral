! examples/Fortran/rutherford_boeing/rb_write.f90
! Example code for SPRAL_RUTHERFORD_BOEING
program rb_write_example
   use spral_rutherford_boeing
   use spral_matrix_util, only: SPRAL_MATRIX_REAL_SYM_INDEF
   implicit none

   ! Parameters
   integer, parameter :: long = selected_int_kind(18)
   integer, parameter :: wp = kind(0d0)

   ! Matrix data
   integer :: n
   integer(long) :: ptr(6)
   integer :: row(8)
   real(wp) :: val(8)

   ! rb_read options and inform
   type(rb_write_options) :: options
   integer :: inform

   ! Data for symmetric matrix
   ! ( 2  1         )
   ! ( 1  4  1    8 )
   ! (    1  3  2   )
   ! (       2      )
   ! (    8       2 )
   n = 5
   ptr(1:n+1)        = (/ 1,        3,             6,      8,8,   9 /)
   row(1:ptr(n+1)-1) = (/ 1,   2,   2,   3,   5,   3,   4,   5   /)
   val(1:ptr(n+1)-1) = (/ 2.0, 1.0, 4.0, 1.0, 8.0, 3.0, 2.0, 2.0 /)

   ! Write matrix
   call rb_write("matrix.rb", SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, &
      options, inform, val=val, title="SPRAL_RUTHERFORD_BOEING test matrix")

end program rb_write_example
