! examples/Fortran/rutherford_boeing/rb_read.f90
! Example code for SPRAL_RUTHERFORD_BOEING
program rb_read_example
   use spral_rutherford_boeing
   use spral_matrix_util, only: print_matrix, &
                                SPRAL_MATRIX_REAL_SYM_INDEF
   implicit none

   ! Parameters
   integer, parameter :: long = selected_int_kind(18)
   integer, parameter :: wp = kind(0d0)

   ! Matrix data
   character(len=72) :: title
   integer :: m, n
   integer(long), dimension(:), allocatable :: ptr
   integer, dimension(:), allocatable :: row
   real(wp), dimension(:), allocatable :: val

   ! rb_read options and inform
   type(rb_read_options) :: options
   integer :: inform

   ! Read matrix
   call rb_read("matrix.rb", m, n, ptr, row, val, options, inform, &
      title=title)

   ! Print matrix
   write(*, "(3a)") "Matrix '", trim(title), "'"
   call print_matrix(6, -1, SPRAL_MATRIX_REAL_SYM_INDEF, m, n, ptr, row, val)

end program rb_read_example
