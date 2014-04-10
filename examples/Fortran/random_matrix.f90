! examples/Fortran/random_matrix.f90 - Example code for SPRAL_RANDOM_MATRIX pkg
program random_matrix_example
   use spral_matrix_util, only : print_matrix
   use spral_random, only : random_state
   use spral_random_matrix
   implicit none

   integer, parameter :: wp = kind(0d0)

   integer, parameter :: m=4, n=5, nnz=8
   integer :: ptr(n+1), row(nnz)
   real(wp) :: val(nnz)
   integer :: flag
   type(random_state) :: state

   ! Generate matrix
   write(*, "(a,i3,a,i3,a,i3,a)") &
      "Generating a ", m, " x", n, " non-singular matrix with ", nnz, &
      " non-zeroes"
   call random_matrix_generate(state, 0, m, n, nnz, ptr, row, flag, val=val, &
      nonsingular=.true.)

   ! Print matrix using utility routine from SPRAL_MATRIX_UTILS package
   write(*, "(a)") "Generated matrix:"
   call print_matrix(6, -1, 0, m, n, ptr, row, val=val)

end program random_matrix_example
