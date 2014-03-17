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

   write(*, "(a,i3,a,i3,a,i3,a)") &
      "Generating a ", m, " x", n, " non-singular matrix with ", nnz, &
      " non-zeroes"
   call generate_random_matrix(state, 0, m, n, nnz, ptr, row, flag, val=val, &
      nonsingular=.true.)

   write(*, "(a)") "Generated matrix:"
   call print_matrix(6, -1, 0, m, n, ptr, row, val=val)

end program random_matrix_example
