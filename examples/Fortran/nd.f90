! examples/Fortran/nd.f90 - Example code for SPRAL_ND package
program example
   use spral_nd
   implicit none

   ! Derived types
   type (nd_options) :: options
   type (nd_inform) :: inform

   ! Matrix data
   integer :: n, row(14), ptr(9), perm(8)

   ! Data for matrix:
   !     1  2  3  4  5  6  7  8
   ! 1 (    X     X  X  X       )
   ! 2 ( X        X     X  X    )
   ! 3 (             X  X       )
   ! 4 ( X  X              X    )
   ! 5 ( X     X        X     X )
   ! 6 ( X  X  X     X     X  X )
   ! 7 (    X     X     X       )
   ! 8 (             X  X       )
   n = 8 
   ptr(1:n+1)        = (/ 1,          5,       8,   10, 11,  13,   15, 15, 15 /)
   row(1:ptr(n+1)-1) = (/ 2, 4, 5, 6, 4, 6, 7, 5, 6, 7, 6, 8, 7, 8 /)

   ! Find nested dissection ordering
   call nd_order(0, n, ptr, row, perm, options, inform)

   ! Print out nested dissection ordering
   write(*,"(a, 8i8)") ' Permutation : ', perm(:)
end program example
