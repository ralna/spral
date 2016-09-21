! examples/Fortran/random.f90 - Example code for SPRAL_RANDOM package
program random_example
   use spral_random
   implicit none

   integer, parameter :: long = selected_int_kind(18)

   type(random_state) :: state
   integer :: seed

   ! Store initial random seed so we can reuse it later
   seed = random_get_seed(state)

   ! Generate some random values
   write(*,"(a)") "Some random values"
   write(*,"(a,f16.12)") "Sample Unif(-1,1)               = ", &
      random_real(state)
   write(*,"(a,f16.12)") "Sample Unif(0,1)                = ", &
      random_real(state, positive=.true.)
   write(*,"(a,i16)") "Sample Unif(1, ..., 20)         = ", &
      random_integer(state, 20)
   write(*,"(a,i16)") "Sample Unif(1, ..., 20*huge(0)) = ", &
      random_integer(state, 20_long*huge(0))
   write(*,"(a,l16)") "Sample B(1,0.5)                 = ", &
      random_logical(state)

   ! Restore initial seed
   call random_set_seed(state, seed)

   ! Generate the same random values
   write(*,"(/a)") "The same random values again"
   write(*,"(a,f16.12)") "Sample Unif(-1,1)               = ", &
      random_real(state)
   write(*,"(a,f16.12)") "Sample Unif(0,1)                = ", &
      random_real(state, positive=.true.)
   write(*,"(a,i16)") "Sample Unif(1, ..., 20)         = ", &
      random_integer(state, 20)
   write(*,"(a,i16)") "Sample Unif(1, ..., 20*huge(0)) = ", &
      random_integer(state, 20_long*huge(0))
   write(*,"(a,l16)") "Sample B(1,0.5)                 = ", &
      random_logical(state)

end program random_example
