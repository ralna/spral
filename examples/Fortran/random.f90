program random_example
   use spral_random
   implicit none

   type(random_state) :: state
   integer :: seed

   ! Store initial random seed so we can reuse it later
   seed = random_get_seed(state)

   ! Generate some random values
   print *, "Sample Unif(-1,1)       = ", random_real(state)
   print *, "Sample Unif(0,1)        = ", random_real(state, positive=.true.)
   print *, "Sample Unif(1, ..., 20) = ", random_integer(state, 20)
   print *, "Sample B(1,0.5)         = ", random_logical(state)

   ! Restore initial seed
   call random_set_seed(state, seed)

   ! Generate the same random values
   print *, "Sample Unif(-1,1)       = ", random_real(state)
   print *, "Sample Unif(0,1)        = ", random_real(state, positive=.true.)
   print *, "Sample Unif(1, ..., 20) = ", random_integer(state, 20)
   print *, "Sample B(1,0.5)         = ", random_logical(state)

end program random_example
