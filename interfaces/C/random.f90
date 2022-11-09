real(C_DOUBLE) function spral_random_real(cstate, cpositive) bind(C)
   use iso_c_binding
   use spral_random
   implicit none

   integer(C_INT), intent(inout) :: cstate
   logical(C_BOOL), value :: cpositive

   type(random_state) :: fstate
   logical :: fpositive

   ! Initialize state
   call random_set_seed(fstate, cstate)

   ! Call Fortran routine
   fpositive = cpositive
   spral_random_real = random_real(fstate, positive=fpositive)

   ! Recover state
   cstate = random_get_seed(fstate)
end function spral_random_real

integer(C_INT) function spral_random_integer(cstate, n) bind(C)
   use iso_c_binding
   use spral_random
   implicit none

   integer(C_INT), intent(inout) :: cstate
   integer(C_INT), value :: n

   type(random_state) :: fstate

   ! Initialize state
   call random_set_seed(fstate, cstate)

   ! Call Fortran routine
   spral_random_integer = random_integer(fstate, n)

   ! Recover state
   cstate = random_get_seed(fstate)
end function spral_random_integer

integer(C_INT64_T) function spral_random_long(cstate, n) bind(C)
   use iso_c_binding
   use spral_random
   implicit none

   integer(C_INT), intent(inout) :: cstate
   integer(C_INT64_T), value :: n

   type(random_state) :: fstate

   ! Initialize state
   call random_set_seed(fstate, cstate)

   ! Call Fortran routine
   spral_random_long = random_integer(fstate, n)

   ! Recover state
   cstate = random_get_seed(fstate)
end function spral_random_long

logical(C_BOOL) function spral_random_logical(cstate) bind(C)
   use iso_c_binding
   use spral_random
   implicit none

   integer(C_INT), intent(inout) :: cstate

   type(random_state) :: fstate

   ! Initialize state
   call random_set_seed(fstate, cstate)

   ! Call Fortran routine
   spral_random_logical = random_logical(fstate)

   ! Recover state
   cstate = random_get_seed(fstate)
end function spral_random_logical
