real(C_DOUBLE) function spral_random_real(cstate, cpositive) bind(C)
   integer(C_INT), intent(inout) :: cstate
   logical(C_BOOL), value :: cpositive

   type(random_state) :: fstate

   call random_set_seed(fstate, cstate)
   spral_random_real = random_real(fstate, positive=cpositive)
   cstate = random_get_seed(fstate)
end function spral_random_real
