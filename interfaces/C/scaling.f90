module spral_scaling_ciface
   use, intrinsic :: iso_c_binding
   use spral_scaling
   implicit none

   type, bind(C) :: spral_scaling_hungarian_options
      integer(C_INT) :: array_base
      logical(C_BOOL) :: scale_if_singular
   end type spral_scaling_hungarian_options

   type, bind(C) :: spral_scaling_hungarian_inform
      integer(C_INT) :: flag
      integer(C_INT) :: stat
   end type spral_scaling_hungarian_inform

contains
   subroutine copy_hungarian_options_in(coptions, foptions, cindexed)
      type(spral_scaling_hungarian_options), intent(in) :: coptions
      type(hungarian_options), intent(out) :: foptions
      logical, intent(out) :: cindexed

      cindexed                   = (coptions%array_base.eq.0)
      foptions%scale_if_singular = coptions%scale_if_singular
   end subroutine copy_hungarian_options_in

   subroutine copy_hungarian_inform_out(finform, cinform)
      type(hungarian_inform), intent(in) :: finform
      type(spral_scaling_hungarian_inform), intent(out) :: cinform

      cinform%flag = finform%flag
      cinform%stat = finform%stat
   end subroutine copy_hungarian_inform_out
end module spral_scaling_ciface

subroutine spral_scaling_hungarian_default_options(coptions) bind(C)
   use spral_scaling_ciface
   implicit none

   type(spral_scaling_hungarian_options), intent(out) :: coptions

   type(hungarian_options) :: default_options

   coptions%array_base        = 0 ! C
   coptions%scale_if_singular = default_options%scale_if_singular
end subroutine spral_scaling_hungarian_default_options

subroutine spral_scaling_hungarian_sym(n, ptr, row, val, cmatch, scaling, &
      coptions, cinform) bind(C)
   use spral_scaling_ciface
   implicit none

   integer(C_INT), value :: n
   integer(C_INT), dimension(*), intent(in) :: ptr
   integer(C_INT), dimension(*), intent(in) :: row
   real(C_DOUBLE), dimension(*), intent(in) :: val
   type(C_PTR), value :: cmatch
   real(C_DOUBLE), dimension(*), intent(out) :: scaling
   type(spral_scaling_hungarian_options), intent(in) :: coptions
   type(spral_scaling_hungarian_inform), intent(out) :: cinform

   logical :: cindexed
   integer, dimension(:), allocatable :: ptr_alloc, row_alloc
   integer(C_INT), dimension(:), pointer :: fmatch
   type(hungarian_options) :: foptions
   type(hungarian_inform) :: finform

   ! Associate arrays and copy control in
   call copy_hungarian_options_in(coptions, foptions, cindexed)
   if(cindexed) then
      allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
      ptr_alloc(1:n+1) = ptr(1:n+1) + 1
      row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
   endif
   if(C_ASSOCIATED(cmatch)) then
      call C_F_POINTER(cmatch, fmatch, shape=(/ n /))
   else
      nullify(fmatch)
   endif

   ! Call Fortran routine
   if(cindexed) then
      if(associated(fmatch)) then
         call hungarian_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
            foptions, finform, match=fmatch)
      else
         call hungarian_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
            foptions, finform)
      endif
   else
      if(associated(fmatch)) then
         call hungarian_scale_sym(n, ptr, row, val, scaling, foptions, &
            finform, match=fmatch)
      else
         call hungarian_scale_sym(n, ptr, row, val, scaling, foptions, finform)
      endif
   endif

   ! Copy info out
   call copy_hungarian_inform_out(finform, cinform)
   if(cindexed .and. associated(fmatch)) &
      fmatch(:) = fmatch(:) - 1

end subroutine spral_scaling_hungarian_sym
