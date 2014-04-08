module spral_ssids_ciface
   use iso_c_binding
   implicit none

   type, bind(C) :: spral_ssids_options
   end type spral_ssids_analyse

   type, bind(C) :: spral_ssids_inform
   end type spral_ssids_inform

contains
   subroutine copy_options_in(coptions, foptions, cindexed)
      type(spral_ssids_options), intent(in) :: coptions
      type(ssids_options), intent(inout) :: foptions
      logical, intent(out) :: cindexed
   end subroutine copy_options_in

   subroutine copy_inform_out(finform, cinform)
      type(ssids_inform), intent(inout) :: finform
      type(spral_ssids_inform), intent(in) :: cinform
   end subroutine copy_inform_out
end module spral_ssids_ciface

subroutine spral_ssids_analyse(ccheck, n, corder, ptr, row, cval, cakeep, &
      coptions, cinform) bind(C)
   use spral_ssids_ciface
   implicit none

   logical(C_BOOL), value :: ccheck
   integer(C_INT), value :: n
   type(C_PTR), value :: corder
   type(C_PTR), value :: ptr
   type(C_PTR), value :: row
   type(C_PTR), value :: cval
   type(C_PTR), intent(inout) :: cakeep
   type(spral_ssids_options), intent(in) :: coptions
   type(spral_ssids_inform), intent(in) :: cinform

   integer(C_INT), dimension(:), pointer :: fptr_alloc
   integer(C_INT), dimension(:), allocatable :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow_alloc
   integer(C_INT), dimension(:), allocatable :: frow_alloc
   logical :: fcheck
   integer(C_INT), dimension(:), pointer :: forder
   real(C_DOUBLE), dimension(:), pointer :: fval
   type(ssids_akeep), pointer :: fakeep
   type(ssids_options) :: foptions
   type(ssids_inform) :: finform

   logical :: cindexed

   ! Copy options in first to find out whether we use Fortran or C indexing
   call copy_options_in(coptions, foptions, cindexed)

   ! Translate arguments
   fcheck = ccheck
   if(C_ASSOCIATED(corder)) then
      call C_F_POINTER(corder, forder, shape=(/ n /))
   else
      nullify(forder)
   endif
   call C_F_POINTER(cptr, fptr, shape=(/ n+1 /))
   if(.not.cindexed) then
      allocate(fptr_alloc(n+1))
      fptr_alloc(:) = fptr(:) + 1
      fptr => fptr_alloc
   endif
   call C_F_POINTER(crow, frow, shape=(/ fptr(n+1)-1 /))
   if(.not.cindexed) then
      allocate(frow_alloc(n+1))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape=(/ fptr(n+1)-1 /))
   else
      nullify(fval)
   endif
   if(C_ASSOCIATED(cakeep)) then
      ! Reuse old pointer
      call C_F_POINTER(cakeep, fakeep)
   else
      ! Create new pointer
      allocate(fakeep)
      cakeep = C_LOC(fakeep)
   endif

   ! Call Fortran routine
   if(ASSOCIATED(forder)) then
      if(ASSOCIATED(fval)) then
         call ssids_analyse(fcheck, n, ptr, row, fakeep, foptions, finform, &
            order=forder, val=fval)
      else
         call ssids_analyse(fcheck, n, ptr, row, fakeep, foptions, finform, &
            order=forder)
      endif
   else
      if(ASSOCIATED(fval)) then
         call ssids_analyse(fcheck, n, ptr, row, fakeep, foptions, finform, &
            val=fval)
      else
         call ssids_analyse(fcheck, n, ptr, row, fakeep, foptions, finform)
      endif
   endif

   ! Copy arguments out
   if(ASSOCIATED(forder) .and. cindexed) forder(:) = forder(:) - 1
   call copy_inform_out(finform, cinform)
end subroutine spral_ssids_analyse

subroutine spral_ssids_analyse_coord() bind(C)
   use spral_ssids_ciface
   implicit none
end subroutine spral_ssids_analyse_coord

subroutine spral_ssids_factor() bind(C)
   use spral_ssids_ciface
   implicit none
end subroutine spral_ssids_factor

subroutine spral_ssids_solve() bind(C)
   use spral_ssids_ciface
   implicit none
end subroutine spral_ssids_solve

subroutine spral_ssids_free() bind(C)
   use spral_ssids_ciface
   implicit none
end subroutine spral_ssids_free

subroutine spral_ssids_enquire_posdef() bind(C)
   use spral_ssids_ciface
   implicit none
end subroutine spral_ssids_enquire_posdef

subroutine spral_ssids_enquire_indef() bind(C)
   use spral_ssids_ciface
   implicit none
end subroutine spral_ssids_enquire_indef

subroutine spral_ssids_alter() bind(C)
   use spral_ssids_ciface
   implicit none
end subroutine spral_ssids_alter
