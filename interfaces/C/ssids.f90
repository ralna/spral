module spral_ssids_ciface
   use iso_c_binding
   use spral_ssids
   implicit none

   type, bind(C) :: spral_ssids_options
      integer(C_INT) :: array_base
      integer(C_INT) :: print_level
      integer(C_INT) :: unit_diagnostics
      integer(C_INT) :: unit_error
      integer(C_INT) :: unit_warning
      integer(C_INT) :: ordering
      integer(C_INT) :: nemin
      integer(C_INT) :: scaling
      logical(C_BOOL) :: action
      real(C_DOUBLE) :: small
      real(C_DOUBLE) :: u
      logical(C_BOOL) :: use_gpu_solve
      integer(C_INT) :: presolve
      character(C_CHAR) :: unused(80)
   end type spral_ssids_options

   type, bind(C) :: spral_ssids_inform
      integer(C_INT) :: flag
      integer(C_INT) :: matrix_dup
      integer(C_INT) :: matrix_missing_diag
      integer(C_INT) :: matrix_outrange
      integer(C_INT) :: matrix_rank
      integer(C_INT) :: maxdepth
      integer(C_INT) :: maxfront
      integer(C_INT) :: num_delay
      integer(C_LONG) :: num_factor
      integer(C_LONG) :: num_flops
      integer(C_INT) :: num_neg
      integer(C_INT) :: num_sup
      integer(C_INT) :: num_two
      integer(C_INT) :: stat
      integer(C_INT) :: cuda_error
      integer(C_INT) :: cublas_error
      character(C_CHAR) :: unused(80)
   end type spral_ssids_inform

contains
   subroutine copy_options_in(coptions, foptions, cindexed)
      type(spral_ssids_options), intent(in) :: coptions
      type(ssids_options), intent(inout) :: foptions ! inherit some defaults!
      logical, intent(out) :: cindexed

      cindexed                   = (coptions%array_base.eq.0)
      foptions%print_level       = coptions%print_level
      foptions%unit_diagnostics  = coptions%unit_diagnostics
      foptions%unit_error        = coptions%unit_error
      foptions%unit_warning      = coptions%unit_warning
      foptions%ordering          = coptions%ordering
      foptions%nemin             = coptions%nemin
      foptions%scaling           = coptions%scaling
      foptions%action            = coptions%action
      foptions%small             = coptions%small
      foptions%u                 = coptions%u
      foptions%use_gpu_solve     = coptions%use_gpu_solve
      foptions%presolve          = coptions%presolve
   end subroutine copy_options_in

   subroutine copy_inform_out(finform, cinform)
      type(ssids_inform), intent(in) :: finform
      type(spral_ssids_inform), intent(out) :: cinform

      cinform%flag                  = finform%flag
      cinform%matrix_dup            = finform%matrix_dup
      cinform%matrix_missing_diag   = finform%matrix_missing_diag
      cinform%matrix_outrange       = finform%matrix_outrange
      cinform%matrix_rank           = finform%matrix_rank
      cinform%maxdepth              = finform%maxdepth
      cinform%maxfront              = finform%maxfront
      cinform%num_delay             = finform%num_delay
      cinform%num_factor            = finform%num_factor
      cinform%num_flops             = finform%num_flops
      cinform%num_neg               = finform%num_neg
      cinform%num_sup               = finform%num_sup
      cinform%num_two               = finform%num_two
      cinform%stat                  = finform%stat
      cinform%cuda_error            = finform%cuda_error
      cinform%cublas_error          = finform%cublas_error
   end subroutine copy_inform_out
end module spral_ssids_ciface

subroutine spral_ssids_default_options(coptions)
   type(spral_ssids_options), intent(out) :: coptions

   type(ssids_options) :: default_options

   coptions%array_base        = 0 ! C
   coptions%print_level       = default_options%print_level
   coptions%unit_diagnostics  = default_options%unit_diagnostics
   coptions%unit_error        = default_options%unit_error
   coptions%unit_warning      = default_options%unit_warning
   coptions%ordering          = default_options%ordering
   coptions%nemin             = default_options%nemin
   coptions%scaling           = default_options%scaling
   coptions%action            = default_options%action
   coptions%small             = default_options%small
   coptions%u                 = default_options%u
   coptions%use_gpu_solve     = default_options%use_gpu_solve
   coptions%presolve          = default_options%presolve
end subroutine spral_ssids_default_options

subroutine spral_ssids_analyse(ccheck, n, corder, cptr, crow, cval, cakeep, &
      coptions, cinform) bind(C)
   use spral_ssids_ciface
   implicit none

   logical(C_BOOL), value :: ccheck
   integer(C_INT), value :: n
   type(C_PTR), value :: corder
   type(C_PTR), value :: cptr
   type(C_PTR), value :: crow
   type(C_PTR), value :: cval
   type(C_PTR), intent(inout) :: cakeep
   type(spral_ssids_options), intent(in) :: coptions
   type(spral_ssids_inform), intent(in) :: cinform

   integer(C_INT), dimension(:), pointer :: fptr
   integer(C_INT), dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer(C_INT), dimension(:), allocatable, target :: frow_alloc
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
         call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
            order=forder, val=fval)
      else
         call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
            order=forder)
      endif
   else
      if(ASSOCIATED(fval)) then
         call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
            val=fval)
      else
         call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform)
      endif
   endif

   ! Copy arguments out
   if(ASSOCIATED(forder) .and. cindexed) forder(:) = forder(:) - 1
   call copy_inform_out(finform, cinform)
end subroutine spral_ssids_analyse

subroutine spral_ssids_analyse_coord(n, corder, ne, crow, ccol, cval, cakeep, &
      coptions, cinform) bind(C)
   use spral_ssids_ciface
   implicit none

   integer(C_INT), value :: n
   type(C_PTR), value :: corder
   integer(C_INT), value :: ne
   type(C_PTR), value :: crow
   type(C_PTR), value :: ccol
   type(C_PTR), value :: cval
   type(C_PTR), intent(inout) :: cakeep
   type(spral_ssids_options), intent(in) :: coptions
   type(spral_ssids_inform), intent(in) :: cinform

   integer(C_INT), dimension(:), pointer :: frow
   integer(C_INT), dimension(:), allocatable, target :: frow_alloc
   integer(C_INT), dimension(:), pointer :: fcol
   integer(C_INT), dimension(:), allocatable, target :: fcol_alloc
   integer(C_INT), dimension(:), pointer :: forder
   real(C_DOUBLE), dimension(:), pointer :: fval
   type(ssids_akeep), pointer :: fakeep
   type(ssids_options) :: foptions
   type(ssids_inform) :: finform

   logical :: cindexed

   ! Copy options in first to find out whether we use Fortran or C indexing
   call copy_options_in(coptions, foptions, cindexed)

   ! Translate arguments
   if(C_ASSOCIATED(corder)) then
      call C_F_POINTER(corder, forder, shape=(/ n /))
   else
      nullify(forder)
   endif
   call C_F_POINTER(crow, frow, shape=(/ ne /))
   if(.not.cindexed) then
      allocate(frow_alloc(n+1))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   call C_F_POINTER(ccol, fcol, shape=(/ ne /))
   if(.not.cindexed) then
      allocate(fcol_alloc(n+1))
      fcol_alloc(:) = fcol(:) + 1
      fcol => fcol_alloc
   endif
   if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape=(/ ne /))
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
         call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform,&
            order=forder, val=fval)
      else
         call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform,&
            order=forder)
      endif
   else
      if(ASSOCIATED(fval)) then
         call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform,&
            val=fval)
      else
         call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform)
      endif
   endif

   ! Copy arguments out
   if(ASSOCIATED(forder) .and. cindexed) forder(:) = forder(:) - 1
   call copy_inform_out(finform, cinform)
end subroutine spral_ssids_analyse_coord

subroutine spral_ssids_factor(cposdef, cptr, crow, val, cscale, cakeep, cfkeep,&
      coptions, cinform) bind(C)
   use spral_ssids_ciface
   implicit none

   logical(C_BOOL), value :: cposdef
   type(C_PTR), value :: cptr
   type(C_PTR), value :: crow
   real(C_DOUBLE), dimension(*), intent(in) :: val
   type(C_PTR), value :: cscale
   type(C_PTR), value :: cakeep
   type(C_PTR), intent(inout) :: cfkeep
   type(spral_ssids_options), intent(in) :: coptions
   type(spral_ssids_inform), intent(out) :: cinform

   logical :: fposdef
   integer(C_INT), dimension(:), pointer :: fptr
   integer(C_INT), dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer(C_INT), dimension(:), allocatable, target :: frow_alloc
   real(C_DOUBLE), dimension(:), pointer :: fscale
   type(ssids_akeep), pointer :: fakeep
   type(ssids_akeep), pointer :: ffkeep
   type(ssids_options) :: foptions
   type(ssids_inform) :: finform

   logical :: cindexed

   ! Copy options in first to find out whether we use Fortran or C indexing
   call copy_options_in(coptions, foptions, cindexed)

   ! Translate arguments
   fposdef = cposdef
   if(C_ASSOCIATED(cptr) .and. C_ASSOCIATED(crow)) then
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
   else
      nullify(fptr)
      nullify(frow)
   endif
   if(C_ASSOCIATED(cscale)) then
      call C_F_POINTER(cscale, fscale, shape=(/ n /))
   else
      nullify(fscale)
   endif
   call C_F_POINTER(cakeep, fakeep)
   if(C_ASSOCIATED(cfkeep)) then
      ! Reuse old pointer
      call C_F_POINTER(cfkeep, ffkeep)
   else
      ! Create new pointer
      allocate(ffkeep)
      cfkeep = C_LOC(ffkeep)
   endif

   ! Call Fortran routine
   if(ASSOCIATED(fptr) .and. ASSOCIATED(frow)) then
      if(ASSOCIATED(fscale)) then
         call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
            ptr=fptr, row=frow, scale=fscale)
      else
         call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
            ptr=fptr, row=frow)
      endif
   else
      if(ASSOCIATED(fscale)) then
         call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
            scale=fscale)
      else
         call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform)
      endif
   endif

   ! Copy arguments out
   call copy_inform_out(finform, cinform)
end subroutine spral_ssids_factor

subroutine spral_ssids_solve1(job, x1, cakeep, cfkeep, coptions, cinform) &
      bind(C)
   use spral_ssids_ciface
   implicit none

   integer(C_INT), value :: job
   real(C_DOUBLE), dimension(*), intent(in) :: x1
   type(C_PTR), value :: cakeep
   type(C_PTR), value :: cfkeep
   type(spral_ssids_options), intent(in) :: coptions
   type(spral_ssids_inform), intent(out) :: cinform

   type(ssids_akeep), pointer :: fakeep
   type(ssids_akeep), pointer :: ffkeep
   type(ssids_options) :: foptions
   type(ssids_inform) :: finform

   ! Copy options in first to find out whether we use Fortran or C indexing
   call copy_options_in(coptions, foptions, cindexed)

   ! Translate arguments
   call C_F_POINTER(cakeep, fakeep)
   call C_F_POINTER(cfkeep, ffkeep)

   ! Call Fortran routine
   call ssids_solve(x1, fakeep, ffkeep, foptions, finform, job=job)

   ! Copy arguments out
   call copy_inform_out(finform, cinform)
end subroutine spral_ssids_solve1

subroutine spral_ssids_solve(job, nrhs, x, ldx, cakeep, cfkeep, coptions, &
      cinform) bind(C)
   use spral_ssids_ciface
   implicit none

   integer(C_INT), value :: job
   integer(C_INT), value :: nrhs
   real(C_DOUBLE), dimension(*), intent(in) :: x
   integer(C_INT), value :: ldx
   type(C_PTR), value :: cakeep
   type(C_PTR), value :: cfkeep
   type(spral_ssids_options), intent(in) :: coptions
   type(spral_ssids_inform), intent(out) :: cinform

   type(ssids_akeep), pointer :: fakeep
   type(ssids_akeep), pointer :: ffkeep
   type(ssids_options) :: foptions
   type(ssids_inform) :: finform

   ! Copy options in first to find out whether we use Fortran or C indexing
   call copy_options_in(coptions, foptions, cindexed)

   ! Translate arguments
   call C_F_POINTER(cakeep, fakeep)
   call C_F_POINTER(cfkeep, ffkeep)

   ! Call Fortran routine
   call ssids_solve(nrhs, x, ldx, fakeep, ffkeep, foptions, finform, job=job)

   ! Copy arguments out
   call copy_inform_out(finform, cinform)
end subroutine spral_ssids_solve

integer(C_INT) function spral_ssids_free_akeep(cakeep) bind(C)
   use spral_ssids_ciface
   implicit none
   
   type(C_PTR), intent(inout) :: cakeep

   type(ssids_akeep), pointer :: fakeep

   if(.not.C_ASSOCIATED(cakeep)) then
      ! Nothing to free
      spral_ssids_free_akeep = 0
      return
   endif

   call C_F_POINTER(cakeep)
   call ssids_free(cakeep, spral_ssids_free_akeep)
   deallocate(cakeep)
   cakeep = C_NULL_PTR
end function spral_ssids_free_akeep

integer(C_INT) function spral_ssids_free_fkeep(cfkeep) bind(C)
   use spral_ssids_ciface
   implicit none
   
   type(C_PTR), intent(inout) :: cfkeep

   type(ssids_fkeep), pointer :: ffkeep

   if(.not.C_ASSOCIATED(cfkeep)) then
      ! Nothing to free
      spral_ssids_free_fkeep = 0
      return
   endif

   call C_F_POINTER(cfkeep)
   call ssids_free(cfkeep, spral_ssids_free_fkeep)
   deallocate(cfkeep)
   cfkeep = C_NULL_PTR
end function spral_ssids_free_fkeep

integer(C_INT) function spral_ssids_free(cakeep, cfkeep) bind(C)
   use spral_ssids_ciface
   implicit none
   
   type(C_PTR), intent(inout) :: cakeep
   type(C_PTR), intent(inout) :: cfkeep

   spral_ssids_free = spral_ssids_free_akeep(cakeep)
   if(spral_ssids_free.ne.0) return
   spral_ssids_free = spral_ssids_free_fkeep(cfkeep)
end function spral_ssids_free_fkeep

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
