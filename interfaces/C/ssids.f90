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
     logical(C_BOOL) :: ignore_numa
     logical(C_BOOL) :: use_gpu
     integer(C_INT64_T) :: min_gpu_work
     real(C_FLOAT) :: max_load_inbalance
     real(C_FLOAT) :: gpu_perf_coeff
     integer(C_INT) :: scaling
     integer(C_INT64_T) :: small_subtree_threshold
     integer(C_INT) :: cpu_block_size
     logical(C_BOOL) :: action
     integer(C_INT) :: pivot_method
     real(C_DOUBLE) :: small
     real(C_DOUBLE) :: u
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
     integer(C_INT64_T) :: num_factor
     integer(C_INT64_T) :: num_flops
     integer(C_INT) :: num_neg
     integer(C_INT) :: num_sup
     integer(C_INT) :: num_two
     integer(C_INT) :: stat
     integer(C_INT) :: cuda_error
     integer(C_INT) :: cublas_error
     integer(C_INT) :: maxsupernode
     character(C_CHAR) :: unused(76)
  end type spral_ssids_inform

contains
  subroutine copy_options_in(coptions, foptions, cindexed)
    implicit none
    type(spral_ssids_options), intent(in) :: coptions
    type(ssids_options), intent(inout) :: foptions ! inherit some defaults!
    logical, intent(out) :: cindexed

    cindexed                   = (coptions%array_base .eq. 0)
    foptions%print_level       = coptions%print_level
    foptions%unit_diagnostics  = coptions%unit_diagnostics
    foptions%unit_error        = coptions%unit_error
    foptions%unit_warning      = coptions%unit_warning
    foptions%ordering          = coptions%ordering
    foptions%nemin             = coptions%nemin
    foptions%ignore_numa       = coptions%ignore_numa
    foptions%use_gpu           = coptions%use_gpu
    foptions%min_gpu_work      = coptions%min_gpu_work
    foptions%max_load_inbalance= coptions%max_load_inbalance
    foptions%gpu_perf_coeff    = coptions%gpu_perf_coeff
    foptions%scaling           = coptions%scaling
    foptions%small_subtree_threshold = coptions%small_subtree_threshold
    foptions%cpu_block_size    = coptions%cpu_block_size
    foptions%action            = coptions%action
    foptions%pivot_method      = coptions%pivot_method
    foptions%small             = coptions%small
    foptions%u                 = coptions%u
  end subroutine copy_options_in

  subroutine copy_inform_out(finform, cinform)
    implicit none
    type(ssids_inform), intent(in) :: finform
    type(spral_ssids_inform), intent(out) :: cinform

    cinform%flag                  = finform%flag
    cinform%matrix_dup            = finform%matrix_dup
    cinform%matrix_missing_diag   = finform%matrix_missing_diag
    cinform%matrix_outrange       = finform%matrix_outrange
    cinform%matrix_rank           = finform%matrix_rank
    cinform%maxdepth              = finform%maxdepth
    cinform%maxfront              = finform%maxfront
    cinform%maxsupernode          = finform%maxsupernode
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

subroutine spral_ssids_default_options(coptions) bind(C)
  use spral_ssids_ciface
  implicit none

  type(spral_ssids_options), intent(out) :: coptions

  type(ssids_options) :: default_options

  coptions%array_base        = 0 ! C
  coptions%print_level       = default_options%print_level
  coptions%unit_diagnostics  = default_options%unit_diagnostics
  coptions%unit_error        = default_options%unit_error
  coptions%unit_warning      = default_options%unit_warning
  coptions%ordering          = default_options%ordering
  coptions%nemin             = default_options%nemin
  coptions%ignore_numa       = default_options%ignore_numa
  coptions%use_gpu           = default_options%use_gpu
  coptions%min_gpu_work      = default_options%min_gpu_work
  coptions%max_load_inbalance= default_options%max_load_inbalance
  coptions%gpu_perf_coeff    = default_options%gpu_perf_coeff
  coptions%scaling           = default_options%scaling
  coptions%small_subtree_threshold = default_options%small_subtree_threshold
  coptions%cpu_block_size    = default_options%cpu_block_size
  coptions%action            = default_options%action
  coptions%pivot_method      = default_options%pivot_method
  coptions%small             = default_options%small
  coptions%u                 = default_options%u
end subroutine spral_ssids_default_options

subroutine spral_ssids_analyse(ccheck, n, corder, cptr, crow, cval, cakeep, &
     coptions, cinform) bind(C)
  use spral_ssids_ciface
  implicit none

  logical(C_BOOL), value :: ccheck
  integer(C_INT), value :: n
  type(C_PTR), value :: corder
  integer(C_INT64_T), target, dimension(n+1) :: cptr
  type(C_PTR), value :: cval
  type(C_PTR), intent(inout) :: cakeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform
  integer(C_INT), target, dimension(cptr(n+1)-coptions%array_base) :: crow

  integer(C_INT64_T), dimension(:), pointer :: fptr
  integer(C_INT64_T), dimension(:), allocatable, target :: fptr_alloc
  integer(C_INT), dimension(:), pointer :: frow
  integer(C_INT), dimension(:), allocatable, target :: frow_alloc
  logical :: fcheck
  integer(C_INT), dimension(:), pointer :: forder
  integer(C_INT), dimension(:), allocatable, target :: forder_alloc
  real(C_DOUBLE), dimension(:), pointer :: fval
  type(ssids_akeep), pointer :: fakeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  fcheck = ccheck
  if (C_ASSOCIATED(corder)) then
     call C_F_POINTER(corder, forder, shape=(/ n /))
  else
     nullify(forder)
  end if
  if (ASSOCIATED(forder) .and. cindexed) then
     allocate(forder_alloc(n))
     forder_alloc(:) = forder(:) + 1
     forder => forder_alloc
  endif
  fptr => cptr
  if (cindexed) then
     allocate(fptr_alloc(n+1))
     fptr_alloc(:) = fptr(:) + 1
     fptr => fptr_alloc
  end if
  frow => crow
  if (cindexed) then
     allocate(frow_alloc(fptr(n+1)-1))
     frow_alloc(:) = frow(:) + 1
     frow => frow_alloc
  end if
  if (C_ASSOCIATED(cval)) then
     call C_F_POINTER(cval, fval, shape=(/ fptr(n+1)-1 /))
  else
     nullify(fval)
  end if
  if (C_ASSOCIATED(cakeep)) then
     ! Reuse old pointer
     call C_F_POINTER(cakeep, fakeep)
  else
     ! Create new pointer
     allocate(fakeep)
     cakeep = C_LOC(fakeep)
  end if

  ! Call Fortran routine
  if (ASSOCIATED(forder)) then
     if (ASSOCIATED(fval)) then
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
             order=forder, val=fval)
     else
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
             order=forder)
     end if
  else
     if (ASSOCIATED(fval)) then
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
             val=fval)
     else
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform)
     end if
  end if

  ! Copy arguments out
  if (ASSOCIATED(forder) .and. cindexed) then
     call C_F_POINTER(corder, forder, shape = (/ n /) )
     forder(:) = forder_alloc(:) - 1
  endif
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_analyse

subroutine spral_ssids_analyse_ptr32(ccheck, n, corder, cptr, crow, cval, &
     cakeep, coptions, cinform) bind(C)
  use spral_ssids_ciface
  implicit none

  logical(C_BOOL), value :: ccheck
  integer(C_INT), value :: n
  type(C_PTR), value :: corder
  integer(C_INT), target, dimension(n+1) :: cptr
  type(C_PTR), value :: cval
  type(C_PTR), intent(inout) :: cakeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform
  integer(C_INT), target, dimension(cptr(n+1)-coptions%array_base) :: crow

  integer(C_INT), dimension(:), pointer :: fptr
  integer(C_INT), dimension(:), allocatable, target :: fptr_alloc
  integer(C_INT), dimension(:), pointer :: frow
  integer(C_INT), dimension(:), allocatable, target :: frow_alloc
  logical :: fcheck
  integer(C_INT), dimension(:), pointer :: forder
  integer(C_INT), dimension(:), allocatable, target :: forder_alloc
  real(C_DOUBLE), dimension(:), pointer :: fval
  type(ssids_akeep), pointer :: fakeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  fcheck = ccheck
  if (C_ASSOCIATED(corder)) then
     call C_F_POINTER(corder, forder, shape=(/ n /))
  else
     nullify(forder)
  end if
  if (ASSOCIATED(forder) .and. cindexed) then
     allocate(forder_alloc(n))
     forder_alloc(:) = forder(:) + 1
     forder => forder_alloc
  endif
  fptr => cptr
  if (cindexed) then
     allocate(fptr_alloc(n+1))
     fptr_alloc(:) = fptr(:) + 1
     fptr => fptr_alloc
  end if
  frow => crow
  if (cindexed) then
     allocate(frow_alloc(fptr(n+1)-1))
     frow_alloc(:) = frow(:) + 1
     frow => frow_alloc
  end if
  if (C_ASSOCIATED(cval)) then
     call C_F_POINTER(cval, fval, shape=(/ fptr(n+1)-1 /))
  else
     nullify(fval)
  end if
  if (C_ASSOCIATED(cakeep)) then
     ! Reuse old pointer
     call C_F_POINTER(cakeep, fakeep)
  else
     ! Create new pointer
     allocate(fakeep)
     cakeep = C_LOC(fakeep)
  end if

  ! Call Fortran routine
  if (ASSOCIATED(forder)) then
     if (ASSOCIATED(fval)) then
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
             order=forder, val=fval)
     else
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
             order=forder)
     end if
  else
     if (ASSOCIATED(fval)) then
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform, &
             val=fval)
     else
        call ssids_analyse(fcheck, n, fptr, frow, fakeep, foptions, finform)
     end if
  end if

  ! Copy arguments out
  if (ASSOCIATED(forder) .and. cindexed) then
     call C_F_POINTER(corder, forder, shape = (/ n /) )
     forder(:) = forder_alloc(:) - 1
  endif
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_analyse_ptr32

subroutine spral_ssids_analyse_coord(n, corder, ne, crow, ccol, cval, cakeep, &
     coptions, cinform) bind(C)
  use spral_ssids_ciface
  implicit none

  integer(C_INT), value :: n
  type(C_PTR), value :: corder
  integer(C_INT64_T), value :: ne
  integer(C_INT), target, dimension(ne) :: crow
  integer(C_INT), target, dimension(ne) :: ccol
  type(C_PTR), value :: cval
  type(C_PTR), intent(inout) :: cakeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform

  integer(C_INT), dimension(:), pointer :: frow
  integer(C_INT), dimension(:), allocatable, target :: frow_alloc
  integer(C_INT), dimension(:), pointer :: fcol
  integer(C_INT), dimension(:), allocatable, target :: fcol_alloc
  integer(C_INT), dimension(:), pointer :: forder
  integer(C_INT), dimension(:), allocatable, target :: forder_alloc
  real(C_DOUBLE), dimension(:), pointer :: fval
  type(ssids_akeep), pointer :: fakeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (C_ASSOCIATED(corder)) then
     call C_F_POINTER(corder, forder, shape=(/ n /))
  else
     nullify(forder)
  end if
  if (ASSOCIATED(forder) .and. cindexed) then
     allocate(forder_alloc(n))
     forder_alloc(:) = forder(:) + 1
     forder => forder_alloc
  endif
  frow => crow
  if (cindexed) then
     allocate(frow_alloc(ne))
     frow_alloc(:) = frow(:) + 1
     frow => frow_alloc
  end if
  fcol => ccol
  if (cindexed) then
     allocate(fcol_alloc(ne))
     fcol_alloc(:) = fcol(:) + 1
     fcol => fcol_alloc
  end if
  if (C_ASSOCIATED(cval)) then
     call C_F_POINTER(cval, fval, shape=(/ ne /))
  else
     nullify(fval)
  end if
  if (C_ASSOCIATED(cakeep)) then
     ! Reuse old pointer
     call C_F_POINTER(cakeep, fakeep)
  else
     ! Create new pointer
     allocate(fakeep)
     cakeep = C_LOC(fakeep)
  end if

  ! Call Fortran routine
  if (ASSOCIATED(forder)) then
     if (ASSOCIATED(fval)) then
        call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform,&
             order=forder, val=fval)
     else
        call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform,&
             order=forder)
     end if
  else
     if (ASSOCIATED(fval)) then
        call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform,&
             val=fval)
     else
        call ssids_analyse_coord(n, ne, frow, fcol, fakeep, foptions, finform)
     end if
  end if

  ! Copy arguments out
  if (ASSOCIATED(forder) .and. cindexed) then
     call C_F_POINTER(corder, forder, shape = (/ n /) )
     forder(:) = forder_alloc(:) - 1
  endif
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
  integer(C_INT64_T), dimension(:), pointer :: fptr
  integer(C_INT64_T), dimension(:), allocatable, target :: fptr_alloc
  integer(C_INT), dimension(:), pointer :: frow
  integer(C_INT), dimension(:), allocatable, target :: frow_alloc
  real(C_DOUBLE), dimension(:), pointer :: fscale
  type(ssids_akeep), pointer :: fakeep
  type(ssids_fkeep), pointer :: ffkeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  fposdef = cposdef
  call C_F_POINTER(cakeep, fakeep) ! Pulled forward so we can use it
  if (C_ASSOCIATED(cptr) .and. C_ASSOCIATED(crow)) then
     call C_F_POINTER(cptr, fptr, shape=(/ fakeep%n+1 /))
     if (cindexed) then
        allocate(fptr_alloc(fakeep%n+1))
        fptr_alloc(:) = fptr(:) + 1
        fptr => fptr_alloc
     end if
     call C_F_POINTER(crow, frow, shape=(/ fptr(fakeep%n+1)-1 /))
     if (cindexed) then
        allocate(frow_alloc(fakeep%n+1))
        frow_alloc(:) = frow(:) + 1
        frow => frow_alloc
     end if
  else
     nullify(fptr)
     nullify(frow)
  end if
  if (C_ASSOCIATED(cscale)) then
     call C_F_POINTER(cscale, fscale, shape=(/ fakeep%n /))
  else
     nullify(fscale)
  end if
  if (C_ASSOCIATED(cfkeep)) then
     ! Reuse old pointer
     call C_F_POINTER(cfkeep, ffkeep)
  else
     ! Create new pointer
     allocate(ffkeep)
     cfkeep = C_LOC(ffkeep)
  end if

  ! Call Fortran routine
  if (ASSOCIATED(fptr) .and. ASSOCIATED(frow)) then
     if (ASSOCIATED(fscale)) then
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
             ptr=fptr, row=frow, scale=fscale)
     else
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
             ptr=fptr, row=frow)
     end if
  else
     if (ASSOCIATED(fscale)) then
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
             scale=fscale)
     else
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform)
     end if
  end if

  ! Copy arguments out
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_factor

subroutine spral_ssids_factor_ptr32(cposdef, cptr, crow, val, cscale, cakeep, &
     cfkeep, coptions, cinform) bind(C)
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
  type(ssids_fkeep), pointer :: ffkeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  fposdef = cposdef
  call C_F_POINTER(cakeep, fakeep) ! Pulled forward so we can use it
  if (C_ASSOCIATED(cptr) .and. C_ASSOCIATED(crow)) then
     call C_F_POINTER(cptr, fptr, shape=(/ fakeep%n+1 /))
     if (cindexed) then
        allocate(fptr_alloc(fakeep%n+1))
        fptr_alloc(:) = fptr(:) + 1
        fptr => fptr_alloc
     end if
     call C_F_POINTER(crow, frow, shape=(/ fptr(fakeep%n+1)-1 /))
     if (cindexed) then
        allocate(frow_alloc(fakeep%n+1))
        frow_alloc(:) = frow(:) + 1
        frow => frow_alloc
     end if
  else
     nullify(fptr)
     nullify(frow)
  end if
  if (C_ASSOCIATED(cscale)) then
     call C_F_POINTER(cscale, fscale, shape=(/ fakeep%n /))
  else
     nullify(fscale)
  end if
  if (C_ASSOCIATED(cfkeep)) then
     ! Reuse old pointer
     call C_F_POINTER(cfkeep, ffkeep)
  else
     ! Create new pointer
     allocate(ffkeep)
     cfkeep = C_LOC(ffkeep)
  end if

  ! Call Fortran routine
  if (ASSOCIATED(fptr) .and. ASSOCIATED(frow)) then
     if (ASSOCIATED(fscale)) then
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
             ptr=fptr, row=frow, scale=fscale)
     else
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
             ptr=fptr, row=frow)
     end if
  else
     if (ASSOCIATED(fscale)) then
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform, &
             scale=fscale)
     else
        call ssids_factor(fposdef, val, fakeep, ffkeep, foptions, finform)
     end if
  end if

  ! Copy arguments out
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_factor_ptr32

subroutine spral_ssids_solve1(job, cx1, cakeep, cfkeep, coptions, cinform) &
     bind(C)
  use spral_ssids_ciface
  implicit none

  integer(C_INT), value :: job
  real(C_DOUBLE), target, dimension(*) :: cx1
  type(C_PTR), value :: cakeep
  type(C_PTR), value :: cfkeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform

  real(C_DOUBLE), dimension(:), pointer :: fx1
  type(ssids_akeep), pointer :: fakeep
  type(ssids_fkeep), pointer :: ffkeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (C_ASSOCIATED(cakeep)) then
     call C_F_POINTER(cakeep, fakeep)
  else
     nullify(fakeep)
  end if
  if (C_ASSOCIATED(cfkeep)) then
     call C_F_POINTER(cfkeep, ffkeep)
  else
     nullify(ffkeep)
  end if
  fx1 => cx1(1:fakeep%n)

  ! Call Fortran routine
  if (job .eq. 0) then
     ! Note: job=0 is an out of range value (but is valid internally!)
     call ssids_solve(fx1, fakeep, ffkeep, foptions, finform)
  else
     call ssids_solve(fx1, fakeep, ffkeep, foptions, finform, job=job)
  end if

  ! Copy arguments out
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_solve1

subroutine spral_ssids_solve(job, nrhs, x, ldx, cakeep, cfkeep, coptions, &
     cinform) bind(C)
  use spral_ssids_ciface
  implicit none

  integer(C_INT), value :: job
  integer(C_INT), value :: nrhs
  real(C_DOUBLE), dimension(ldx,nrhs) :: x
  integer(C_INT), value :: ldx
  type(C_PTR), value :: cakeep
  type(C_PTR), value :: cfkeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform

  type(ssids_akeep), pointer :: fakeep
  type(ssids_fkeep), pointer :: ffkeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (C_ASSOCIATED(cakeep)) then
     call C_F_POINTER(cakeep, fakeep)
  else
     nullify(fakeep)
  end if
  if (C_ASSOCIATED(cfkeep)) then
     call C_F_POINTER(cfkeep, ffkeep)
  else
     nullify(ffkeep)
  end if

  ! Call Fortran routine
  if (job .eq. 0) then
     call ssids_solve(nrhs, x, ldx, fakeep, ffkeep, foptions, finform)
  else
     call ssids_solve(nrhs, x, ldx, fakeep, ffkeep, foptions, finform, job=job)
  end if

  ! Copy arguments out
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_solve

integer(C_INT) function spral_ssids_free_akeep(cakeep) bind(C)
  use spral_ssids_ciface
  implicit none

  type(C_PTR), intent(inout) :: cakeep

  type(ssids_akeep), pointer :: fakeep

  if (.not. C_ASSOCIATED(cakeep)) then
     ! Nothing to free
     spral_ssids_free_akeep = 0
     return
  end if

  call C_F_POINTER(cakeep, fakeep)
  call ssids_free(fakeep, spral_ssids_free_akeep)
  deallocate(fakeep)
  cakeep = C_NULL_PTR
end function spral_ssids_free_akeep

integer(C_INT) function spral_ssids_free_fkeep(cfkeep) bind(C)
  use spral_ssids_ciface
  implicit none

  type(C_PTR), intent(inout) :: cfkeep

  type(ssids_fkeep), pointer :: ffkeep

  if (.not. C_ASSOCIATED(cfkeep)) then
     ! Nothing to free
     spral_ssids_free_fkeep = 0
     return
  end if

  call C_F_POINTER(cfkeep, ffkeep)
  call ssids_free(ffkeep, spral_ssids_free_fkeep)
  deallocate(ffkeep)
  cfkeep = C_NULL_PTR
end function spral_ssids_free_fkeep

integer(C_INT) function spral_ssids_free(cakeep, cfkeep) bind(C)
  use spral_ssids_ciface
  implicit none

  type(C_PTR), intent(inout) :: cakeep
  type(C_PTR), intent(inout) :: cfkeep

  interface
     integer(C_INT) function spral_ssids_free_akeep(cakeep) bind(C)
       use iso_c_binding
       implicit none
       type(C_PTR), intent(inout) :: cakeep
     end function spral_ssids_free_akeep
     integer(C_INT) function spral_ssids_free_fkeep(cfkeep) bind(C)
       use iso_c_binding
       implicit none
       type(C_PTR), intent(inout) :: cfkeep
     end function spral_ssids_free_fkeep
  end interface

  spral_ssids_free = spral_ssids_free_akeep(cakeep)
  if (spral_ssids_free .ne. 0) return
  spral_ssids_free = spral_ssids_free_fkeep(cfkeep)
end function spral_ssids_free

subroutine spral_ssids_enquire_posdef(cakeep, cfkeep, coptions, cinform, d) &
     bind(C)
  use spral_ssids_ciface
  implicit none

  type(C_PTR), value :: cakeep
  type(C_PTR), value :: cfkeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform
  real(C_DOUBLE), dimension(*), intent(out) :: d

  type(ssids_akeep), pointer :: fakeep
  type(ssids_fkeep), pointer :: ffkeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (C_ASSOCIATED(cakeep)) then
     call C_F_POINTER(cakeep, fakeep)
  else
     nullify(fakeep)
  end if
  if (C_ASSOCIATED(cfkeep)) then
     call C_F_POINTER(cfkeep, ffkeep)
  else
     nullify(ffkeep)
  end if

  ! Call Fortran routine
  call ssids_enquire_posdef(fakeep, ffkeep, foptions, finform, d)

  ! Copy arguments out
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_enquire_posdef

subroutine spral_ssids_enquire_indef(cakeep, cfkeep, coptions, cinform, &
     cpiv_order, cd) bind(C)
  use spral_ssids_ciface
  implicit none

  type(C_PTR), value :: cakeep
  type(C_PTR), value :: cfkeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform
  type(C_PTR), value :: cpiv_order
  type(C_PTR), value :: cd

  type(ssids_akeep), pointer :: fakeep
  type(ssids_fkeep), pointer :: ffkeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform
  integer(C_INT), dimension(:), pointer :: fpiv_order
  real(C_DOUBLE), dimension(:,:), pointer :: fd

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (C_ASSOCIATED(cakeep)) then
     call C_F_POINTER(cakeep, fakeep)
  else
     nullify(fakeep)
  end if
  if (C_ASSOCIATED(cfkeep)) then
     call C_F_POINTER(cfkeep, ffkeep)
  else
     nullify(ffkeep)
  end if
  if (C_ASSOCIATED(cpiv_order)) then
     call C_F_POINTER(cpiv_order, fpiv_order, shape=(/ fakeep%n /))
  else
     nullify(fpiv_order)
  end if
  if (C_ASSOCIATED(cd)) then
     call C_F_POINTER(cd, fd, shape=(/ 2,fakeep%n /))
  else
     nullify(fd)
  end if

  ! Call Fortran routine
  if (ASSOCIATED(fpiv_order)) then
     if (ASSOCIATED(fd)) then
        call ssids_enquire_indef(fakeep, ffkeep, foptions, finform, &
             piv_order=fpiv_order, d=fd)
     else
        call ssids_enquire_indef(fakeep, ffkeep, foptions, finform, &
             piv_order=fpiv_order)
     end if
  else
     if (ASSOCIATED(fd)) then
        call ssids_enquire_indef(fakeep, ffkeep, foptions, finform, d=fd)
     else
        call ssids_enquire_indef(fakeep, ffkeep, foptions, finform)
     end if
  end if

  ! Copy arguments out
  ! Note: we use abs value of piv_order in C indexing, as 0 and -0 are the same
  if (ASSOCIATED(fpiv_order) .and. cindexed) &
       fpiv_order(:) = abs(fpiv_order(:)) - 1
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_enquire_indef

subroutine spral_ssids_alter(d, cakeep, cfkeep, coptions, cinform) bind(C)
  use spral_ssids_ciface
  implicit none

  real(C_DOUBLE), dimension(2,*), intent(in) :: d
  type(C_PTR), value :: cakeep
  type(C_PTR), value :: cfkeep
  type(spral_ssids_options), intent(in) :: coptions
  type(spral_ssids_inform), intent(out) :: cinform

  type(ssids_akeep), pointer :: fakeep
  type(ssids_fkeep), pointer :: ffkeep
  type(ssids_options) :: foptions
  type(ssids_inform) :: finform

  logical :: cindexed

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (C_ASSOCIATED(cakeep)) then
     call C_F_POINTER(cakeep, fakeep)
  else
     nullify(fakeep)
  end if
  if (C_ASSOCIATED(cfkeep)) then
     call C_F_POINTER(cfkeep, ffkeep)
  else
     nullify(ffkeep)
  end if

  ! Call Fortran routine
  call ssids_alter(d, fakeep, ffkeep, foptions, finform)

  ! Copy arguments out
  call copy_inform_out(finform, cinform)
end subroutine spral_ssids_alter
