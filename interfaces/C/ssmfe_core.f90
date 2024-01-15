module spral_ssmfe_core_ciface
  use iso_c_binding
  use spral_ssmfe_core
  implicit none

  type, bind(C) :: spral_ssmfe_rcid
     integer(C_INT) :: job
     integer(C_INT) :: nx
     integer(C_INT) :: jx
     integer(C_INT) :: kx
     integer(C_INT) :: ny
     integer(C_INT) :: jy
     integer(C_INT) :: ky
     integer(C_INT) :: i
     integer(C_INT) :: j
     integer(C_INT) :: k
     real(C_DOUBLE) :: alpha
     real(C_DOUBLE) :: beta
     type(C_PTR) :: x
     type(C_PTR) :: y
     character(C_CHAR) :: unused(80)
  end type spral_ssmfe_rcid

  type, bind(C) :: spral_ssmfe_rciz
     integer(C_INT) :: job
     integer(C_INT) :: nx
     integer(C_INT) :: jx
     integer(C_INT) :: kx
     integer(C_INT) :: ny
     integer(C_INT) :: jy
     integer(C_INT) :: ky
     integer(C_INT) :: i
     integer(C_INT) :: j
     integer(C_INT) :: k
     complex(C_DOUBLE_COMPLEX) :: alpha
     complex(C_DOUBLE_COMPLEX) :: beta
     type(C_PTR) :: x
     type(C_PTR) :: y
     character(C_CHAR) :: unused(80)
  end type spral_ssmfe_rciz

  type, bind(C) :: spral_ssmfe_core_options
     integer(C_INT) :: array_base
     real(C_DOUBLE) :: cf_max
     integer(C_INT) :: err_est
     integer(C_INT) :: extra_left
     integer(C_INT) :: extra_right
     real(C_DOUBLE) :: min_gap
     logical(C_BOOL) :: minAprod
     logical(C_BOOL) :: minBprod
     character(C_CHAR) :: unused(80)
  end type spral_ssmfe_core_options

  type, bind(C) :: spral_ssmfe_inform
     integer(C_INT) :: flag
     integer(C_INT) :: stat
     integer(C_INT) :: non_converged
     integer(C_INT) :: iteration
     integer(C_INT) :: left
     integer(C_INT) :: right
     type(C_PTR) :: converged
     real(C_DOUBLE) :: next_left
     real(C_DOUBLE) :: next_right
     type(C_PTR) :: residual_norms
     type(C_PTR) :: err_lambda
     type(C_PTR) :: err_X
     character(C_CHAR) :: unused(80)
  end type spral_ssmfe_inform

  type ssmfe_core_ciface_keep
     type(ssmfe_core_keep) :: keep
     type(ssmfe_rcid) :: rcid
     type(ssmfe_rciz) :: rciz
     type(ssmfe_inform) :: inform
  end type ssmfe_core_ciface_keep

  interface copy_rci_out
     module procedure copy_rci_out_double, copy_rci_out_double_complex
  end interface copy_rci_out

contains
  subroutine copy_core_options_in(coptions, foptions, cindexed)
    implicit none
    type(spral_ssmfe_core_options), intent(in) :: coptions
    type(ssmfe_core_options), intent(inout) :: foptions ! inherits defaults!
    logical, intent(out) :: cindexed

    cindexed             = (coptions%array_base.eq.0)
    foptions%cf_max      = coptions%cf_max
    foptions%err_est     = coptions%err_est
    foptions%extra_left  = coptions%extra_left
    foptions%extra_right = coptions%extra_right
    foptions%min_gap     = coptions%min_gap
    foptions%minAprod    = coptions%minAprod
    foptions%minBprod    = coptions%minBprod
  end subroutine copy_core_options_in

  subroutine copy_rci_out_double_complex(frci, crci, cindexed)
    implicit none
    type(ssmfe_rciz), target, intent(in) :: frci
    type(spral_ssmfe_rciz), intent(inout) :: crci
    logical, intent(in) :: cindexed

    integer :: coffset

    coffset = 0
    if (cindexed) coffset = -1

    crci%job = frci%job ! enum, not an array index
    crci%nx = frci%nx ! count, not an array index
    crci%jx = frci%jx + coffset
    crci%kx = frci%kx ! is an array index, but expected to start at 0 anyway
    crci%ny = frci%ny ! count, not an array index
    crci%jy = frci%jy + coffset
    crci%ky = frci%ky ! is an array index, but expected to start at 0 anyway
    select case (crci%job)
    case(5,11,999)
       ! i, j, k are NOT array indices
       crci%i = frci%i
       crci%j = frci%j
       crci%k = frci%k
    case default ! defintely 12, 14, 15, 16, 17
       ! i, j, k are array indices
       crci%i = frci%i + coffset
       crci%j = frci%j + coffset
       crci%k = frci%k + coffset
    end select
    crci%alpha = frci%alpha ! floating point, not an array index
    crci%beta = frci%beta ! floating point, not an array index
    if (associated(frci%x)) &
         crci%x = c_loc(frci%x(1,1))
    if (associated(frci%y)) &
         crci%y = c_loc(frci%y(1,1))
  end subroutine copy_rci_out_double_complex

  subroutine copy_rci_out_double(frci, crci, cindexed)
    implicit none
    type(ssmfe_rcid), target, intent(in) :: frci
    type(spral_ssmfe_rcid), intent(inout) :: crci
    logical, intent(in) :: cindexed

    integer :: coffset

    coffset = 0
    if (cindexed) coffset = -1

    crci%job = frci%job ! enum, not an array index
    crci%nx = frci%nx ! count, not an array index
    crci%jx = frci%jx + coffset
    crci%kx = frci%kx ! is an array index, but expected to start at 0 anyway
    crci%ny = frci%ny ! count, not an array index
    crci%jy = frci%jy + coffset
    crci%ky = frci%ky ! is an array index, but expected to start at 0 anyway
    select case (crci%job)
    case(5,11,999)
       ! i, j, k are NOT array indices
       crci%i = frci%i
       crci%j = frci%j
       crci%k = frci%k
    case default ! defintely 12, 14, 15, 16, 17
       ! i, j, k are array indices
       crci%i = frci%i + coffset
       crci%j = frci%j + coffset
       crci%k = frci%k + coffset
    end select
    crci%alpha = frci%alpha ! floating point, not an array index
    crci%beta = frci%beta ! floating point, not an array index
    if (associated(frci%x)) &
         crci%x = c_loc(frci%x(1,1))
    if (associated(frci%y)) &
         crci%y = c_loc(frci%y(1,1))
  end subroutine copy_rci_out_double

  ! NB: Note that as we take address of components of finform, finform must
  ! persist in memory - achieve this by keeping it as part of ckeep data
  ! structure
  subroutine copy_inform_out(finform, cinform)
    implicit none
    type(ssmfe_inform), target, intent(in) :: finform
    type(spral_ssmfe_inform), intent(out) :: cinform

    cinform%flag            = finform%flag
    cinform%stat            = finform%stat
    cinform%non_converged   = finform%non_converged
    cinform%iteration       = finform%iteration
    cinform%left            = finform%left
    cinform%right           = finform%right
    if (allocated(finform%converged)) &
         cinform%converged = c_loc(finform%converged(1))
    cinform%next_left       = finform%next_left
    cinform%next_right      = finform%next_right
    if (allocated(finform%residual_norms))  &
         cinform%residual_norms = c_loc(finform%residual_norms(1))
    if (allocated(finform%err_lambda)) &
         cinform%err_lambda = c_loc(finform%err_lambda(1))
    if (allocated(finform%err_x)) &
         cinform%err_x = c_loc(finform%err_x(1))
  end subroutine copy_inform_out
end module spral_ssmfe_core_ciface

subroutine spral_ssmfe_core_default_options(coptions) bind(C)
  use spral_ssmfe_core_ciface
  implicit none

  type(spral_ssmfe_core_options), intent(out) :: coptions

  type(ssmfe_core_options) :: default_options

  coptions%array_base  = 0 ! C
  coptions%cf_max      = default_options%cf_max
  coptions%err_est     = default_options%err_est
  coptions%extra_left  = default_options%extra_left
  coptions%extra_right = default_options%extra_right
  coptions%min_gap     = default_options%min_gap
  coptions%minAprod    = default_options%minAprod
  coptions%minBprod    = default_options%minBprod
end subroutine spral_ssmfe_core_default_options

subroutine spral_ssmfe_double(crci, problem, left, right, m, lambda, rr, &
     ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_core_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  integer(C_INT), value :: problem
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(m), intent(inout) :: lambda
  real(C_DOUBLE), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_core_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_core_ciface_keep), pointer :: fcikeep
  type(ssmfe_core_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_core_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (c_associated(ckeep)) then
     call c_f_pointer(ckeep, fcikeep)
  else
     allocate(fcikeep)
     ckeep = c_loc(fcikeep)
  end if
  if (crci%job .eq. 0) fcikeep%rcid%job = 0
  if ((fcikeep%rcid%job .eq. 999) .and. (fcikeep%rcid%k .gt. 0)) then
     fcikeep%rcid%i = crci%i
     fcikeep%rcid%j = crci%j
  end if

  ! Call Fortran routine
  call ssmfe(fcikeep%rcid, problem, left, right, m, lambda, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rcid, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_double

subroutine spral_ssmfe_double_complex(crci, problem, left, right, m, &
     lambda, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_core_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  integer(C_INT), value :: problem
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(m), intent(inout) :: lambda
  complex(C_DOUBLE_COMPLEX), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_core_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_core_ciface_keep), pointer :: fcikeep
  type(ssmfe_core_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_core_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (c_associated(ckeep)) then
     call c_f_pointer(ckeep, fcikeep)
  else
     allocate(fcikeep)
     ckeep = c_loc(fcikeep)
  end if
  if (crci%job .eq. 0) fcikeep%rciz%job = 0
  if ((fcikeep%rciz%job .eq. 999) .and. (fcikeep%rciz%k .gt. 0)) then
     fcikeep%rciz%i = crci%i
     fcikeep%rciz%j = crci%j
  end if

  ! Call Fortran routine
  call ssmfe(fcikeep%rciz, problem, left, right, m, lambda, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rciz, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_double_complex

subroutine spral_ssmfe_largest_double(crci, problem, nep, m, lambda, rr, &
     ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_core_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  integer(C_INT), value :: problem
  integer(C_INT), value :: nep
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(m), intent(inout) :: lambda
  real(C_DOUBLE), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_core_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_core_ciface_keep), pointer :: fcikeep
  type(ssmfe_core_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_core_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (c_associated(ckeep)) then
     call c_f_pointer(ckeep, fcikeep)
  else
     allocate(fcikeep)
     ckeep = c_loc(fcikeep)
  end if
  if (crci%job .eq. 0) fcikeep%rcid%job = 0
  if ((fcikeep%rcid%job .eq. 999) .and. (fcikeep%rcid%k .gt. 0)) then
     fcikeep%rcid%i = crci%i
     fcikeep%rcid%j = crci%j
  end if

  ! Call Fortran routine
  call ssmfe_largest(fcikeep%rcid, problem, nep, m, lambda, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rcid, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_largest_double

subroutine spral_ssmfe_largest_double_complex(crci, problem, nep, m, lambda, &
     rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_core_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  integer(C_INT), value :: problem
  integer(C_INT), value :: nep
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(m), intent(inout) :: lambda
  complex(C_DOUBLE_COMPLEX), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_core_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_core_ciface_keep), pointer :: fcikeep
  type(ssmfe_core_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_core_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (c_associated(ckeep)) then
     call c_f_pointer(ckeep, fcikeep)
  else
     allocate(fcikeep)
     ckeep = c_loc(fcikeep)
  end if
  if (crci%job .eq. 0) fcikeep%rcid%job = 0
  if ((fcikeep%rcid%job .eq. 999) .and. (fcikeep%rcid%k .gt. 0)) then
     fcikeep%rcid%i = crci%i
     fcikeep%rcid%j = crci%j
  end if

  ! Call Fortran routine
  call ssmfe_largest(fcikeep%rciz, problem, nep, m, lambda, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rciz, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_largest_double_complex

subroutine spral_ssmfe_core_free(ckeep, cinform) bind(C)
  use spral_ssmfe_core_ciface
  implicit none

  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_inform), intent(inout) :: cinform

  type(ssmfe_core_ciface_keep), pointer :: fcikeep

  ! Nullify pointer components of cinform
  cinform%converged       = C_NULL_PTR
  cinform%residual_norms  = C_NULL_PTR
  cinform%err_lambda      = C_NULL_PTR
  cinform%err_x           = C_NULL_PTR

  ! Check ckeep is not null
  if (.not. c_associated(ckeep)) return

  ! Associate fkeep
  call c_f_pointer(ckeep, fcikeep)

  ! Call Fortran cleanup
  call ssmfe_free(fcikeep%keep, fcikeep%inform)

  ! Free fcikeep
  deallocate(fcikeep)
  ckeep = C_NULL_PTR
end subroutine spral_ssmfe_core_free
