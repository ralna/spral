module spral_ssmfe_expert_ciface
  use iso_c_binding
  use spral_ssmfe_core_ciface
  use spral_ssmfe_expert
  implicit none

  type, bind(C) :: spral_ssmfe_options
     integer(C_INT) :: array_base
     integer(C_INT) :: print_level
     integer(C_INT) :: unit_error
     integer(C_INT) :: unit_warning
     integer(C_INT) :: unit_diagnostic
     integer(C_INT) :: max_iterations
     integer(C_INT) :: user_x
     integer(C_INT) :: err_est
     real(C_DOUBLE) :: abs_tol_lambda
     real(C_DOUBLE) :: rel_tol_lambda
     real(C_DOUBLE) :: abs_tol_residual
     real(C_DOUBLE) :: rel_tol_residual
     real(C_DOUBLE) :: tol_x
     real(C_DOUBLE) :: left_gap
     real(C_DOUBLE) :: right_gap
     integer(C_INT) :: extra_left
     integer(C_INT) :: extra_right
     integer(C_INT) :: max_left
     integer(C_INT) :: max_right
     logical(C_BOOL) :: minAprod
     logical(C_BOOL) :: minBprod
  end type spral_ssmfe_options

  type ssmfe_expert_ciface_keep
     type(ssmfe_expert_keep) :: keep
     type(ssmfe_rcid) :: rcid
     type(ssmfe_rciz) :: rciz
     type(ssmfe_inform) :: inform
  end type ssmfe_expert_ciface_keep

contains
  subroutine copy_options_in(coptions, foptions, cindexed)
    implicit none
    type(spral_ssmfe_options), intent(in) :: coptions
    type(ssmfe_options), intent(inout) :: foptions ! inherits defaults!
    logical, intent(out) :: cindexed

    cindexed                   = (coptions%array_base .eq. 0)
    foptions%print_level       = coptions%print_level
    foptions%unit_error        = coptions%unit_error
    foptions%unit_warning      = coptions%unit_warning
    foptions%unit_diagnostic   = coptions%unit_diagnostic
    foptions%max_iterations    = coptions%max_iterations
    foptions%user_x            = coptions%user_x
    foptions%err_est           = coptions%err_est
    foptions%abs_tol_lambda    = coptions%abs_tol_lambda
    foptions%rel_tol_lambda    = coptions%rel_tol_lambda
    foptions%abs_tol_residual  = coptions%abs_tol_residual
    foptions%rel_tol_residual  = coptions%rel_tol_residual
    foptions%tol_x             = coptions%tol_x
    foptions%left_gap          = coptions%left_gap
    foptions%right_gap         = coptions%right_gap
    foptions%extra_left        = coptions%extra_left
    foptions%extra_right       = coptions%extra_right
    foptions%max_left          = coptions%max_left
    foptions%max_right         = coptions%max_right
    foptions%minAprod          = coptions%minAprod
    foptions%minBprod          = coptions%minBprod
  end subroutine copy_options_in
end module spral_ssmfe_expert_ciface

subroutine spral_ssmfe_default_options(coptions) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_options), intent(out) :: coptions

  type(ssmfe_options) :: default_options

  coptions%array_base  = 0 ! C
  coptions%print_level       = default_options%print_level
  coptions%unit_error        = default_options%unit_error
  coptions%unit_warning      = default_options%unit_warning
  coptions%unit_diagnostic   = default_options%unit_diagnostic
  coptions%max_iterations    = default_options%max_iterations
  coptions%user_x            = default_options%user_x
  coptions%err_est           = default_options%err_est
  coptions%abs_tol_lambda    = default_options%abs_tol_lambda
  coptions%rel_tol_lambda    = default_options%rel_tol_lambda
  coptions%abs_tol_residual  = default_options%abs_tol_residual
  coptions%rel_tol_residual  = default_options%rel_tol_residual
  coptions%tol_x             = default_options%tol_x
  coptions%left_gap          = default_options%left_gap
  coptions%right_gap         = default_options%right_gap
  coptions%extra_left        = default_options%extra_left
  coptions%extra_right       = default_options%extra_right
  coptions%max_left          = default_options%max_left
  coptions%max_right         = default_options%max_right
  coptions%minAprod          = default_options%minAprod
  coptions%minBprod          = default_options%minBprod
end subroutine spral_ssmfe_default_options

subroutine spral_ssmfe_expert_standard_double(crci, left, mep, lambda, m, rr, &
     ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_standard(fcikeep%rcid, left, mep, lambda, m, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rcid, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq.11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_standard_double

subroutine spral_ssmfe_expert_standard_double_complex(crci, left, mep, lambda, &
     m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  complex(C_DOUBLE_COMPLEX), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_standard(fcikeep%rciz, left, mep, lambda, m, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rciz, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_standard_double_complex

subroutine spral_ssmfe_expert_standard_shift_double(crci, sigma, left, right, &
     mep, lambda, m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_standard_shift(fcikeep%rcid, sigma, left, right, mep, lambda, &
       m, rr, ind, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rcid, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_standard_shift_double

subroutine spral_ssmfe_expert_standard_shift_double_complex(crci, sigma, left, &
     right, mep, lambda, m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  complex(C_DOUBLE_COMPLEX), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_standard_shift(fcikeep%rciz, sigma, left, right, mep, lambda, &
       m, rr, ind, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rciz, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_standard_shift_double_complex

subroutine spral_ssmfe_expert_generalized_double(crci, left, mep, lambda, m, &
     rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_generalized(fcikeep%rcid, left, mep, lambda, m, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rcid, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_generalized_double

subroutine spral_ssmfe_expert_generalized_double_complex(crci, left, mep, &
     lambda, m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  complex(C_DOUBLE_COMPLEX), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_generalized(fcikeep%rciz, left, mep, lambda, m, rr, ind, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rciz, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_generalized_double_complex

subroutine spral_ssmfe_expert_generalized_shift_double(crci, sigma, left, &
     right, mep, lambda, m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_generalized_shift(fcikeep%rcid, sigma, left, right, mep, lambda, &
       m, rr, ind, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rcid, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_generalized_shift_double

subroutine spral_ssmfe_expert_generalized_shift_double_complex(crci, sigma, &
     left, right, mep, lambda, m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  complex(C_DOUBLE_COMPLEX), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

  ! Translate arguments
  if (c_associated(ckeep)) then
     call c_f_pointer(ckeep, fcikeep)
  else
     allocate(fcikeep)
     ckeep = c_loc(fcikeep)
  endif
  if (crci%job .eq. 0) fcikeep%rciz%job = 0
  if ((fcikeep%rciz%job .eq. 999) .and. (fcikeep%rciz%k .gt. 0)) then
     fcikeep%rciz%i = crci%i
     fcikeep%rciz%j = crci%j
  end if

  ! Call Fortran routine
  call ssmfe_generalized_shift(fcikeep%rciz, sigma, left, right, mep, lambda, &
       m, rr, ind, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rciz, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_generalized_shift_double_complex

subroutine spral_ssmfe_expert_buckling_double(crci, sigma, left, right, mep, &
     lambda, m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  real(C_DOUBLE), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_buckling(fcikeep%rcid, sigma, left, right, mep, lambda, &
       m, rr, ind, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rcid, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_buckling_double

subroutine spral_ssmfe_expert_buckling_double_complex(crci, sigma, left, &
     right, mep, lambda, m, rr, ind, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: m
  complex(C_DOUBLE_COMPLEX), dimension(2*m, 2*m, 3), intent(inout) :: rr
  integer(C_INT), dimension(m), intent(out) :: ind
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_expert_ciface_keep), pointer :: fcikeep
  type(ssmfe_options) :: foptions

  ! Copy options in first to find out whether we use Fortran or C indexing
  call copy_options_in(coptions, foptions, cindexed)

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
  call ssmfe_buckling(fcikeep%rciz, sigma, left, right, mep, lambda, &
       m, rr, ind, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rciz, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
  if ((crci%job .eq. 11) .and. cindexed) &
       ind(1:crci%nx) = ind(1:crci%nx) - 1
end subroutine spral_ssmfe_expert_buckling_double_complex

subroutine spral_ssmfe_expert_free(ckeep, cinform) bind(C)
  use spral_ssmfe_expert_ciface
  implicit none

  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_inform), intent(inout) :: cinform

  type(ssmfe_expert_ciface_keep), pointer :: fcikeep

  ! Nullify pointer components of cinform
  cinform%converged       = C_NULL_PTR
  cinform%residual_norms  = C_NULL_PTR
  cinform%err_lambda      = C_NULL_PTR
  cinform%err_x           = C_NULL_PTR

  ! Check ckeep is not null
  if(.not. c_associated(ckeep)) return

  ! Associate fkeep
  call c_f_pointer(ckeep, fcikeep)

  ! Call Fortran cleanup
  call ssmfe_free(fcikeep%keep, fcikeep%inform)

  ! Free fcikeep
  deallocate(fcikeep)
  ckeep = C_NULL_PTR
end subroutine spral_ssmfe_expert_free
