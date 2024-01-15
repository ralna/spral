module spral_ssmfe_ciface
  use iso_c_binding
  use spral_ssmfe_expert_ciface
  use spral_ssmfe
  implicit none

  type ssmfe_ciface_keepd
     type(ssmfe_rcid) :: rci
     type(ssmfe_keepd) :: keep
     type(ssmfe_inform) :: inform
  end type ssmfe_ciface_keepd

  type ssmfe_ciface_keepz
     type(ssmfe_rciz) :: rci
     type(ssmfe_keepz) :: keep
     type(ssmfe_inform) :: inform
  end type ssmfe_ciface_keepz
end module spral_ssmfe_ciface

subroutine spral_ssmfe_standard_double(crci, left, mep, lambda, n, x, ldx, &
     ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  real(C_DOUBLE), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepd), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_standard(fcikeep%rci, left, mep, lambda, n, x, ldx, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_standard_double

subroutine spral_ssmfe_standard_double_complex(crci, left, mep, lambda, n, &
     x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  complex(C_DOUBLE_COMPLEX), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepz), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_standard(fcikeep%rci, left, mep, lambda, n, x, ldx, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_standard_double_complex

subroutine spral_ssmfe_standard_shift_double(crci, sigma, left, right, mep, &
     lambda, n, x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  real(C_DOUBLE), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepd), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_standard_shift(fcikeep%rci, sigma, left, right, mep, lambda, n, &
       x, ldx, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_standard_shift_double

subroutine spral_ssmfe_standard_shift_double_complex(crci, sigma, left, right, &
     mep, lambda, n, x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  complex(C_DOUBLE_COMPLEX), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepz), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_standard_shift(fcikeep%rci, sigma, left, right, mep, lambda, n, &
       x, ldx, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_standard_shift_double_complex

subroutine spral_ssmfe_generalized_double(crci, left, mep, lambda, n, x, ldx, &
     ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  real(C_DOUBLE), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepd), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_generalized(fcikeep%rci, left, mep, lambda, n, x, ldx, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_generalized_double

subroutine spral_ssmfe_generalized_double_complex(crci, left, mep, lambda, n, &
     x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  integer(C_INT), value :: left
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  real(C_DOUBLE), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepd), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_generalized(fcikeep%rci, left, mep, lambda, n, x, ldx, &
       fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_generalized_double_complex

subroutine spral_ssmfe_generalized_shift_double(crci, sigma, left, right, mep, &
     lambda, n, x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  real(C_DOUBLE), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepd), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_generalized_shift(fcikeep%rci, sigma, left, right, mep, lambda, &
       n, x, ldx, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_generalized_shift_double

subroutine spral_ssmfe_generalized_shift_double_complex(crci, sigma, left, &
     right, mep, lambda, n, x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  complex(C_DOUBLE_COMPLEX), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepz), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_generalized_shift(fcikeep%rci, sigma, left, right, mep, lambda, &
       n, x, ldx, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_generalized_shift_double_complex

subroutine spral_ssmfe_buckling_double(crci, sigma, left, right, mep, &
     lambda, n, x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rcid), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  real(C_DOUBLE), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepd), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_buckling(fcikeep%rci, sigma, left, right, mep, lambda, &
       n, x, ldx, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_buckling_double

subroutine spral_ssmfe_buckling_double_complex(crci, sigma, left, right, mep, &
     lambda, n, x, ldx, ckeep, coptions, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(spral_ssmfe_rciz), intent(inout) :: crci
  real(C_DOUBLE), value :: sigma
  integer(C_INT), value :: left
  integer(C_INT), value :: right
  integer(C_INT), value :: mep
  real(C_DOUBLE), dimension(mep), intent(inout) :: lambda
  integer(C_INT), value :: n
  integer(C_INT), value :: ldx
  complex(C_DOUBLE_COMPLEX), dimension(ldx, mep), intent(inout) :: x
  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_options), intent(in) :: coptions
  type(spral_ssmfe_inform), intent(inout) :: cinform

  logical :: cindexed
  type(ssmfe_ciface_keepz), pointer :: fcikeep
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
  if (crci%job .eq. 0) fcikeep%rci%job = 0

  ! Call Fortran routine
  call ssmfe_buckling(fcikeep%rci, sigma, left, right, mep, lambda, &
       n, x, ldx, fcikeep%keep, foptions, fcikeep%inform)

  ! Copy arguments out
  call copy_rci_out(fcikeep%rci, crci, cindexed)
  call copy_inform_out(fcikeep%inform, cinform)
end subroutine spral_ssmfe_buckling_double_complex

subroutine spral_ssmfe_free_double(ckeep, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_inform), intent(inout) :: cinform

  type(ssmfe_ciface_keepd), pointer :: fcikeep

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
end subroutine spral_ssmfe_free_double

subroutine spral_ssmfe_free_double_complex(ckeep, cinform) bind(C)
  use spral_ssmfe_ciface
  implicit none

  type(C_PTR), intent(inout) :: ckeep
  type(spral_ssmfe_inform), intent(inout) :: cinform

  type(ssmfe_ciface_keepz), pointer :: fcikeep

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
end subroutine spral_ssmfe_free_double_complex
