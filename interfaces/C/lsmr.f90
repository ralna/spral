module spral_lsmr_ciface
  use iso_c_binding
  use spral_lsmr
  implicit none

  type, bind(C) :: spral_lsmr_options
     real(C_DOUBLE) ::  atol
     real(C_DOUBLE) ::  btol
     real(C_DOUBLE) ::  conlim
     integer(C_INT) ::  ctest
     integer(C_INT) ::  itnlim
     integer(C_INT) ::  itn_test
     integer(C_INT) ::  localSize
     integer(C_INT) ::  print_freq_head
     integer(C_INT) ::  print_freq_itn
     integer(C_INT) ::  unit_diagnostics
     integer(C_INT) ::  unit_error
  end type spral_lsmr_options

  type, bind(C) :: spral_lsmr_inform
     integer(C_INT) ::  flag
     integer(C_INT) ::  itn
     integer(C_INT) ::  stat
     real(C_DOUBLE) ::  normb
     real(C_DOUBLE) ::  normAP
     real(C_DOUBLE) ::  condAP
     real(C_DOUBLE) ::  normr
     real(C_DOUBLE) ::  normAPr
     real(C_DOUBLE) ::  normy
  end type spral_lsmr_inform

contains
  subroutine copy_options_in(coptions, foptions)
    implicit none
    type(spral_lsmr_options), intent(in) :: coptions
    type(lsmr_options), intent(inout) :: foptions ! may inherit some defaults

    foptions%atol              = coptions%atol
    foptions%btol              = coptions%btol
    foptions%conlim            = coptions%conlim
    foptions%ctest             = coptions%ctest
    foptions%itnlim            = coptions%itnlim
    foptions%itn_test          = coptions%itn_test
    foptions%localSize         = coptions%localSize
    foptions%print_freq_head   = coptions%print_freq_head
    foptions%print_freq_itn    = coptions%print_freq_itn
    foptions%unit_diagnostics  = coptions%unit_diagnostics
    foptions%unit_error        = coptions%unit_error
  end subroutine copy_options_in

  subroutine copy_inform_in(cinform, finform)
    implicit none
    type(spral_lsmr_inform), intent(in) :: cinform
    type(lsmr_inform), intent(out) :: finform

    finform%flag      = cinform%flag
    finform%itn       = cinform%itn
    finform%stat      = cinform%stat
    finform%normb     = cinform%normb
    finform%normAP    = cinform%normAP
    finform%condAP    = cinform%condAP
    finform%normr     = cinform%normr
    finform%normAPr   = cinform%normAPr
    finform%normy     = cinform%normy
  end subroutine copy_inform_in

  subroutine copy_inform_out(finform, cinform)
    implicit none
    type(lsmr_inform), intent(in) :: finform
    type(spral_lsmr_inform), intent(out) :: cinform

    cinform%flag      = finform%flag
    cinform%itn       = finform%itn
    cinform%stat      = finform%stat
    cinform%normb     = finform%normb
    cinform%normAP    = finform%normAP
    cinform%condAP    = finform%condAP
    cinform%normr     = finform%normr
    cinform%normAPr   = finform%normAPr
    cinform%normy     = finform%normy
  end subroutine copy_inform_out
end module spral_lsmr_ciface

subroutine spral_lsmr_default_options(coptions) bind(C)
  use spral_lsmr_ciface
  implicit none

  type(spral_lsmr_options), intent(out) :: coptions

  type(lsmr_options) :: default_options

  coptions%atol              = default_options%atol
  coptions%btol              = default_options%btol
  coptions%conlim            = default_options%conlim
  coptions%ctest             = default_options%ctest
  coptions%itnlim            = default_options%itnlim
  coptions%itn_test          = default_options%itn_test
  coptions%localSize         = default_options%localSize
  coptions%print_freq_head   = default_options%print_freq_head
  coptions%print_freq_itn    = default_options%print_freq_itn
  coptions%unit_diagnostics  = default_options%unit_diagnostics
  coptions%unit_error        = default_options%unit_error
end subroutine spral_lsmr_default_options

integer(C_INT) function spral_lsmr_solve(action, m, n, u, v, y, ckeep, &
     coptions, cinform, cdamp) bind(C)
  use spral_lsmr_ciface
  implicit none

  integer(C_INT) :: action
  integer(C_INT), value :: m
  integer(C_INT), value :: n
  real(C_DOUBLE), dimension(*), intent(inout) :: u
  real(C_DOUBLE), dimension(*), intent(inout) :: v
  real(C_DOUBLE), dimension(*), intent(inout) :: y
  type(C_PTR), intent(inout) :: ckeep
  type(spral_lsmr_options), intent(in) :: coptions
  type(spral_lsmr_inform), intent(inout) :: cinform
  type(C_PTR), value :: cdamp

  type(lsmr_keep), pointer :: fkeep
  type(lsmr_options) :: foptions
  type(lsmr_inform) :: finform
  real(C_DOUBLE), pointer :: fdamp
  integer :: st

  ! Copy options in
  call copy_options_in(coptions, foptions)

  ! We must also copy inform in because lsmr_solve assumes that these are
  ! persistent between calls.
  call copy_inform_in(cinform, finform)

  ! Translate arguments
  if (c_associated(ckeep)) then
     call c_f_pointer(ckeep, fkeep)
  else
     allocate(fkeep, stat=st)
     if (st .ne. 0) then
        call copy_inform_out(finform, cinform)
        cinform%flag = 8 ! Allocation failed
        cinform%stat = st
        spral_lsmr_solve = cinform%flag
        return
     end if
     ckeep = c_loc(fkeep)
  end if
  nullify(fdamp)
  if (c_associated(cdamp)) call c_f_pointer(cdamp, fdamp)

  ! Call Fortran routine
  if (associated(fdamp)) then
     call lsmr_solve(action, m, n, u, v, y, fkeep, foptions, finform, &
          damp=fdamp)
  else
     call lsmr_solve(action, m, n, u, v, y, fkeep, foptions, finform)
  end if

  ! Copy arguments out
  call copy_inform_out(finform, cinform)

  ! Return status
  spral_lsmr_solve = cinform%flag
end function spral_lsmr_solve

integer(C_INT) function spral_lsmr_free(ckeep) bind(C)
  use spral_lsmr_ciface
  implicit none

  type(C_PTR), intent(inout) :: ckeep

  type(lsmr_keep), pointer :: fkeep

  if (.not. c_associated(ckeep)) then
     ! Nothing to do, already NULL
     spral_lsmr_free = 0 ! Success
     return
  end if

  ! Otherwise, ckeep is not NULL
  call c_f_pointer(ckeep, fkeep)

  ! Call Fortran routine, capture stat as return parameter
  call lsmr_free(fkeep, spral_lsmr_free)

  ! Set ckeep to NULL if deallocation succeeded
  if (spral_lsmr_free .eq. 0) ckeep = C_NULL_PTR
end function spral_lsmr_free
