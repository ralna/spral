module spral_scaling_ciface
  use, intrinsic :: iso_c_binding
  use spral_scaling
  implicit none

  integer, parameter :: long = selected_int_kind(18)

  type, bind(C) :: spral_scaling_auction_options
     integer(C_INT) :: array_base
     integer(C_INT) :: max_iterations
     integer(C_INT) :: max_unchanged(3)
     real(C_FLOAT) :: min_proportion(3)
     real(C_FLOAT) :: eps_initial
     character(C_CHAR) :: unused(80)
  end type spral_scaling_auction_options

  type, bind(C) :: spral_scaling_auction_inform
     integer(C_INT) :: flag
     integer(C_INT) :: stat
     integer(C_INT) :: matched
     integer(C_INT) :: iterations
     integer(C_INT) :: unmatchable
     character(C_CHAR) :: unused(80)
  end type spral_scaling_auction_inform

  type, bind(C) :: spral_scaling_equilib_options
     integer(C_INT) :: array_base
     integer(C_INT) :: max_iterations
     real(C_FLOAT) :: tol
     character(C_CHAR) :: unused(80)
  end type spral_scaling_equilib_options

  type, bind(C) :: spral_scaling_equilib_inform
     integer(C_INT) :: flag
     integer(C_INT) :: stat
     integer(C_INT) :: iterations
     character(C_CHAR) :: unused(80)
  end type spral_scaling_equilib_inform

  type, bind(C) :: spral_scaling_hungarian_options
     integer(C_INT) :: array_base
     logical(C_BOOL) :: scale_if_singular
     character(C_CHAR) :: unused(80)
  end type spral_scaling_hungarian_options

  type, bind(C) :: spral_scaling_hungarian_inform
     integer(C_INT) :: flag
     integer(C_INT) :: stat
     integer(C_INT) :: matched
     character(C_CHAR) :: unused(80)
  end type spral_scaling_hungarian_inform

contains
   ! Auction-related data types
  subroutine copy_auction_options_in(coptions, foptions, cindexed)
    implicit none
    type(spral_scaling_auction_options), intent(in) :: coptions
    type(auction_options), intent(out) :: foptions
    logical, intent(out) :: cindexed

    cindexed                   = (coptions%array_base.eq.0)
    foptions%max_iterations    = coptions%max_iterations
    foptions%max_unchanged(:)  = coptions%max_unchanged(:)
    foptions%min_proportion(:) = coptions%min_proportion(:)
    foptions%eps_initial       = coptions%eps_initial
  end subroutine copy_auction_options_in
  subroutine copy_auction_inform_out(finform, cinform)
    implicit none
    type(auction_inform), intent(in) :: finform
    type(spral_scaling_auction_inform), intent(out) :: cinform

    cinform%flag        = finform%flag
    cinform%stat        = finform%stat
    cinform%matched     = finform%matched
    cinform%iterations  = finform%iterations
    cinform%unmatchable = finform%unmatchable
  end subroutine copy_auction_inform_out

  ! Equilib-related data types
  subroutine copy_equilib_options_in(coptions, foptions, cindexed)
    implicit none
    type(spral_scaling_equilib_options), intent(in) :: coptions
    type(equilib_options), intent(out) :: foptions
    logical, intent(out) :: cindexed

    cindexed                = (coptions%array_base.eq.0)
    foptions%max_iterations = coptions%max_iterations
    foptions%tol            = coptions%tol
  end subroutine copy_equilib_options_in
  subroutine copy_equilib_inform_out(finform, cinform)
    implicit none
    type(equilib_inform), intent(in) :: finform
    type(spral_scaling_equilib_inform), intent(out) :: cinform

    cinform%flag       = finform%flag
    cinform%stat       = finform%stat
    cinform%iterations = finform%iterations
  end subroutine copy_equilib_inform_out

  ! Hungarian-related data types
  subroutine copy_hungarian_options_in(coptions, foptions, cindexed)
    implicit none
    type(spral_scaling_hungarian_options), intent(in) :: coptions
    type(hungarian_options), intent(out) :: foptions
    logical, intent(out) :: cindexed

    cindexed                   = (coptions%array_base.eq.0)
    foptions%scale_if_singular = coptions%scale_if_singular
  end subroutine copy_hungarian_options_in
  subroutine copy_hungarian_inform_out(finform, cinform)
    implicit none
    type(hungarian_inform), intent(in) :: finform
    type(spral_scaling_hungarian_inform), intent(out) :: cinform

    cinform%flag    = finform%flag
    cinform%stat    = finform%stat
    cinform%matched = finform%matched
  end subroutine copy_hungarian_inform_out
end module spral_scaling_ciface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Default option routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spral_scaling_auction_default_options(coptions) bind(C)
  use spral_scaling_ciface
  implicit none

  type(spral_scaling_auction_options), intent(out) :: coptions

  type(auction_options) :: default_options

  coptions%array_base        = 0 ! C
  coptions%max_iterations    = default_options%max_iterations
  coptions%max_unchanged(:)  = default_options%max_unchanged(:)
  coptions%min_proportion(:) = default_options%min_proportion(:)
  coptions%eps_initial       = default_options%eps_initial
end subroutine spral_scaling_auction_default_options

subroutine spral_scaling_equilib_default_options(coptions) bind(C)
  use spral_scaling_ciface
  implicit none

  type(spral_scaling_equilib_options), intent(out) :: coptions

  type(equilib_options) :: default_options

  coptions%array_base     = 0 ! C
  coptions%max_iterations = default_options%max_iterations
  coptions%tol            = default_options%tol
end subroutine spral_scaling_equilib_default_options

subroutine spral_scaling_hungarian_default_options(coptions) bind(C)
  use spral_scaling_ciface
  implicit none

  type(spral_scaling_hungarian_options), intent(out) :: coptions

  type(hungarian_options) :: default_options

  coptions%array_base        = 0 ! C
  coptions%scale_if_singular = default_options%scale_if_singular
end subroutine spral_scaling_hungarian_default_options

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Symmetric scaling routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spral_scaling_auction_sym(n, ptr, row, val, scaling, cmatch, &
     coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: n
  integer(C_INT), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: scaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_auction_options), intent(in) :: coptions
  type(spral_scaling_auction_inform), intent(out) :: cinform

  logical :: cindexed
  integer, dimension(:), allocatable :: ptr_alloc, row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(auction_options) :: foptions
  type(auction_inform) :: finform

  ! Associate arrays and copy control in
  call copy_auction_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ n /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call auction_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform, match=fmatch)
     else
        call auction_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call auction_scale_sym(n, ptr, row, val, scaling, foptions, &
             finform, match=fmatch)
     else
        call auction_scale_sym(n, ptr, row, val, scaling, foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_auction_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
       fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_auction_sym

subroutine spral_scaling_auction_sym_long(n, ptr, row, val, scaling, cmatch, &
     coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: scaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_auction_options), intent(in) :: coptions
  type(spral_scaling_auction_inform), intent(out) :: cinform

  logical :: cindexed
  integer(long), dimension(:), allocatable :: ptr_alloc
  integer, dimension(:), allocatable :: row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(auction_options) :: foptions
  type(auction_inform) :: finform

  ! Associate arrays and copy control in
  call copy_auction_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ n /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call auction_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform, match=fmatch)
     else
        call auction_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call auction_scale_sym(n, ptr, row, val, scaling, foptions, &
             finform, match=fmatch)
     else
        call auction_scale_sym(n, ptr, row, val, scaling, foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_auction_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
      fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_auction_sym_long

subroutine spral_scaling_equilib_sym(n, ptr, row, val, scaling, &
     coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: n
  integer(C_INT), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: scaling
  type(spral_scaling_equilib_options), intent(in) :: coptions
  type(spral_scaling_equilib_inform), intent(out) :: cinform

  logical :: cindexed
  integer, dimension(:), allocatable :: ptr_alloc, row_alloc
  type(equilib_options) :: foptions
  type(equilib_inform) :: finform

  ! Associate arrays and copy control in
  call copy_equilib_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if

  ! Call Fortran routine
  if (cindexed) then
     call equilib_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
          foptions, finform)
  else
     call equilib_scale_sym(n, ptr, row, val, scaling, foptions, finform)
  end if

  ! Copy info out
  call copy_equilib_inform_out(finform, cinform)
end subroutine spral_scaling_equilib_sym

subroutine spral_scaling_equilib_sym_long(n, ptr, row, val, scaling, &
     coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: scaling
  type(spral_scaling_equilib_options), intent(in) :: coptions
  type(spral_scaling_equilib_inform), intent(out) :: cinform

  logical :: cindexed
  integer(long), dimension(:), allocatable :: ptr_alloc
  integer, dimension(:), allocatable :: row_alloc
  type(equilib_options) :: foptions
  type(equilib_inform) :: finform

  ! Associate arrays and copy control in
  call copy_equilib_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if

  ! Call Fortran routine
  if (cindexed) then
     call equilib_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
          foptions, finform)
  else
     call equilib_scale_sym(n, ptr, row, val, scaling, foptions, finform)
  end if

  ! Copy info out
  call copy_equilib_inform_out(finform, cinform)
end subroutine spral_scaling_equilib_sym_long

subroutine spral_scaling_hungarian_sym(n, ptr, row, val, scaling, cmatch, &
     coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: n
  integer(C_INT), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: scaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_hungarian_options), intent(in) :: coptions
  type(spral_scaling_hungarian_inform), intent(out) :: cinform

  logical :: cindexed
  integer, dimension(:), allocatable :: ptr_alloc, row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(hungarian_options) :: foptions
  type(hungarian_inform) :: finform

  ! Associate arrays and copy control in
  call copy_hungarian_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ n /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call hungarian_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform, match=fmatch)
     else
        call hungarian_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call hungarian_scale_sym(n, ptr, row, val, scaling, foptions, &
             finform, match=fmatch)
     else
        call hungarian_scale_sym(n, ptr, row, val, scaling, foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_hungarian_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
       fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_hungarian_sym

subroutine spral_scaling_hungarian_sym_long(n, ptr, row, val, scaling, cmatch, &
     coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: scaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_hungarian_options), intent(in) :: coptions
  type(spral_scaling_hungarian_inform), intent(out) :: cinform

  logical :: cindexed
  integer(long), dimension(:), allocatable :: ptr_alloc
  integer, dimension(:), allocatable :: row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(hungarian_options) :: foptions
  type(hungarian_inform) :: finform

  ! Associate arrays and copy control in
  call copy_hungarian_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ n /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call hungarian_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform, match=fmatch)
     else
        call hungarian_scale_sym(n, ptr_alloc, row_alloc, val, scaling, &
             foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call hungarian_scale_sym(n, ptr, row, val, scaling, foptions, &
             finform, match=fmatch)
     else
        call hungarian_scale_sym(n, ptr, row, val, scaling, foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_hungarian_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
       fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_hungarian_sym_long

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unsymmetric scaling routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spral_scaling_auction_unsym(m, n, ptr, row, val, rscaling, cscaling,&
     cmatch, coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: rscaling
  real(C_DOUBLE), dimension(*), intent(out) :: cscaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_auction_options), intent(in) :: coptions
  type(spral_scaling_auction_inform), intent(out) :: cinform

  logical :: cindexed
  integer, dimension(:), allocatable :: ptr_alloc, row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(auction_options) :: foptions
  type(auction_inform) :: finform

  ! Associate arrays and copy control in
  call copy_auction_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ m /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call auction_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform, match=fmatch)
     else
        call auction_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call auction_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform, match=fmatch)
     else
        call auction_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_auction_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
      fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_auction_unsym

subroutine spral_scaling_auction_unsym_long(m, n, ptr, row, val, rscaling, &
     cscaling, cmatch, coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: rscaling
  real(C_DOUBLE), dimension(*), intent(out) :: cscaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_auction_options), intent(in) :: coptions
  type(spral_scaling_auction_inform), intent(out) :: cinform

  logical :: cindexed
  integer(long), dimension(:), allocatable :: ptr_alloc
  integer, dimension(:), allocatable :: row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(auction_options) :: foptions
  type(auction_inform) :: finform

  ! Associate arrays and copy control in
  call copy_auction_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ m /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call auction_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform, match=fmatch)
     else
        call auction_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call auction_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform, match=fmatch)
     else
        call auction_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_auction_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
       fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_auction_unsym_long

subroutine spral_scaling_equilib_unsym(m, n, ptr, row, val, rscaling, &
     cscaling, coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: rscaling
  real(C_DOUBLE), dimension(*), intent(out) :: cscaling
  type(spral_scaling_equilib_options), intent(in) :: coptions
  type(spral_scaling_equilib_inform), intent(out) :: cinform

  logical :: cindexed
  integer, dimension(:), allocatable :: ptr_alloc, row_alloc
  type(equilib_options) :: foptions
  type(equilib_inform) :: finform

  ! Associate arrays and copy control in
  call copy_equilib_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if

  ! Call Fortran routine
  if (cindexed) then
     call equilib_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
          cscaling, foptions, finform)
  else
     call equilib_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
          foptions, finform)
  end if

  ! Copy info out
  call copy_equilib_inform_out(finform, cinform)
end subroutine spral_scaling_equilib_unsym

subroutine spral_scaling_equilib_unsym_long(m, n, ptr, row, val, rscaling, &
     cscaling, coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: rscaling
  real(C_DOUBLE), dimension(*), intent(out) :: cscaling
  type(spral_scaling_equilib_options), intent(in) :: coptions
  type(spral_scaling_equilib_inform), intent(out) :: cinform

  logical :: cindexed
  integer(long), dimension(:), allocatable :: ptr_alloc
  integer, dimension(:), allocatable :: row_alloc
  type(equilib_options) :: foptions
  type(equilib_inform) :: finform

  ! Associate arrays and copy control in
  call copy_equilib_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if

  ! Call Fortran routine
  if (cindexed) then
     call equilib_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
          cscaling, foptions, finform)
  else
     call equilib_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
          foptions, finform)
  end if

  ! Copy info out
  call copy_equilib_inform_out(finform, cinform)
end subroutine spral_scaling_equilib_unsym_long

subroutine spral_scaling_hungarian_unsym(m, n, ptr, row, val, rscaling, &
     cscaling, cmatch, coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: rscaling
  real(C_DOUBLE), dimension(*), intent(out) :: cscaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_hungarian_options), intent(in) :: coptions
  type(spral_scaling_hungarian_inform), intent(out) :: cinform

  logical :: cindexed
  integer, dimension(:), allocatable :: ptr_alloc, row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(hungarian_options) :: foptions
  type(hungarian_inform) :: finform

  ! Associate arrays and copy control in
  call copy_hungarian_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ m /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call hungarian_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform, match=fmatch)
     else
        call hungarian_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call hungarian_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform, match=fmatch)
     else
        call hungarian_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_hungarian_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
       fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_hungarian_unsym

subroutine spral_scaling_hungarian_unsym_long(m, n, ptr, row, val, rscaling, &
     cscaling, cmatch, coptions, cinform) bind(C)
  use spral_scaling_ciface
  implicit none

  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  real(C_DOUBLE), dimension(*), intent(in) :: val
  real(C_DOUBLE), dimension(*), intent(out) :: rscaling
  real(C_DOUBLE), dimension(*), intent(out) :: cscaling
  type(C_PTR), value :: cmatch
  type(spral_scaling_hungarian_options), intent(in) :: coptions
  type(spral_scaling_hungarian_inform), intent(out) :: cinform

  logical :: cindexed
  integer(long), dimension(:), allocatable :: ptr_alloc
  integer, dimension(:), allocatable :: row_alloc
  integer(C_INT), dimension(:), pointer, contiguous :: fmatch
  type(hungarian_options) :: foptions
  type(hungarian_inform) :: finform

  ! Associate arrays and copy control in
  call copy_hungarian_options_in(coptions, foptions, cindexed)
  if (cindexed) then
     allocate(ptr_alloc(n+1), row_alloc(ptr(n+1)))
     ptr_alloc(1:n+1) = ptr(1:n+1) + 1
     row_alloc(1:ptr(n+1)) = row(1:ptr(n+1)) + 1
  end if
  if (C_ASSOCIATED(cmatch)) then
     call C_F_POINTER(cmatch, fmatch, shape=(/ m /))
  else
     nullify(fmatch)
  end if

  ! Call Fortran routine
  if (cindexed) then
     if (associated(fmatch)) then
        call hungarian_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform, match=fmatch)
     else
        call hungarian_scale_unsym(m, n, ptr_alloc, row_alloc, val, rscaling, &
             cscaling, foptions, finform)
     end if
  else
     if (associated(fmatch)) then
        call hungarian_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform, match=fmatch)
     else
        call hungarian_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
             foptions, finform)
     end if
  end if

  ! Copy info out
  call copy_hungarian_inform_out(finform, cinform)
  if (cindexed .and. associated(fmatch)) &
       fmatch(:) = fmatch(:) - 1
end subroutine spral_scaling_hungarian_unsym_long
