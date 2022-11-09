subroutine spral_print_matrix(lines, matrix_type, m, n, ptr, row, cval, base) &
     bind(C)
  use iso_c_binding
  use spral_matrix_util
  implicit none

  integer(C_INT), value :: lines
  integer(C_INT), value :: matrix_type
  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  type(C_PTR), value :: cval
  integer(C_INT), value :: base

  real(C_DOUBLE), dimension(:), pointer, contiguous :: fval

  if (C_ASSOCIATED(cval)) then
     call C_F_POINTER(cval, fval, shape = (/ ptr(n+1)-1 /) )
  else
     nullify(fval)
  end if

  if (ASSOCIATED(fval)) then
     call print_matrix(6, lines, matrix_type, m, n, ptr, row, &
          cbase=(base .eq. 0), val=fval)
  else
     call print_matrix(6, lines, matrix_type, m, n, ptr, row, &
          cbase=(base .eq. 0))
  end if
end subroutine spral_print_matrix

subroutine spral_print_matrix_i64d(lines, matrix_type, m, n, ptr, row, cval, &
     base) bind(C)
  use iso_c_binding
  use spral_matrix_util
  implicit none

  integer(C_INT), value :: lines
  integer(C_INT), value :: matrix_type
  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(in) :: ptr
  integer(C_INT), dimension(*), intent(in) :: row
  type(C_PTR), value :: cval
  integer(C_INT), value :: base

  real(C_DOUBLE), dimension(:), pointer, contiguous :: fval

  if (C_ASSOCIATED(cval)) then
     call C_F_POINTER(cval, fval, shape = (/ ptr(n+1)-1 /) )
  else
     nullify(fval)
  end if

  if (ASSOCIATED(fval)) then
     call print_matrix(6, lines, matrix_type, m, n, ptr, row, &
          cbase=(base .eq. 0), val=fval)
  else
     call print_matrix(6, lines, matrix_type, m, n, ptr, row, &
          cbase=(base .eq. 0))
  end if
end subroutine spral_print_matrix_i64d

subroutine spral_half_to_full_i64d(n, ptr, row, cval, base) bind(C)
  use spral_matrix_util
  use, intrinsic :: iso_c_binding
  implicit none

  integer(C_INT), value :: n
  integer(C_INT64_T), dimension(*), intent(inout) :: ptr
  integer(C_INT),  dimension(*),intent(inout) :: row
  type(C_PTR), value :: cval
  integer(C_INT), value :: base

  integer, dimension(:), allocatable :: iw
  real(C_DOUBLE), dimension(:), pointer, contiguous :: fval

  ! Import arguments
  if (C_ASSOCIATED(cval)) then
     call c_f_pointer(cval, fval, shape=(/ ptr(n+1)-1 /))
  else
     nullify(fval)
  end if

  ! Make main call
  allocate(iw(n))
  if (associated(fval)) then
     call half_to_full(n, row, ptr, iw, cbase=(base .eq. 0), a=fval)
  else
     call half_to_full(n, row, ptr, iw, cbase=(base .eq. 0))
  end if
end subroutine spral_half_to_full_i64d
