integer(C_INT) function spral_random_matrix_generate(cstate, matrix_type, m, &
     n, nnz, ptr, row, cval, flags) bind(C)
  use iso_c_binding
  use spral_random, only: random_state, random_get_seed, random_set_seed
  use spral_random_matrix, only: random_matrix_generate
  implicit none

  integer(C_INT), intent(inout) :: cstate
  integer(C_INT), value :: matrix_type
  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT), value :: nnz
  integer(C_INT), dimension(n+1), intent(out) :: ptr
  integer(C_INT), dimension(nnz), intent(out) :: row
  type(C_PTR), value :: cval
  integer(C_INT), value :: flags

  integer, parameter :: wp = C_DOUBLE
  integer, parameter :: SPRAL_RANDOM_MATRIX_FINDEX       = 1
  integer, parameter :: SPRAL_RANDOM_MATRIX_NONSINGULAR  = 2
  integer, parameter :: SPRAL_RANDOM_MATRIX_SORT         = 4

  type(random_state) :: fstate
  real(wp), dimension(:), pointer, contiguous :: fval
  logical :: findex, nonsingular, sort

  ! Set random generator state
  call random_set_seed(fstate, cstate)

  ! Decipher flags
  findex      = (iand(flags, SPRAL_RANDOM_MATRIX_FINDEX)      .ne. 0)
  nonsingular = (iand(flags, SPRAL_RANDOM_MATRIX_NONSINGULAR) .ne. 0)
  sort        = (iand(flags, SPRAL_RANDOM_MATRIX_SORT)        .ne. 0)

  ! Check if we have a val vector
  if (C_ASSOCIATED(cval)) then
     call C_F_POINTER(cval, fval, shape = (/ nnz /))
  else
     nullify(fval)
  end if

  if (ASSOCIATED(fval)) then
     call random_matrix_generate(fstate, matrix_type, m, n, nnz, ptr, row, &
          spral_random_matrix_generate, nonsingular=nonsingular, sort=sort, &
          val=fval)
  else
     call random_matrix_generate(fstate, matrix_type, m, n, nnz, ptr, row, &
          spral_random_matrix_generate, nonsingular=nonsingular, sort=sort)
  end if

  ! Convert to C indexing if required
  if (.not. findex) then
     ptr(:) = ptr(:) - 1
     row(:) = row(:) - 1
  end if

  ! Recover new random genenerator state
  cstate = random_get_seed(fstate)
end function spral_random_matrix_generate

integer(C_INT) function spral_random_matrix_generate_long(cstate, matrix_type, &
     m, n, nnz, ptr, row, cval, flags) bind(C)
  use iso_c_binding
  use spral_random, only: random_state, random_get_seed, random_set_seed
  use spral_random_matrix, only: random_matrix_generate
  implicit none

  integer(C_INT), intent(inout) :: cstate
  integer(C_INT), value :: matrix_type
  integer(C_INT), value :: m
  integer(C_INT), value :: n
  integer(C_INT64_T), value :: nnz
  integer(C_INT64_T), dimension(n+1), intent(out) :: ptr
  integer(C_INT), dimension(nnz), intent(out) :: row
  type(C_PTR), value :: cval
  integer(C_INT), value :: flags

  integer, parameter :: wp = C_DOUBLE
  integer, parameter :: SPRAL_RANDOM_MATRIX_FINDEX       = 1
  integer, parameter :: SPRAL_RANDOM_MATRIX_NONSINGULAR  = 2
  integer, parameter :: SPRAL_RANDOM_MATRIX_SORT         = 4

  type(random_state) :: fstate
  real(wp), dimension(:), pointer, contiguous :: fval
  logical :: findex, nonsingular, sort

  ! Set random generator state
  call random_set_seed(fstate, cstate)

  ! Decipher flags
  findex      = (iand(flags, SPRAL_RANDOM_MATRIX_FINDEX)      .ne. 0)
  nonsingular = (iand(flags, SPRAL_RANDOM_MATRIX_NONSINGULAR) .ne. 0)
  sort        = (iand(flags, SPRAL_RANDOM_MATRIX_SORT)        .ne. 0)

  ! Check if we have a val vector
  if (C_ASSOCIATED(cval)) then
     call C_F_POINTER(cval, fval, shape = (/ nnz /))
  else
     nullify(fval)
  end if

  if (ASSOCIATED(fval)) then
     call random_matrix_generate(fstate, matrix_type, m, n, nnz, ptr, row, &
          spral_random_matrix_generate_long, nonsingular=nonsingular, sort=sort,&
          val=fval)
  else
     call random_matrix_generate(fstate, matrix_type, m, n, nnz, ptr, row, &
          spral_random_matrix_generate_long, nonsingular=nonsingular, sort=sort)
  end if

  ! Convert to C indexing if required
  if (.not. findex) then
     ptr(:) = ptr(:) - 1
     row(:) = row(:) - 1
  end if

  ! Recover new random genenerator state
  cstate = random_get_seed(fstate)
end function spral_random_matrix_generate_long
