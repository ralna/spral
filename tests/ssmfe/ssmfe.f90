program main

  implicit none

  integer, parameter :: wp = kind(0d0)
  real(wp), parameter :: ONE = 1.0_wp
  real(wp), parameter :: ZERO = 0.0_wp
  complex(wp), parameter :: ZZERO = (ZERO, ZERO)
  complex(wp), parameter :: ZONE = (ONE, ZERO)

  integer, parameter :: WRONG_RCI_JOB      = -1
  integer, parameter :: WRONG_BLOCK_SIZE   = -2
  integer, parameter :: WRONG_ERR_EST      = -3
  integer, parameter :: WRONG_MINPROD      = -4
  integer, parameter :: WRONG_EXTRAS       = -5
  integer, parameter :: WRONG_MIN_GAP      = -6
  integer, parameter :: WRONG_CF_MAX       = -7
  integer, parameter :: WRONG_PROBLEM_SIZE = -9
  integer, parameter :: WRONG_LDX          = -10
  integer, parameter :: WRONG_LEFT         = -11
  integer, parameter :: WRONG_RIGHT        = -12
  integer, parameter :: WRONG_STORAGE_SIZE = -13
  integer, parameter :: WRONG_SIGMA        = -14

  integer, parameter :: OUT_OF_MEMORY       = -100
  integer, parameter :: INDEFINITE_B_OR_XBX = -200

  integer, parameter :: NO_SEARCH_DIRECTIONS_LEFT = 1
  integer, parameter :: MAX_ITERATIONS_EXCEEDED   = 2
  integer, parameter :: OUT_OF_STORAGE            = 3

  integer, parameter :: we_unit = 11
  integer, parameter :: dl_unit = 12

  character(len=6) :: we_file = "we.out"
  character(len=6) :: dl_file = "dl.out"

  integer :: errors

  if(we_unit.gt.6) open(unit=we_unit,file=we_file,status="replace")
  if(dl_unit.gt.6) open(unit=dl_unit,file=dl_file,status="replace")

  errors = 0

  call test_expert
  call test_ssmfe
  call test_core

  write(*, "(/a)") "============================="
  write(*, "(a,i4)") "Total number of errors = ", errors

  if(we_unit.gt.6) close(we_unit)
  if(dl_unit.gt.6) close(dl_unit)

  if(errors.ne.0) stop 1 ! ERROR CODE for make check script

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe

  write(*,"(/a)") "=============="
  write(*,"(a)")  "Testing ssmfe:"
  write(*,"(a)")  "=============="

  call test_ssmfe_errors_d
  call test_ssmfe_errors_z
  call test_ssmfe_warnings_d
  call test_ssmfe_warnings_z
  call test_ssmfe_options_d
  call test_ssmfe_options_z
  call test_ssmfe_misc_d
  call test_ssmfe_misc_z

end subroutine test_ssmfe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_core

  write(*,"(/a)") "==================="
  write(*,"(a)")  "Testing ssmfe_core:"
  write(*,"(a)")  "==================="

  call test_core_errors_d
  call test_core_errors_z
  call test_core_misc_d
  call test_core_misc_z

end subroutine test_core

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_core_errors_d

  use spral_ssmfe_core !, ssmfe_work => ssmfe_keep !, ssmfe_opts => ssmfe_options

  integer :: n
  integer :: m
  integer :: left, right

  integer, allocatable :: ind(:)

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: rr(:,:,:)

  type(ssmfe_rcid        ) :: rci
  type(ssmfe_core_options) :: options
  type(ssmfe_core_keep   ) :: keep
  type(ssmfe_inform      ) :: inform

  write(*,"(/a)") "======================"
  write(*,"(a)")  "Testing errors (real):"
  write(*,"(a)")  "======================"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n = 6
  m = 2
  left = 1
  right = 1
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )
  allocate ( rr(m + m, m + m, 3), ind(m) )

  write(*,"(a)",advance="no") " * Testing extra_left < 0...................."
  rci%job = 0
  options%extra_left = -1
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_EXTRAS )
  call ssmfe_free( keep, inform )
  options%extra_left = 0

  write(*,"(a)",advance="no") " * Testing min_gap < 0......................."
  rci%job = 0
  options%min_gap = -1
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_MIN_GAP )
  call ssmfe_free( keep, inform )
  options%min_gap = 0

  write(*,"(a)",advance="no") " * Testing min_gap > 1......................."
  rci%job = 0
  options%min_gap = 2
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_MIN_GAP )
  call ssmfe_free( keep, inform )
  options%min_gap = 0

  write(*,"(a)",advance="no") " * Testing cf_max < 0.5......................"
  rci%job = 0
  options%cf_max = 0.49999
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_CF_MAX )
  call ssmfe_free( keep, inform )
  options%cf_max = ONE

  write(*,"(a)",advance="no") " * Testing cf_max > 1........................"
  rci%job = 0
  options%cf_max = 1.00001
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_CF_MAX )
  call ssmfe_free( keep, inform )
  options%cf_max = ONE

  write(*,"(a)",advance="no") " * Testing wrong block size.................."
  rci%job = 0
  call ssmfe &
    ( rci, 0, left, 0, 0, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE, .true. )
  call ssmfe_free( keep, inform )
  rci%job = 0
  call ssmfe &
    ( rci, 0, left, right, 1, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE, .true. )
  call ssmfe_free( keep, inform )
  rci%job = 0
  call ssmfe &
    ( rci, 0, -1, right, 1, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing wrong error estimation flag......."
  rci%job = 0
  options%err_est = 0
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_ERR_EST )
  call ssmfe_free( keep, inform )
  options%err_est = 2

  write(*,"(a)",advance="no") " * Testing wrong minAprod...................."
  rci%job = 0
  options%minAprod = .false.
  call ssmfe &
    ( rci, -1, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_MINPROD )
  call ssmfe_free( keep, inform )
  options%minAprod = .true.

  deallocate ( a, b, t, x, lambda, rr, ind )

end subroutine test_core_errors_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_core_errors_z

  use spral_ssmfe_core !, ssmfe_work => ssmfe_keep !, ssmfe_opts => ssmfe_options

  integer :: n
  integer :: m
  integer :: left, right

  integer, allocatable :: ind(:)

  real(wp), allocatable :: lambda(:)
  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)
  complex(wp), allocatable :: rr(:,:,:)

  type(ssmfe_rciz        ) :: rci
  type(ssmfe_core_options) :: options
  type(ssmfe_core_keep   ) :: keep
  type(ssmfe_inform      ) :: inform

  write(*,"(/a)") "========================="
  write(*,"(a)")  "Testing errors (complex):"
  write(*,"(a)")  "========================="

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n = 6
  m = 2
  left = 1
  right = 1
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )
  allocate ( rr(m + m, m + m, 3), ind(m) )

  write(*,"(a)",advance="no") " * Testing extra_left < 0...................."
  rci%job = 0
  options%extra_left = -1
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_EXTRAS )
  call ssmfe_free( keep, inform )
  options%extra_left = 0

  write(*,"(a)",advance="no") " * Testing min_gap < 0......................."
  rci%job = 0
  options%min_gap = -1
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_MIN_GAP )
  call ssmfe_free( keep, inform )
  options%min_gap = 0

  write(*,"(a)",advance="no") " * Testing min_gap > 1......................."
  rci%job = 0
  options%min_gap = 2
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_MIN_GAP )
  call ssmfe_free( keep, inform )
  options%min_gap = 0

  write(*,"(a)",advance="no") " * Testing cf_max < 0.5......................"
  rci%job = 0
  options%cf_max = 0.49999
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_CF_MAX )
  call ssmfe_free( keep, inform )
  options%cf_max = ONE

  write(*,"(a)",advance="no") " * Testing cf_max > 1........................"
  rci%job = 0
  options%cf_max = 1.00001
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_CF_MAX )
  call ssmfe_free( keep, inform )
  options%cf_max = ONE

  write(*,"(a)",advance="no") " * Testing wrong block size.................."
  rci%job = 0
  call ssmfe &
    ( rci, 0, left, 0, 0, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE, .true. )
  call ssmfe_free( keep, inform )
  rci%job = 0
  call ssmfe &
    ( rci, 0, left, right, 1, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE, .true. )
  call ssmfe_free( keep, inform )
  rci%job = 0
  call ssmfe &
    ( rci, 0, -1, right, 1, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing wrong error estimation flag......."
  rci%job = 0
  options%err_est = 0
  call ssmfe &
    ( rci, 0, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_ERR_EST )
  call ssmfe_free( keep, inform )
  options%err_est = 2

  write(*,"(a)",advance="no") " * Testing wrong minAprod...................."
  rci%job = 0
  options%minAprod = .false.
  call ssmfe &
    ( rci, -1, left, right, m, lambda, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_MINPROD )
  call ssmfe_free( keep, inform )
  options%minAprod = .true.

  deallocate ( a, b, t, x, lambda, rr, ind )

end subroutine test_core_errors_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_core_misc_d

  use spral_ssmfe_core !, ssmfe_work => ssmfe_keep !, ssmfe_opts => ssmfe_options

  integer :: n
  integer :: m
  integer :: maxit
  integer :: verb
  integer :: i, j

  integer, allocatable :: ind(:)

  real(wp) :: tol

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: rr(:,:,:)

  type(ssmfe_core_options) :: options
  type(ssmfe_inform      ) :: inform

  write(*,"(/a)") "==========================="
  write(*,"(a)")  "Miscellaneous tests (real):"
  write(*,"(a)")  "==========================="

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n = 50
  m = 2
  maxit = 300
  verb = 0
  tol = 1e-3

  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )
  allocate ( rr(m + m, m + m, 3), ind(m) )

  a = ZERO
  forall ( i = 1 : n ) a(i, i) = i*10 - n*5

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing left = 1, right = 0, m = 1........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_d &
    ( options, 0, n, a, b, t, 1, 0, tol, maxit, verb, 1, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing left = 2, right = 0, m = 1........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_d &
    ( options, 0, n, a, b, t, 2, 0, tol, maxit, verb, 1, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing left = 0, right = 1, m = 1........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_d &
    ( options, 0, n, a, b, t, 0, 1, tol, maxit, verb, 1, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing left = 1, right = 1, m = 2........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_d &
    ( options, 0, n, a, b, t, 1, 1, tol, maxit, verb, 2, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing left = 2, right = 1, m = 2........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_d &
    ( options, 0, n, a, b, t, 2, 1, tol, maxit, verb, 2, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing left = 10, right = 10, m = 10....."
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  do i = 2, n
    x(i, 1) = ZERO
    x(1, i ) = ZERO
  end do
  call run_ssmfe_d &
    ( options, 0, n, a, b, t, 10, 10, tol, maxit, verb, 10, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing largest = 1, m = 2................"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_largest_d &
    ( options, 0, n, a, b, t, 1, tol, maxit, verb, 2, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing largest = 10, m = 5..............."
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_largest_d &
    ( options, 0, n, a, b, t, 10, tol, maxit, verb, 5, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing left = 10, right = 10, m = 10....."
  forall ( i = 1 : n ) b(i, i) = 2.0D0**i
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_d &
    ( options, 0, n, a, b, t, 10, 10, tol, maxit, verb, 10, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  deallocate ( a, b, t, x, lambda, rr, ind )

end subroutine test_core_misc_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_core_misc_z

  use spral_ssmfe_core !, ssmfe_work => ssmfe_keep !, ssmfe_opts => ssmfe_options

  integer :: n
  integer :: m
  integer :: maxit
  integer :: verb
  integer :: i, j

  integer, allocatable :: ind(:)

  real(wp) :: tol

  real(wp), allocatable :: lambda(:)
  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)
  complex(wp), allocatable :: rr(:,:,:)

  type(ssmfe_core_options) :: options
  type(ssmfe_core_keep   ) :: keep
  type(ssmfe_inform      ) :: inform

  write(*,"(/a)") "=============================="
  write(*,"(a)")  "Miscellaneous tests (complex):"
  write(*,"(a)")  "=============================="

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n = 100
  m = 2
  maxit = 300
  verb = 0
  tol = 1e-3

  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )
  allocate ( rr(m + m, m + m, 3), ind(m) )

  a = ZERO
  forall ( i = 1 : n ) a(i, i) = i*10 - n*5

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing left = 1, right = 0, m = 1........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_z &
    ( options, 0, n, a, b, t, 1, 0, tol, maxit, verb, 1, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left = 2, right = 0, m = 1........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_z &
    ( options, 0, n, a, b, t, 2, 0, tol, maxit, verb, 1, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left = 0, right = 1, m = 1........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_z &
    ( options, 0, n, a, b, t, 0, 1, tol, maxit, verb, 1, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left = 1, right = 1, m = 2........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_z &
    ( options, 0, n, a, b, t, 1, 1, tol, maxit, verb, 2, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left = 2, right = 1, m = 2........"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_z &
    ( options, 0, n, a, b, t, 2, 1, tol, maxit, verb, 2, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left = 10, right = 10, m = 10....."
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  do i = 2, n
    x(i, 1) = ZERO
    x(1, i ) = ZERO
  end do
  call run_ssmfe_z &
    ( options, 0, n, a, b, t, 10, 10, tol, maxit, verb, 10, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing largest = 1, m = 2................"
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_largest_z &
    ( options, 0, n, a, b, t, 1, tol, maxit, verb, 2, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing largest = 10, m = 5..............."
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_largest_z &
    ( options, 0, n, a, b, t, 10, tol, maxit, verb, 5, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left = 10, right = 10, m = 10....."
  forall ( i = 1 : n ) b(i, i) = 2.0D0**i
  do i = 1, n
    do j = 1, n
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_ssmfe_z &
    ( options, 0, n, a, b, t, 10, 10, tol, maxit, verb, 10, &
      n, lambda, x, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  deallocate ( a, b, t, x, lambda, rr, ind )

end subroutine test_core_misc_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_expert

  write(*,"(/a)") "====================="
  write(*,"(a)")  "Testing ssmfe_expert:"
  write(*,"(a)")  "====================="

  call test_expert_errors_d
  call test_expert_errors_z
  call test_expert_options_d
  call test_expert_options_z

end subroutine test_expert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_expert_errors_d

  use spral_ssmfe_expert

  integer :: n
  integer :: m
  integer :: nep
  integer :: mep
  integer :: left, right

  integer, allocatable :: ind(:)

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: rr(:,:,:)

  type(ssmfe_rcid       ) :: rci
  type(ssmfe_options    ) :: options
  type(ssmfe_expert_keep) :: keep
  type(ssmfe_inform     ) :: inform

  write(*,"(/a)") "======================"
  write(*,"(a)")  "Testing errors (real):"
  write(*,"(a)")  "======================"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n = 6
  m = 2
  nep = 1
  mep = n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )
  allocate ( rr(m + m, m + m, 3), ind(m) )

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing block_size < 1...................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, -1, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left > max_left..................."
  options%max_left = 0
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )
  options%max_left = -1

  write(*,"(a)",advance="no") " * Testing mep < left........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, 0, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing wrong err_est....................."
  options%err_est = 0
  rci%job = 0
  call ssmfe_generalized &
    ( rci, nep, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_ERR_EST )
  call ssmfe_free( keep, inform )
  options%err_est = 2

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = ZERO
  left = 1
  right = 1

  write(*,"(a)",advance="no") " * Testing block_size < 2...................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, -1, left, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing right < 0........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, -1, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_RIGHT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing right > max_right................."
  options%max_right = 0
  rci%job = 0
  call ssmfe_generalized_shift &
    ( rci, sigma, left, right, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_RIGHT )
  call ssmfe_free( keep, inform )
  options%max_right = -1

  write(*,"(a)",advance="no") " * Testing mep < left + right................"
  rci%job = 0
  call ssmfe_generalized_shift &
    ( rci, sigma, left, right, 0, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

  write ( *, '(a)' ) ' * Testing buckling...'

  write(*,"(a)",advance="no") " * Testing zero sigma........................"
  sigma = ZERO
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result(inform%flag, WRONG_SIGMA)
  call ssmfe_free( keep, inform )
  sigma = ONE

  write(*,"(a)",advance="no") " * Testing wrong minAprod...................."
  options%minAprod = .false.
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result(inform%flag, WRONG_MINPROD)
  call ssmfe_free( keep, inform )
  options%minAprod = .true.

  deallocate ( a, b, t, x, lambda, rr, ind )

end subroutine test_expert_errors_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_expert_errors_z

  use spral_ssmfe_expert

  integer :: n
  integer :: m
  integer :: nep
  integer :: mep
  integer :: left, right

  integer, allocatable :: ind(:)

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)
  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)
  complex(wp), allocatable :: rr(:,:,:)

  type(ssmfe_rciz       ) :: rci
  type(ssmfe_options    ) :: options
  type(ssmfe_expert_keep) :: keep
  type(ssmfe_inform     ) :: inform

  write(*,"(/a)") "========================="
  write(*,"(a)")  "Testing errors (complex):"
  write(*,"(a)")  "========================="

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n = 6
  m = 2
  nep = 1
  mep = n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )
  allocate ( rr(m + m, m + m, 3), ind(m) )

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing block_size < 1...................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, -1, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left > max_left..................."
  options%max_left = 0
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )
  options%max_left = -1

  write(*,"(a)",advance="no") " * Testing mep < left........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, 0, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing wrong err_est....................."
  options%err_est = 0
  rci%job = 0
  call ssmfe_generalized &
    ( rci, nep, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_ERR_EST )
  call ssmfe_free( keep, inform )
  options%err_est = 2

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = ZERO
  left = 1
  right = 1

  write(*,"(a)",advance="no") " * Testing block_size < 2...................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_BLOCK_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, -1, right, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing right < 0........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, -1, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_RIGHT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing right > max_right................."
  options%max_right = 0
  rci%job = 0
  call ssmfe_generalized_shift &
    ( rci, sigma, left, right, mep, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_RIGHT )
  call ssmfe_free( keep, inform )
  options%max_right = -1

  write(*,"(a)",advance="no") " * Testing mep < left + right................"
  rci%job = 0
  call ssmfe_generalized_shift &
    ( rci, sigma, left, right, 0, lambda, 0, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

  write ( *, '(a)' ) ' * Testing buckling...'

  write(*,"(a)",advance="no") " * Testing zero sigma........................"
  sigma = ZERO
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_SIGMA )
  call ssmfe_free( keep, inform )
  sigma = ONE

  write(*,"(a)",advance="no") " * Testing wrong minAprod...................."
  options%minAprod = .false.
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, m, rr, ind, keep, options, inform )
  call print_result( inform%flag, WRONG_MINPROD )
  call ssmfe_free( keep, inform )
  options%minAprod = .true.

  deallocate ( a, b, t, x, lambda, rr, ind )

end subroutine test_expert_errors_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_expert_options_d

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "======================"
  write(*,"(a)") "Testing options (real)"
  write(*,"(a)") "======================"

  n = 50
  nep = 5
  mep = n
  lwork = n*n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  forall ( i = 1 : n ) a(i, i) = i*10

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing residual-based error estimates...."
  options%err_est = 1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%err_est = 2

  write(*,"(a)",advance="no") " * Testing minAprod = .false. ..............."
  options%minAprod = .false.
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%minAprod = .true.

  write(*,"(a)",advance="no") " * Testing minBprod = .false. ..............."
  options%minBprod = .false.
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%minBprod = .true.

  write(*,"(a)",advance="no") " * Testing minAprod and minBprod = .false. .."
  options%minAprod = .false.
  options%minBprod = .false.
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%minAprod = .true.
  options%minBprod = .true.

  write(*,"(a)",advance="no") " * Testing zero extra_left..................."
  options%extra_left = 0
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1

  write(*,"(a)",advance="no") " * Testing nonzero extra_left................"
  options%extra_left = 5
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1

  write ( *, '(a)' ) ' * Testing shift-invert...'
  left = 1
  right = 1

  write(*,"(a)",advance="no") " * Testing zero extras......................."
  sigma = 355
  options%extra_left = 0
  options%extra_right = 0
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1
  options%extra_right = -1

  write(*,"(a)",advance="no") " * Testing nonzero extras...................."
  sigma = 355
  options%extra_left = 5
  options%extra_right = 5
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1
  options%extra_right = -1

  write(*,"(a)",advance="no") " * Testing residual-based error estimates...."
  sigma = 355
  options%err_est = 1
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%err_est = 2

  deallocate ( a, b, t, x, w, lambda, ipiv ) !, rr, ind )

end subroutine test_expert_options_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_expert_options_z

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)
  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)
  complex(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "========================="
  write(*,"(a)") "Testing options (complex)"
  write(*,"(a)") "========================="

  n = 50
  nep = 5
  mep = n
  lwork = n*n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  forall ( i = 1 : n ) a(i, i) = i*10

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing zero extra_left..................."
  options%extra_left = 0
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1

  write(*,"(a)",advance="no") " * Testing nonzero extra_left................"
  options%extra_left = 5
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1

  write(*,"(a)",advance="no") " * Testing residual-based error estimates...."
  options%err_est = 1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%err_est = 2

  write(*,"(a)",advance="no") " * Testing minAprod = .false. ..............."
  options%minAprod = .false.
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%minAprod = .true.

  write(*,"(a)",advance="no") " * Testing minBprod = .false. ..............."
  options%minBprod = .false.
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%minBprod = .true.

  write(*,"(a)",advance="no") " * Testing minAprod and minBprod = .false. .."
  options%minAprod = .false.
  options%minBprod = .false.
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%minAprod = .true.
  options%minBprod = .true.

  write ( *, '(a)' ) ' * Testing shift-invert...'
  left = 1
  right = 1

  write(*,"(a)",advance="no") " * Testing zero extras......................."
  sigma = 355
  options%extra_left = 0
  options%extra_right = 0
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1
  options%extra_right = -1

  write(*,"(a)",advance="no") " * Testing nonzero extras...................."
  sigma = 355
  options%extra_left = 5
  options%extra_right = 5
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%extra_left = -1
  options%extra_right = -1

  write(*,"(a)",advance="no") " * Testing residual-based error estimates...."
  sigma = 355
  options%err_est = 1
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%err_est = 2

  deallocate ( a, b, t, x, w, lambda, ipiv ) !, rr, ind )

end subroutine test_expert_options_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe_errors_d

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: ldx
  integer :: i

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)

  type(ssmfe_rcid   ) :: rci
  type(ssmfe_options) :: options
  type(ssmfe_keepd  ) :: keep
  type(ssmfe_inform ) :: inform

  write(*,"(/a)") "======================"
  write(*,"(a)")  "Testing errors (real):"
  write(*,"(a)")  "======================"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  n = 6
  nep = 1
  mep = n
  ldx = n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_RCI_JOB )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_PROBLEM_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, 1, keep, options, inform )
  call print_result( inform%flag, WRONG_LDX )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, -1, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left > n/2........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, n, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing mep < left........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, 0, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

  a = ZERO
  t = ZERO
  forall ( i = 1 : n ) a(i, i) = ONE*i
  forall ( i = 1 : n ) t(i, i) = ONE

!if ( .false. ) then
  write(*,"(a)",advance="no") " * Testing bad initial x....................."
  x = ZERO
  options%user_x = 1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, INDEFINITE_B_OR_XBX )
  options%user_x = 0
  call ssmfe_free( inform )
!end if

!if ( .false. ) then
  write(*,"(a)",advance="no") " * Testing indefinite b......................"
  b = ZERO
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, INDEFINITE_B_OR_XBX )
  call ssmfe_free( inform )
!end if

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = ZERO
  left = 1
  right = 1

  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_RCI_JOB )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_PROBLEM_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, 1, keep, options, inform )
  call print_result( inform%flag, WRONG_LDX )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, -1, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing right < 0........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, -1, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_RIGHT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > n/2................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, 2, 2, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > mep................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, 1, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

  write ( *, '(a)' ) ' * Testing buckling...'
  sigma = ZERO
  write(*,"(a)",advance="no") " * Testing zero sigma........................"
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_SIGMA )
  call ssmfe_free( keep, inform )

  deallocate ( a, b, t, x, lambda )

end subroutine test_ssmfe_errors_d

subroutine test_ssmfe_errors_z

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: ldx
  integer :: i

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)
  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)

  type(ssmfe_rciz   ) :: rci
  type(ssmfe_options) :: options
  type(ssmfe_keepz  ) :: keep
  type(ssmfe_inform ) :: inform

  write(*,"(/a)") "========================="
  write(*,"(a)")  "Testing errors (complex):"
  write(*,"(a)")  "========================="

  n = 6
  nep = 1
  mep = n
  ldx = n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), lambda(n) )

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_RCI_JOB )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_PROBLEM_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, 1, keep, options, inform )
  call print_result( inform%flag, WRONG_LDX )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, -1, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left > n/2........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, n, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing mep < left........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, 0, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

  a = ZERO
  t = ZERO
  forall ( i = 1 : n ) a(i, i) = ONE*i
  forall ( i = 1 : n ) t(i, i) = ONE
!if ( .false. ) then
  write(*,"(a)",advance="no") " * Testing bad initial x....................."
  x = ZERO
  options%user_x = 1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, INDEFINITE_B_OR_XBX )
  options%user_x = 0
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing indefinite b......................"
  b = ZERO
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, INDEFINITE_B_OR_XBX )
  call ssmfe_free( inform )
!end if

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = ZERO
  left = 1
  right = 1

  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_RCI_JOB )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_PROBLEM_SIZE )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, 1, keep, options, inform )
  call print_result( inform%flag, WRONG_LDX )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, -1, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing right < 0........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, -1, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_RIGHT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > n/2................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, 2, 2, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_free( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > mep................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, 1, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_free( keep, inform )

!if ( .false. ) then
  write ( *, '(a)' ) ' * Testing buckling...'
  sigma = ZERO
  write(*,"(a)",advance="no") " * Testing zero sigma........................"
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_SIGMA )
  call ssmfe_free( keep, inform )
!end if

  deallocate ( a, b, t, x, lambda )

end subroutine test_ssmfe_errors_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe_warnings_d

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i, j

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "======================="
  write(*,"(a)") "Testing warnings (real)"
  write(*,"(a)") "======================="

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  n = 50
  nep = 2
  mep = n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  forall ( i = 1 : n/2 ) a(i, i) = ONE
  forall ( i = n/2 + 1 : n ) a(i, i) = i

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  write(*,"(a)",advance="no") " * Testing warning flag 1...................."
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, NO_SEARCH_DIRECTIONS_LEFT )
  call ssmfe_free( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, MAX_ITERATIONS_EXCEEDED )
  call ssmfe_free( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = -2
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, OUT_OF_STORAGE )
  call ssmfe_free( inform )
  mep = n
  options%left_gap = 0

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = 2
  left = 1
  right = 1
  lwork = n*n

  write(*,"(a)",advance="no") " * Testing warning flag 1...................."
  b(1,1) = -0.01
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, NO_SEARCH_DIRECTIONS_LEFT )
  call ssmfe_free( inform )
  b(1,1) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, MAX_ITERATIONS_EXCEEDED )
  call ssmfe_free( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = -2
  options%right_gap = -2
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, OUT_OF_STORAGE )
  call ssmfe_free( inform )

  deallocate ( a, b, t, x, w, lambda, ipiv )

end subroutine test_ssmfe_warnings_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe_warnings_z

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i, j

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma

  real(wp), allocatable :: lambda(:)

  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)
  complex(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "=========================="
  write(*,"(a)") "Testing warnings (complex)"
  write(*,"(a)") "=========================="

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  n = 50
  nep = 2
  mep = n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  forall ( i = 1 : n/2 ) a(i, i) = ONE
  forall ( i = n/2 + 1 : n ) a(i, i) = i

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 1...................."
  t = ZERO
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, NO_SEARCH_DIRECTIONS_LEFT )
  call ssmfe_free( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, MAX_ITERATIONS_EXCEEDED )
  call ssmfe_free( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = -2
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, OUT_OF_STORAGE )
  call ssmfe_free( inform )
  mep = n
  options%left_gap = 0

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = 2
  left = 1
  right = 1
  lwork = n*n

  write(*,"(a)",advance="no") " * Testing warning flag 1...................."
  b(1,1) = -0.01
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, NO_SEARCH_DIRECTIONS_LEFT )
  call ssmfe_free( inform )
  b(1,1) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, MAX_ITERATIONS_EXCEEDED )
  call ssmfe_free( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = -2
  options%right_gap = -2
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, OUT_OF_STORAGE )
  call ssmfe_free( inform )

  deallocate ( a, b, t, x, w, lambda, ipiv )

end subroutine test_ssmfe_warnings_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe_options_d

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i, j

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma
  real(wp) :: eps

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "======================"
  write(*,"(a)") "Testing options (real)"
  write(*,"(a)") "======================"

  n = 50
  nep = 5
  mep = n
  lwork = n*n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  forall ( i = 1 : n ) a(i, i) = i*10

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing printing suppresssion............."
  options%print_level = -1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing units suppression........"
  options%print_level = 3
  options%unit_error = -1
  options%unit_warning = -1
  options%unit_diagnostic = -1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing error printing suppresssion......."
  options%print_level = -1
  n = 0
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, WRONG_PROBLEM_SIZE )
  call ssmfe_free( inform )
  n = 50

  write(*,"(a)",advance="no") " * Testing printing level 0.................."
  options%print_level = 0
  write ( dl_unit, '(/1x, a)' ) 'print level 0'
  t = ZERO
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, NO_SEARCH_DIRECTIONS_LEFT )
  call ssmfe_free( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing printing level 1.................."
  options%print_level = 1
  write ( dl_unit, '(/1x, a)' ) 'print level 1'
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing level 2.................."
  options%print_level = 2
  write ( dl_unit, '(/1x, a)' ) 'print level 2'
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing level 3.................."
  options%print_level = 3
  write ( dl_unit, '(/1x, a)' ) 'print level 3'
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%print_level = 1

  write(*,"(a)",advance="no") " * Testing zero tolerances..................."
  options%tol_x = 0
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1

  write(*,"(a)",advance="no") " * Testing all non-zero tolerances..........."
  options%tol_x = 1
  options%abs_tol_lambda = 1
  options%rel_tol_lambda = 1
  options%abs_tol_residual = 1
  options%rel_tol_residual = 1
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%abs_tol_lambda = 0
  options%rel_tol_lambda = 0
  options%abs_tol_residual = 0
  options%rel_tol_residual = 0
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  nep = 3
  options%left_gap = 1
  options%tol_x = 1e-12
  forall ( i = 5 : n ) a(i, i) = i*10 + 2*eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%left, 3, .true. )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%tol_x = -1
  options%user_x = 0

  eps = 1D-3
  nep = 3
  options%left_gap = 10 + eps
  forall ( i = 5 : n ) a(i, i) = i*10 + 2*eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%left, 4 )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%user_x = 0

  write(*,"(a)",advance="no") " * Testing max left eigenvalues.............."
  options%max_left = n
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%max_left = -1

  write ( *, '(a)' ) ' * Testing shift-invert...'

  write(*,"(a)",advance="no") " * Testing zero tolerances..................."
  left = 5
  right = 5
  sigma = 255
  options%print_level = 3
  options%tol_x = 0
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%user_x = 0
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing absolute residual tolerance......."
  options%tol_x = 0
  options%abs_tol_residual = 1
  sigma = 255
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%abs_tol_residual = 0

  write(*,"(a)",advance="no") " * Testing all non-zero tolerances..........."
  options%tol_x = 1
  options%abs_tol_lambda = 1
  options%rel_tol_lambda = 1
  options%abs_tol_residual = 1
  options%rel_tol_residual = 1
  sigma = 345
  left = 3
  right = 3
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%abs_tol_lambda = 0
  options%rel_tol_lambda = 0
  options%abs_tol_residual = 0
  options%rel_tol_residual = 0

  write(*,"(a)",advance="no") " * Testing printing level 2.................."
  options%print_level = 2
  write ( dl_unit, '(/1x, a)' ) 'print level 2'
  sigma = 345
  left = 3
  right = 3
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing level 3.................."
  options%print_level = 3
  write ( dl_unit, '(/1x, a)' ) 'print level 3'
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  sigma = 440
  left = 3
  right = 0
  options%tol_x = 1e-3
  options%left_gap = 1
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 11 ) a(i, i) = i*10
  forall ( i = n - 10 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%left, 3, .true. )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%tol_x = -1

  eps = 1D-3
  sigma = 440
  left = 3
  right = 0
  options%left_gap = 10 + eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 11 ) a(i, i) = i*10
  forall ( i = n - 10 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%left, 4 )
  call ssmfe_free( inform )
  options%left_gap = 0

  write(*,"(a)",advance="no") " * Testing minimal right eigenvalue gap......"
  eps = 1D-3
  sigma = 439
  left = 0
  right = 3
  options%tol_x = 1e-12
  options%right_gap = 1
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 3 ) a(i, i) = i*10
  forall ( i = n - 2 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%right, 3, .true. )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%right_gap = 0
  options%print_level = -1

  eps = 1D-3
  sigma = 439
  left = 0
  right = 3
  options%right_gap = 10 + eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 3 ) a(i, i) = i*10
  forall ( i = n - 2 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%right, 4 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%right_gap = 0
  options%print_level = -1

  deallocate ( a, b, t, x, w, lambda, ipiv )

end subroutine test_ssmfe_options_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe_options_z

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i, j

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma
  real(wp) :: eps

  real(wp), allocatable :: lambda(:)
  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)
  complex(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "========================="
  write(*,"(a)") "Testing options (complex)"
  write(*,"(a)") "========================="

  n = 50
  nep = 5
  mep = n
  lwork = n*n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  forall ( i = 1 : n ) a(i, i) = i*10

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing error printing suppresssion......."
  options%print_level = -1
  n = 0
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, WRONG_PROBLEM_SIZE )
  call ssmfe_free( inform )
  n = 50

  write(*,"(a)",advance="no") " * Testing printing suppresssion............."
  options%print_level = 0
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing units suppression........"
  options%print_level = 3
  options%unit_error = -1
  options%unit_warning = -1
  options%unit_diagnostic = -1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing printing level 0.................."
  options%print_level = 0
  write ( dl_unit, '(/1x, a)' ) 'print level 0'
  t = ZERO
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, NO_SEARCH_DIRECTIONS_LEFT )
  call ssmfe_free( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing printing level 1.................."
  options%print_level = 1
  write ( dl_unit, '(/1x, a)' ) 'print level 1'
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing level 2.................."
  options%print_level = 2
  write ( dl_unit, '(/1x, a)' ) 'print level 2'
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing level 3.................."
  options%print_level = 3
  write ( dl_unit, '(/1x, a)' ) 'print level 3'
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing zero tolerances..................."
  options%tol_x = 0
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1

  write(*,"(a)",advance="no") " * Testing all non-zero tolerances..........."
  options%tol_x = 1
  options%abs_tol_lambda = 1
  options%rel_tol_lambda = 1
  options%abs_tol_residual = 1
  options%rel_tol_residual = 1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%abs_tol_lambda = 0
  options%rel_tol_lambda = 0
  options%abs_tol_residual = 0
  options%rel_tol_residual = 0

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  nep = 3
  options%left_gap = 1
  options%tol_x = 1e-12
  forall ( i = 5 : n ) a(i, i) = i*10 + 2*eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%left, 3, .true. )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%tol_x = -1
  options%user_x = 0

  eps = 1D-3
  nep = 3
  options%left_gap = 10 + eps
  forall ( i = 5 : n ) a(i, i) = i*10 + 2*eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%left, 4 )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%user_x = 0

  write(*,"(a)",advance="no") " * Testing max left eigenvalues.............."
  options%max_left = n
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%max_left = -1

  write ( *, '(a)' ) ' * Testing shift-invert...'

  write(*,"(a)",advance="no") " * Testing zero tolerances..................."
  sigma = 255
  left = 2
  right = 7
  options%tol_x = 0
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%user_x = 0
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing absolute residual tolerance......."
  options%tol_x = 0
  options%abs_tol_residual = 1
  sigma = 255
  left = 3
  right = 3
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%abs_tol_residual = 0

  write(*,"(a)",advance="no") " * Testing all non-zero tolerances..........."
  options%tol_x = 1
  options%abs_tol_lambda = 1
  options%rel_tol_lambda = 1
  options%abs_tol_residual = 1
  options%rel_tol_residual = 1
  sigma = 345
  left = 3
  right = 3
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%abs_tol_lambda = 0
  options%rel_tol_lambda = 0
  options%abs_tol_residual = 0
  options%rel_tol_residual = 0

  write(*,"(a)",advance="no") " * Testing printing level 2.................."
  options%print_level = 2
  write ( dl_unit, '(/1x, a)' ) 'print level 2'
  sigma = 345
  left = 3
  right = 3
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing printing level 3.................."
  options%print_level = 3
  write ( dl_unit, '(/1x, a)' ) 'print level 3'
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  sigma = 440
  left = 3
  right = 0
  options%tol_x = 1e-3
  options%left_gap = 1
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 11 ) a(i, i) = i*10
  forall ( i = n - 10 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%left, 3, .true. )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%tol_x = -1

  eps = 1D-3
  sigma = 440
  left = 3
  right = 0
  options%left_gap = 10 + eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 11 ) a(i, i) = i*10
  forall ( i = n - 10 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%left, 4 )
  call ssmfe_free( inform )
  options%left_gap = 0

  write(*,"(a)",advance="no") " * Testing minimal right eigenvalue gap......"
  eps = 1D-3
  sigma = 439
  left = 0
  right = 3
  options%tol_x = 1e-12
  options%right_gap = 1
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 3 ) a(i, i) = i*10
  forall ( i = n - 2 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%right, 3,  .true. )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%right_gap = 0
  options%print_level = -1

  eps = 1D-3
  sigma = 439
  left = 0
  right = 3
  options%right_gap = 10 + eps
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n - 3 ) a(i, i) = i*10
  forall ( i = n - 2 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%right, 4 )
  call ssmfe_free( inform )
  options%tol_x = -1
  options%right_gap = 0
  options%print_level = -1

  deallocate ( a, b, t, x, w, lambda, ipiv )

end subroutine test_ssmfe_options_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe_misc_d

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i, j

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma
  real(wp) :: eps

  real(wp), allocatable :: lambda(:)
  real(wp), allocatable :: a(:,:)
  real(wp), allocatable :: b(:,:)
  real(wp), allocatable :: t(:,:)
  real(wp), allocatable :: x(:,:)
  real(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "=========================="
  write(*,"(a)") "Miscellaneous tests (real)"
  write(*,"(a)") "=========================="

  n = 50
  nep = 5
  mep = n
  lwork = n*n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  eps = 1D-3
  forall ( i = 1 : 10 ) a(i, i) = i*10
  forall ( i = 11 : n ) a(i, i) = i*10 + 2*eps

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing zero nep.........................."
  call run_std_d( n, a, t, 0, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing restart..........................."
!  eps = 1D-3
  nep = 3
  options%left_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 4 )
  call ssmfe_free( inform )
  options%max_left = -1
  options%left_gap = 0
  options%user_X = 0

!if ( .false. ) then
  write(*,"(a)",advance="no") " * Testing restart..........................."
  nep = 3
  options%max_left = 6
  options%left_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 4 )
  call ssmfe_free( inform )
  options%max_left = -1
  options%left_gap = 0
  options%user_X = 0
!end if

  write ( *, '(a)' ) ' * Testing shift-invert...'

  write(*,"(a)",advance="no") " * Testing zero nep.........................."
  sigma = 255
  call run_std_si_d &
    ( options, n, a, sigma, 0, 0, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, 0)
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing standard.........................."
  sigma = 255
  left = 5
  right = 5
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing generalized......................."
  sigma = 255
  left = 5
  right = 5
  options%user_x = 1
  x = ONE
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%user_x = 0

!if ( .false. ) then
  write(*,"(a)",advance="no") " * Testing restart..........................."
  sigma = 101
  left = 5
  right = 0
  options%left_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 4, .true. )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%user_x = 0

  sigma = 409
  left = 0
  right = 5
  options%right_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 4, .true. )
  call ssmfe_free( inform )
  options%right_gap = 0
  options%user_x = 0
!end if

  mep = 4
  sigma = 2
  options%left_gap = -2
  options%right_gap = -2
  options%user_x = mep
  forall ( i = 1 : n/2 ) a(i, i) = ONE
  forall ( i = n/2 + 1 : n ) a(i, i) = i
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  a = -a
  sigma = -sigma
  call run_std_si_d &
    ( options, n, a, sigma, 1, 1, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, OUT_OF_STORAGE )
  call ssmfe_free( inform )
  mep = n
  a = -a
  sigma = -sigma
  options%left_gap = 0
  options%right_gap = 0
  options%user_x = 0

  write ( *, '(a)' ) ' * Testing buckling...'

  write(*,"(a)",advance="no") " * Testing restart..........................."
  forall ( i = 1 : 30 ) a(i, i) = i*10
  forall ( i = 31 : n ) a(i, i) = i*10 + 2*eps
  a(n, n) = n*10 + 4*eps
  sigma = 399
  left = 6
  right = 5
  options%right_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_buckling_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 4 )
  call ssmfe_free( inform )
  options%right_gap = 0
  options%user_x = 0

  deallocate ( a, b, t, x, w, lambda, ipiv )

end subroutine test_ssmfe_misc_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_ssmfe_misc_z

  use spral_ssmfe

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: lwork
  integer :: i, j

  integer, allocatable :: ipiv(:)

  real(wp) :: sigma
  real(wp) :: eps

  real(wp), allocatable :: lambda(:)
  complex(wp), allocatable :: a(:,:)
  complex(wp), allocatable :: b(:,:)
  complex(wp), allocatable :: t(:,:)
  complex(wp), allocatable :: x(:,:)
  complex(wp), allocatable :: w(:,:)

  type(ssmfe_options) :: options
  type(ssmfe_inform ) :: inform

  write(*,"(a)")
  write(*,"(a)") "============================="
  write(*,"(a)") "Miscellaneous tests (complex)"
  write(*,"(a)") "============================="

  n = 50
  nep = 5
  mep = n
  lwork = n*n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )

  a = ZERO
  eps = 1D-3
  forall ( i = 1 : 10 ) a(i, i) = i*10
  forall ( i = 11 : n ) a(i, i) = i*10 + 2*eps

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  options%print_level = 0
  options%unit_error = we_unit
  options%unit_warning = we_unit
  options%unit_diagnostic = dl_unit

  write(*,"(a)",advance="no") " * Testing zero nep.........................."
  call run_std_z( n, a, t, 0, mep, lambda, x, options, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing restart..........................."
  eps = 1D-3
  nep = 3
  options%left_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 4 )
  call ssmfe_free( inform )
  options%max_left = -1
  options%left_gap = 0
  options%user_X = 0

!if ( .false. ) then
  write(*,"(a)",advance="no") " * Testing restart..........................."
  eps = 1D-3
  nep = 3
  options%max_left = 6
  options%left_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, 4 )
  call ssmfe_free( inform )
  options%max_left = -1
  options%left_gap = 0
  options%user_X = 0
!end if

  write ( *, '(a)' ) ' * Testing shift-invert...'

  write(*,"(a)",advance="no") " * Testing zero nep.........................."
  sigma = 255
  call run_std_si_z &
    ( options, n, a, sigma, 0, 0, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing standard.........................."
  sigma = 255
  left = 5
  right = 5
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )

  write(*,"(a)",advance="no") " * Testing generalized......................."
  sigma = 255
  left = 5
  right = 5
  options%user_x = 1
  x = ONE
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 0 )
  call ssmfe_free( inform )
  options%user_x = 0

  write(*,"(a)",advance="no") " * Testing restart..........................."
  sigma = 101
  left = 5
  right = 0
  options%left_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 4, .true. )
  call ssmfe_free( inform )
  options%left_gap = 0
  options%user_x = 0

  sigma = 409
  left = 0
  right = 5
  options%right_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 4, .true. )
  call ssmfe_free( inform )
  options%right_gap = 0
  options%user_x = 0

  mep = 4
  options%left_gap = -2
  options%right_gap = -2
  options%user_x = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  forall ( i = 1 : n/2 ) a(i, i) = ONE
  forall ( i = n/2 + 1 : n ) a(i, i) = i
  a = -a
  call run_std_si_z &
    ( options, n, a, -2.0D0, 1, 1, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, OUT_OF_STORAGE )
  call ssmfe_free( inform )
  mep = n
  a = -a
  options%left_gap = 0
  options%right_gap = 0
  options%user_x = 0

  write ( *, '(a)' ) ' * Testing buckling...'

  write(*,"(a)",advance="no") " * Testing restart..........................."
  forall ( i = 1 : 30 ) a(i, i) = i*10
  forall ( i = 31 : n ) a(i, i) = i*10 + 2*eps
  a(n, n) = n*10 + 4*eps
  sigma = 399
  left = 6
  right = 5
  options%left_gap = 10 + eps
  options%right_gap = 10 + eps
  options%user_X = mep
  do i = 1, n
    do j = 1, mep
      x(i, j) = sin(i*j*ONE)
    end do
  end do
  call run_buckling_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result( inform%flag, 4 )
  call ssmfe_free( inform )
  options%right_gap = 0
  options%user_x = 0

  deallocate ( a, b, t, x, w, lambda, ipiv )

end subroutine test_ssmfe_misc_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_d( n, a, t, nep, mep, lambda, x, options, inform )

  use spral_ssmfe

  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(inout) :: x(n, mep)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform

  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep

  logical :: restarted

  restarted = .false.
  rci%job = 0
  do
    call ssmfe_standard &
      ( rci, nep, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm( 'N', 'N', n, rci%nx, n, ONE, a, n, rci%x, n, ZERO, rci%y, n )
    case ( 2 )
      call dgemm( 'N', 'N', n, rci%nx, n, ONE, t, n, rci%x, n, ZERO, rci%y, n )
    case ( 999 )
      restarted = .true.
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_std_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_z( n, a, t, nep, mep, lambda, x, options, inform )

  use spral_ssmfe

  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(inout) :: x(n, mep)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform

  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep

  logical :: restarted

  restarted = .false.
  rci%job = 0
  do
    call ssmfe_standard &
      ( rci, nep, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, a, n, rci%x, n, ZZERO, rci%y, n )
    case ( 2 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, t, n, rci%x, n, ZZERO, rci%y, n )
    case ( 999 )
      restarted = .true.
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_std_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  use spral_ssmfe

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(inout) :: x(n, mep)
  real(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  integer, intent(in) :: lwork
  real(wp), intent(out) :: work(lwork)
  type(ssmfe_inform), intent(out) :: inform

  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep

  integer :: i

  ldlt = a
  forall ( i = 1 : n ) ldlt(i, i) = a(i, i) - sigma
  call dsytrf( 'L', n, ldlt, n, ipiv, work, lwork, i )
  if ( i /= 0 ) print *, 'dsytrf info', i
  options%max_left = num_neg_D(n, ldlt, n, ipiv)
  options%max_right = n - options%max_left
  rci%job = 0
  do
    call ssmfe_standard_shift &
      ( rci, sigma, left, right, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, a, n, rci%x, n, ZERO, rci%y, n )
    case ( 9 )
      call dcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call dsytrs( 'L', n, rci%nx, ldlt, n, ipiv, rci%y, n, i )
      if ( i /= 0 ) print *, 'dsytrs info', i
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  options%max_left = -1
  options%max_right = -1

end subroutine run_std_si_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  use spral_ssmfe

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(inout) :: x(n, mep)
  complex(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  integer, intent(in) :: lwork
  complex(wp), intent(out) :: work(lwork)
  type(ssmfe_inform), intent(out) :: inform

  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep

  integer :: i

  ldlt = a
  forall ( i = 1 : n ) ldlt(i, i) = a(i, i) - sigma
  call zhetrf( 'L', n, ldlt, n, ipiv, work, lwork, i )
  if ( i /= 0 ) print *, 'zhetrf info', i
  options%max_left = num_neg_Dz(n, ldlt, n, ipiv)
  options%max_right = n - options%max_left
  rci%job = 0
  do
    call ssmfe_standard_shift &
      ( rci, sigma, left, right, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, a, n, rci%x, n, ZZERO, rci%y, n )
    case ( 9 )
      call zcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call zhetrs( 'L', n, rci%nx, ldlt, n, ipiv, rci%y, n, i )
      if ( i /= 0 ) print *, 'zhetrs info', i
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  options%max_left = -1
  options%max_right = -1

end subroutine run_std_si_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )

  use spral_ssmfe

  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(inout) :: x(n, mep)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform

  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep

  rci%job = 0
  do
    call ssmfe_generalized &
      ( rci, nep, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm( 'N', 'N', n, rci%nx, n, ONE, a, n, rci%x, n, ZERO, rci%y, n )
    case ( 2 )
      call dgemm( 'N', 'N', n, rci%nx, n, ONE, t, n, rci%x, n, ZERO, rci%y, n )
    case ( 3 )
      call dgemm( 'N', 'N', n, rci%nx, n, ONE, b, n, rci%x, n, ZERO, rci%y, n )
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )

end subroutine run_gen_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )

  use spral_ssmfe

  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  complex(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(inout) :: x(n, mep)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform

  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep

  rci%job = 0
  do
    call ssmfe_generalized &
      ( rci, nep, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, a, n, rci%x, n, ZZERO, rci%y, n )
    case ( 2 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, t, n, rci%x, n, ZZERO, rci%y, n )
    case ( 3 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, b, n, rci%x, n, ZZERO, rci%y, n )
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )

end subroutine run_gen_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  use spral_ssmfe

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(inout) :: x(n, mep)
  real(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  integer, intent(in) :: lwork
  real(wp), intent(out) :: work(lwork)
  type(ssmfe_inform), intent(out) :: inform

  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep

  logical :: restarted

  integer :: i

  ldlt = a
  forall ( i = 1 : n ) ldlt(i, i) = a(i, i) - sigma*b(i, i)
  call dsytrf( 'L', n, ldlt, n, ipiv, work, lwork, i )
  if ( i /= 0 ) print *, 'dsytrf info', i
  options%max_left = num_neg_D(n, ldlt, n, ipiv)
  options%max_right = n - options%max_left
  restarted = .false.
  rci%job = 0
  do
    call ssmfe_generalized_shift &
      ( rci, sigma, left, right, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, a, n, rci%x, n, ZERO, rci%y, n )
    case ( 3 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, b, n, rci%x, n, ZERO, rci%y, n )
    case ( 9 )
      call dcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call dsytrs( 'L', n, rci%nx, ldlt, n, ipiv, rci%y, n, i )
      if ( i /= 0 ) print *, 'dsytrs info', i
    case ( 999 )
      restarted = .true.
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_gen_si_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  use spral_ssmfe

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(inout) :: x(n, mep)
  complex(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  integer, intent(in) :: lwork
  complex(wp), intent(out) :: work(lwork)
  type(ssmfe_inform), intent(out) :: inform

  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep

  logical :: restarted

  integer :: i

  ldlt = a
  forall ( i = 1 : n ) ldlt(i, i) = a(i, i) - sigma*b(i, i)
  call zhetrf( 'L', n, ldlt, n, ipiv, work, lwork, i )
  if ( i /= 0 ) print *, 'zhetrf info', i
  options%max_left = num_neg_Dz(n, ldlt, n, ipiv)
  options%max_right = n - options%max_left
  restarted = .false.
  rci%job = 0
  do
    call ssmfe_generalized_shift &
      ( rci, sigma, left, right, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, a, n, rci%x, n, ZZERO, rci%y, n )
    case ( 3 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, b, n, rci%x, n, ZZERO, rci%y, n )
    case ( 9 )
      call zcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call zhetrs( 'L', n, rci%nx, ldlt, n, ipiv, rci%y, n, i )
      if ( i /= 0 ) print *, 'dsytrs info', i
    case ( 999 )
      restarted = .true.
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_gen_si_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_buckling_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  use spral_ssmfe

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(inout) :: x(n, mep)
  real(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  integer, intent(in) :: lwork
  real(wp), intent(out) :: work(lwork)
  type(ssmfe_inform), intent(out) :: inform

  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep

  logical :: restarted

  integer :: i

  ldlt = a
  forall ( i = 1 : n ) ldlt(i, i) = a(i, i) - sigma*b(i, i)
  call dsytrf( 'L', n, ldlt, n, ipiv, work, lwork, i )
  if ( i /= 0 ) print *, 'dsytrf info', i
  options%max_left = num_neg_D(n, ldlt, n, ipiv)
  options%max_right = n - options%max_left
  restarted = .false.
  rci%job = 0
  do
    call ssmfe_buckling &
      ( rci, sigma, left, right, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, b, n, rci%x, n, ZERO, rci%y, n )
    case ( 3 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, a, n, rci%x, n, ZERO, rci%y, n )
    case ( 9 )
      call dcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call dsytrs( 'L', n, rci%nx, ldlt, n, ipiv, rci%y, n, i )
      if ( i /= 0 ) print *, 'dsytrs info', i
    case ( 999 )
      restarted = .true.
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag >= 0 ) inform%flag = 4

end subroutine run_buckling_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_buckling_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  use spral_ssmfe

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(inout) :: x(n, mep)
  complex(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  integer, intent(in) :: lwork
  complex(wp), intent(out) :: work(lwork)
  type(ssmfe_inform), intent(out) :: inform

  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep

  logical :: restarted

  integer :: i

  ldlt = a
  forall ( i = 1 : n ) ldlt(i, i) = a(i, i) - sigma*b(i, i)
  call zhetrf( 'L', n, ldlt, n, ipiv, work, lwork, i )
  if ( i /= 0 ) print *, 'zhetrf info', i
  options%max_left = num_neg_Dz(n, ldlt, n, ipiv)
  options%max_right = n - options%max_left
  restarted = .false.
  rci%job = 0
  do
    call ssmfe_buckling &
      ( rci, sigma, left, right, mep, lambda, n, x, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, b, n, rci%x, n, ZZERO, rci%y, n )
    case ( 3 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, a, n, rci%x, n, ZZERO, rci%y, n )
    case ( 9 )
      call zcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call zhetrs( 'L', n, rci%nx, ldlt, n, ipiv, rci%y, n, i )
      if ( i /= 0 ) print *, 'zhetrs info', i
    case ( 999 )
      restarted = .true.
    case ( :-1 )
      exit
    end select
  end do
  call ssmfe_free ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag >= 0 ) inform%flag = 4

end subroutine run_buckling_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_ssmfe_d &
    ( options, problem, n, a, b, t, left, right, tol, maxit, verb, m, &
      mep, lambda, x, inform )

  use spral_ssmfe_core

  type(ssmfe_core_options), intent(in) :: options
  integer, intent(in) :: problem
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: t(n, n)
  integer, intent(in) :: left
  integer, intent(in) :: right
  real(wp), intent(in) :: tol
  integer, intent(in) :: maxit
  integer, intent(in) :: verb
  integer, intent(in) :: m
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(inout) :: x(n, mep)
  type(ssmfe_inform), intent(out) :: inform

  character(7) :: word
  character(80) :: head, neck, line, form

  integer :: lcon
  integer :: rcon
  integer :: u_diag
  integer :: i, j, k

  integer, allocatable :: ind(:)

  real(wp) :: s
  real(wp) :: dnrm2, ddot

  real(wp), allocatable :: lmd(:)
  real(wp), allocatable :: w(:, :, :)
  real(wp), allocatable :: rr(:, :, :)
  real(wp), allocatable :: u(:, :)
  real(wp), allocatable :: bx(:, :)

  type(ssmfe_rcid     ) :: rci
  type(ssmfe_core_keep) :: keep

  head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
  neck = &
 '                       |           |           |   errors   |    errors'
  line = &
 '-------------------------------------------------------------------------'
  form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
  u_diag = dl_unit

!  print *, n, m, mep

  allocate ( ind(m), lmd(m), w(n, m, 0:7), rr(m + m, m + m, 3), u(mep, m) )
!  print *, 'ok'
  if ( problem /= 0 ) allocate ( bx(n, mep) )
!  print *, 'ok'

  call dcopy( n*m, x, 1, w, 1 )
!  print *, 'ok'

  lcon = 0
  rcon = 0
  rci%job = 0
  do
    call ssmfe &
      ( rci, problem, left, right, m, lmd, rr, ind, keep, options, inform )
!    print *, 'job', rci%job
    select case ( rci%job )
    case ( 1 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, a, n, w(1, rci%jx, rci%kx), n, &
          ZERO, w(1, rci%jy, rci%ky), n )
    case ( 2 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, t, n, w(1, rci%jx, rci%kx), n, &
          ZERO, w(1, rci%jy, rci%ky), n )
    case ( 3 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, b, n, w(1, rci%jx, rci%kx), n, &
          ZERO, w(1, rci%jy, rci%ky), n )
    case ( 4 )
      do j = 1, m
        if ( inform%converged(j) /= 0 ) then
          cycle
        else if ( inform%err_x(j) >= 0 .and. inform%err_x(j) <= tol ) then
          inform%converged(j) = 1
        end if
      end do
    case ( 5 )
      do i = 1, rci%nx
        k = (i - 1)*rci%i
        if ( rci%i > 0 ) then
          j = lcon + i
        else
          j = mep - rcon - i + 1
        end if
        lambda(j) = lmd(rci%jx + k)
        call dcopy( n, w(1, rci%jx + k, 0), 1, x(1, j), 1 )
        if ( problem /= 0 ) &
          call dcopy( n, w(1, rci%jy + k, rci%ky), 1, bx(1, j), 1 )
      end do
      if ( rci%i > 0 ) then
        lcon = lcon + rci%nx
      else
        rcon = rcon + rci%nx
      end if
      if ( rci%i < 0 ) then
        if ( verb > 0 ) then
          write( u_diag, '(/a, i5, a, i4, i4)' ) &
            'Iteration: ', inform%iteration, ', not converged: ', &
            max(0, left - lcon), max(0, right - rcon)
          write( u_diag, '(a/a)' ) &
            trim(line), trim(head)
          write( u_diag, '(a)' ) trim(neck)
          write( u_diag, '(a)' ) trim(line)
          do i = 1, m
            if ( inform%converged(i) /= 0 ) then
              word = '   yes'
            else
              word = '    no'
            end if
            write( u_diag, form ) &
              lmd(i), ' |', word, ' |', inform%residual_norms(i), '  |', &
              inform%err_lambda(m + i), '  |', inform%err_X(m + i)
          end do ! i = 1, m
          write( u_diag, '(a)' ) trim(line)
        end if
        if ( lcon >= left .and. rcon >= right .or. inform%iteration > maxit ) &
          exit
      end if
    case ( 11 )
      if ( rci%i == 0 ) then
        call dcopy &
          ( n * rci%nx, w(1, rci%jx, rci%kx), 1, w(1, rci%jy, rci%ky), 1 )
      else
        do i = 1, n
          do j = 1, rci%nx
            u(i, j) = w(i, ind(j), rci%kx)
          end do
          do j = 1, rci%nx
            w(i, j, rci%kx) = u(i, j)
          end do
          if ( rci%ky /= rci%kx ) then
            do j = 1, rci%nx
              u(i, j) = w(i, ind(j), rci%ky)
            end do
            do j = 1, rci%nx
              w(i, j, rci%ky) = u(i, j)
            end do
          end if
        end do
      end if
    case ( 12 )
      do i = 0, rci%nx - 1
        rr(rci%i + i, rci%j + i, rci%k) = &
          ddot(n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)
      end do
    case ( 13 )
      do i = 0, rci%nx - 1
        if ( rci%kx == rci%ky ) then
          s = dnrm2(n, w(1, rci%jx + i, rci%kx), 1)
          if ( s > 0 ) &
            call dscal( n, 1/s, w(1, rci%jx + i, rci%kx), 1 )
        else
          s = sqrt(abs(ddot &
            (n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)))
          if ( s > 0 ) then
            call dscal( n, 1/s, w(1, rci%jx + i, rci%kx), 1 )
            call dscal( n, 1/s, w(1, rci%jy + i, rci%ky), 1 )
          end if
        end if
      end do
    case ( 14 )
      do i = 0, rci%nx - 1
        s = -rr(rci%i + i, rci%j + i, rci%k)
        call daxpy&
          ( n, s, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1 )
      end do
    case ( 15 )
      if ( rci%nx > 0 .and. rci%ny > 0 ) &
        call dgemm &
          ( 'T', 'N', rci%nx, rci%ny, n, &
            rci%alpha, w(1, rci%jx, rci%kx), n, w(1, rci%jy, rci%ky), n, &
            rci%beta, rr(rci%i, rci%j, rci%k), 2*m )
    case ( 16, 17 )
      if ( rci%ny < 1 ) cycle
      if ( rci%job == 17 ) then
        call dgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            1.0D0, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            0.0D0, w(1, rci%jy, rci%ky), n )
        call dcopy &
          ( n * rci%ny, w(1, rci%jy, rci%ky), 1, w(1, rci%jx, rci%kx), 1 )
      else
        call dgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            rci%alpha, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            rci%beta, w(1, rci%jy, rci%ky), n )
      end if
    case ( 21 )
      if ( lcon > 0 ) then
        call dgemm &
          ( 'T', 'N', lcon, rci%nx, n, &
            ONE, x, n, w(1, rci%jy, rci%ky), n, ZERO, u, mep )
        call dgemm &
          ( 'N', 'N', n, rci%nx, lcon, -ONE, x, n, u, mep, &
            ONE, w(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call dgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ONE, bx, n, u, mep, ONE, w(1, rci%jy, rci%ky), n )
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call dgemm &
          ( 'T', 'N', rcon, rci%nx, n, &
            ONE, x(1, j), n, w(1, rci%jy, rci%ky), n, ZERO, u(j, 1), mep )
        call dgemm &
          ( 'N', 'N', n, rci%nx, rcon, &
            -ONE, x(1, j), n, u(j, 1), mep, ONE, w(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call dgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ONE, bx(1, j), n, u(j, 1), mep, ONE, w(1, rci%jy, rci%ky), n )
      end if
    case ( 22 )
      if ( lcon > 0 ) then
        call dgemm &
          ( 'T', 'N', lcon, rci%nx, n, &
            ONE, x, n, w(1, rci%jx, rci%kx), n, ZERO, u, mep )
        if ( problem /= 0 ) then
          call dgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ONE, bx, n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        else
          call dgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ONE, x, n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call dgemm &
          ( 'T', 'N', rcon, rci%nx, n, &
            ONE, x(1, j), n, w(1, rci%jx, rci%kx), n, ZERO, u, mep )
        if ( problem /= 0 ) then
          call dgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ONE, bx(1, j), n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        else
          call dgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ONE, x(1, j), n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
    case ( 999 )
      if ( rci%k == 0 ) then
        if ( rci%jx > 1 ) call random_number( w(:, 1 : rci%jx - 1, 0) )
        if ( rci%jx + rci%nx - 1 < m ) &
          call random_number( w(:, rci%jx + rci%nx : m, 0) )
      end if
    case ( :-1 )
      exit
    end select
  end do

  deallocate ( ind, lmd, w, rr, u )
  if ( problem /= 0 ) deallocate ( bx )

end subroutine run_ssmfe_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_ssmfe_z &
    ( options, problem, n, a, b, t, left, right, tol, maxit, verb, m, &
      mep, lambda, x, inform )

  use spral_ssmfe_core

  type(ssmfe_core_options), intent(in) :: options
  integer, intent(in) :: problem
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  complex(wp), intent(in) :: t(n, n)
  integer, intent(in) :: left
  integer, intent(in) :: right
  real(wp), intent(in) :: tol
  integer, intent(in) :: maxit
  integer, intent(in) :: verb
  integer, intent(in) :: m
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(inout) :: x(n, mep)
  type(ssmfe_inform), intent(out) :: inform

  character(7) :: word
  character(80) :: head, neck, line, form

  integer :: lcon
  integer :: rcon
  integer :: u_diag
  integer :: i, j, k

  integer, allocatable :: ind(:)

  real(wp) :: s
  real(wp) :: dznrm2

  real(wp), allocatable :: lmd(:)
  real(wp), allocatable :: dwork(:,:)

  complex(wp) :: z
  complex(wp) :: zdotc

  complex(wp), allocatable :: w(:, :, :)
  complex(wp), allocatable :: rr(:, :, :)
  complex(wp), allocatable :: u(:, :)
  complex(wp), allocatable :: bx(:, :)

  type(ssmfe_rciz     ) :: rci
  type(ssmfe_core_keep) :: keep

  head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
  neck = &
 '                       |           |           |   errors   |    errors'
  line = &
 '-------------------------------------------------------------------------'
  form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
  u_diag = dl_unit

  allocate ( ind(m), lmd(m), w(n, m, 0:7), rr(m + m, m + m, 3), u(mep, m) )
  if ( problem /= 0 ) allocate ( bx(n, mep) )

  call zcopy( n*m, x, 1, w, 1 )

  lcon = 0
  rcon = 0
  rci%job = 0
  do
    call ssmfe &
      ( rci, problem, left, right, m, lmd, rr, ind, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, a, n, w(1, rci%jx, rci%kx), n, &
          ZZERO, w(1, rci%jy, rci%ky), n )
    case ( 2 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, t, n, w(1, rci%jx, rci%kx), n, &
          ZZERO, w(1, rci%jy, rci%ky), n )
    case ( 3 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, b, n, w(1, rci%jx, rci%kx), n, &
          ZZERO, w(1, rci%jy, rci%ky), n )
    case ( 4 )
      do j = 1, m
        if ( inform%converged(j) /= 0 ) then
          cycle
        else if ( inform%err_x(j) >= 0 .and. inform%err_x(j) <= tol ) then
          inform%converged(j) = 1
        end if
      end do
    case ( 5 )
      do i = 1, rci%nx
        k = (i - 1)*rci%i
        if ( rci%i > 0 ) then
          j = lcon + i
        else
          j = mep - rcon - i + 1
        end if
        lambda(j) = lmd(rci%jx + k)
        call zcopy( n, w(1, rci%jx + k, 0), 1, x(1, j), 1 )
        if ( problem /= 0 ) &
          call zcopy( n, w(1, rci%jy + k, rci%ky), 1, bx(1, j), 1 )
      end do
      if ( rci%i > 0 ) then
        lcon = lcon + rci%nx
      else
        rcon = rcon + rci%nx
      end if
      if ( rci%i < 0 ) then
        if ( verb > 0 ) then
          write( u_diag, '(/a, i5, a, i4, i4)' ) &
            'Iteration: ', inform%iteration, ', not converged: ', &
            max(0, left - lcon), max(0, right - rcon)
          write( u_diag, '(a/a)' ) &
            trim(line), trim(head)
          write( u_diag, '(a)' ) trim(neck)
          write( u_diag, '(a)' ) trim(line)
          do i = 1, m
            if ( inform%converged(i) /= 0 ) then
              word = '   yes'
            else
              word = '    no'
            end if
            write( u_diag, form ) &
              lmd(i), ' |', word, ' |', inform%residual_norms(i), '  |', &
              inform%err_lambda(m + i), '  |', inform%err_X(m + i)
          end do ! i = 1, m
          write( u_diag, '(a)' ) trim(line)
        end if
        if ( lcon >= left .and. rcon >= right .or. inform%iteration > maxit ) &
          exit
      end if
    case ( 11 )
      if ( rci%i == 0 ) then
        call zcopy &
          ( n * rci%nx, w(1, rci%jx, rci%kx), 1, w(1, rci%jy, rci%ky), 1 )
      else
        do i = 1, n
          do j = 1, rci%nx
            u(i, j) = w(i, ind(j), rci%kx)
          end do
          do j = 1, rci%nx
            w(i, j, rci%kx) = u(i, j)
          end do
          if ( rci%ky /= rci%kx ) then
            do j = 1, rci%nx
              u(i, j) = w(i, ind(j), rci%ky)
            end do
            do j = 1, rci%nx
              w(i, j, rci%ky) = u(i, j)
            end do
          end if
        end do
      end if
    case ( 12 )
      do i = 0, rci%nx - 1
        rr(rci%i + i, rci%j + i, rci%k) = &
          zdotc(n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)
      end do
    case ( 13 )
      do i = 0, rci%nx - 1
        if ( rci%kx == rci%ky ) then
          s = dznrm2(n, w(1, rci%jx + i, rci%kx), 1)
          if ( s > 0 ) &
            call zscal( n, ZONE/s, w(1, rci%jx + i, rci%kx), 1 )
        else
          s = sqrt(abs(zdotc &
            (n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)))
          if ( s > 0 ) then
            call zscal( n, ZONE/s, w(1, rci%jx + i, rci%kx), 1 )
            call zscal( n, ZONE/s, w(1, rci%jy + i, rci%ky), 1 )
          end if
        end if
      end do
    case ( 14 )
      do i = 0, rci%nx - 1
        z = -rr(rci%i + i, rci%j + i, rci%k)
        call zaxpy&
          ( n, z, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1 )
      end do
    case ( 15 )
      if ( rci%nx > 0 .and. rci%ny > 0 ) &
        call zgemm &
          ( 'C', 'N', rci%nx, rci%ny, n, &
            rci%alpha, w(1, rci%jx, rci%kx), n, w(1, rci%jy, rci%ky), n, &
            rci%beta, rr(rci%i, rci%j, rci%k), 2*m )
    case ( 16, 17 )
      if ( rci%ny < 1 ) cycle
      if ( rci%job == 17 ) then
        call zgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            ZONE, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            ZZERO, w(1, rci%jy, rci%ky), n )
        call zcopy &
          ( n * rci%ny, w(1, rci%jy, rci%ky), 1, w(1, rci%jx, rci%kx), 1 )
      else
        call zgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            rci%alpha, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            rci%beta, w(1, rci%jy, rci%ky), n )
      end if
    case ( 21 )
      if ( lcon > 0 ) then
        call zgemm &
          ( 'C', 'N', lcon, rci%nx, n, &
            ZONE, x, n, w(1, rci%jy, rci%ky), n, ZZERO, u, mep )
        call zgemm &
          ( 'N', 'N', n, rci%nx, lcon, -ZONE, x, n, u, mep, &
            ZONE, w(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call zgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ZONE, bx, n, u, mep, ZONE, w(1, rci%jy, rci%ky), n )
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call zgemm &
          ( 'C', 'N', rcon, rci%nx, n, &
            ZONE, x(1, j), n, w(1, rci%jy, rci%ky), n, ZZERO, u(j, 1), mep )
        call zgemm &
          ( 'N', 'N', n, rci%nx, rcon, &
            -ZONE, x(1, j), n, u(j, 1), mep, ZONE, W(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call zgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ZONE, bx(1, j), n, u(j, 1), mep, ZONE, w(1, rci%jy, rci%ky), n )
      end if
    case ( 22 )
      if ( lcon > 0 ) then
        call zgemm &
          ( 'C', 'N', lcon, rci%nx, n, &
            ZONE, x, n, w(1, rci%jx, rci%kx), n, ZZERO, u, mep )
        if ( problem /= 0 ) then
          call zgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ZONE, bx, n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        else
          call zgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ZONE, x, n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call zgemm &
          ( 'C', 'N', rcon, rci%nx, n, &
            ZONE, x(1, j), n, w(1, rci%jx, rci%kx), n, ZZERO, u, mep )
        if ( problem /= 0 ) then
          call zgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ZONE, bx(1, j), n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        else
          call zgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ZONE, x(1, j), n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
    case ( 999 )
      if ( rci%k == 0 ) then
        allocate ( dwork(n, m), stat = inform%stat )
        if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
        if ( inform%stat /= 0 ) rci%job = -3
        if ( inform%stat /= 0 ) exit
        call random_number( dwork )
        if ( rci%jx > 1 ) &
          w(:, 1 : rci%jx - 1, 0) = 2*dwork(:, 1 : rci%jx - 1) - ONE
        if ( rci%jx + rci%nx - 1 < m ) &
          w(:, rci%jx + rci%nx : m, 0) = &
            2*dwork(:, rci%jx + rci%nx : m) - ONE
        deallocate ( dwork )
      end if
    case ( :-1 )
      exit
    end select
  end do

  deallocate ( ind, lmd, w, rr, u )
  if ( problem /= 0 ) deallocate ( bx )

end subroutine run_ssmfe_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_ssmfe_largest_d &
    ( options, problem, n, a, b, t, nep, tol, maxit, verb, m, &
      mep, lambda, x, inform )

  use spral_ssmfe_core

  type(ssmfe_core_options), intent(in) :: options
  integer, intent(in) :: problem
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  real(wp), intent(in) :: tol
  integer, intent(in) :: maxit
  integer, intent(in) :: verb
  integer, intent(in) :: m
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(inout) :: x(n, mep)
  type(ssmfe_inform), intent(out) :: inform

  character(7) :: word
  character(80) :: head, neck, line, form

  integer :: lcon
  integer :: rcon
  integer :: u_diag
  integer :: i, j, k

  integer, allocatable :: ind(:)

  real(wp) :: s
  real(wp) :: dnrm2, ddot

  real(wp), allocatable :: lmd(:)
  real(wp), allocatable :: w(:, :, :)
  real(wp), allocatable :: rr(:, :, :)
  real(wp), allocatable :: u(:, :)
  real(wp), allocatable :: bx(:, :)

  type(ssmfe_rcid     ) :: rci
  type(ssmfe_core_keep) :: keep

  head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
  neck = &
 '                       |           |           |   errors   |    errors'
  line = &
 '-------------------------------------------------------------------------'
  form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
  u_diag = dl_unit

  allocate ( ind(m), lmd(m), w(n, m, 0:7), rr(m + m, m + m, 3), u(mep, m) )
  if ( problem /= 0 ) allocate ( bx(n, mep) )

  call dcopy( n*m, x, 1, w, 1 )

  lcon = 0
  rcon = 0
  rci%job = 0
  do
    call ssmfe_largest &
      ( rci, problem, nep, m, lmd, rr, ind, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, a, n, w(1, rci%jx, rci%kx), n, &
          ZERO, w(1, rci%jy, rci%ky), n )
    case ( 2 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, t, n, w(1, rci%jx, rci%kx), n, &
          ZERO, w(1, rci%jy, rci%ky), n )
    case ( 3 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, b, n, w(1, rci%jx, rci%kx), n, &
          ZERO, w(1, rci%jy, rci%ky), n )
    case ( 4 )
      do j = 1, m
        if ( inform%converged(j) /= 0 ) then
          cycle
        else if ( inform%err_x(j) >= 0 .and. inform%err_x(j) <= tol ) then
          inform%converged(j) = 1
        end if
      end do
    case ( 5 )
      do i = 1, rci%nx
        k = (i - 1)*rci%i
        if ( rci%i > 0 ) then
          j = lcon + i
        else
          j = mep - rcon - i + 1
        end if
        lambda(j) = lmd(rci%jx + k)
        call dcopy( n, w(1, rci%jx + k, 0), 1, x(1, j), 1 )
        if ( problem /= 0 ) &
          call dcopy( n, w(1, rci%jy + k, rci%ky), 1, bx(1, j), 1 )
      end do
      if ( rci%i > 0 ) then
        lcon = lcon + rci%nx
      else
        rcon = rcon + rci%nx
      end if
      if ( rci%i < 0 ) then
        if ( verb > 0 ) then
          write( u_diag, '(/a, i5, a, i4, i4)' ) &
            'Iteration: ', inform%iteration, ', not converged: ', &
            max(0, nep - lcon - rcon)
          write( u_diag, '(a/a)' ) &
            trim(line), trim(head)
          write( u_diag, '(a)' ) trim(neck)
          write( u_diag, '(a)' ) trim(line)
          do i = 1, m
            if ( inform%converged(i) /= 0 ) then
              word = '   yes'
            else
              word = '    no'
            end if
            write( u_diag, form ) &
              lmd(i), ' |', word, ' |', inform%residual_norms(i), '  |', &
              inform%err_lambda(m + i), '  |', inform%err_X(m + i)
          end do ! i = 1, m
          write( u_diag, '(a)' ) trim(line)
        end if
        if ( lcon + rcon >= nep .or. inform%iteration > maxit ) &
          exit
      end if
    case ( 11 )
      if ( rci%i == 0 ) then
        call dcopy &
          ( n * rci%nx, w(1, rci%jx, rci%kx), 1, w(1, rci%jy, rci%ky), 1 )
      else
        do i = 1, n
          do j = 1, rci%nx
            u(i, j) = w(i, ind(j), rci%kx)
          end do
          do j = 1, rci%nx
            w(i, j, rci%kx) = u(i, j)
          end do
          if ( rci%ky /= rci%kx ) then
            do j = 1, rci%nx
              u(i, j) = w(i, ind(j), rci%ky)
            end do
            do j = 1, rci%nx
              w(i, j, rci%ky) = u(i, j)
            end do
          end if
        end do
      end if
    case ( 12 )
      do i = 0, rci%nx - 1
        rr(rci%i + i, rci%j + i, rci%k) = &
          ddot(n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)
      end do
    case ( 13 )
      do i = 0, rci%nx - 1
        if ( rci%kx == rci%ky ) then
          s = dnrm2(n, w(1, rci%jx + i, rci%kx), 1)
          if ( s > 0 ) &
            call dscal( n, 1/s, w(1, rci%jx + i, rci%kx), 1 )
        else
          s = sqrt(abs(ddot &
            (n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)))
          if ( s > 0 ) then
            call dscal( n, 1/s, w(1, rci%jx + i, rci%kx), 1 )
            call dscal( n, 1/s, w(1, rci%jy + i, rci%ky), 1 )
          end if
        end if
      end do
    case ( 14 )
      do i = 0, rci%nx - 1
        s = -rr(rci%i + i, rci%j + i, rci%k)
        call daxpy&
          ( n, s, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1 )
      end do
    case ( 15 )
      if ( rci%nx > 0 .and. rci%ny > 0 ) &
        call dgemm &
          ( 'T', 'N', rci%nx, rci%ny, n, &
            rci%alpha, w(1, rci%jx, rci%kx), n, w(1, rci%jy, rci%ky), n, &
            rci%beta, rr(rci%i, rci%j, rci%k), 2*m )
    case ( 16, 17 )
      if ( rci%ny < 1 ) cycle
      if ( rci%job == 17 ) then
        call dgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            1.0D0, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            0.0D0, w(1, rci%jy, rci%ky), n )
        call dcopy &
          ( n * rci%ny, w(1, rci%jy, rci%ky), 1, w(1, rci%jx, rci%kx), 1 )
      else
        call dgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            rci%alpha, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            rci%beta, w(1, rci%jy, rci%ky), n )
      end if
    case ( 21 )
      if ( lcon > 0 ) then
        call dgemm &
          ( 'T', 'N', lcon, rci%nx, n, &
            ONE, x, n, w(1, rci%jy, rci%ky), n, ZERO, u, mep )
        call dgemm &
          ( 'N', 'N', n, rci%nx, lcon, -ONE, x, n, u, mep, &
            ONE, w(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call dgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ONE, bx, n, u, mep, ONE, w(1, rci%jy, rci%ky), n )
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call dgemm &
          ( 'T', 'N', rcon, rci%nx, n, &
            ONE, x(1, j), n, w(1, rci%jy, rci%ky), n, ZERO, u(j, 1), mep )
        call dgemm &
          ( 'N', 'N', n, rci%nx, rcon, &
            -ONE, x(1, j), n, u(j, 1), mep, ONE, W(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call dgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ONE, bx(1, j), n, u(j, 1), mep, ONE, w(1, rci%jy, rci%ky), n )
      end if
    case ( 22 )
      if ( lcon > 0 ) then
        call dgemm &
          ( 'T', 'N', lcon, rci%nx, n, &
            ONE, x, n, w(1, rci%jx, rci%kx), n, ZERO, u, mep )
        if ( problem /= 0 ) then
          call dgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ONE, bx, n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        else
          call dgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ONE, x, n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call dgemm &
          ( 'T', 'N', rcon, rci%nx, n, &
            ONE, x(1, j), n, w(1, rci%jx, rci%kx), n, ZERO, u, mep )
        if ( problem /= 0 ) then
          call dgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ONE, bx(1, j), n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        else
          call dgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ONE, x(1, j), n, u, mep, ONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
    case ( 999 )
      if ( rci%k == 0 ) then
        if ( rci%jx > 1 ) call random_number( w(:, 1 : rci%jx - 1, 0) )
        if ( rci%jx + rci%nx - 1 < m ) &
          call random_number( w(:, rci%jx + rci%nx : m, 0) )
      end if
    case ( :-1 )
      exit
    end select
  end do

  deallocate ( ind, lmd, w, rr, u )
  if ( problem /= 0 ) deallocate ( bx )

end subroutine run_ssmfe_largest_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_ssmfe_largest_z &
    ( options, problem, n, a, b, t, nep, tol, maxit, verb, m, &
      mep, lambda, x, inform )

  use spral_ssmfe_core

  type(ssmfe_core_options), intent(in) :: options
  integer, intent(in) :: problem
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  complex(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  real(wp), intent(in) :: tol
  integer, intent(in) :: maxit
  integer, intent(in) :: verb
  integer, intent(in) :: m
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(inout) :: x(n, mep)
  type(ssmfe_inform), intent(out) :: inform

  character(7) :: word
  character(80) :: head, neck, line, form

  integer :: lcon
  integer :: rcon
  integer :: u_diag
  integer :: i, j, k

  integer, allocatable :: ind(:)

  real(wp) :: s
  real(wp) :: dznrm2

  real(wp), allocatable :: lmd(:)
  real(wp), allocatable :: dwork(:,:)

  complex(wp) :: z
  complex(wp) :: zdotc

  complex(wp), allocatable :: w(:, :, :)
  complex(wp), allocatable :: rr(:, :, :)
  complex(wp), allocatable :: u(:, :)
  complex(wp), allocatable :: bx(:, :)

  type(ssmfe_rciz     ) :: rci
  type(ssmfe_core_keep) :: keep

  head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
  neck = &
 '                       |           |           |   errors   |    errors'
  line = &
 '-------------------------------------------------------------------------'
  form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
  u_diag = dl_unit

  allocate ( ind(m), lmd(m), w(n, m, 0:7), rr(m + m, m + m, 3), u(mep, m) )
  if ( problem /= 0 ) allocate ( bx(n, mep) )

  call zcopy( n*m, x, 1, w, 1 )

  lcon = 0
  rcon = 0
  rci%job = 0
  do
    call ssmfe_largest &
      ( rci, problem, nep, m, lmd, rr, ind, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, a, n, w(1, rci%jx, rci%kx), n, &
          ZZERO, w(1, rci%jy, rci%ky), n )
    case ( 2 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, t, n, w(1, rci%jx, rci%kx), n, &
          ZZERO, w(1, rci%jy, rci%ky), n )
    case ( 3 )
      call zgemm &
        ( 'N', 'N', n, rci%nx, n, ZONE, b, n, w(1, rci%jx, rci%kx), n, &
          ZZERO, w(1, rci%jy, rci%ky), n )
    case ( 4 )
      do j = 1, m
        if ( inform%converged(j) /= 0 ) then
          cycle
        else if ( inform%err_x(j) >= 0 .and. inform%err_x(j) <= tol ) then
          inform%converged(j) = 1
        end if
      end do
    case ( 5 )
      do i = 1, rci%nx
        k = (i - 1)*rci%i
        if ( rci%i > 0 ) then
          j = lcon + i
        else
          j = mep - rcon - i + 1
        end if
        lambda(j) = lmd(rci%jx + k)
        call zcopy( n, w(1, rci%jx + k, 0), 1, x(1, j), 1 )
        if ( problem /= 0 ) &
          call zcopy( n, w(1, rci%jy + k, rci%ky), 1, bx(1, j), 1 )
      end do
      if ( rci%i > 0 ) then
        lcon = lcon + rci%nx
      else
        rcon = rcon + rci%nx
      end if
      if ( rci%i < 0 ) then
        if ( verb > 0 ) then
          write( u_diag, '(/a, i5, a, i4, i4)' ) &
            'Iteration: ', inform%iteration, ', not converged: ', &
            max(0, nep - lcon - rcon)
          write( u_diag, '(a/a)' ) &
            trim(line), trim(head)
          write( u_diag, '(a)' ) trim(neck)
          write( u_diag, '(a)' ) trim(line)
          do i = 1, m
            if ( inform%converged(i) /= 0 ) then
              word = '   yes'
            else
              word = '    no'
            end if
            write( u_diag, form ) &
              lmd(i), ' |', word, ' |', inform%residual_norms(i), '  |', &
              inform%err_lambda(m + i), '  |', inform%err_X(m + i)
          end do ! i = 1, m
          write( u_diag, '(a)' ) trim(line)
        end if
        if ( lcon + rcon >= nep .or. inform%iteration > maxit ) &
          exit
      end if
    case ( 11 )
      if ( rci%i == 0 ) then
        call zcopy &
          ( n * rci%nx, w(1, rci%jx, rci%kx), 1, w(1, rci%jy, rci%ky), 1 )
      else
        do i = 1, n
          do j = 1, rci%nx
            u(i, j) = w(i, ind(j), rci%kx)
          end do
          do j = 1, rci%nx
            w(i, j, rci%kx) = u(i, j)
          end do
          if ( rci%ky /= rci%kx ) then
            do j = 1, rci%nx
              u(i, j) = w(i, ind(j), rci%ky)
            end do
            do j = 1, rci%nx
              w(i, j, rci%ky) = u(i, j)
            end do
          end if
        end do
      end if
    case ( 12 )
      do i = 0, rci%nx - 1
        rr(rci%i + i, rci%j + i, rci%k) = &
          zdotc(n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)
      end do
    case ( 13 )
      do i = 0, rci%nx - 1
        if ( rci%kx == rci%ky ) then
          s = dznrm2(n, w(1, rci%jx + i, rci%kx), 1)
          if ( s > 0 ) &
            call zscal( n, ZONE/s, w(1, rci%jx + i, rci%kx), 1 )
        else
          s = sqrt(abs(zdotc &
            (n, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1)))
          if ( s > 0 ) then
            call zscal( n, ZONE/s, w(1, rci%jx + i, rci%kx), 1 )
            call zscal( n, ZONE/s, w(1, rci%jy + i, rci%ky), 1 )
          end if
        end if
      end do
    case ( 14 )
      do i = 0, rci%nx - 1
        z = -rr(rci%i + i, rci%j + i, rci%k)
        call zaxpy&
          ( n, z, w(1, rci%jx + i, rci%kx), 1, w(1, rci%jy + i, rci%ky), 1 )
      end do
    case ( 15 )
      if ( rci%nx > 0 .and. rci%ny > 0 ) &
        call zgemm &
          ( 'C', 'N', rci%nx, rci%ny, n, &
            rci%alpha, w(1, rci%jx, rci%kx), n, w(1, rci%jy, rci%ky), n, &
            rci%beta, rr(rci%i, rci%j, rci%k), 2*m )
    case ( 16, 17 )
      if ( rci%ny < 1 ) cycle
      if ( rci%job == 17 ) then
        call zgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            ZONE, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            ZZERO, w(1, rci%jy, rci%ky), n )
        call zcopy &
          ( n * rci%ny, w(1, rci%jy, rci%ky), 1, w(1, rci%jx, rci%kx), 1 )
      else
        call zgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            rci%alpha, w(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            rci%beta, w(1, rci%jy, rci%ky), n )
      end if
    case ( 21 )
      if ( lcon > 0 ) then
        call zgemm &
          ( 'C', 'N', lcon, rci%nx, n, &
            ZONE, x, n, w(1, rci%jy, rci%ky), n, ZZERO, u, mep )
        call zgemm &
          ( 'N', 'N', n, rci%nx, lcon, -ZONE, x, n, u, mep, &
            ZONE, w(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call zgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ZONE, bx, n, u, mep, ZONE, w(1, rci%jy, rci%ky), n )
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call zgemm &
          ( 'C', 'N', rcon, rci%nx, n, &
            ZONE, x(1, j), n, w(1, rci%jy, rci%ky), n, ZZERO, u(j, 1), mep )
        call zgemm &
          ( 'N', 'N', n, rci%nx, rcon, &
            -ZONE, x(1, j), n, u(j, 1), mep, ZONE, W(1, rci%jx, rci%kx), n )
        if ( problem /= 0 ) &
          call zgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ZONE, bx(1, j), n, u(j, 1), mep, ZONE, w(1, rci%jy, rci%ky), n )
      end if
    case ( 22 )
      if ( lcon > 0 ) then
        call zgemm &
          ( 'C', 'N', lcon, rci%nx, n, &
            ZONE, x, n, w(1, rci%jx, rci%kx), n, ZZERO, u, mep )
        if ( problem /= 0 ) then
          call zgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ZONE, bx, n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        else
          call zgemm &
            ( 'N', 'N', n, rci%nx, lcon, &
              -ZONE, x, n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
      if ( rcon > 0 ) then
        j = mep - rcon + 1
        call zgemm &
          ( 'C', 'N', rcon, rci%nx, n, &
            ZONE, x(1, j), n, w(1, rci%jx, rci%kx), n, ZZERO, u, mep )
        if ( problem /= 0 ) then
          call zgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ZONE, bx(1, j), n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        else
          call zgemm &
            ( 'N', 'N', n, rci%nx, rcon, &
              -ZONE, x(1, j), n, u, mep, ZONE, w(1, rci%jx, rci%kx), n )
        end if
      end if
    case ( 999 )
      if ( rci%k == 0 ) then
        allocate ( dwork(n, m), stat = inform%stat )
        if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
        if ( inform%stat /= 0 ) rci%job = -3
        if ( inform%stat /= 0 ) exit
        call random_number( dwork )
        if ( rci%jx > 1 ) &
          w(:, 1 : rci%jx - 1, 0) = 2*dwork(:, 1 : rci%jx - 1) - ONE
        if ( rci%jx + rci%nx - 1 < m ) &
          w(:, rci%jx + rci%nx : m, 0) = &
            2*dwork(:, rci%jx + rci%nx : m) - ONE
        deallocate ( dwork )
      end if
    case ( :-1 )
      exit
    end select
  end do

  deallocate ( ind, lmd, w, rr, u )
  if ( problem /= 0 ) deallocate ( bx )

end subroutine run_ssmfe_largest_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_result(actual, expected, continued)
   integer :: actual
   integer :: expected
   logical, optional :: continued

   logical :: mycontinued

   mycontinued = .false.
   if(present(continued)) mycontinued = continued

   if(actual.eq.expected) then
      if(mycontinued) then
         write(*,"(a)", advance="no") "ok..."
      else
         write(*,"(a)") "ok"
      endif
      return
   endif

   write(*,"(a)") "fail"
   write(*,"(2(a,i4))") "returned ", actual, ", expected ", expected
   errors = errors + 1
end subroutine print_result

integer function num_neg_D( n, LDLT, ld, ipiv )
! Counts negative eigenvalues of factor D

  double precision, parameter :: ZERO = 0.0D0

  integer, intent(in) :: n, ld
  real(wp), intent(in) :: LDLT(ld, n)
  integer, intent(in) :: ipiv(n)
  integer :: nneg
  integer :: i
  double precision :: r, s, t

  i = 1
  nneg = 0
  do while ( i <= n )
    s = LDLT(i, i)
    if ( ipiv(i) < 0 ) then
      t = LDLT(i + 1, i)
      r = LDLT(i + 1, i + 1)
      if ( s*r - t*t < ZERO ) then
        nneg = nneg + 1
      else if ( s*r - t*t > ZERO .and. s + r < ZERO ) then
        nneg = nneg + 2
      end if
      i = i + 2
    else
      if ( s < ZERO ) nneg = nneg + 1
      i = i + 1
    end if
  end do
  num_neg_D = nneg

end function num_neg_D

integer function num_neg_Dz( n, LDLT, ld, ipiv )
! Counts negative eigenvalues of factor D

  double precision, parameter :: ZERO = 0.0D0

  integer, intent(in) :: n, ld
  complex(wp), intent(in) :: LDLT(ld, n)
  integer, intent(in) :: ipiv(n)
  integer :: nneg
  integer :: i
  double precision :: r, s, t

  i = 1
  nneg = 0
  do while ( i <= n )
    s = real(LDLT(i, i))
    if ( ipiv(i) < 0 ) then
      t = abs(LDLT(i + 1, i))
      r = real(LDLT(i + 1, i + 1))
      if ( s*r - t*t < ZERO ) then
        nneg = nneg + 1
      else if ( s*r - t*t > ZERO .and. s + r < ZERO ) then
        nneg = nneg + 2
      end if
      i = i + 2
    else
      if ( s < ZERO ) nneg = nneg + 1
      i = i + 1
    end if
  end do
  num_neg_Dz = nneg

end function num_neg_Dz

end program


