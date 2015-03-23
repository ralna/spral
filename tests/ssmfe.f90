program main

   use spral_ssmfe

   implicit none

   integer, parameter :: wp = kind(0d0)
   real(wp), parameter :: ONE = 1.0_wp
   real(wp), parameter :: ZERO = 0.0_wp
   complex(wp), parameter :: ZZERO = (ZERO, ZERO)
   complex(wp), parameter :: ZONE = (ONE, ZERO)

   real(wp), parameter :: err_tol = 5e-11

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
   
   open ( 1, file = 'ssmfe_msg.txt' )
   call test_errors_d
   call test_errors_z
   call test_warnings_d
   call test_warnings_z
   close ( 1 )
   open ( 1, file = 'ssmfe_output.txt' )
   call test_options_d
   call test_options_z
   call test_misc_d
   close ( 1 )

   write(*, "(/a)") "============================="
   write(*, "(a,i4)") "Total number of errors = ", errors

   if(we_unit.gt.6) close(we_unit)
   if(dl_unit.gt.6) close(dl_unit)

   if(errors.ne.0) stop 1 ! ERROR CODE for make check script

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_errors_d

!  logical :: continued = .true.

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
  options%unit_error = 1
  options%unit_warning = 1
  options%unit_diagnostic = 1
  
  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_RCI_JOB)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_PROBLEM_SIZE)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, 1, keep, options, inform )
  call print_result(inform%flag, WRONG_LDX)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 1.........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, 0, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_LEFT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left > n/2........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, n, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_LEFT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing mep < left........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, 0, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_STORAGE_SIZE)
  call ssmfe_terminate( keep, inform )
  
  a = ZERO
  t = ZERO
  forall ( i = 1 : n ) a(i, i) = ONE*i
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing bad initial x....................."
  x = ZERO
  options%user_x = 1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, INDEFINITE_B_OR_XBX)
  options%user_x = 0
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing indefinite b......................"
  b = ZERO
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, INDEFINITE_B_OR_XBX)
  call ssmfe_terminate( inform )

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = ZERO
  left = 1
  right = 1
  
  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_RCI_JOB)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_PROBLEM_SIZE)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, 1, keep, options, inform )
  call print_result(inform%flag, WRONG_LDX)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, -1, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_LEFT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing right < 0........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, -1, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_RIGHT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > n/2................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, 2, 2, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_LEFT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > mep................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, 1, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_STORAGE_SIZE)
  call ssmfe_terminate( keep, inform )

  write ( *, '(a)' ) ' * Testing buckling...'
  sigma = ZERO
  write(*,"(a)",advance="no") " * Testing zero sigma........................"
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_SIGMA)
  call ssmfe_terminate( keep, inform )

  deallocate ( a, b, t, x, lambda )

end subroutine test_errors_d

subroutine test_errors_z

!  logical :: continued = .true.

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
  options%unit_error = 1
  options%unit_warning = 1
  options%unit_diagnostic = 1
  
  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_RCI_JOB )
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_PROBLEM_SIZE )
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, mep, lambda, n, x, 1, keep, options, inform )
  call print_result( inform%flag, WRONG_LDX )
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 1.........................."
  rci%job = 0
  call ssmfe_standard &
    ( rci, 0, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left > n/2........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, n, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_LEFT )
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing mep < left........................"
  rci%job = 0
  call ssmfe_standard &
    ( rci, nep, 0, lambda, n, x, ldx, keep, options, inform )
  call print_result( inform%flag, WRONG_STORAGE_SIZE )
  call ssmfe_terminate( keep, inform )
  
  a = ZERO
  t = ZERO
  forall ( i = 1 : n ) a(i, i) = ONE*i
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing bad initial x....................."
  x = ZERO
  options%user_x = 1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, INDEFINITE_B_OR_XBX )
  options%user_x = 0
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing indefinite b......................"
  b = ZERO
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result( inform%flag, INDEFINITE_B_OR_XBX )
  call ssmfe_terminate( inform )

  write ( *, '(a)' ) ' * Testing shift-invert...'
  sigma = ZERO
  left = 1
  right = 1
  
  write(*,"(a)",advance="no") " * Testing bad rci%job......................."
  rci%job = 1
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_RCI_JOB)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing n < 1............................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, 0, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_PROBLEM_SIZE)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing ldx < n..........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, mep, lambda, n, x, 1, keep, options, inform )
  call print_result(inform%flag, WRONG_LDX)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left < 0.........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, -1, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_LEFT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing right < 0........................."
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, -1, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_RIGHT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > n/2................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, 2, 2, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_LEFT)
  call ssmfe_terminate( keep, inform )

  write(*,"(a)",advance="no") " * Testing left + right > mep................"
  rci%job = 0
  call ssmfe_standard_shift &
    ( rci, sigma, left, right, 1, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_STORAGE_SIZE)
  call ssmfe_terminate( keep, inform )

  write ( *, '(a)' ) ' * Testing buckling...'
  sigma = ZERO
  write(*,"(a)",advance="no") " * Testing zero sigma........................"
  rci%job = 0
  call ssmfe_buckling &
    ( rci, sigma, left, right, mep, lambda, n, x, ldx, keep, options, inform )
  call print_result(inform%flag, WRONG_SIGMA)
  call ssmfe_terminate( keep, inform )

  deallocate ( a, b, t, x, lambda )

end subroutine test_errors_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_warnings_d

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: ldx
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
   write(*,"(a)") "======================="
   write(*,"(a)") "Testing warnings (real)"
   write(*,"(a)") "======================="

  options%print_level = -1
  
  n = 50
  nep = 2
  mep = n
  ldx = n
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
  call print_result(inform%flag, NO_SEARCH_DIRECTIONS_LEFT)
  call ssmfe_terminate( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, MAX_ITERATIONS_EXCEEDED)
  call ssmfe_terminate( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = 2
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, OUT_OF_STORAGE)
  call ssmfe_terminate( inform )
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
  call print_result(inform%flag, NO_SEARCH_DIRECTIONS_LEFT)
  call ssmfe_terminate( inform )
  b(1,1) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, MAX_ITERATIONS_EXCEEDED)
  call ssmfe_terminate( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = 2
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, OUT_OF_STORAGE)
  call ssmfe_terminate( inform )
  mep = n
  options%left_gap = 0

  deallocate ( a, b, t, x, lambda )

end subroutine test_warnings_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_warnings_z

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: ldx
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
   write(*,"(a)") "=========================="
   write(*,"(a)") "Testing warnings (complex)"
   write(*,"(a)") "=========================="

  options%print_level = -1
  
  n = 50
  nep = 2
  mep = n
  ldx = n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )
  
  a = ZERO
  forall ( i = 1 : n/2 ) a(i, i) = ONE
  forall ( i = n/2 + 1 : n ) a(i, i) = i

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  write(*,"(a)",advance="no") " * Testing warning flag 1...................."
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, NO_SEARCH_DIRECTIONS_LEFT)
  call ssmfe_terminate( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, MAX_ITERATIONS_EXCEEDED)
  call ssmfe_terminate( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = 2
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, OUT_OF_STORAGE)
  call ssmfe_terminate( inform )
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
  call print_result(inform%flag, NO_SEARCH_DIRECTIONS_LEFT)
  call ssmfe_terminate( inform )
  b(1,1) = ONE

  write(*,"(a)",advance="no") " * Testing warning flag 2...................."
  options%max_iterations = 1
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, MAX_ITERATIONS_EXCEEDED)
  call ssmfe_terminate( inform )
  options%max_iterations = 100

  write(*,"(a)",advance="no") " * Testing warning flag 3...................."
  mep = 4
  options%left_gap = 2
  call run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, OUT_OF_STORAGE)
  call ssmfe_terminate( inform )
  mep = n
  options%left_gap = 0

  deallocate ( a, b, t, x, lambda )

end subroutine test_warnings_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_options_d

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: ldx
  integer :: lwork
  integer :: i
  
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
  ldx = n
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
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, -9)
  call ssmfe_terminate( inform )
  n = 50

  write(*,"(a)",advance="no") " * Testing printing suppresssion............."
  options%print_level = 0
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing printing units suppression........"
  options%print_level = 3
  options%unit_error = -1
  options%unit_warning = -1
  options%unit_diagnostic = -1
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  options%unit_error = 1
  options%unit_warning = 1
  options%unit_diagnostic = 1

  write(*,"(a)",advance="no") " * Testing printing level 0.................."
  options%print_level = 0
  write ( 1, '(/x, a)' ) 'print level 0'
  t = ZERO
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, NO_SEARCH_DIRECTIONS_LEFT)
  call ssmfe_terminate( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing printing level 1.................."
  options%print_level = 1
  write ( 1, '(/x, a)' ) 'print level 1'
  call run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing printing level 2.................."
  options%print_level = 2
  write ( 1, '(/x, a)' ) 'print level 2'
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing printing level 3.................."
  options%print_level = 3
  write ( 1, '(/x, a)' ) 'print level 3'
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  options%unit_error = 6
  options%unit_warning = 6
  options%unit_diagnostic = 6
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing zero tolerances..................."
  options%tol_x = 0
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )
  options%tol_x = -1

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  nep = 3
  options%tol_x = 1e-3
  options%left_gap = 10 + eps
  forall ( i = 5 : n ) a(i, i) = i*10 + 2*eps
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%left, 4)
  call ssmfe_terminate( inform )
  options%left_gap = 0

  write ( *, '(a)' ) ' * Testing shift-invert...'

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  sigma = 440
  left = 3
  right = 0
  options%tol_x = 1e-3
  options%left_gap = 10 + eps
  forall ( i = 1 : n - 11 ) a(i, i) = i*10
  forall ( i = n - 10 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%left, 4)
  call ssmfe_terminate( inform )
  options%left_gap = 0

  write(*,"(a)",advance="no") " * Testing minimal right eigenvalue gap......"
  eps = 1D-3
  sigma = 439
  left = 0
  right = 3
  options%right_gap = 10 + eps
  forall ( i = 1 : n - 3 ) a(i, i) = i*10
  forall ( i = n - 2 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%right, 4)
  call ssmfe_terminate( inform )
  options%tol_x = -1
  options%right_gap = 0
  options%print_level = -1

  deallocate ( a, b, t, x, lambda )

end subroutine test_options_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_options_z

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: ldx
  integer :: lwork
  integer :: i
  
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
  ldx = n
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
  call print_result(inform%flag, -9)
  call ssmfe_terminate( inform )
  n = 50

  write(*,"(a)",advance="no") " * Testing printing units suppression........"
  options%print_level = 3
  options%unit_error = -1
  options%unit_warning = -1
  options%unit_diagnostic = -1
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  options%unit_error = 1
  options%unit_warning = 1
  options%unit_diagnostic = 1

  write(*,"(a)",advance="no") " * Testing printing level 0.................."
  options%print_level = 0
  write ( 1, '(/x, a)' ) 'print level 0'
  t = ZERO
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, NO_SEARCH_DIRECTIONS_LEFT)
  call ssmfe_terminate( inform )
  forall ( i = 1 : n ) t(i, i) = ONE

  write(*,"(a)",advance="no") " * Testing printing level 1.................."
  options%print_level = 1
  write ( 1, '(/x, a)' ) 'print level 1'
  call run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing printing level 2.................."
  options%print_level = 2
  write ( 1, '(/x, a)' ) 'print level 2'
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing printing level 3.................."
  options%print_level = 3
  write ( 1, '(/x, a)' ) 'print level 3'
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  options%unit_error = 6
  options%unit_warning = 6
  options%unit_diagnostic = 6
  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing zero tolerances..................."
  options%tol_x = 0
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )
  options%tol_x = -1

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  nep = 3
  options%tol_x = 1e-3
  options%left_gap = 10 + eps
  forall ( i = 5 : n ) a(i, i) = i*10 + 2*eps
  call run_std_z( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%left, 4)
  call ssmfe_terminate( inform )
  options%left_gap = 0

  write ( *, '(a)' ) ' * Testing shift-invert...'

  write(*,"(a)",advance="no") " * Testing minimal left eigenvalue gap......."
  eps = 1D-3
  sigma = 440
  left = 3
  right = 0
  options%tol_x = 1e-3
  options%left_gap = 10 + eps
  forall ( i = 1 : n - 11 ) a(i, i) = i*10
  forall ( i = n - 10 : n ) a(i, i) = i*10 + 2*eps
  call run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%left, 4)
  call ssmfe_terminate( inform )
  options%left_gap = 0

  deallocate ( a, b, t, x, lambda )

end subroutine test_options_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_misc_d

  integer :: n
  integer :: nep
  integer :: mep
  integer :: left, right
  integer :: ldx
  integer :: lwork
  integer :: i
  
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
   write(*,"(a)") "==================="
   write(*,"(a)") "Miscellaneous tests"
   write(*,"(a)") "==================="

  n = 50
  nep = 5
  mep = n
  ldx = n
  lwork = n*n
  allocate ( a(n, n), b(n, n), t(n, n), x(n, n), w(n, n) )
  allocate ( lambda(n), ipiv(n) )
  
  a = ZERO
  forall ( i = 1 : n ) a(i, i) = i*10

  b = ZERO
  forall ( i = 1 : n ) b(i, i) = ONE

  t = ZERO
  forall ( i = 1 : n ) t(i, i) = ONE

  eps = 1D-3
  forall ( i = 1 : 10 ) a(i, i) = i*10
  forall ( i = 11 : n ) a(i, i) = i*10 + 2*eps

  write(*,"(a)",advance="no") " * Testing restart..........................."
  eps = 1D-3
  nep = 3
!  options%tol_x = 1e-3
  options%left_gap = 10 + eps
  call run_std_d( n, a, t, nep, mep, lambda, x, options, inform )
  call print_result(inform%flag, 4)
  call ssmfe_terminate( inform )
  options%left_gap = 0

  options%print_level = -1

  write ( *, '(a)' ) ' * Testing shift-invert...'

  write(*,"(a)",advance="no") " * Testing standard.........................."
  sigma = 255
  left = 5
  right = 5
  call run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )

  write(*,"(a)",advance="no") " * Testing generalized......................."
  sigma = 255
  left = 5
  right = 5
  options%user_x = 1
  x = ONE
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, 0)
  call ssmfe_terminate( inform )
  options%user_x = 0

  options%print_level = -1

  write(*,"(a)",advance="no") " * Testing restart..........................."
  sigma = 409
  left = 0
  right = 5
  options%right_gap = 10 + eps
  call run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, 4)
  call ssmfe_terminate( inform )

  write ( *, '(a)' ) ' * Testing buckling...'

!  options%print_level = 3

  write(*,"(a)",advance="no") " * Testing restart..........................."
  a(n, n) = n*10 + 4*eps
  sigma = 399
  left = 0
  right = 5
  options%right_gap = 10 + eps
  call run_buckling_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, t, ipiv, w, lwork, inform )
  call print_result(inform%flag, 4)
  call ssmfe_terminate( inform )

  deallocate ( a, b, t, x, lambda )

end subroutine test_misc_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_d( n, a, t, nep, mep, lambda, x, options, inform )

  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(n)
  real(wp), intent(out) :: x(n, n)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform
  
  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep
  
  logical :: restarted

  restarted = .false.
  rci%job = 0
  do ! reverse communication loop
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
  call ssmfe_terminate ( keep )
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_std_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_z( n, a, t, nep, mep, lambda, x, options, inform )

  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(n)
  complex(wp), intent(out) :: x(n, n)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform
  
  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep
  
  logical :: restarted

  restarted = .false.
  rci%job = 0
  do ! reverse communication loop
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
  call ssmfe_terminate ( keep )
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_std_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_si_d &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(out) :: x(n, mep)
  real(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  real(wp), intent(out) :: work(lwork)
  integer, intent(in) :: lwork
  type(ssmfe_inform), intent(out) :: inform
  
  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep
  
  integer :: i

  ! perform LDLT factorization of the shifted matrix
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
!  print *, 'ssmfe status', inform%flag
  call ssmfe_terminate ( keep )

end subroutine run_std_si_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_std_si_z &
    ( options, n, a, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(out) :: x(n, mep)
  complex(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  complex(wp), intent(out) :: work(lwork)
  integer, intent(in) :: lwork
  type(ssmfe_inform), intent(out) :: inform
  
  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep
  
  integer :: i

  ! perform LDLT factorization of the shifted matrix
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
  call ssmfe_terminate ( keep )

end subroutine run_std_si_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_d( n, a, b, t, nep, mep, lambda, x, options, inform )

  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(n)
  real(wp), intent(out) :: x(n, n)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform
  
  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep

  rci%job = 0
  do ! reverse communication loop
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
  call ssmfe_terminate ( keep )

end subroutine run_gen_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_z( n, a, b, t, nep, mep, lambda, x, options, inform )

  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  complex(wp), intent(in) :: t(n, n)
  integer, intent(in) :: nep
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(n)
  complex(wp), intent(out) :: x(n, n)
  type(ssmfe_options), intent(in ) :: options
  type(ssmfe_inform ), intent(out) :: inform
  
  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep

  rci%job = 0
  do ! reverse communication loop
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
  call ssmfe_terminate ( keep )

end subroutine run_gen_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_si_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(out) :: x(n, mep)
  real(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  real(wp), intent(out) :: work(lwork)
  integer, intent(in) :: lwork
  type(ssmfe_inform), intent(out) :: inform
  
  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep
  
  logical :: restarted
  
  integer :: i

  ! perform LDLT factorization of the shifted matrix
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
  call ssmfe_terminate ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_gen_si_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_gen_si_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(out) :: x(n, mep)
  complex(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  complex(wp), intent(out) :: work(lwork)
  integer, intent(in) :: lwork
  type(ssmfe_inform), intent(out) :: inform
  
  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep
  
  logical :: restarted
  
  integer :: i

  ! perform LDLT factorization of the shifted matrix
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
  call ssmfe_terminate ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_gen_si_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_buckling_d &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  real(wp), intent(in) :: a(n, n)
  real(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  real(wp), intent(out) :: x(n, mep)
  real(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  real(wp), intent(out) :: work(lwork)
  integer, intent(in) :: lwork
  type(ssmfe_inform), intent(out) :: inform
  
  type(ssmfe_rcid) :: rci
  type(ssmfe_keepd) :: keep
  
  logical :: restarted
  
  integer :: i

  ! perform LDLT factorization of the shifted matrix
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
  call ssmfe_terminate ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_buckling_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine run_buckling_z &
    ( options, n, a, b, sigma, left, right, &
      mep, lambda, x, ldlt, ipiv, work, lwork, inform )

  type(ssmfe_options), intent(inout) :: options
  integer, intent(in) :: n
  complex(wp), intent(in) :: a(n, n)
  complex(wp), intent(in) :: b(n, n)
  real(wp), intent(in) :: sigma
  integer, intent(in) :: left
  integer, intent(in) :: right
  integer, intent(in) :: mep
  real(wp), intent(out) :: lambda(mep)
  complex(wp), intent(out) :: x(n, mep)
  complex(wp), intent(out) :: ldlt(n, n)
  integer, intent(out) :: ipiv(n)
  complex(wp), intent(out) :: work(lwork)
  integer, intent(in) :: lwork
  type(ssmfe_inform), intent(out) :: inform
  
  type(ssmfe_rciz) :: rci
  type(ssmfe_keepz) :: keep
  
  logical :: restarted
  
  integer :: i

  ! perform LDLT factorization of the shifted matrix
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
  call ssmfe_terminate ( keep )
  options%max_left = -1
  options%max_right = -1
  if ( restarted .and. inform%flag == 0 ) inform%flag = 4

end subroutine run_buckling_z

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

    implicit none

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

    implicit none

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


