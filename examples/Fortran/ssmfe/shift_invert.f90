! examples/Fortran/ssmfe/shift_invert.f90
! Laplacian on a rectangular grid by shift-invert via LDLT factorization
program ssmfe_shift_invert_example
  use spral_ssmfe
  use laplace2d ! implement Lapalacian and preconditioners
  use ldltf     ! implements LDLT support routines
  implicit none

  integer, parameter :: wp = kind(0d0) ! Working precision

  integer, parameter :: nx = 8    ! grid points along x
  integer, parameter :: ny = 8    ! grid points along y
  integer, parameter :: n = nx*ny ! problem size
  real(wp), parameter :: sigma = 1.0 ! shift

  integer :: ipiv(n)               ! LDLT pivot index
  real(wp) :: lambda(n)            ! eigenvalues
  real(wp) :: X(n, n)              ! eigenvectors
  real(wp) :: A(n, n)              ! matrix
  real(wp) :: LDLT(n, n)           ! factors
  real(wp) :: work(n*n)            ! work array for dsytrf
  integer :: lwork = n*n           ! size of work
  integer :: left, right           ! wanted eigenvalues left and right of sigma
  integer :: i                     ! index
  type(ssmfe_options) :: options   ! eigensolver options
  type(ssmfe_inform ) :: inform    ! information
  type(ssmfe_rcid   ) :: rci       ! reverse communication data
  type(ssmfe_keepd  ) :: keep      ! private data

  call set_laplacian_matrix( nx, ny, A, n )

  ! perform LDLT factorization of the shifted matrix
  LDLT = A
  forall ( i = 1 : n ) LDLT(i, i) = A(i, i) - sigma
  lwork = n*n
  call dsytrf( 'L', n, LDLT, n, ipiv, work, lwork, i )

  left = num_neg_D(n, LDLT, n, ipiv) ! all eigenvalues to the left from sigma
  right = 5                          !   5 eigenvalues to the right from sigma
  rci%job = 0
  do
    call ssmfe_standard_shift &
      ( rci, sigma, left, right, n, lambda, n, X, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, 1.0_wp, A, n, rci%x, n, 0.0_wp, rci%y, n )
    case ( 9 )
      call dcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call dsytrs( 'L', n, rci%nx, LDLT, n, ipiv, rci%y, n, i )
    case ( :-1 )
      exit
    end select
  end do
  print '(1x, a, es10.2, 1x, a, i3, 1x, a)', 'Eigenvalues near', sigma, &
     '(took', inform%iteration, 'iterations)'
  print '(1x, a, i2, a, es13.7)', &
    ('lambda(', i, ') = ', lambda(i), i = 1, inform%left + inform%right)
  call ssmfe_free( keep, inform )
end program ssmfe_shift_invert_example
