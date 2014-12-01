program spec_test 
! Laplacian on a rectangular grid by shift-invert via LDLT factorization
  use spral_ssmfe
  use laplace2d
  use ldltf
  implicit none
  double precision, parameter :: ZERO = 0.0D0
  double precision, parameter :: ONE = 1.0D0
  integer, parameter :: nx = 8    ! grid points along x
  integer, parameter :: ny = 8    ! grid points along y
  integer, parameter :: n = nx*ny ! problem size
  integer :: ipiv(n)               ! LDLT pivot index
  double precision :: sigma = 1.0 ! shift
  double precision :: lambda(n)   ! eigenvalues
  double precision :: X(n, n)     ! eigenvectors
  double precision :: A(n, n)     ! matrix
  double precision :: LDLT(n, n)  ! factors
  double precision :: work(n*n)   ! work array for dsytrf
  integer :: lwork = n*n           ! size of work
  integer :: left, right           ! wanted eigenvalues left and right of sigma
  integer :: i                     ! index
  type(ssmfe_options) :: options ! eigensolver options
  type(ssmfe_inform ) :: inform  ! information
  type(ssmfe_rci    ) :: rci     ! reverse communication data
  type(ssmfe_keep   ) :: keep    ! private data
  call set_laplacian_matrix( nx, ny, A, n )
  ! perform LDLT factorization of the shifted matrix
  LDLT = A
  forall ( i = 1 : n ) LDLT(i, i) = A(i, i) - sigma
  lwork = n*n
  call dsytrf( 'L', n, LDLT, n, ipiv, work, lwork, i )
  ! all eigenvalues to the left from sigma are wanted
  left = num_neg_D(n, LDLT, n, ipiv)
  ! and 5 eigenvalues to the right from sigma are wanted
  right = 5
  rci%job = 0
  do
    call ssmfe_shift &
      ( rci, sigma, left, right, n, lambda, n, X, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call dgemm &
        ( 'N', 'N', n, rci%nx, n, ONE, A, n, rci%x, n, ZERO, rci%y, n )
    case ( 9 )
      call dcopy( n * rci%nx, rci%x, 1, rci%y, 1 )
      call dsytrs( 'L', n, rci%nx, LDLT, n, ipiv, rci%y, n, i )
    case ( :-1 )
      exit
    end select
  end do
  print '(1x, a, es10.2)', 'Eigenvalues near', sigma
  print '(1x, a, i2, a, es13.7)', &
    ('lambda(', i, ') = ', lambda(i), i = 1, inform%left + inform%right)
  call ssmfe_terminate( keep, inform )
end program spec_test ! Laplacian on a square grid

