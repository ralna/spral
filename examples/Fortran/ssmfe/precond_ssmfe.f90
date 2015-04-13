! examples/Fortran/ssmfe/precond_ssmfe.f90
! Laplacian on a square grid (using SPRAL_SSMFE routines)
program ssmfe_precond_example
  use spral_ssmfe
  use laplace2d ! implement Lapalacian and preconditioners
  implicit none

  integer, parameter :: wp = kind(0d0) ! Working precision is double

  integer, parameter :: m   = 20    ! grid points along each side
  integer, parameter :: n   = m*m   ! problem size
  integer, parameter :: nep = 5     ! eigenpairs wanted

  real(wp) :: lambda(2*nep)         ! eigenvalues
  real(wp) :: X(n, 2*nep)           ! eigenvectors
  type(ssmfe_rcid   ) :: rci        ! reverse communication data
  type(ssmfe_options) :: options    ! options
  type(ssmfe_keepd  ) :: keep       ! private data
  type(ssmfe_inform ) :: inform     ! information
  integer :: i                      ! loop index

  ! the gap between the last converged eigenvalue and the rest of the spectrum
  ! must be at least 0.1 times average gap between computed eigenvalues
  options%left_gap = -0.1
  rci%job = 0
  do ! reverse communication loop
    call ssmfe_standard &
      ( rci, nep, 2*nep, lambda, n, X, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call apply_laplacian( m, m, rci%nx, rci%x, rci%y )
    case ( 2 )
      call apply_gauss_seidel_step( m, m, rci%nx, rci%x, rci%y )
    case ( :-1 )
      exit
    end select
  end do
  print '(i3, 1x, a, i3, 1x, a)', inform%left, 'eigenpairs converged in', &
     inform%iteration, 'iterations'
  print '(1x, a, i2, a, es13.7)', &
    ('lambda(', i, ') = ', lambda(i), i = 1, inform%left)
  call ssmfe_free( keep, inform )
end program ssmfe_precond_example
