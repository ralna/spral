! examples/Fortran/ssmfe/hermitian.f90 - Example code for SPRAL_SSMFE package
! Hermitian operator example
program ssmfe_hermitian_example
  use spral_ssmfe
  implicit none

  integer, parameter :: wp = kind(0d0) ! working precision

  integer, parameter :: n   = 80     ! problem size
  integer, parameter :: nep = 5      ! eigenpairs wanted

  real(wp) :: lambda(nep)            ! eigenvalues
  complex(wp) :: X(n, nep)           ! eigenvectors
  type(ssmfe_rciz   ) :: rci         ! reverse communication data
  type(ssmfe_options) :: options     ! options
  type(ssmfe_keepz  ) :: keep        ! private data
  type(ssmfe_inform ) :: inform      ! information
  integer :: i                       ! loop index

  rci%job = 0
  do ! reverse communication loop
    call ssmfe_standard( rci, nep, nep, lambda, n, X, n, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call apply_idx( n, rci%nx, rci%x, rci%y )
    case ( :-1 )
      exit
    end select
  end do
  print '(i3, 1x, a, i3, 1x, a)', inform%left, 'eigenpairs converged in', &
     inform%iteration, 'iterations'
  print '(1x, a, i2, a, es14.7)', &
    ('lambda(', i, ') = ', lambda(i), i = 1, inform%left)
  call ssmfe_free( keep, inform )

contains

  subroutine apply_idx( n, m, x, y ) ! central differences for i d/dx
    implicit none
    complex(wp), parameter :: IM_ONE = (0.0D0, 1.0D0)
    integer, intent(in) :: n, m
    complex(wp), intent(in) :: x(n, m)
    complex(wp), intent(out) :: y(n, m)
    integer :: i, j, il, ir
    do j = 1, m
      do i = 1, n
        if ( i == 1 ) then
          il = n
        else
          il = i - 1
        end if
        if ( i == n ) then
          ir = 1
        else
          ir = i + 1
        end if
        y(i, j) = IM_ONE*(x(ir, j) - x(il, j))
      end do
    end do
  end subroutine apply_idx

end program ssmfe_hermitian_example
