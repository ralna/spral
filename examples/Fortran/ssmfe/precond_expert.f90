! examples/Fortran/ssmfe/precond_expert.f90
! Laplacian on a square grid (using SPRAL_SSMFE_EXPERT routines)
program ssmfe_expert_precond_example
  use spral_random
  use spral_ssmfe_expert
  use laplace2d ! implement Lapalacian and preconditioners
  implicit none

  integer, parameter :: wp = kind(0d0) ! Working precision is double

  integer, parameter :: ngrid = 20           ! grid points along each side
  integer, parameter :: n     = ngrid*ngrid  ! problem size
  integer, parameter :: nep   = 5            ! eigenpairs wanted
  integer, parameter :: m     = 3            ! dimn of the iterated subspace

  type(random_state) :: state                ! PRNG state

  real(wp), external :: dnrm2, ddot          ! BLAS functions

  integer :: ncon                 ! number of converged eigenpairs
  integer :: i, j, k
  integer :: ind(m)               ! permutation index
  real(wp) :: s
  real(wp) :: lambda(n)           ! eigenvalues
  real(wp) :: X(n, n)             ! eigenvectors
  real(wp) :: rr(2*m, 2*m, 3)     ! work array
  real(wp) :: W(n, m, 0:7)        ! work array
  real(wp) :: U(n, m)             ! work array
  type(ssmfe_rcid       ) :: rci      ! reverse communication data
  type(ssmfe_options    ) :: options  ! options
  type(ssmfe_expert_keep) :: keep     ! private data
  type(ssmfe_inform     ) :: inform   ! information

  ! the gap between the last converged eigenvalue and the rest of the spectrum
  ! must be at least 0.1 times average gap between computed eigenvalues
  options%left_gap = -0.1
  ! block size is small, convergence may be slow, allow more iterations
  options%max_iterations = 200

  ! Initialize W to lin indep vectors by randomizing
  do i=1,n
    do j=1,m
      W(i,j,0) = random_real(state, positive=.true.)
    end do
  end do

  ncon = 0
  rci%job = 0
  do ! reverse communication loop
!    call ssmfe_standard &
    call ssmfe_generalized &
      ( rci, nep, n, lambda, m, rr, ind, keep, options, inform )
    select case ( rci%job )
    case ( 1 )
      call apply_laplacian &
        ( ngrid, ngrid, rci%nx, W(1 : n, rci%jx : rci%jx + rci%nx - 1, rci%kx),&
          W(1 : n, rci%jy : rci%jy + rci%nx - 1, rci%ky ) )
    case ( 2 )
      call apply_gauss_seidel_step &
        ( ngrid, ngrid, rci%nx, W(1 : n, rci%jx : rci%jx + rci%nx - 1, rci%kx),&
          W(1 : n, rci%jy : rci%jy + rci%nx - 1, rci%ky ) )
    case ( 3 )
      call dcopy( n*rci%nx, W(1, rci%jx, rci%kx), 1, W(1, rci%jy, rci%ky), 1 )
    case ( 5 )
      if ( rci%i < 0 ) cycle
      do k = 0, rci%nx - 1
        j = ncon + k + 1
        call dcopy( n, W(1, rci%jx + k, 0), 1, X(1, j), 1 )
      end do
      ncon = ncon + rci%nx
    case ( 11 )
      if ( rci%i == 0 ) then
        if ( rci%kx /= rci%ky .or. rci%jx > rci%jy ) then
          call dcopy &
            ( n * rci%nx, W(1, rci%jx, rci%kx), 1, W(1, rci%jy, rci%ky), 1 )
        else if ( rci%jx < rci%jy ) then
          do j = rci%nx - 1, 0, -1
            call dcopy &
              ( n, W(1, rci%jx + j, rci%kx), 1, W(1, rci%jy + j, rci%ky), 1 )
          end do
        end if
      else
        do i = 1, n
          do j = 1, rci%nx
            U(i, j) = W(i, ind(j), rci%kx)
          end do
          do j = 1, rci%nx
            W(i, j, rci%kx) = U(i, j)
          end do
          if ( rci%ky /= rci%kx ) then
            do j = 1, rci%nx
              U(i, j) = W(i, ind(j), rci%ky)
            end do
            do j = 1, rci%nx
              W(i, j, rci%ky) = U(i, j)
            end do
          end if
        end do
      end if
    case ( 12 )
      do i = 0, rci%nx - 1
        rr(rci%i + i, rci%j + i, rci%k) = &
          ddot(n, W(1, rci%jx + i, rci%kx), 1, W(1, rci%jy + i, rci%ky), 1)
      end do
    case ( 13 )
      do i = 0, rci%nx - 1
        if ( rci%kx == rci%ky ) then
          s = dnrm2(n, W(1, rci%jx + i, rci%kx), 1)
          if ( s > 0 ) &
            call dscal( n, 1/s, W(1, rci%jx + i, rci%kx), 1 )
        else
          s = sqrt(abs(ddot &
            (n, W(1, rci%jx + i, rci%kx), 1, W(1, rci%jy + i, rci%ky), 1)))
          if ( s > 0 ) then
            call dscal( n, 1/s, W(1, rci%jx + i, rci%kx), 1 )
            call dscal( n, 1/s, W(1, rci%jy + i, rci%ky), 1 )
          else
            call dcopy( n, 0.0D0, 0, W(1, rci%jy + i, rci%ky), 1 )
          end if
        end if
      end do
    case ( 14 )
      do i = 0, rci%nx - 1
        s = -rr(rci%i + i, rci%j + i, rci%k)
        call daxpy&
          ( n, s, W(1, rci%jx + i, rci%kx), 1, W(1, rci%jy + i, rci%ky), 1 )
      end do
    case ( 15 )
      if ( rci%nx > 0 .and. rci%ny > 0 ) &
        call dgemm &
          ( 'T', 'N', rci%nx, rci%ny, n, &
            rci%alpha, W(1, rci%jx, rci%kx), n, W(1, rci%jy, rci%ky), n, &
            rci%beta, rr(rci%i, rci%j, rci%k), 2*m )
    case ( 16, 17 )
      if ( rci%ny < 1 ) cycle
      if ( rci%nx < 1 ) then
        if ( rci%job == 17 ) cycle
        if ( rci%beta == 1.0D0 ) cycle
        do j = rci%jy, rci%jy + rci%ny - 1
          call dscal( n, rci%beta, W(1, j, rci%ky), 1 )
        end do
        cycle
      end if
      if ( rci%job == 17 ) then
        call dgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            1.0D0, W(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            0.0D0, W(1, rci%jy, rci%ky), n )
        call dcopy &
          ( n * rci%ny, W(1, rci%jy, rci%ky), 1, W(1, rci%jx, rci%kx), 1 )
      else
        call dgemm &
          ( 'N', 'N', n, rci%ny, rci%nx, &
            rci%alpha, W(1, rci%jx, rci%kx), n, rr(rci%i, rci%j, rci%k), 2*m, &
            rci%beta, W(1, rci%jy, rci%ky), n )
      end if
    case ( 21, 22 )
      if ( ncon > 0 ) then
        call dgemm &
          ( 'T', 'N', ncon, rci%nx, n, &
            1.0D0, X, n, W(1, rci%jy, rci%ky), n, 0.0D0, U, n )
        call dgemm &
          ( 'N', 'N', n, rci%nx, ncon, &
            -1.0D0, X, n, U, n, 1.0D0, W(1, rci%jx, rci%kx), n )
        call dgemm &
          ( 'N', 'N', n, rci%nx, ncon, &
            -1.0D0, X, n, U, n, 1.0D0, W(1, rci%jy, rci%ky), n )
      end if
    case ( 999 )
      if ( rci%k == 0 ) then
        if ( rci%jx > 1 ) then
          do j = 1, rci%jx-1
            do i = 1, n
              W(i, j, 0) = random_real(state, positive=.true.)
            end do
          end do
        endif
      end if
    case ( :-1 )
      exit
    end select
  end do
  print '(i3, 1x, a, i4, 1x, a)', &
    inform%left, 'eigenpairs converged in', inform%iteration, 'iterations'
  print '(1x, a, i2, a, es13.7)', &
    ('lambda(', i, ') = ', lambda(i), i = 1, inform%left)
  call ssmfe_free( keep, inform )
end program ssmfe_expert_precond_example
