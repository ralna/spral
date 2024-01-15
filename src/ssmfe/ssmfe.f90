! COPYRIGHT (c) 2013      Evgueni Ovtchinnikov
! COPYRIGHT (c) 2014,2015 The Science and Technology Facilities Council (STFC)
! Version 1.0.0
!
! Written by: Evgueni Ovtchinnikov
!
module SPRAL_ssmfe
  use SPRAL_ssmfe_expert, only: &
    ssmfe_errmsg, ssmfe_solve, ssmfe_free, &
    ssmfe_inform, ssmfe_expert_keep, ssmfe_options, ssmfe_rcid, ssmfe_rciz
  implicit none

  private

  ! double precision to be used
  integer, parameter :: PRECISION = kind(1.0D0)

  ! input error flags
  integer, parameter :: WRONG_RCI_JOB      =  -1
  integer, parameter :: WRONG_PROBLEM_SIZE =  -9
  integer, parameter :: WRONG_LDX          = -10
  integer, parameter :: WRONG_LEFT         = -11
  integer, parameter :: WRONG_RIGHT        = -12
  integer, parameter :: WRONG_STORAGE_SIZE = -13

  ! fatal execution error flags
  integer, parameter :: OUT_OF_MEMORY = -100

  ! termination status
  integer, parameter :: SSMFE_DONE  = -1 ! success
  integer, parameter :: SSMFE_QUIT  = -2 ! non-fatal error
  integer, parameter :: SSMFE_ABORT = -3 ! fatal error

  !
  ! data type to hold information between calls to solver procedures
  ! in RCI loop, real version
  !
  type ssmfe_keepd

    private

    integer :: block_size = 0 ! BJCG block size
    ! no. converged leftmost/left-of-the-shift eigenvalues
    integer :: lcon = 0
    integer :: rcon = 0 ! no. converged eigenvalues right of the shift
    integer :: step = 0 ! main computational step number

    ! integer work array for reordering work vectors
    integer, dimension(:), allocatable :: ind

    ! array for B*X
    real(PRECISION), dimension(:,:  ), allocatable :: BX
    ! work array for orthogonalization to converged eigenvectors
    real(PRECISION), dimension(:,:  ), allocatable :: U
    ! work array for the Rayleig-Ritz matrices
    real(PRECISION), dimension(:,:,:), allocatable :: V
    ! array of approximate eigenvectors and work vectors
    real(PRECISION), dimension(:,:,:), allocatable :: W

    ! expert interface keep (see expert.f90)
    type(ssmfe_expert_keep) :: keep

  end type ssmfe_keepd

  !
  ! data type to hold information between calls to solver procedures
  ! in RCI loop, complex version
  !
  type ssmfe_keepz

    private

    integer :: block_size = 0
    integer :: lcon = 0
    integer :: rcon = 0
    integer :: step = 0

    integer, dimension(:), allocatable :: ind

    complex(PRECISION), dimension(:,:  ), allocatable :: BX
    complex(PRECISION), dimension(:,:  ), allocatable :: U
    complex(PRECISION), dimension(:,:,:), allocatable :: V
    complex(PRECISION), dimension(:,:,:), allocatable :: W

    type(ssmfe_expert_keep) :: keep

  end type ssmfe_keepz

  interface ssmfe_standard
    module procedure ssmfe_std_double, ssmfe_std_double_complex
  end interface

  interface ssmfe_generalized
    module procedure ssmfe_gen_double, ssmfe_gen_double_complex
  end interface

  interface ssmfe_standard_shift
    module procedure ssmfe_shift_double, ssmfe_shift_double_complex
  end interface

  interface ssmfe_generalized_shift
    module procedure ssmfe_gen_shift_double,  ssmfe_gen_shift_double_complex
  end interface

  interface ssmfe_buckling
    module procedure ssmfe_buckling_double, ssmfe_buckling_double_complex
  end interface

  interface ssmfe_vector_operations
    module procedure &
      ssmfe_vector_operations_double, &
      ssmfe_vector_operations_double_complex
  end interface

  interface ssmfe_free
    module procedure &
      ssmfe_free_double, &
      ssmfe_free_double_complex, &
      ssmfe_free_keep_double, &
      ssmfe_free_keep_double_complex
  end interface

  public :: ssmfe_standard, ssmfe_standard_shift
  public :: ssmfe_generalized, ssmfe_generalized_shift
  public :: ssmfe_buckling
  public :: ssmfe_free
  public :: ssmfe_options, ssmfe_keepd, ssmfe_keepz
  public :: ssmfe_rcid, ssmfe_rciz, ssmfe_inform

contains

!*************************************************************************
! real RCI for computing leftmost eigenpairs of A x = lambda x
! see ssmfe spec for the arguments
!
  subroutine ssmfe_std_double &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    real(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keepd  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_direct_srci_double &
      ( 0, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_std_double

!*************************************************************************
! real RCI for computing leftmost eigenpairs of A x = lambda B x, B > 0,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_gen_double &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    real(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keepd  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_direct_srci_double &
      ( 1, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_gen_double

!*************************************************************************
! real RCI for computing eigenpairs {lambda, x} of A x = lambda x
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_shift_double &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    real(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keepd  ), intent(inout), target :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_inverse_srci_double &
      ( 0, sigma, left, right, mep, lambda, n, X, ldX, &
        rci, keep, options, inform )

  end subroutine ssmfe_shift_double

!**************************************************************************
! real RCI for computing eigenpairs {lambda, x} of A x = lambda B x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_gen_shift_double &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    real(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keepd  ), intent(inout), target :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_inverse_srci_double &
      ( 1, sigma, left, right, mep, lambda, n, X, ldX, &
        rci, keep, options, inform )

  end subroutine ssmfe_gen_shift_double

!**************************************************************************
! real RCI for computing eigenpairs {lambda, x} of B x = lambda A x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_buckling_double &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    real(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keepd  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_inverse_srci_double &
      ( 2, sigma, left, right, mep, lambda, n, X, ldX, &
        rci, keep, options, inform )

  end subroutine ssmfe_buckling_double

!*************************************************************************
! complex RCI for computing leftmost eigenpairs of A x = lambda x,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_std_double_complex &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    complex(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keepz  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_direct_srci_double_complex &
      ( 0, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_std_double_complex

!*************************************************************************
! complex RCI for computing leftmost eigenpairs of A x = lambda B x, B > 0,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_gen_double_complex &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    complex(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keepz  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_direct_srci_double_complex &
      ( 1, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_gen_double_complex

!*************************************************************************
! complex RCI for computing eigenpairs {lambda, x} of A x = lambda x
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_shift_double_complex &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    complex(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keepz  ), intent(inout), target :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_inverse_srci_double_complex &
      ( 0, sigma, left, right, mep, lambda, n, X, ldX, &
        rci, keep, options, inform )

  end subroutine ssmfe_shift_double_complex

!*****************************************************************************
! complex RCI for computing eigenpairs {lambda, x} of A x = lambda B x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_gen_shift_double_complex &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    complex(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keepz  ), intent(inout), target :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_inverse_srci_double_complex &
      ( 1, sigma, left, right, mep, lambda, n, X, ldX, &
        rci, keep, options, inform )

  end subroutine ssmfe_gen_shift_double_complex

!*****************************************************************************
! complex RCI for computing eigenpairs {lambda, x} of B x = lambda A x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_buckling_double_complex &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    complex(PRECISION), intent(inout) :: X(ldX, mep)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keepz  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_inverse_srci_double_complex &
      ( 2, sigma, left, right, mep, lambda, n, X, ldX, &
        rci, keep, options, inform )

  end subroutine ssmfe_buckling_double_complex

!*************************************************************************
! real RCI for computing leftmost eigenpairs
!
  subroutine ssmfe_direct_srci_double &
      ( problem, left, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, copy => dcopy, gemm => dgemm
    use SPRAL_lapack_iface, mxcopy => dlacpy

    character, parameter :: TR = 'T'

    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

    integer, parameter :: NONE = -1

    real(PRECISION), parameter :: ZERO = 0.0D0
    real(PRECISION), parameter :: ONE = 1.0D0

    ! problem type
    ! 0: A x = lambda x
    ! 1: A x = lambda B x
    ! A = A', B = B' > 0
    integer, intent(in) :: problem

    ! number of wanted leftmost eigenpairs
    integer, intent(in) :: left

    ! eigenpairs storage size
    integer, intent(in) :: max_nep

    ! eigenvalue storage
    real(PRECISION), intent(inout) :: lambda(max_nep)

    ! problem size
    integer, intent(in) :: n

    ! leading dimension of X
    integer, intent(in) :: ldX

    ! eigenvector storage
    real(PRECISION), intent(inout) :: X(ldX, max_nep)

    ! ssmfe types - see the spec
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    type(ssmfe_keepd  ), intent(inout), target :: keep

    integer :: extra
    integer :: total
    integer :: m
    integer :: kw
    integer :: ldBX
    integer :: i, j, k

    ! check for the input data errors

    if ( n < 1 ) then
      inform%flag = WRONG_PROBLEM_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong problem size', n
      return
    end if

    if ( ldx < n ) then
      inform%flag = WRONG_LDX
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong leading dimension of x', ldx
      return
    end if

    if ( left < 0 .or. left >= n/2 ) then
      inform%flag = WRONG_LEFT
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong number of eigenpairs', left
      return
    end if

    if ( left > max_nep ) then
      inform%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong eigenpairs storage size', max_nep
      return
    end if

    if ( rci%job == 0 ) then

      ! set the number of extra vectors iterated for better convergence
      ! and reliability
      if ( left == 1 ) then
        extra = 7
      else
        extra = max(left/10, 10)
      end if
      total = min(max(1, left + extra), max(1, n/2 - 1))
      if ( options%max_left >= 0 ) &
        total = min(total, options%max_left)

      keep%block_size = total ! BJCG block size
      keep%lcon = 0 ! number of converged eigenpairs

      ! (re)allocate work arrays
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate &
        ( keep%ind(keep%block_size), &
          keep%U(max_nep, 2*keep%block_size + max_nep), &
          keep%V(2*keep%block_size, 2*keep%block_size, 3), stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return

      if ( problem == 0 ) then
        ! the number of blocks in the main work array keep%W
        kw = 5
      else
        ! (re)allocate storage for B*X
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep), stat = inform%stat )
        if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
        if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
        if ( inform%stat /= 0 ) return
        ! the number of blocks in the main work array keep%W
        kw = 7
      end if
      allocate ( keep%W(n, keep%block_size, 0:kw), stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return
      ! fill the first block of keep%W with random numbers in (-1,1),
      ! its columns to be used as the initial eigenvector approximations
      call random_number( keep%W(:,:,0) )
      keep%W(:,:,0) = 2*keep%W(:,:,0) - ONE
      if ( options%user_X > 0 ) &
        ! overwrite with user-supplied initial vectors
        call mxcopy &
          ( 'A', n, min(keep%block_size, options%user_X), X, ldx, keep%W, n )

    else

      if ( .not. allocated(keep%ind) ) then
        ! missing the call with rci%job = 0 detected
        inform%flag = WRONG_RCI_JOB
        rci%job = SSMFE_QUIT
        if ( options%unit_error > NONE .and. options%print_level > NONE ) &
          write( options%unit_error, '(/a, i3/)' ) &
            '??? Wrong rci%job', rci%job
        return
      end if

    end if

    ldBX = n

    do
      ! call ssmfe_expert solver
      call ssmfe_solve &
        ( problem, left, max_nep, lambda, keep%block_size, keep%V, keep%ind, &
          rci, keep%keep, options, inform )

      select case ( rci%job )

      case ( 11:19 )

        ! apply some standard operations to vectors in W
        k = max(1, rci%k) ! cover for the case rci%k = 0
        call ssmfe_vector_operations_double &
          ( rci, n, keep%block_size, keep%W, n, keep%V(1, 1, k), &
            keep%ind, keep%U )

      case ( SSMFE_SAVE_CONVERGED )

        if ( rci%nx < 1 ) cycle

        ! only leftmost eigenpairs computed (cf. ssmfe_expert solver)
        if ( rci%i < 0 ) cycle

        ! save the converged eigenvectors in X
        do i = 1, rci%nx
          k = (i - 1)*rci%i
          j = keep%lcon + i
          call copy( n, keep%W(1, rci%jx + k, 0), 1, X(1,j), 1 )
          if ( problem /= 0 ) &
            ! save their products with the matrix B in keep%BX
            call copy( n, keep%W(1, rci%jy + k, rci%ky), 1, keep%BX(1,j), 1 )
        end do

        ! update the Gram matrix for the converged eigenvectors

        ! newly converged eigenvectors are in columns
        ! k, ..., k + rci%nx - 1 of X
        k = keep%lcon + 1

        m = keep%block_size
        if ( keep%lcon > 0 ) then
          if ( problem == 0 ) then
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                ONE, X, ldX, X(1, k), ldX, &
                ZERO, keep%U(1, m + k), max_nep )
          else
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                ONE, X, ldX, keep%BX(1, k), ldX, &
                ZERO, keep%U(1, m + k), max_nep )
          end if
          do j = 1, rci%nx
            do i = 1, keep%lcon
              keep%U(k + j - 1, m + i) = keep%U(i, m + k + j - 1)
            end do
          end do
        end if
        if ( problem == 0 ) then
          call gemm &
            ( TR, 'N', rci%nx, rci%nx, n, &
              ONE, X(1, k), ldX, X(1, k), ldX, &
              ZERO, keep%U(k, m + k), max_nep)
        else
          call gemm &
            ( TR, 'N', rci%nx, rci%nx, n, &
              ONE, X(1, k), ldX, keep%BX(1, k), ldX, &
              ZERO, keep%U(k, m + k), max_nep)
        end if

        ! update the number of converged eigenpairs
        keep%lcon = keep%lcon + rci%nx

      case ( SSMFE_APPLY_ADJ_CONSTRS )

        ! apply I - B X X' to columns of keep%W specified by rci

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
        if ( keep%lcon > 0 ) then
          call gemm &
            ( TR, 'N', keep%lcon, rci%nx, n, &
              ONE, X, ldX, rci%x, n, ZERO, keep%U, max_nep )
          if ( problem == 0 ) then
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -ONE, X, ldX, keep%U, max_nep, ONE, rci%x, n )
          else
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -ONE, keep%BX, ldBX, keep%U, max_nep, ONE, rci%x, n )
          end if
        end if

      case ( SSMFE_APPLY_CONSTRAINTS )

        ! apply I - Y H Y'B to W, where W holds the columns of keep%W
        ! specified by rci%jx, rci%kx and rci%nx,
        ! Y is the left part of X containing the converged eigenvectors
        ! and H is the inverse of the gram matrix G = Y' B Y

        ! rci%x points to W
        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        ! rci%y points to B W
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

        if ( keep%lcon > 0 ) then

          ! compute Y'B W
          call gemm &
            ( TR, 'N', keep%lcon, rci%nx, n, &
              ONE, X, ldX, rci%y, n, ZERO, keep%U, max_nep )

          ! since the Gram matrix G is close to identity, its inverse H
          ! is close to 2I - G, the difference between the two matrices
          ! being of the order (I - G)^2
          m = keep%block_size
          k = m + max_nep + 1
          call mxcopy &
            ( 'A', keep%lcon, rci%nx, keep%U, max_nep, keep%U(1, k), max_nep )
          call gemm &
            ( 'N', 'N', keep%lcon, rci%nx, keep%lcon, &
              -ONE, keep%U(1, m + 1), max_nep, keep%U(1, k), max_nep, &
              2*ONE, keep%U, max_nep )

          ! update W
          call gemm &
            ( 'N', 'N', n, rci%nx, keep%lcon, &
              -ONE, X, ldX, keep%U, max_nep, ONE, rci%x, n )
          if ( problem /= 0 ) &
            ! update B W
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -ONE, keep%BX, ldBX, keep%U, max_nep, ONE, rci%y, n )

        end if

      case ( SSMFE_APPLY_A, SSMFE_APPLY_B )

        ! pass pointers to columns of keep%W to the user for
        ! multiplying by A or B

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)

        return

      case ( SSMFE_APPLY_PREC )

        ! pass pointers to columns of keep%W to the user for
        ! applying the preconditioner

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)

        ! no action needs to be taken by the user if no preconditioning is used
        call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

        return

      case ( SSMFE_RESTART )

        if ( rci%k == 0 .and. rci%jx > 1 ) then
          ! fill the columns specified by rci%jx with random numbers
          call random_number( keep%W(:, 1 : rci%jx - 1, 0) )
          keep%W(:, 1 : rci%jx - 1, 0) = 2*keep%W(:, 1 : rci%jx - 1, 0) - ONE
        end if

        return

      case ( : SSMFE_DONE )

        return

      end select

    end do

  end subroutine ssmfe_direct_srci_double

  subroutine ssmfe_inverse_srci_double &
      ( problem, sigma, left, right, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, copy => dcopy, gemm => dgemm
    use SPRAL_lapack_iface, mxcopy => dlacpy

    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE   = 9
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

    integer, parameter :: NONE = -1

    integer, intent(in) :: problem
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: max_nep
    real(PRECISION), intent(inout) :: lambda(max_nep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    real(PRECISION), intent(inout) :: X(ldX, max_nep)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    type(ssmfe_keepd  ), intent(inout), target :: keep

    character, parameter :: TR = 'T'
    real(PRECISION), parameter :: ZERO = 0.0D0
    real(PRECISION), parameter :: ONE = 1.0D0

    integer :: nep
    integer :: extra_left, extra_right
    integer :: total_left, total_right
    integer :: kw
    integer :: ldBX
    integer :: i, j, k, m

    integer :: ilaenv

    real(PRECISION) :: s

    ! check for the input data errors

    if ( n < 1 ) then
      inform%flag = WRONG_PROBLEM_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong problem size', n
      return
    end if

    if ( ldx < n ) then
      inform%flag = WRONG_LDX
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong leading dimension of x', ldx
      return
    end if

    if ( left < 0 .or. left >= right .and. left + right > n/2 ) then
      inform%flag = WRONG_LEFT
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong number of eigenpairs on the left', left
      return
    end if

    if ( right < 0 .or. right > left .and. left + right > n/2 ) then
      inform%flag = WRONG_RIGHT
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong number of eigenpairs on the right', right
      return
    end if

    if ( left + right > max_nep ) then
      inform%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong eigenpairs storage size', max_nep
      return
    end if

    if ( rci%job == 0 ) then

      keep%step = 0 ! solve step
      keep%lcon = 0 ! number of converged eigenvalues left of sigma
      keep%rcon = 0 ! number of converged eigenvalues right of sigma
      ! set the number of extra vectors iterated for better convergence
      ! and reliability
      if ( left > 0 ) then
        extra_left = max(10, left/10)
      else
        extra_left = 0
      end if
      if ( right > 0 ) then
        extra_right = max(10, right/10)
      else
        extra_right = 0
      end if
      total_left = left + extra_left
      total_right = right + extra_right
      if ( options%max_left >= 0 ) &
        total_left = min(total_left, options%max_left)
      if ( options%max_right >= 0 ) &
        total_right = min(total_right, options%max_right)

      ! BJCG block size
      keep%block_size = max(2, total_left + total_right)
      keep%block_size = min(keep%block_size, max(1, n/2 - 1))

      ! (re)allocate work arrays
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate &
        ( keep%ind(keep%block_size), &
          keep%U(max_nep, 2*keep%block_size + max_nep), &
          keep%V(2*keep%block_size, 2*keep%block_size, 3), stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return

      if ( problem == 0 ) then
        ! the number of blocks in the main work array keep%W
        kw = 5
      else
        ! (re)allocate storage for B*X
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep), stat = inform%stat )
        if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
        if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
        if ( inform%stat /= 0 ) return
        ! the number of blocks in the main work array keep%W
        kw = 7
      end if
      allocate ( keep%W(n, keep%block_size, 0:kw), stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return
      ! fill the first block of keep%W with random numbers in (-1,1),
      ! its columns to be used as the initial eigenvector approximations
      call random_number( keep%W(:,:,0) )
      keep%W(:,:,0) = 2*keep%W(:,:,0) - ONE
      if ( options%user_X > 0 ) &
        ! overwrite with user-supplied initial vectors
        call mxcopy &
          ( 'A', n, min(keep%block_size, options%user_X), X, ldx, keep%W, n )

      ! initialize the Gram matrix for the converged eigenvectors to the
      ! identity matrix of size max_nep
      keep%U = ZERO
      do i = 1, max_nep
        keep%U(i, keep%block_size + i) = ONE
      end do

    else

      if ( .not. allocated(keep%ind) ) then
        ! missing the call with rci%job = 0 detected
        inform%flag = WRONG_RCI_JOB
        rci%job = SSMFE_QUIT
        if ( options%unit_error > NONE .and. options%print_level > NONE ) &
          write( options%unit_error, '(/a/)' ) &
            '??? Wrong rci%job'
        return
      end if

    end if

    ldBX = n
    nep = inform%left + inform%right

    select case ( keep%step )

    case ( 0 )

      do
        ! call ssmfe_expert solver
        call ssmfe_solve &
          ( problem, sigma, left, right, &
            max_nep, lambda, keep%block_size, keep%V, keep%ind, &
            rci, keep%keep, options, inform )

        select case ( rci%job )

        case ( 11:19 )

          ! apply some standard operations to vectors in W
          k = max(1, rci%k) ! cover for the case rci%k = 0
          call ssmfe_vector_operations_double &
            ( rci, n, keep%block_size, keep%W, n, keep%V(1, 1, k), &
              keep%ind, keep%U )

        case ( SSMFE_SAVE_CONVERGED )

          if ( rci%nx < 1 ) cycle

          ! save the converged eigenvectors in X
          do i = 1, rci%nx
            k = (i - 1)*rci%i
            if ( rci%i > 0 ) then
              j = keep%lcon + i ! left eigenvectors are saved on the left
            else
              j = max_nep - keep%rcon - i + 1 ! and right on the right
            end if
            call copy( n, keep%W(1, rci%jx + k, 0), 1, X(1, j), 1 )
            if ( problem /= 0 ) &
              ! save their products with the matrix B in keep%BX
              call copy( n, keep%W(1, rci%jy + k, rci%ky), 1, keep%BX(1, j), 1 )
          end do

          ! update the Gram matrix for the converged eigenvectors
          ! stored as four blocks (two diagonal and two off-diagonal)
          ! of a matrix G_X of size max_nep in the array keep%U:

          !       | G_l  0   F  |
          ! G_X = |  0   I   0  |
          !       !  F'  0  G_r |

          ! where

          ! G_l = X(:, 1:keep%lcon)' B X(:, 1:keep%lcon)
          ! G_r = X(:, 1:keep%rcon)' B X(:, 1:keep%rcon)
          ! F   = X(:, 1:keep%lcon)' B X(:, 1:keep%rcon)

          ! newly converged eigenvectors are in columns
          ! k, ..., k + rci%nx - 1 of X
          if ( rci%i > 0 ) then
            k = keep%lcon + 1
          else
            k = max_nep - keep%rcon - rci%nx + 1
          end if

          m = keep%block_size
          if ( keep%lcon > 0 ) then
            if ( problem == 0 ) then
              call gemm &
                ( TR, 'N', keep%lcon, rci%nx, n, &
                  ONE, X, ldX, X(1, k), ldX, &
                  ZERO, keep%U(1, m + k), max_nep )
            else
              call gemm &
                ( TR, 'N', keep%lcon, rci%nx, n, &
                  ONE, X, ldX, keep%BX(1, k), ldX, &
                  ZERO, keep%U(1, m + k), max_nep )
            end if
            do j = 1, rci%nx
              do i = 1, keep%lcon
                keep%U(k + j - 1, m + i) = keep%U(i, m + k + j - 1)
              end do
            end do
          end if
          if ( keep%rcon > 0 ) then
            if ( problem == 0 ) then
              call gemm &
                ( TR, 'N', keep%rcon, rci%nx, n, &
                  ONE, X(1, max_nep - keep%rcon + 1), ldX, X(1, k), ldX, &
                  ZERO, keep%U(max_nep - keep%rcon + 1, m + k), max_nep )
            else
              call gemm &
                ( TR, 'N', keep%rcon, rci%nx, n, &
                  ONE, X(1, max_nep - keep%rcon + 1), ldX, keep%BX(1, k), ldX, &
                  ZERO, keep%U(max_nep - keep%rcon + 1, m + k), max_nep )
            end if
            do j = 1, rci%nx
              do i = 1, keep%rcon
                keep%U(k + j - 1, m + max_nep - keep%rcon + i) = &
                  keep%U(max_nep - keep%rcon + i, m + k + j - 1)
              end do
            end do
          end if
          if ( problem == 0 ) then
            call gemm &
              ( TR, 'N', rci%nx, rci%nx, n, &
                ONE, X(1, k), ldX, X(1, k), ldX, &
                ZERO, keep%U(k, m + k), max_nep)
          else
            call gemm &
              ( TR, 'N', rci%nx, rci%nx, n, &
                ONE, X(1, k), ldX, keep%BX(1, k), ldX, &
                ZERO, keep%U(k, m + k), max_nep)
          end if

          ! update the numbers of converged eigenpairs
          if ( rci%i > 0 ) then
            keep%lcon = keep%lcon + rci%nx
          else
            keep%rcon = keep%rcon + rci%nx
          end if

        case ( SSMFE_APPLY_ADJ_CONSTRS )

          ! apply I - B Y Y' to W, where W holds the  columns of keep%W
          ! specified by rci%jx, rci%kx and rci%nx, and Y is holds
          ! the columns of X containing the converged eigenvectors

          ! rci%x points to W
          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          ! rci%y points to B W
          rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

          if ( keep%lcon > 0 ) then
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                ONE, X, ldX, rci%x, n, ZERO, keep%U, max_nep )
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -ONE, X, ldX, keep%U, max_nep, ONE, rci%x, n )
          end if

          if ( keep%rcon > 0 ) then
            j = max_nep - keep%rcon + 1
            call gemm &
              ( TR, 'N', keep%rcon, rci%nx, n, &
                ONE, X(1, j), ldX, rci%x, n, ZERO, keep%U, max_nep )
              call gemm &
                ( 'N', 'N', n, rci%nx, keep%rcon, &
                  -ONE, X(1,j), ldX, keep%U, max_nep, ONE, rci%x, n )
          end if

        case ( SSMFE_APPLY_CONSTRAINTS )

          if ( keep%lcon == 0 .and. keep%rcon == 0 ) cycle

          ! apply I - Y H Y'B to W, where W holds the  columns of keep%W
          ! specified by rci%jx, rci%kx and rci%nx, Y holds columns on the
          ! left and right margins of X containing the converged eigenvectors,
          ! and H is the inverse of the Gram matrix G = Y' B Y

          ! rci%x points to W
          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          ! rci%y points to B W
          rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

          ! compute Y'B W
          if ( keep%lcon > 0 ) then
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                ONE, X, ldX, rci%y, n, ZERO, keep%U, max_nep )
          end if
          if ( keep%rcon > 0 ) then
            j = max_nep - keep%rcon + 1
            call gemm &
              ( TR, 'N', keep%rcon, rci%nx, n, &
                ONE, X(1, j), ldX, rci%y, n, ZERO, keep%U(j, 1), max_nep )
          end if

          ! since the Gram matrix G is close to identity, its inverse H
          ! is close to 2I - G, the difference between the two matrices
          ! being of the order (I - G)^2

          ! H is applied to

          ! | U_l |
          ! | U_r |

          ! where U_l is keep%lcon by rci%nx and U_r is keep%rcon by rci%nx,
          ! which is equivalent to applying H_X = 2I - G_X to

          ! | U_l |
          ! !  0  |
          ! | U_r |

          ! and collecting the first keep%lcon and last keep%rcon rows
          ! of the result

          m = keep%block_size
          k = m + max_nep + 1
          call mxcopy &
            ( 'A', max_nep, rci%nx, keep%U, max_nep, keep%U(1, k), max_nep )
          call gemm &
            ( 'N', 'N', max_nep, rci%nx, max_nep, &
              -ONE, keep%U(1, m + 1), max_nep, keep%U(1, k), max_nep, &
              2*ONE, keep%U, max_nep )

          ! update W and, if B is not the identtity, B W
          if ( keep%lcon > 0 ) then
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -ONE, X, ldX, keep%U, max_nep, ONE, rci%x, n )
            if ( problem /= 0 ) &
              call gemm &
                ( 'N', 'N', n, rci%nx, keep%lcon, &
                  -ONE, keep%BX, ldBX, keep%U, max_nep, ONE, rci%y, n )
          end if
          if ( keep%rcon > 0 ) then
            j = max_nep - keep%rcon + 1
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -ONE, X(1, j), ldX, keep%U(j, 1), max_nep, ONE, rci%x, n )
            if ( problem /= 0 ) &
              call gemm &
                ( 'N', 'N', n, rci%nx, keep%rcon, &
                  -ONE, keep%BX(1, j), ldBX, keep%U(j, 1), max_nep, &
                  ONE, rci%y, n )
          end if

        case ( SSMFE_DO_SHIFTED_SOLVE, SSMFE_APPLY_B )

          ! pass pointers to columns of keep%W to the user for
          ! multiplying by A or B

          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)

          return

        case ( SSMFE_APPLY_PREC )

          ! pass pointers to columns of keep%W to the user for
          ! applying the preconditioner

          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)

          ! no action is to be taken by the user if no preconditioning is used
          call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

          return

        case ( SSMFE_RESTART )

          if ( rci%k == 0 ) then
            ! fill the columns specified by rci with random numbers
            if ( rci%jx > 1 ) then
              call random_number( keep%W(:, 1 : rci%jx - 1, 0) )
              keep%W(:, 1 : rci%jx - 1, 0) = &
                2*keep%W(:, 1 : rci%jx - 1, 0) - ONE
            end if
            if ( rci%jx + rci%nx - 1 < keep%block_size ) then
              call random_number &
                ( keep%W(:, rci%jx + rci%nx : keep%block_size, 0) )
              keep%W(:, rci%jx + rci%nx : keep%block_size, 0) = &
                2*keep%W(:, rci%jx + rci%nx : keep%block_size, 0) - ONE
            end if
          end if

          return

        case ( SSMFE_QUIT, SSMFE_DONE )

          ! total number of converged eigenpairs
          nep = inform%left + inform%right

          if ( nep < 1 ) return

          rci%k = rci%job

          if ( inform%right > 0 ) then
            ! bring the eigenpairs together
            do j = 1, inform%right
              lambda(inform%left + j) = lambda(max_nep - j + 1)
            end do
            call mxcopy &
              ( 'A', n, inform%right, X(1, max_nep - inform%right + 1), ldX, &
                X(1, inform%left + 1), ldX )
          end if

          ! apply one more block shifted solve to the converged eigenvectors
          ! and the Rayleigh-Ritz procedure in the trial subspace spanned by
          ! them

          k = ilaenv(1, 'DSYTRD', 'U', nep, -1, -1, -1)
          k = max(3*nep - 1, (k + 2)*nep)

          ! reallocate the work arrays
          deallocate ( keep%W, keep%V, keep%U )
          allocate &
            ( keep%W(n, nep, 3), keep%V(nep, nep, 3), keep%U(k, 1), &
              stat = inform%stat )
          if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
          if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
          if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
          if ( inform%stat /= 0 ) return

          if ( problem == 0 ) then
            ! compute Y = (A - sigma I)^{-1} X
            call mxcopy( 'A', n, nep, X, ldX, keep%W, n )
            rci%job = SSMFE_DO_SHIFTED_SOLVE
            rci%nx = nep
            rci%jx = 1
            rci%kx = 0
            rci%jy = 1
            rci%ky = 2
            rci%x => keep%W(:,:,1)
            rci%y => keep%W(:,:,2)
            keep%step = 2
          else
            ! compute B X
            call mxcopy( 'A', n, nep, X, ldX, keep%W(1,1,2), n )
            rci%job = SSMFE_APPLY_B
            rci%nx = nep
            rci%jx = 1
            rci%kx = 2
            rci%jy = 1
            rci%ky = 0
            rci%x => keep%W(:,:,2)
            rci%y => keep%W(:,:,1)
            keep%step = 1
          end if

          return

        case ( SSMFE_ABORT )

          return

        end select

      end do

    case ( 1 )

      ! compute Y = (A - sigma B)^{-1} B X
      ! (or Y = (B - sigma A)^{-1} B X in the buckling case)
      rci%job = SSMFE_DO_SHIFTED_SOLVE
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 2
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,2)
      keep%step = 2

    case ( 2 )

      ! compute A Y
      call mxcopy( 'A', n, nep, keep%W(1,1,2), n, keep%W, n )
      if ( problem == 0 ) &
        call mxcopy( 'A', n, nep, keep%W, n, keep%W(1,1,3), n )
      rci%job = SSMFE_APPLY_A
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 2
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,2)
      if ( problem == 0 ) then
        keep%step = 4
      else
        keep%step = 3
      end if

    case ( 3 )

      ! compute B Y
      rci%job = SSMFE_APPLY_B
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 4
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,3)
      keep%step = 4

    case ( 4 )

      ! perform the Rayleigh-Ritz procedure in the subspace spanned by
      ! the columns of Y to obtain better approximations to eigenvectors

      if ( problem == 2 ) then
        ! in the buckling case, A and B are swapped, so reverse the order
        ! of eigenvalues by multiplying Y'A Y by -1
        s = -ONE
      else
        s = ONE
      end if
      call gemm &
        ( TR, 'N', nep, nep, n, s, keep%W, n, keep%W(1,1,2), n, &
          ZERO, keep%V, nep )
      call gemm &
        ( TR, 'N', nep, nep, n, ONE, keep%W, n, keep%W(1,1,3), n, &
          ZERO, keep%V(1,1,2), nep )
      k = ilaenv(1, 'DSYTRD', 'U', nep, -1, -1, -1)
      k = max(3*nep - 1, (k + 2)*nep)
      call dsygv &
        ( 1, 'V', 'U', nep, keep%V, nep, keep%V(1,1,2), nep, lambda, &
          keep%U, k, i )
      call gemm &
        ( 'N', 'N', n, nep, nep, ONE, keep%W, n, keep%V, nep, ZERO, X, ldX )

      if ( problem == 2 ) then
        ! recover buckling eigenvalues in the correct (ascending) order
        do i = 1, nep
          lambda(i) = -ONE/lambda(i)
        end do
      end if

      rci%job = rci%k

    end select

  end subroutine ssmfe_inverse_srci_double

  subroutine ssmfe_free_double( keep, inform )
    type(ssmfe_keepd), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: inform

    call ssmfe_free( keep )
    call ssmfe_free( inform )

  end subroutine ssmfe_free_double

  subroutine ssmfe_free_keep_double( keep )
    type(ssmfe_keepd), intent(inout) :: keep

    if ( allocated(keep%ind) ) deallocate ( keep%ind )
    if ( allocated(keep%V  ) ) deallocate ( keep%V   )
    if ( allocated(keep%W  ) ) deallocate ( keep%W   )
    if ( allocated(keep%U  ) ) deallocate ( keep%U   )
    if ( allocated(keep%BX ) ) deallocate ( keep%BX  )
    call ssmfe_free( keep%keep )

  end subroutine ssmfe_free_keep_double

  subroutine ssmfe_direct_srci_double_complex &
      ( problem, left, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, copy => zcopy, gemm => zgemm
    use SPRAL_lapack_iface, mxcopy => zlacpy

    character, parameter :: TR = 'C'

    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

    integer, parameter :: NONE = -1

    real(PRECISION), parameter :: ZERO = 0.0D0
    real(PRECISION), parameter :: ONE = 1.0D0

    complex(PRECISION), parameter :: NIL = ZERO
    complex(PRECISION), parameter :: UNIT = ONE

    integer, intent(in) :: problem
    integer, intent(in) :: left
    integer, intent(in) :: max_nep
    real(PRECISION), intent(inout) :: lambda(max_nep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    complex(PRECISION), intent(inout) :: X(ldX, max_nep)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    type(ssmfe_keepz  ), intent(inout), target :: keep

    integer :: extra
    integer :: total
    integer :: m
    integer :: kw
    integer :: ldBX
    integer :: i, j, k

    real(PRECISION), allocatable :: dwork(:,:)

    if ( n < 1 ) then
      inform%flag = WRONG_PROBLEM_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong problem size', n
      return
    end if

    if ( ldx < n ) then
      inform%flag = WRONG_LDX
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong leading dimension of x', ldx
      return
    end if

    if ( left < 0 .or. left > n/2 ) then
      inform%flag = WRONG_LEFT
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong number of eigenpairs', left
      return
    end if

    if ( left > max_nep ) then
      inform%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong eigenpairs storage size', max_nep
      return
    end if

    if ( rci%job == 0 ) then
      if ( left == 1 ) then
        extra = 7
      else
        extra = max(left/10, 10)
      end if
      total = min(max(1, left + extra), max(1, n/2 - 1))
      if ( options%max_left >= 0 ) &
        total = min(total, options%max_left)
      keep%block_size = total
      keep%lcon = 0
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate &
        ( keep%ind(keep%block_size), &
          keep%U(max_nep, 2*keep%block_size + max_nep), &
          keep%V(2*keep%block_size, 2*keep%block_size, 3), stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return
      if ( problem == 0 ) then
        kw = 5
      else
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep), stat = inform%stat )
        if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
        if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
        if ( inform%stat /= 0 ) return
        kw = 7
      end if
      allocate &
        ( keep%W(n, keep%block_size, 0:kw), dwork(n, keep%block_size), &
          stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return
      call random_number( dwork )
      keep%W(:,:,0) = 2*dwork - ONE
      deallocate ( dwork )
      if ( options%user_X > 0 ) &
        call mxcopy &
          ( 'A', n, min(keep%block_size, options%user_X), X, ldx, keep%W, n )
    else
      if ( .not. allocated(keep%ind) ) then
        inform%flag = WRONG_RCI_JOB
        rci%job = SSMFE_QUIT
        if ( options%unit_error > NONE .and. options%print_level > NONE ) &
          write( options%unit_error, '(/a, i3/)' ) &
            '??? Wrong rci%job', rci%job
        return
      end if
    end if

    ldBX = n

    do

      call ssmfe_solve &
        ( problem, left, max_nep, lambda, keep%block_size, keep%V, keep%ind, &
          rci, keep%keep, options, inform )

      select case ( rci%job )

      case ( 11:19 )

        ! apply some standard vector operations to vectors in W
        k = max(1, rci%k)
        call ssmfe_vector_operations_double_complex &
          ( rci, n, keep%block_size, keep%W, n, keep%V(1, 1, k), &
            keep%ind, keep%U )

      case ( SSMFE_SAVE_CONVERGED )

        if ( rci%nx < 1 ) cycle
        if ( rci%i < 0 ) cycle
        do i = 1, rci%nx
          k = (i - 1)*rci%i
          j = keep%lcon + i
          call copy( n, keep%W(1, rci%jx + k, 0), 1, X(1,j), 1 )
          if ( problem /= 0 ) &
            call copy( n, keep%W(1, rci%jy + k, rci%ky), 1, keep%BX(1,j), 1 )
        end do
        k = keep%lcon + 1
        m = keep%block_size
        if ( keep%lcon > 0 ) then
          if ( problem == 0 ) then
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                UNIT, X, ldX, X(1, k), ldX, &
                NIL, keep%U(1, m + k), max_nep )
          else
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                UNIT, X, ldX, keep%BX(1, k), ldX, &
                NIL, keep%U(1, m + k), max_nep )
          end if
          do j = 1, rci%nx
            do i = 1, keep%lcon
              keep%U(k + j - 1, m + i) = keep%U(i, m + k + j - 1)
            end do
          end do
        end if
        if ( problem == 0 ) then
          call gemm &
            ( TR, 'N', rci%nx, rci%nx, n, &
              UNIT, X(1, k), ldX, X(1, k), ldX, &
              NIL, keep%U(k, m + k), max_nep)
        else
          call gemm &
            ( TR, 'N', rci%nx, rci%nx, n, &
              UNIT, X(1, k), ldX, keep%BX(1, k), ldX, &
              NIL, keep%U(k, m + k), max_nep)
        end if
        keep%lcon = keep%lcon + rci%nx

      case ( SSMFE_APPLY_ADJ_CONSTRS )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
        if ( keep%lcon > 0 ) then
          call gemm &
            ( TR, 'N', keep%lcon, rci%nx, n, &
              UNIT, X, ldX, rci%x, n, NIL, keep%U, max_nep )
          if ( problem == 0 ) then
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -UNIT, X, ldX, keep%U, max_nep, UNIT, rci%x, n )
          else
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -UNIT, keep%BX, ldBX, keep%U, max_nep, UNIT, rci%x, n )
          end if
        end if

      case ( SSMFE_APPLY_CONSTRAINTS )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
        if ( keep%lcon > 0 ) then
          call gemm &
            ( TR, 'N', keep%lcon, rci%nx, n, &
              UNIT, X, ldX, rci%y, n, NIL, keep%U, max_nep )
          m = keep%block_size
          k = m + max_nep + 1
          call mxcopy &
            ( 'A', keep%lcon, rci%nx, keep%U, max_nep, keep%U(1, k), max_nep )
          call gemm &
            ( 'N', 'N', keep%lcon, rci%nx, keep%lcon, &
              -UNIT, keep%U(1, m + 1), max_nep, keep%U(1, k), max_nep, &
              2*UNIT, keep%U, max_nep )
          call gemm &
            ( 'N', 'N', n, rci%nx, keep%lcon, &
              -UNIT, X, ldX, keep%U, max_nep, UNIT, rci%x, n )
          if ( problem /= 0 ) &
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -UNIT, keep%BX, ldBX, keep%U, max_nep, UNIT, rci%y, n )
        end if

      case ( SSMFE_APPLY_A, SSMFE_APPLY_B )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)

        return

      case ( SSMFE_APPLY_PREC )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)
        call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

        return

      case ( SSMFE_RESTART )

        if ( rci%k == 0 ) then
          allocate ( dwork(n, keep%block_size), stat = inform%stat )
          if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
          if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
          if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
          if ( inform%stat /= 0 ) return
          call random_number( dwork )
          if ( rci%jx > 1 ) &
            keep%W(:, 1 : rci%jx - 1, 0) = 2*dwork(:, 1 : rci%jx - 1) - ONE
          deallocate ( dwork )
        end if

        return

      case ( : SSMFE_DONE )

        return

      end select

    end do

  end subroutine ssmfe_direct_srci_double_complex

  subroutine ssmfe_inverse_srci_double_complex &
      ( problem, sigma, left, right, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, copy => zcopy, gemm => zgemm
    use SPRAL_lapack_iface, mxcopy => zlacpy

    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE   = 9
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

    integer, parameter :: NONE = -1

    integer, intent(in) :: problem
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: max_nep
    real(PRECISION), intent(inout) :: lambda(max_nep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    complex(PRECISION), intent(inout) :: X(ldX, max_nep)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    type(ssmfe_keepz  ), intent(inout), target :: keep

    character, parameter :: TR = 'C'
    real(PRECISION), parameter :: ZERO = 0.0D0
    real(PRECISION), parameter :: ONE = 1.0D0
    complex(PRECISION), parameter :: NIL = ZERO
    complex(PRECISION), parameter :: UNIT = ONE

    integer :: nep
    integer :: extra_left, extra_right
    integer :: total_left, total_right
    integer :: m
    integer :: kw
    integer :: ldBX
    integer :: i, j, k

    integer :: ilaenv

    real(PRECISION) :: s

    real(PRECISION), allocatable :: dwork(:,:)

    complex(PRECISION) :: z

    if ( n < 1 ) then
      inform%flag = WRONG_PROBLEM_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong problem size', n
      return
    end if

    if ( ldx < n ) then
      inform%flag = WRONG_LDX
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong leading dimension of x', ldx
      return
    end if

    if ( left < 0 .or. left >= right .and. left + right > n/2 ) then
      inform%flag = WRONG_LEFT
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong number of eigenpairs on the left', left
      return
    end if

    if ( right < 0 .or. right > left .and. left + right > n/2 ) then
      inform%flag = WRONG_RIGHT
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong number of eigenpairs on the right', right
      return
    end if

    if ( left + right > max_nep ) then
      inform%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Wrong eigenpairs storage size', max_nep
      return
    end if

    if ( rci%job == 0 ) then
      keep%step = 0
      keep%lcon = 0
      keep%rcon = 0
      if ( left > 0 ) then
        extra_left = max(10, left/10)
      else
        extra_left = 0
      end if
      if ( right > 0 ) then
        extra_right = max(10, right/10)
      else
        extra_right = 0
      end if
      total_left = left + extra_left
      total_right = right + extra_right
      if ( options%max_left >= 0 ) &
        total_left = min(total_left, options%max_left)
      if ( options%max_right >= 0 ) &
        total_right = min(total_right, options%max_right)
      keep%block_size = max(2, total_left + total_right)
      keep%block_size = min(keep%block_size, max(1, n/2 - 1))
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate &
        ( keep%ind(keep%block_size), &
          keep%U(max_nep, 2*keep%block_size + max_nep), &
          keep%V(2*keep%block_size, 2*keep%block_size, 3), stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return
      keep%U = ZERO
      do i = 1, max_nep
        keep%U(i, keep%block_size + i) = ONE
      end do
      if ( problem == 0 ) then
        kw = 5
      else
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep), stat = inform%stat )
        if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
        if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
        if ( inform%stat /= 0 ) return
        kw = 7
      end if
      allocate &
        ( keep%W(n, keep%block_size, 0:kw), dwork(n, keep%block_size), &
          stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return
      call random_number( dwork )
      keep%W(:,:,0) = 2*dwork - ONE
      deallocate ( dwork )
      if ( options%user_X > 0 ) &
        call mxcopy &
          ( 'A', n, min(keep%block_size, options%user_X), X, ldx, keep%W, n )
    else
      if ( .not. allocated(keep%ind) ) then
        inform%flag = WRONG_RCI_JOB
        rci%job = SSMFE_QUIT
        if ( options%unit_error > NONE .and. options%print_level > NONE ) &
          write( options%unit_error, '(/a/)' ) &
            '??? Wrong rci%job'
        return
      end if
    end if

    ldBX = n
    nep = inform%left + inform%right

    select case ( keep%step )

    case ( 0 )

      do

        call ssmfe_solve &
          ( problem, sigma, left, right, &
            max_nep, lambda, keep%block_size, keep%V, keep%ind, &
            rci, keep%keep, options, inform )

        select case ( rci%job )

        case ( 11:19 )

          ! apply some standard vector operations to vectors in W
          k = max(1, rci%k)
          call ssmfe_vector_operations_double_complex &
            ( rci, n, keep%block_size, keep%W, n, keep%V(1, 1, k), &
              keep%ind, keep%U )

        case ( SSMFE_SAVE_CONVERGED )

          if ( rci%nx < 1 ) cycle

          do i = 1, rci%nx
            k = (i - 1)*rci%i
            if ( rci%i > 0 ) then
              j = keep%lcon + i
            else
              j = max_nep - keep%rcon - i + 1
            end if
            call copy( n, keep%W(1, rci%jx + k, 0), 1, X(1,j), 1 )
            if ( problem /= 0 ) &
              call copy( n, keep%W(1, rci%jy + k, rci%ky), 1, keep%BX(1,j), 1 )
          end do

          if ( rci%i > 0 ) then
            k = keep%lcon + 1
          else
            k = max_nep - keep%rcon - rci%nx + 1
          end if
          m = keep%block_size
          if ( keep%lcon > 0 ) then
            if ( problem == 0 ) then
              call gemm &
                ( TR, 'N', keep%lcon, rci%nx, n, &
                  UNIT, X, ldX, X(1, k), ldX, &
                  NIL, keep%U(1, m + k), max_nep )
            else
              call gemm &
                ( TR, 'N', keep%lcon, rci%nx, n, &
                  UNIT, X, ldX, keep%BX(1, k), ldX, &
                  NIL, keep%U(1, m + k), max_nep )
            end if
            do j = 1, rci%nx
              do i = 1, keep%lcon
                keep%U(k + j - 1, m + i) = keep%U(i, m + k + j - 1)
              end do
            end do
          end if
          if ( keep%rcon > 0 ) then
            if ( problem == 0 ) then
              call gemm &
                ( TR, 'N', keep%rcon, rci%nx, n, &
                  UNIT, X(1, max_nep - keep%rcon + 1), ldX, X(1, k), ldX, &
                  NIL, keep%U(max_nep - keep%rcon + 1, m + k), max_nep )
            else
              call gemm &
                ( TR, 'N', keep%rcon, rci%nx, n, &
                  UNIT, X(1, max_nep - keep%rcon + 1), ldX, &
                  keep%BX(1, k), ldX, &
                  NIL, keep%U(max_nep - keep%rcon + 1, m + k), max_nep )
            end if
            do j = 1, rci%nx
              do i = 1, keep%rcon
                keep%U(k + j - 1, m + max_nep - keep%rcon + i) = &
                  keep%U(max_nep - keep%rcon + i, m + k + j - 1)
              end do
            end do
          end if
          if ( problem == 0 ) then
            call gemm &
              ( TR, 'N', rci%nx, rci%nx, n, &
                UNIT, X(1, k), ldX, X(1, k), ldX, &
                NIL, keep%U(k, m + k), max_nep)
          else
            call gemm &
              ( TR, 'N', rci%nx, rci%nx, n, &
                UNIT, X(1, k), ldX, keep%BX(1, k), ldX, &
                NIL, keep%U(k, m + k), max_nep)
          end if

          if ( rci%i > 0 ) then
            keep%lcon = keep%lcon + rci%nx
          else
            keep%rcon = keep%rcon + rci%nx
          end if

        case ( SSMFE_APPLY_ADJ_CONSTRS )

          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

          if ( keep%lcon > 0 ) then
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                UNIT, X, ldX, rci%x, n, NIL, keep%U, max_nep )
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -UNIT, X, ldX, keep%U, max_nep, UNIT, rci%x, n )
          end if

          if ( keep%rcon > 0 ) then
            j = max_nep - keep%rcon + 1
            call gemm &
              ( TR, 'N', keep%rcon, rci%nx, n, &
                UNIT, X(1, j), ldX, rci%x, n, NIL, keep%U, max_nep )
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -UNIT, X(1,j), ldX, keep%U, max_nep, UNIT, rci%x, n )
          end if

        case ( SSMFE_APPLY_CONSTRAINTS )

          if ( keep%lcon == 0 .and. keep%rcon == 0 ) cycle

          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

          if ( keep%lcon > 0 ) then
            call gemm &
              ( TR, 'N', keep%lcon, rci%nx, n, &
                UNIT, X, ldX, rci%y, n, NIL, keep%U, max_nep )
          end if

          if ( keep%rcon > 0 ) then
            j = max_nep - keep%rcon + 1
            call gemm &
              ( TR, 'N', keep%rcon, rci%nx, n, &
                UNIT, X(1, j), ldX, rci%y, n, NIL, keep%U(j, 1), max_nep )
          end if

          m = keep%block_size
          k = m + max_nep + 1
          call mxcopy &
            ( 'A', max_nep, rci%nx, keep%U, max_nep, keep%U(1, k), max_nep )
          call gemm &
            ( 'N', 'N', max_nep, rci%nx, max_nep, &
              -UNIT, keep%U(1, m + 1), max_nep, keep%U(1, k), max_nep, &
              2*UNIT, keep%U, max_nep )

          if ( keep%lcon > 0 ) then
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%lcon, &
                -UNIT, X, ldX, keep%U, max_nep, UNIT, rci%x, n )
            if ( problem /= 0 ) &
              call gemm &
                ( 'N', 'N', n, rci%nx, keep%lcon, &
                  -UNIT, keep%BX, ldBX, keep%U, max_nep, UNIT, rci%y, n )
          end if

          if ( keep%rcon > 0 ) then
            j = max_nep - keep%rcon + 1
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -UNIT, X(1, j), ldX, keep%U(j, 1), max_nep, UNIT, rci%x, n )
            if ( problem /= 0 ) &
              call gemm &
                ( 'N', 'N', n, rci%nx, keep%rcon, &
                  -UNIT, keep%BX(1, j), ldBX, keep%U(j, 1), max_nep, &
                  UNIT, rci%y, n )
          end if

        case ( SSMFE_DO_SHIFTED_SOLVE, SSMFE_APPLY_B )

          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)

          return

        case ( SSMFE_APPLY_PREC )

          rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
          rci%y => keep%W(:, rci%jy : rci%jy + rci%nx - 1, rci%ky)

          call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

          return

        case ( SSMFE_RESTART )

          if ( rci%k == 0 ) then
            allocate ( dwork(n, keep%block_size), stat = inform%stat )
            if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
            if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
            if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
            if ( inform%stat /= 0 ) return
            call random_number( dwork )
            if ( rci%jx > 1 ) &
              keep%W(:, 1 : rci%jx - 1, 0) = 2*dwork(:, 1 : rci%jx - 1) - ONE
            if ( rci%jx + rci%nx - 1 < keep%block_size ) &
              keep%W(:, rci%jx + rci%nx : keep%block_size, 0) = &
                2*dwork(:, rci%jx + rci%nx : keep%block_size) - ONE
            deallocate ( dwork )
          end if

          return

        case ( SSMFE_DONE, SSMFE_QUIT )

          nep = inform%left + inform%right
          if ( nep < 1 ) return

          rci%k = rci%job

          if ( inform%left > 0 ) then
            do j = 1, inform%left/2
              s = lambda(j)
              lambda(j) = lambda(inform%left - j + 1)
              lambda(inform%left - j + 1) = s
              do i = 1, n
                z = X(i, j)
                X(i, j) = X(i, inform%left - j + 1)
                X(i, inform%left - j + 1) = z
              end do
            end do
          end if
          if ( inform%right > 0 ) then
            do j = 1, inform%right
              lambda(inform%left + j) = lambda(max_nep - j + 1)
            end do
            call mxcopy &
              ( 'A', n, inform%right, X(1, max_nep - inform%right + 1), ldX, &
                X(1, inform%left + 1), ldX )
            do j = 1, inform%right/2
              do i = 1, n
                z = X(i, inform%left + j)
                X(i, inform%left + j) = X(i, nep - j + 1)
                X(i, nep - j + 1) = z
              end do
            end do
          end if

          k = ilaenv(1, 'ZHETRD', 'U', nep, -1, -1, -1)
          k = max(2*nep - 1, (k + 2)*nep) + 3*nep - 1

          deallocate ( keep%W, keep%V, keep%U )
          allocate &
            ( keep%W(n, nep, 3), keep%V(nep, nep, 3), keep%U(k, 1), &
              stat = inform%stat )
          if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
          if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
          if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
          if ( inform%stat /= 0 ) return

          if ( problem == 0 ) then
            call mxcopy( 'A', n, nep, X, ldX, keep%W, n )
            rci%job = SSMFE_DO_SHIFTED_SOLVE
            rci%nx = nep
            rci%jx = 1
            rci%kx = 0
            rci%jy = 1
            rci%ky = 2
            rci%x => keep%W(:,:,1)
            rci%y => keep%W(:,:,2)
            keep%step = 2
          else
            call mxcopy( 'A', n, nep, X, ldX, keep%W(1,1,2), n )
            rci%job = SSMFE_APPLY_B
            rci%nx = nep
            rci%jx = 1
            rci%kx = 2
            rci%jy = 1
            rci%ky = 0
            rci%x => keep%W(:,:,2)
            rci%y => keep%W(:,:,1)
            keep%step = 1
          end if

          return

        case ( SSMFE_ABORT )

          return

        end select

      end do

    case ( 1 )

      rci%job = SSMFE_DO_SHIFTED_SOLVE
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 2
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,2)
      keep%step = 2

    case ( 2 )

      call mxcopy( 'A', n, nep, keep%W(1,1,2), n, keep%W, n )
      if ( problem == 0 ) &
        call mxcopy( 'A', n, nep, keep%W, n, keep%W(1,1,3), n )
      rci%job = SSMFE_APPLY_A
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 2
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,2)
      if ( problem == 0 ) then
        keep%step = 4
      else
        keep%step = 3
      end if

    case ( 3 )

      rci%job = SSMFE_APPLY_B
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 4
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,3)
      keep%step = 4

    case ( 4 )

      if ( problem == 2 ) then
        z = -UNIT
      else
        z = UNIT
      end if
      call gemm &
        ( TR, 'N', nep, nep, n, z, keep%W, n, keep%W(1,1,2), n, &
          NIL, keep%V, nep )
      call gemm &
        ( TR, 'N', nep, nep, n, UNIT, keep%W, n, keep%W(1,1,3), n, &
          NIL, keep%V(1,1,2), nep )

      k = ilaenv(1, 'ZHETRD', 'U', nep, -1, -1, -1)
      k = max(2*nep - 1, (k + 2)*nep) + 3*nep - 1
      allocate ( dwork(max(3*nep - 2, 1), 1), stat = inform%stat )
      if ( inform%stat /= 0 ) inform%flag = OUT_OF_MEMORY
      if ( inform%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( inform%stat /= 0 ) call ssmfe_errmsg( options, inform )
      if ( inform%stat /= 0 ) return
      call zhegv &
        ( 1, 'V', 'U', nep, keep%V, nep, keep%V(1,1,2), nep, lambda, &
          keep%U(3*nep, 1), k, dwork, i )
      deallocate ( dwork )
      call gemm &
        ( 'N', 'N', n, nep, nep, UNIT, keep%W, n, keep%V, nep, NIL, X, ldX )

      if ( problem == 2 ) then
        do i = 1, nep
          lambda(i) = -ONE/lambda(i)
        end do
      end if

      rci%job = rci%k

    end select

  end subroutine ssmfe_inverse_srci_double_complex

  subroutine ssmfe_vector_operations_double( rci, n, m, W, ldW, V, ind, U )

    use SPRAL_blas_iface, &
      norm => dnrm2, &
      dot  => ddot,  &
      scal => dscal, &
      axpy => daxpy, &
      gemm => dgemm

    use SPRAL_lapack_iface, mxcopy => dlacpy
    integer :: n, m, ldW
    real(PRECISION) :: W(ldW, m, 0:*), V(m + m, m + m), U(m)
    integer :: ind(m)
    type(ssmfe_rcid) :: rci
    intent(in) :: n, m, ldW, rci
    intent(inout) :: W, V, U, ind

    integer, parameter :: SSMFE_COPY_VECTORS       = 11
    integer, parameter :: SSMFE_COMPUTE_DOTS       = 12
    integer, parameter :: SSMFE_SCALE_VECTORS      = 13
    integer, parameter :: SSMFE_COMPUTE_YMXD       = 14
    integer, parameter :: SSMFE_COMPUTE_XY         = 15
    integer, parameter :: SSMFE_COMPUTE_XQ         = 16
    integer, parameter :: SSMFE_TRANSFORM_X        = 17

    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0
    real(PRECISION), parameter :: NIL = ZERO
    real(PRECISION), parameter :: UNIT = ONE

    integer :: i, j, mm

    real(PRECISION) :: alpha, beta, s
    character, parameter :: TRANS = 'T'

    mm = m + m

    if ( rci%job == SSMFE_TRANSFORM_X ) then
      alpha = UNIT
      beta = NIL
    else
      alpha = rci%alpha
      beta = rci%beta
    end if

    select case (rci%job)

    case (SSMFE_COPY_VECTORS)

      if ( rci%nx < 1 ) return

      if ( rci%i == 0 ) then

        call mxcopy &
          ( 'A', n, rci%nx, &
            W(1, rci%jx, rci%kx), ldW, &
            W(1, rci%jy, rci%ky), ldW )

      else

        do i = 1, n
          do j = 1, rci%nx
            U(j) = W(i, ind(j), rci%kx)
          end do
          do j = 1, rci%nx
            W(i, j, rci%kx) = U(j)
          end do
          if ( rci%ky /= rci%kx ) then
            do j = 1, rci%nx
              U(j) = W(i, ind(j), rci%ky)
            end do
            do j = 1, rci%nx
              W(i, j, rci%ky) = U(j)
            end do
          end if
        end do

      end if

    case (SSMFE_COMPUTE_DOTS)

      do i = 0, rci%nx - 1
        V(rci%i + i, rci%j + i) = &
          dot( n, W(1, rci%jx + i, rci%kx), 1, &
                  W(1, rci%jy + i, rci%ky), 1 )
      end do

    case (SSMFE_SCALE_VECTORS)

      do i = 0, rci%nx - 1
        if ( rci%kx == rci%ky ) then
          s = norm(n, W(1, rci%jx + i, rci%kx), 1)
          if ( s > 0 ) &
            call scal( n, ONE/s, W(1, rci%jx + i, rci%kx), 1 )
        else
          s = sqrt(abs(dot &
            (n, W(1, rci%jx + i, rci%kx), 1, W(1, rci%jy + i, rci%ky), 1)))
          if ( s > 0 ) then
            call scal( n, ONE/s, W(1, rci%jx + i, rci%kx), 1 )
            call scal( n, ONE/s, W(1, rci%jy + i, rci%ky), 1 )
          end if
        end if
      end do

    case (SSMFE_COMPUTE_YMXD)

      do i = 0, rci%nx - 1
        s = -V(rci%i + i, rci%j + i)
        call axpy( n, s, W(1, rci%jx + i, rci%kx), 1, &
                         W(1, rci%jy + i, rci%ky), 1 )
      end do ! i

    case (SSMFE_COMPUTE_XY)

      if ( rci%nx < 1 .or. rci%ny < 1 ) return
      call gemm( TRANS, 'N', rci%nx, rci%ny, n, &
                 alpha, W(1, rci%jx, rci%kx), ldW, W(1, rci%jy, rci%ky), ldW, &
                 beta, V(rci%i, rci%j), mm )

    case (SSMFE_COMPUTE_XQ, SSMFE_TRANSFORM_X)

      if ( rci%ny < 1 ) return
      call gemm &
        ( 'N', 'N', n, rci%ny, rci%nx, &
          alpha, W(1, rci%jx, rci%kx), ldW, V(rci%i, rci%j), mm, &
          beta, W(1, rci%jy, rci%ky), ldW )
      if ( rci%job == SSMFE_TRANSFORM_X ) &
        call mxcopy &
          ( 'A', n, rci%ny, W(1, rci%jy, rci%ky), ldW, &
                            W(1, rci%jx, rci%kx), ldW )

    end select

  end subroutine ssmfe_vector_operations_double

  subroutine ssmfe_vector_operations_double_complex &
      ( rci, n, m, W, ldW, V, ind, U )

    use SPRAL_blas_iface, &
      norm => dznrm2, &
      dot  => zdotc, &
      scal => zscal, &
      axpy => zaxpy, &
      gemm => zgemm

    use SPRAL_lapack_iface, mxcopy => zlacpy
    integer :: n, m, ldW
    complex(PRECISION) :: W(ldW, m, 0:*), V(m + m, m + m), U(m)
    integer :: ind(m)
    type(ssmfe_rciz) :: rci
    intent(in) :: n, m, ldW, rci
    intent(inout) :: W, V, U, ind

    character, parameter :: TRANS = 'C'
    integer, parameter :: SSMFE_COPY_VECTORS       = 11
    integer, parameter :: SSMFE_COMPUTE_DOTS       = 12
    integer, parameter :: SSMFE_SCALE_VECTORS      = 13
    integer, parameter :: SSMFE_COMPUTE_YMXD       = 14
    integer, parameter :: SSMFE_COMPUTE_XY         = 15
    integer, parameter :: SSMFE_COMPUTE_XQ         = 16
    integer, parameter :: SSMFE_TRANSFORM_X        = 17

    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0
    complex(PRECISION), parameter :: NIL(1) = ZERO
    complex(PRECISION), parameter :: UNIT = ONE

    integer :: i, j, mm

    real(PRECISION) ::  s
    complex(PRECISION) :: alpha, beta

    mm = m + m

    if ( rci%job == SSMFE_TRANSFORM_X ) then
      alpha = UNIT
      beta = NIL(1)
    else
      alpha = rci%alpha
      beta = rci%beta
    end if

    select case (rci%job)

    case (SSMFE_COPY_VECTORS)

      if ( rci%nx < 1 ) return

      if ( rci%i == 0 ) then

        call mxcopy &
          ( 'A', n, rci%nx, &
            W(1, rci%jx, rci%kx), ldW, &
            W(1, rci%jy, rci%ky), ldW )

      else

        do i = 1, n
          do j = 1, rci%nx
            U(j) = W(i, ind(j), rci%kx)
          end do
          do j = 1, rci%nx
            W(i, j, rci%kx) = U(j)
          end do
          if ( rci%ky /= rci%kx ) then
            do j = 1, rci%nx
              U(j) = W(i, ind(j), rci%ky)
            end do
            do j = 1, rci%nx
              W(i, j, rci%ky) = U(j)
            end do
          end if
        end do

      end if

    case (SSMFE_COMPUTE_DOTS)

      do i = 0, rci%nx - 1
        V(rci%i + i, rci%j + i) = &
          dot( n, W(1, rci%jx + i, rci%kx), 1, &
                  W(1, rci%jy + i, rci%ky), 1 )
      end do

    case (SSMFE_SCALE_VECTORS)

      do i = 0, rci%nx - 1
        if ( rci%kx == rci%ky ) then
          s = norm(n, W(1, rci%jx + i, rci%kx), 1)
          if ( s > 0 ) &
            call scal( n, UNIT/s, W(1, rci%jx + i, rci%kx), 1 )
        else
          s = sqrt(abs(dot &
            (n, W(1, rci%jx + i, rci%kx), 1, W(1, rci%jy + i, rci%ky), 1)))
          if ( s > 0 ) then
            call scal( n, UNIT/s, W(1, rci%jx + i, rci%kx), 1 )
            call scal( n, UNIT/s, W(1, rci%jy + i, rci%ky), 1 )
          end if
        end if
      end do

    case (SSMFE_COMPUTE_YMXD)

      do i = 0, rci%nx - 1
        call axpy &
          ( n, -V(rci%i + i, rci%j + i), &
            W(1, rci%jx + i, rci%kx), 1, W(1, rci%jy + i, rci%ky), 1 )
      end do ! i

    case (SSMFE_COMPUTE_XY)

      if ( rci%nx < 1 .or. rci%ny < 1 ) return
      call gemm &
        ( TRANS, 'N', rci%nx, rci%ny, n, &
          alpha, W(1, rci%jx, rci%kx), ldW, W(1, rci%jy, rci%ky), ldW, &
          beta, V(rci%i, rci%j), mm )

    case (SSMFE_COMPUTE_XQ, SSMFE_TRANSFORM_X)

      if ( rci%ny < 1 ) return
      call gemm &
        ( 'N', 'N', n, rci%ny, rci%nx, &
          alpha, W(1, rci%jx, rci%kx), ldW, V(rci%i, rci%j), mm, &
          beta, W(1, rci%jy, rci%ky), ldW )
      if ( rci%job == SSMFE_TRANSFORM_X ) &
        call mxcopy &
          ( 'A', n, rci%ny, W(1, rci%jy, rci%ky), ldW, &
                            W(1, rci%jx, rci%kx), ldW )

    end select

  end subroutine ssmfe_vector_operations_double_complex

  subroutine ssmfe_free_double_complex( keep, inform )
    type(ssmfe_keepz), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: inform

    call ssmfe_free( keep )
    call ssmfe_free( inform )

  end subroutine ssmfe_free_double_complex

  subroutine ssmfe_free_keep_double_complex( keep )
    type(ssmfe_keepz), intent(inout) :: keep

    if ( allocated(keep%ind) ) deallocate ( keep%ind )
    if ( allocated(keep%BX ) ) deallocate ( keep%BX  )
    if ( allocated(keep%U  ) ) deallocate ( keep%U   )
    if ( allocated(keep%V  ) ) deallocate ( keep%V   )
    if ( allocated(keep%W  ) ) deallocate ( keep%W   )
    call ssmfe_free( keep%keep )

  end subroutine ssmfe_free_keep_double_complex

end module SPRAL_ssmfe


