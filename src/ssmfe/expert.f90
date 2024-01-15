! COPYRIGHT (c) 2013      Evgueni Ovtchinnikov
! COPYRIGHT (c) 2014,2015 The Science and Technology Facilities Council (STFC)
! Version 1.0.0
!
! Written by: Evgueni Ovtchinnikov
!
module SPRAL_ssmfe_expert
  use SPRAL_ssmfe_core , only: &
    ssmfe, ssmfe_free, &
    ssmfe_core_keep, ssmfe_core_options, ssmfe_inform, ssmfe_rcid, ssmfe_rciz
  implicit none

  private

  ! double precision to be used
  integer, parameter :: PRECISION = kind(1.0D0)

  ! input error flags
  integer, parameter :: WRONG_BLOCK_SIZE   =  -2
  integer, parameter :: WRONG_ERR_EST      =  -3
  integer, parameter :: WRONG_MINPROD      =  -4
  integer, parameter :: WRONG_LEFT         = -11
  integer, parameter :: WRONG_RIGHT        = -12
  integer, parameter :: WRONG_STORAGE_SIZE = -13
  integer, parameter :: WRONG_SIGMA        = -14

  ! fatal execution error flags
  integer, parameter :: OUT_OF_MEMORY           = -100
  integer, parameter :: B_NOT_POSITIVE_DEFINITE = -200

  ! warning flags
  integer, parameter :: NO_SEARCH_DIRECTIONS_LEFT   = 1
  integer, parameter :: MAX_NUM_ITERATIONS_EXCEEDED = 2
  integer, parameter :: OUT_OF_STORAGE              = 3

  ! termination status
  integer, parameter :: SSMFE_DONE  = -1 ! success
  integer, parameter :: SSMFE_QUIT  = -2 ! non-fatal error
  integer, parameter :: SSMFE_ABORT = -3 ! fatal error

  ! default error estimation scheme, see the spec
  integer, parameter :: SSMFE_KINEMATIC = 2

  public :: ssmfe_standard, ssmfe_standard_shift
  public :: ssmfe_generalized, ssmfe_generalized_shift
  public :: ssmfe_buckling
  public :: ssmfe_solve
  public :: ssmfe_free
  public :: ssmfe_errmsg
  public :: ssmfe_options, ssmfe_expert_keep
  public :: ssmfe_rcid, ssmfe_rciz, ssmfe_inform

  interface ssmfe_standard
    module procedure &
      ssmfe_expert_double, &
      ssmfe_expert_double_complex
  end interface

  interface ssmfe_generalized
    module procedure &
      ssmfe_expert_gen_double, &
      ssmfe_expert_gen_double_complex
  end interface

  interface ssmfe_standard_shift
    module procedure &
      ssmfe_expert_shift_double, &
      ssmfe_expert_shift_double_complex
  end interface

  interface ssmfe_generalized_shift
    module procedure &
      ssmfe_expert_gen_shift_double, &
      ssmfe_expert_gen_shift_double_complex
  end interface

  interface ssmfe_buckling
    module procedure &
      ssmfe_expert_buckling_double, &
      ssmfe_expert_buckling_double_complex
  end interface

  interface ssmfe_solve
    module procedure &
      ssmfe_inverse_rci_double, &
      ssmfe_direct_rci_double, &
      ssmfe_inverse_rci_double_complex, &
      ssmfe_direct_rci_double_complex
  end interface

  interface ssmfe_free
    module procedure &
      ssmfe_expert_free_double, &
      ssmfe_expert_free_keep_double
  end interface

  type ssmfe_options ! see the spec

    integer :: print_level = 0

    integer :: unit_error      = 6
    integer :: unit_warning    = 6
    integer :: unit_diagnostic = 6

    integer :: max_iterations = 100

    integer :: user_x = 0

    integer :: err_est = SSMFE_KINEMATIC

    real(PRECISION) :: abs_tol_lambda   =  0.0
    real(PRECISION) :: rel_tol_lambda   =  0.0
    real(PRECISION) :: abs_tol_residual =  0.0
    real(PRECISION) :: rel_tol_residual =  0.0
    real(PRECISION) :: tol_x            = -1.0

    real(PRECISION) :: left_gap  = 0.0
    real(PRECISION) :: right_gap = 0.0

    integer :: extra_left  = -1
    integer :: extra_right = -1

    integer :: max_left  = -1
    integer :: max_right = -1

    logical :: minAprod = .true.
    logical :: minBprod = .true.

  end type ssmfe_options

  ! private subroutines and types

  type ssmfe_expert_keep

    private

    ! convergence flags for eigenvalues
    logical :: left_converged  = .false. ! leftmost/left-of-the-shift
    logical :: right_converged = .false. ! right of the shift

    ! problem type (0: standard, 1: generalized, 2: buckling)
    integer :: problem = 0

    ! numbers of wanted eigenvalues
    integer :: left  = 0 ! leftmost/left-of-the-shift
    integer :: right = 0 ! right of the shift

    ! total actual numbers of eigenvalues
    integer :: max_left  = -1 ! left of the shift
    integer :: max_right = -1 ! right of the shift

    ! numbers of converged eigenvalues
    integer :: lcon = 0 ! leftmost/left-of-the-shift
    integer :: rcon = 0 ! right of the shift

    ! positions of the first and last non-converged eigenvector
    ! in block 0 of the vector workspace
    integer :: first = 1
    integer :: last  = 0

    ! estimated average distance between eigenvalues
    real(PRECISION) :: av_dist = 0

    ! work array for eigenvalues: converged eigenvalues are copied from
    ! this array to the eigenvalue storage
    real(PRECISION), allocatable :: lmd(:)

    type(ssmfe_inform      ) :: info ! information about execution - see spec
    type(ssmfe_core_keep   ) :: keep ! core interface keep - see core.f90
    type(ssmfe_core_options) :: options ! core interface options - see core spec

  end type ssmfe_expert_keep

contains

!*************************************************************************
! real RCI for computing leftmost eigenpairs of A x = lambda x
! see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_double &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep

    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    real(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_expert_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform

    call ssmfe_direct_rci_double &
      ( 0, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_double

!*************************************************************************
! complex RCI for computing leftmost eigenpairs of A x = lambda x,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_double_complex &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    complex(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rciz       ), intent(inout) :: rci
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_direct_rci_double_complex &
      ( 0, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_double_complex

!*************************************************************************
! real RCI for computing leftmost eigenpairs of A x = lambda B x, B > 0,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_gen_double &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep

    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    real(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rcid       ), intent(inout) :: rci
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_direct_rci_double &
      ( 1, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_double

!*************************************************************************
! complex RCI for computing leftmost eigenpairs of A x = lambda B x, B > 0,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_gen_double_complex &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    complex(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rciz       ), intent(inout) :: rci
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_direct_rci_double_complex &
      ( 1, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_double_complex

!*************************************************************************
! real RCI for computing eigenpairs {lambda, x} of A x = lambda x
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_shift_double &
      ( rci, sigma, left, right, mep, lambda, m, rr, ind, &
        keep, options, inform )

    type(ssmfe_rcid), intent(inout) :: rci
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    real(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_inverse_rci_double &
      ( 0, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_shift_double

!*************************************************************************
! complex RCI for computing eigenpairs {lambda, x} of A x = lambda x
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_shift_double_complex &
      ( rci, sigma, left, right, mep, lambda, m, rr, ind, &
        keep, options, inform )

    type(ssmfe_rciz), intent(inout) :: rci
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    complex(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_inverse_rci_double_complex &
      ( 0, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_shift_double_complex

!**************************************************************************
! real RCI for computing eigenpairs {lambda, x} of A x = lambda B x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_gen_shift_double &
      ( rci, sigma, left, right, mep, lambda, m, rr, ind, &
        keep, options, inform )

    type(ssmfe_rcid), intent(inout) :: rci
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    real(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_inverse_rci_double &
      ( 1, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_shift_double

!*****************************************************************************
! complex RCI for computing eigenpairs {lambda, x} of A x = lambda B x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_gen_shift_double_complex &
      ( rci, sigma, left, right, mep, lambda, m, rr, ind, &
        keep, options, inform )

    type(ssmfe_rciz), intent(inout) :: rci
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    complex(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_inverse_rci_double_complex &
      ( 1, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_shift_double_complex

!**************************************************************************
! real RCI for computing eigenpairs {lambda, x} of B x = lambda A x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_buckling_double &
      ( rci, sigma, left, right, mep, lambda, m, rr, ind, &
        keep, options, inform )

    type(ssmfe_rcid), intent(inout) :: rci
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    real(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_inverse_rci_double &
      ( 2, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_buckling_double

!*****************************************************************************
! complex RCI for computing eigenpairs {lambda, x} of B x = lambda A x, B > 0,
! with lambda near the shift sigma, see ssmfe spec for the arguments
!
  subroutine ssmfe_expert_buckling_double_complex &
      ( rci, sigma, left, right, mep, lambda, m, rr, ind, &
        keep, options, inform )

    type(ssmfe_rciz), intent(inout) :: rci
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    complex(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: inform

    call ssmfe_inverse_rci_double_complex &
      ( 2, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_buckling_double_complex

!*************************************************************************
! real RCI solver procedure for eigenvalues near the shift sigma
!
  subroutine ssmfe_inverse_rci_double &
    ( problem, sigma, left, right, &
      max_nep, lambda, block_size, rr_matrices, ind, &
      rci, keep, options, info )
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_APPLY_A           = 1
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE  = 9
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0

    real(kind = PRECISION), parameter :: MIN_REL_GAP = 1E-3

    ! problem type
    ! 0: A x = lambda x
    ! 1: A x = lambda B x
    ! 2: B x = lambda A x
    ! A = A', B = B' > 0
    integer, intent(in) :: problem

    ! shift
    real(PRECISION), intent(in) :: sigma

    ! number of wanted eigenvalues left of sigma
    integer, intent(in) :: left

    ! number of wanted eigenvalues right of sigma
    integer, intent(in) :: right

    ! eigenpairs storage size
    integer, intent(in) :: max_nep

    ! storage for eigenvalues
    real(PRECISION), intent(inout) :: lambda(max_nep)

    ! BJCG block size
    integer, intent(in) :: block_size

    ! matrices for the Rayleigh-Ritz procedure
    real(PRECISION), intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)

    ! work matrix column reordering index
    integer, intent(inout) :: ind(block_size)

    ! ssmfe_expert types - see the spec
    type(ssmfe_rcid       ), intent(inout) :: rci
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, frmt, line
    character(7) :: word

    integer :: u_diag
    integer :: m, mep, new, ncon, first, last, step
    integer :: i, j, k, l

    real(PRECISION) :: left_gap
    real(PRECISION) :: right_gap
    real(PRECISION) :: delta
    real(PRECISION) :: q, r, s, t

    ! short names for options

    u_diag    = options%unit_diagnostic
    left_gap  = options%left_gap
    right_gap = options%right_gap

    info%flag = 0

    ! check for the input data errors

    if ( left < 0 &
        .or. options%max_left >= 0 .and. options%max_left < left ) then
      info%flag = WRONG_LEFT
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( right < 0 &
        .or. options%max_right >= 0 .and. options%max_right < right ) then
      info%flag = WRONG_RIGHT
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( left + right > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( block_size < 2 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( problem == 2 .and. sigma == ZERO ) then
      info%flag = WRONG_SIGMA
      rci%job = SSMFE_QUIT
      call ssmfe_errmsg( options, info )
      return
    end if

    ! things to do on (re)start
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        ! reallocate the information arrays
        mep = max_nep
        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda    ) ) deallocate ( info%err_lambda     )
        if ( allocated(info%err_X         ) ) deallocate ( info%err_X          )
        if ( allocated(info%converged     ) ) deallocate ( info%converged      )
        allocate &
          ( info%residual_norms(mep), info%err_lambda(mep), &
            info%err_X(mep), info%converged(mep), stat = info%stat )
        if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
        if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
        if ( info%stat /= 0 ) return

        ! initialize the execution information
        info%iteration = 0
        info%residual_norms = ZERO
        info%err_lambda = -ONE
        info%err_X = -ONE
        info%converged = 0
        info%left = 0
        info%right = 0
        info%next_left = sigma
        info%next_right = sigma

        ! nothing to do
        if ( left <= 0 .and. right <= 0 ) then
          rci%job = SSMFE_DONE
          return
        end if

        ! save the input data
        keep%problem = -abs(problem)
        keep%options%err_est = options%err_est
        keep%left = left
        keep%right = right
        keep%options%minAprod = options%minAprod
        keep%options%minBprod = options%minBprod

        ! compute the number of extra vectors iterated for the sake
        ! of better convergence
        if ( left < 1 ) then
          keep%options%extra_left = 0
        else if ( options%extra_left < 0 ) then
          keep%options%extra_left = max(10, left/10) ! use default
        else
          keep%options%extra_left = options%extra_left ! take user's value
        end if
        if ( right < 1 ) then
          keep%options%extra_right = 0
        else if ( options%extra_right < 0 ) then
          keep%options%extra_right = max(10, right/10) ! use default
        else
          keep%options%extra_right = options%extra_right ! take user's value
        end if

        if ( options%max_left >= 0 ) keep%options%extra_left = &
          min(keep%options%extra_left, options%max_left - left)
        if ( options%max_right >= 0 ) keep%options%extra_right = &
          min(keep%options%extra_right, options%max_right - right)

      end if

      m = block_size
      ! if the total number of eigenvalues to be computed (wanted + extra)
      ! exceeds the block size m, they are computed in portions
      ! of sizes keep%left (those left of sigma) and keep%right (right of it);
      ! once a portion is computed, iterations are restarted (see below)
      ! (note that more than wanted + extra eigenpairs are needed if gaps in
      ! the spectrum separating wanted eigenvalues from the rest are sought)

      if ( keep%left > 0 ) then
        ! the total numers of eigenvalues left of sigma to be computed
        i = left + keep%options%extra_left
        if ( keep%right > 0 ) then
          ! the total numers of eigenvalues right of sigma to be computed
          j = right + keep%options%extra_right
          ! workspace for eigenvectors corresponding to left eigenvalues
          l = m*i/(i + j)
          keep%left = left ! assumingly can be computed in one go
          if ( l < i .or. left_gap /= ZERO ) then ! not enough space
            keep%options%extra_left = min(keep%options%extra_left, l/2)
            keep%left = l - keep%options%extra_left ! portion size
            ! limiting extras number to l/2 ensures that the portion size
            ! is substantial thus reducing the number of restarts
          end if
          keep%right = right ! assumingly can be computed in one go
          if ( m - l < j .or. right_gap /= ZERO ) then ! not enough space
            keep%options%extra_right = min(keep%options%extra_right, (m - l)/2)
            keep%right = m - l - keep%options%extra_right ! portion size
          end if
        else
          keep%left = left
          if ( m < i .or. left_gap /= ZERO ) then
            keep%options%extra_left = min(keep%options%extra_left, m/2)
            keep%left = m - keep%options%extra_left
          end if
          keep%options%extra_right = 0
          keep%right = 0
        end if
      else
        keep%options%extra_left = 0
        if ( m < right + keep%options%extra_right .or. right_gap /= ZERO ) then
          keep%options%extra_right = min(keep%options%extra_right, m/2)
          keep%right = m - keep%options%extra_right
        end if
      end if
      ! if gaps in the spectrum are sought, the number of computed
      ! eigenpairs must be larger than requested
      if ( left > 0 .and. keep%left == left .and. left_gap /= ZERO &
          .and. keep%options%extra_left > 0 ) then
        keep%left = keep%left + 1
        keep%options%extra_left = keep%options%extra_left - 1
      end if
      if ( right > 0 .and. keep%right == right .and. right_gap /= ZERO &
          .and. keep%options%extra_right > 0 ) then
        keep%right = keep%right + 1
        keep%options%extra_right = keep%options%extra_right - 1
      end if

      if ( rci%job == SSMFE_START ) then

        keep%options%min_gap = 0.05 ! see ssmfe_core spec

        keep%lcon = 0 ! the number of converged eigenpairs left of sigma
        keep%rcon = 0 ! the number of converged eigenpairs right of sigma

        ! initialize the convergence flags
        keep%left_converged = left == 0
        keep%right_converged = right == 0

        if ( options%max_left < 0 ) then ! the number of eigenvalues < sigma
          keep%max_left = max_nep + 1     ! not known, set to no-impact value
        else
          keep%max_left = options%max_left
        end if
        if ( options%max_right < 0 ) then ! ditto on the right
          keep%max_right = max_nep + 1
        else
          keep%max_right = options%max_right
        end if

      else

        if ( allocated(keep%lmd) ) deallocate ( keep%lmd )

      end if

      ! allocate eigenvalues array of ssmfe_core solver
      m = block_size
      allocate ( keep%lmd(m), stat = info%stat )
      if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
      if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
      if ( info%stat /= 0 ) return

      keep%first = 1 ! first non-converged eigenpair index
      keep%last  = m ! last non-converged eigenpair index

      ! print problem summary message
      call ssmfe_msg( problem, options, left, right, m )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) &
        rci%job = SSMFE_START

    else if ( rci%job == SSMFE_SAVE_CONVERGED .and. rci%i < 0 ) then

      ! some converged eigenpairs just saved to storage,
      ! decide on whether to go for more

      mep = max_nep
      m = block_size
      ncon = keep%lcon + keep%rcon

      if ( .not. (keep%left_converged .and. keep%right_converged) ) then

        ! check the number of iterations performed
        if ( info%iteration >= options%max_iterations ) then
          rci%job = SSMFE_QUIT
          info%flag = MAX_NUM_ITERATIONS_EXCEEDED

        ! check the remaining storage
        else if ( ncon == mep &
          .or. .not. keep%left_converged &
               .and. keep%lcon >= mep - max(right, keep%rcon) &
          .or. .not. keep%right_converged &
               .and. keep%rcon >= mep - max(left, keep%lcon) &
            ) then
          rci%job = SSMFE_QUIT
          info%flag = OUT_OF_STORAGE
        end if

      end if

      if ( rci%job /= SSMFE_QUIT ) then
        if ( keep%left_converged ) keep%left = 0 ! no more left eigenvalues
        if ( keep%right_converged ) keep%right = 0 ! ditto right
        if ( keep%left_converged .and. keep%right_converged ) then
          rci%job = SSMFE_DONE
        else
          if ( (keep%left > 0 .and. keep%first > keep%left .or. &
                keep%right > 0 .and. block_size - keep%last >= keep%right) &
              ) then
            ! the current portion of eigenpairs converged,
            ! but not all wanted eigenpairs have been computed,
            ! restart to compute the next portion
            rci%job = SSMFE_RESTART
            rci%jx = keep%first
            rci%nx = keep%last - keep%first + 1
            rci%k = 0
            return
          end if
        end if
      end if

      if ( rci%job < 0 ) then ! computation terminated
        info%non_converged = &
          max(0, left - keep%lcon) + max(0, right - keep%rcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
          info%right = keep%rcon
        end if
        if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )
        return
      end if

    end if

    ! call the core solver
1   call ssmfe &
      ( rci, keep%problem, keep%left, keep%right, block_size, &
        keep%lmd, rr_matrices, ind, keep%keep, keep%options, keep%info )

    select case (rci%job)

    case (SSMFE_APPLY_A)

      rci%job = SSMFE_DO_SHIFTED_SOLVE

    case (SSMFE_SAVE_CONVERGED)

      ! the number of newly converged eigenpairs
      new = rci%nx
      ! number of eigenpairs converged previously
      ncon = keep%lcon + keep%rcon
      mep = max_nep

      ! rci%i > 0: save left converged eigenvalues data
      ! rci%i < 0: save right converged eigenvalues data

      ! if the number of converged eigenpairs exceeds the storage size
      ! reduce the number of newly converged eigenpairs to be saved
      if ( rci%i > 0 ) then
        if ( left_gap == ZERO ) new = min(new, left - keep%lcon)
        k = max(right, keep%rcon) ! space reserved for the right eigenpairs
        if ( new > mep - k - keep%lcon ) new = mep - k - keep%lcon
      else
        if ( right_gap == ZERO ) new = min(new, right - keep%rcon)
        k = max(left, keep%lcon) ! space reserved for the left eigenpairs
        if ( new > mep - k - keep%rcon ) new = mep - k - keep%rcon
      end if
      rci%nx = new

      ! save the newly converged eigenpairs and related data
      ! to the eigenvector storage and information arrays
      do i = 1, new
        k = (i - 1)*rci%i
        if ( rci%i > 0 ) then
          j = keep%lcon + i ! left eigenvalues data are saved on the left
        else
          j = mep - keep%rcon - i + 1 ! right data sre saved on the right
        end if
        s = keep%lmd(rci%jx + k)
        t = keep%info%err_lambda(rci%jx + k)
        info%residual_norms(j) = keep%info%residual_norms(rci%jx + k)
        if ( t < 0 ) then ! eigenvalue error not available yet
          info%err_lambda(j) = t
        else ! apply reverse shift-invert mapping
          info%err_lambda(j) = si_map_err(problem, sigma, s, t)
        end if
        info%err_x(j) = keep%info%err_x(rci%jx + k)
        l = keep%info%converged(rci%jx + k)
        if ( l > 0 ) then ! converged eigenpair
          info%converged(j) = max(1, info%iteration)
        else ! stagnated eigenpair
          info%converged(j) = -max(1, info%iteration)
        end if
        lambda(j) = si_map(problem, sigma, keep%lmd(rci%jx + k))
      end do
      if ( rci%i > 0 ) then
        j = rci%jx + new
        m = block_size
        ! find, if possible, where is the next left eigenvalue
        if ( .not. keep%left_converged ) then
          info%next_left = sigma ! assume not available yet
          if ( j <= m ) then
            if ( keep%lcon + new > 0 .and. keep%lmd(j)*keep%lmd(1) > 0 ) &
              info%next_left = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%first = j ! update first non-converged eigenpair index
      else ! same for the right converged
        j = rci%jx - new
        m = block_size
        if ( .not. keep%right_converged ) then
          info%next_right = sigma
          if ( j > 0 ) then
            if ( keep%rcon + new > 0 .and. keep%lmd(j)*keep%lmd(m) > 0 ) &
              info%next_right = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%last = j
      end if

      ! report the convergence situation
      if ( options%print_level == 2 .and. u_diag > NONE ) then

        ! print data up to the first non-converged eigenpair
        ! on each side of sigma only

        frmt = '(i7, i8, es23.14, a, es10.1, 2es11.1)'

        ! the first and last converged eigenpair
        first = rci%jx
        last = rci%jx + (rci%nx - 1)*rci%i
        ! print one non-convergent if needed and possible
        if ( rci%i > 0 .and. (.not. keep%left_converged) &
          .and. keep%lcon + new < keep%max_left .and. last < block_size ) then
          if ( keep%lmd(last + 1 ) < ZERO ) last = last + 1
        end if
        if ( rci%i < 0 .and. (.not. keep%right_converged) &
          .and. keep%rcon + new < keep%max_right .and. last > 1 ) then
          if ( keep%lmd(last - 1) > ZERO ) last = last - 1
        end if

        step = rci%i

        do i = first, last, step

          if ( keep%info%converged(i) .ne. 0 ) then
            word = '  yes'
          else
            word = '   no'
          end if
          s = keep%lmd(i)
          ! t gets the eigenvalue error estimate when available,
          ! otherwise a rough guess, see the core solver keep%err_lambda meaning
          t = keep%info%err_lambda(block_size + i)
          ! recomputing the estimate for wanted eigenvalues
          ! from that for shifted-inverted
          if ( t >= ZERO ) t = si_map_err(problem, sigma, s, t)
          ! same for eigenvalues
          s = si_map(problem, sigma, keep%lmd(i))

          ! use negative indices for eigenvalues left of sigma
          ! and positive for those on the right
          if ( rci%i > 0 ) then
            ! eigenvalues left of sigma are marked by negative indices:
            ! ... <= lambda_{-2} <= \lambda_{-1} < sigma
            j = -(keep%lcon + i - first + 1)
          else
            ! right of sigma positive indices are used
            ! sigma < lambda_1 <= lambda_2 <= ...
            j = (keep%rcon + first - i + 1)
          end if
          write( u_diag, frmt ) &
            info%iteration, j, s, word, &
            keep%info%residual_norms(i), t, keep%info%err_x(block_size + i)

          ! quit printing after the first non-convergent eigenpair
          if ( keep%info%converged(i) == 0 ) exit

        end do

      end if

      ! update the numbers of converged eigenpairs
      if ( rci%i > 0 ) then
        keep%lcon = keep%lcon + new
      else
        keep%rcon = keep%rcon + new
      end if
      ncon = keep%lcon + keep%rcon

      if ( options%print_level > 2 .and. rci%i < 0 .and. u_diag > NONE ) then

        ! print data for all recently computed eigenpairs

        write( u_diag, '(/a, i5, a, i4, i4)' ) &
          'Iteration: ', info%iteration, ', not converged: ', &
          max(0, left - keep%lcon), max(0, right - keep%rcon)

        head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
        neck = &
 '                       |           |           |   errors   |    errors'
        line = &
 '-------------------------------------------------------------------------'

        frmt = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'

        write( u_diag, '(a/a)' ) &
          trim(line), trim(head)
        write( u_diag, '(a)' ) trim(neck)
        write( u_diag, '(a)' ) trim(line)

        do i = 1, block_size

          if ( keep%info%converged(i) /= 0 ) then
            word = '   yes'
          else
            word = '    no'
          end if

          s = keep%lmd(i)
          t = keep%info%err_lambda(block_size + i)
          if ( t >= ZERO ) t = si_map_err(problem, sigma, s, t)
          s = si_map(problem, sigma, keep%lmd(i))

          write( u_diag, frmt ) &
            s, ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            t, '  |', keep%info%err_X(block_size + i)

        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line)

      end if ! print_level > 2

      ! find whether all wanted eigenpairs have converged

      select case ( rci%i )

      case ( 1 )

        if ( left_gap == ZERO ) then

          ! no gap in the spectrum is required, just look at the
          ! number of converged eigenpairs
          keep%left_converged = keep%lcon >= left
          if ( keep%left_converged ) then
            info%left = keep%lcon
            info%next_left = sigma
            if ( keep%first <= block_size ) then
              if ( keep%lmd(keep%first) < ZERO ) &
                info%next_left = si_map(problem, sigma, keep%lmd(keep%first))
            end if
          end if

        else

          ! check if there is a sufficient gap between the last
          ! converged eigenvalue and the rest of the spectrum

          if ( left_gap > 0 ) then ! user set the absolute value of the gap
            t = left_gap
            q = ZERO
          else ! user set the relative value of the gap
            t = abs(left_gap) * keep%av_dist
            q = MIN_REL_GAP
          end if

          keep%left_converged = keep%left == 0 .or. keep%lcon >= keep%max_left
          ! look for a sufficient gap between the converged eigenvalues
          do i = keep%lcon, left + 1, -1
            j = minloc(lambda(1 : i - 1), 1)
            s = lambda(j) - lambda(i)
            if ( s >= t .and. s > q*(sigma - lambda(i)) ) then
              ! sufficient gap found
              info%next_left = lambda(i)
              info%left = i - 1
              keep%left_converged = .true.
              exit
            end if
          end do
          if ( keep%lcon >= left .and. &
              .not. keep%left_converged .and. keep%first <= block_size ) then
            ! if a sufficient gap is not found between
            ! the converged eigenvalues, look at the gap
            ! between the last converged and the first non-converged
            ! left of sigma
            if ( si_map(problem, sigma, keep%lmd(keep%first)) < sigma ) then
              ! make sure the error in the first non-converged eigenvalue
              ! is small enough for the gap to be trustworthy
              r = si_map_err &
                (problem, sigma, keep%lmd(keep%first), &
                  keep%info%err_lambda(keep%first))
              s = lambda(keep%lcon) &
                - si_map(problem, sigma, keep%lmd(keep%first))
              if ( r > ZERO .and. r < s*sqrt(epsilon(ONE)) .and. s >= t ) then
                info%left = keep%lcon
                i = keep%first
                info%next_left = si_map(problem, sigma, keep%lmd(keep%first))
                keep%left_converged = .true.
              end if
            end if
          end if

        end if

      case ( -1 ) ! same on the right of sigma

        if ( right_gap == ZERO ) then

          keep%right_converged = keep%rcon >= right
          if ( keep%right_converged ) then
            info%right = keep%rcon
            info%next_right = sigma
            if ( keep%last > 0 ) then
              if ( keep%lmd(keep%last) > ZERO ) &
                info%next_right = si_map(problem, sigma, keep%lmd(keep%last))
            end if
          end if


        else

          if ( right_gap > 0 ) then
            t = right_gap
            q = ZERO
          else
            t = abs(right_gap) * keep%av_dist
            q = MIN_REL_GAP
          end if

          keep%right_converged = &
            keep%right == 0 .or. keep%rcon >= keep%max_right
          do i = mep - keep%rcon + 1, mep - right
            j = i + maxloc(lambda(i + 1 : mep), 1)
            s = lambda(i) - lambda(j)
            if ( s >= t .and. s > q*(lambda(i) - sigma) ) then
              info%right = mep - i
              info%next_right = lambda(i)
              keep%right_converged = .true.
              exit
            end if
          end do

          if ( keep%rcon >= right .and. &
              .not. keep%right_converged .and. keep%last >= 1 ) then
            if ( si_map(problem, sigma, keep%lmd(keep%last)) > sigma ) then
              r = si_map_err &
                (problem, sigma, keep%lmd(keep%last), &
                  keep%info%err_lambda(keep%last))
              s = si_map(problem, sigma, keep%lmd(keep%last)) &
                - lambda(mep - keep%rcon + 1)
              if ( r > ZERO .and. r < s*sqrt(epsilon(ONE)) .and. s >= t ) then
                info%right = keep%rcon
                info%next_right = si_map(problem, sigma, keep%lmd(keep%last))
                keep%right_converged = .true.
              end if
            end if
          end if

        end if

      end select ! case ( rci%i )

    case (SSMFE_CHECK_CONVERGENCE)

      mep = max_nep
      m = block_size

      ! estimate the average distance between eigenvalues

      ! maximal eigenvector error for a 'trusted' eigenvalue
      delta = max(options%tol_x, sqrt(options%rel_tol_lambda))
      delta = sqrt(max(delta, epsilon(ONE)))

      r = ZERO
      k = 0

      ! initial 'trusted' eigenvalues count
      l = keep%lcon

      ! t is the last computed eigenvalue on the left or its approximation
      if ( keep%lcon > 0 ) then
        t = lambda(keep%lcon)
      else
        t = si_map(problem, sigma, keep%lmd(1))
      end if

      ! establish the search range for the last 'trusted' eigenvalue
      ! left of sigma
      first = keep%first
      if ( keep%left /= 0 ) then
        last = keep%last
      else
        last = first - 1 ! no eigenvalues available
      end if

      do step = -1, 1, 2
        s = t ! the current 'trusted' eigenvalue
        do i = first, last, -step
          ! stop if we moved to the other side of sigma
          if ( step*(si_map(problem, sigma, keep%lmd(i)) - sigma) < ZERO ) exit
          ! stop if the eigenvector error is too large
          if ( keep%info%err_x(m + i) > delta .and. i /= first ) exit
          ! increment the number of 'trusted' eigenvalues
          l = l + 1
          ! update the current 'trusted' eigenvalue
          s = si_map(problem, sigma, keep%lmd(i))
        end do
        if ( step > 0 ) exit ! both sides done
        ! save the number of 'trusted' eigenvalues left of sigma
        k = l
        ! save the last 'trusted' eigenvalue left of sigma
        r = s
        ! move to the other side of sigma and do the same
        l = keep%rcon
        if ( keep%rcon > 0 ) then
          t = lambda(mep - keep%rcon + 1)
        else
          t = si_map(problem, sigma, keep%lmd(m))
        end if
        first = keep%last
        if ( keep%right > 0 ) then
          last = keep%first
        else
          last = first + 1
        end if
      end do

      if ( k == 0 .and. l == 0 ) then
        keep%av_dist = ZERO ! assume not available yet
      else if ( k == 0 ) then
        ! no 'trusted' eigenvalues left of sigma, use those to the right
        ! t is the 'trusted' eigenvalue closest to sigma
        if ( keep%rcon > 0 ) then
          t = lambda(mep)
        else ! or its approximation
          t = si_map(problem, sigma, keep%lmd(m))
        end if
        ! 'safe' average (l can be 1)
        keep%av_dist = (s - t)/l
      else if ( l == 0 ) then
        ! no 'trusted' eigenvalues left of sigma, use those to the left
        if ( keep%lcon > 0 ) then
          t = lambda(1)
        else
          t = si_map(problem, sigma, keep%lmd(1))
        end if
        keep%av_dist = (t - r)/k
      else
        ! the estimated average distance is the 'trusted' eigenvalues range
        ! divided by their number
        keep%av_dist = (s - r)/(k + l)
      end if

      ! are residual norms to be checked?
      check_res = options%abs_tol_residual > 0 &
             .or. options%rel_tol_residual > 0
      ! are eigenvalue errors to be checked?
      check_lmd = options%abs_tol_lambda > 0 &
             .or. options%rel_tol_lambda > 0

      ! round-off errors level
      t = 10*epsilon(ONE)

      do i = 1, m

        ! skip previously converged
        if ( keep%info%converged(i) /= 0 ) cycle

        converged = check_res .or. check_lmd .or. options%tol_x /= 0

        if ( check_res ) then
          s = abs(keep%lmd(i)) * max(options%rel_tol_residual, t)
          s = max(s, options%abs_tol_residual)
          converged = keep%info%residual_norms(i) <= s
        end if

        if ( check_lmd ) then
          if ( keep%info%err_lambda(i) >= ZERO ) then
            s = keep%av_dist * max(options%rel_tol_lambda, t*t)
            s = max(s, options%abs_tol_lambda)
            r = si_map_err(problem, sigma, keep%lmd(i), keep%info%err_lambda(i))
            converged = converged .and. r <= s
          else
            converged = .false.
          end if
        end if

        if ( options%tol_x > ZERO ) then
          ! check the eigenvector error against user-defined tolerance
          converged = converged &
          .and. keep%info%err_x(i) >= ZERO & ! error estimate available
          .and. keep%info%err_x(i) <= max(options%tol_x, t)
        else if ( options%tol_x < ZERO ) then
          ! use default tolerance
          converged = converged &
          .and. keep%info%err_x(i) >= ZERO &
          .and. keep%info%err_x(i) <= sqrt(epsilon(ONE))
        end if

        ! for extra reliability make sure the eigenvector error is small
        converged = converged .and. abs(keep%info%err_x(m + i)) < 0.05

        if ( converged ) keep%info%converged(i) = 1 ! mark as converged

      end do ! i = 1, m

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1

      go to 1

    end select

    if ( rci%job < 0 ) then
      info%non_converged = 0 ! assume success
      if ( rci%job /= SSMFE_DONE ) then
        ! if the return is not successful, report the numbers of
        ! converged and non-converged eigenpairs
        info%non_converged = &
          max(0, left - keep%lcon) + max(0, right - keep%rcon)
        info%left = keep%lcon
        info%right = keep%rcon
      end if
      info%flag = keep%info%flag
      if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )
    end if

  contains

    ! reverse shift-invert mapping: for a given eigenvalue mu
    ! of (A - sigma B)^{-1} B (problem = 1)
    ! or (B - sigma A)^{-1} B (problem = 2)
    ! returns the corresponding eigenvalue lambda of
    ! A x = lambda B x
    real(PRECISION) function si_map( problem, sigma, mu )
      integer, intent(in) :: problem
      real(PRECISION), intent(in) :: sigma, mu
      if ( problem == 2 ) then
        si_map = sigma*mu/(mu - 1)
      else
        si_map = sigma + 1/mu
      end if
    end function si_map

    ! same for the eigenvalue error
    real(PRECISION) function si_map_err( problem, sigma, mu, delta )
      integer, intent(in) :: problem
      real(PRECISION), intent(in) :: sigma, mu, delta
      if ( problem == 2 ) then
        si_map_err = sigma*delta/(mu - 1)**2
      else
        si_map_err = delta/mu**2
      end if
    end function si_map_err

  end subroutine ssmfe_inverse_rci_double

!*************************************************************************
! real RCI for leftmost eigenvalues
!
  subroutine ssmfe_direct_rci_double &
    ( problem, nep, max_nep, lambda, block_size, rr_matrices, ind, &
      rci, keep, options, info )
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0

    real(kind = PRECISION), parameter :: MIN_REL_GAP = 1E-3

    ! problem type
    ! 0: A x = lambda x
    ! 1: A x = lambda B x
    ! A = A', B = B' > 0
    integer, intent(in) :: problem

    ! number of wanted eigenpairs
    integer, intent(in) :: nep

    ! eigenpairs storage size
    integer, intent(in) :: max_nep

    ! eigenvalue storage
    real(PRECISION), intent(inout) :: lambda(max_nep)

    ! BJCG block size
    integer, intent(in) :: block_size

    ! matrices for the Rayleigh-Ritz procedure
    real(PRECISION), intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)

    ! work matrix column reordering index
    integer, intent(inout) :: ind(block_size)

    ! ssmfe_expert types - see the spec
    type(ssmfe_rcid       ), intent(inout) :: rci
    type(ssmfe_options    ), intent(in) :: options
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_inform     ), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, frmt, line
    character(7) :: word

    integer :: u_diag
    integer :: first, last
    integer :: i, j, k, l

    real(PRECISION) :: gap
    real(PRECISION) :: delta
    real(PRECISION) :: q, r, s, t

    ! short names for options

    u_diag = options%unit_diagnostic
    gap    = options%left_gap

    info%flag = 0

    ! check for the input data errors

    if ( nep < 0 &
        .or. options%max_left >= 0 .and. options%max_left < nep ) then
      info%flag = WRONG_LEFT
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( nep > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( block_size < 1 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    ! things to do on (re)start
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        ! reallocate the information arrays
        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda    ) ) deallocate ( info%err_lambda     )
        if ( allocated(info%err_X         ) ) deallocate ( info%err_X          )
        if ( allocated(info%converged     ) ) deallocate ( info%converged      )
        allocate &
          ( info%residual_norms(max_nep), info%err_lambda(max_nep), &
            info%err_X(max_nep), info%converged(max_nep), stat = info%stat )
        if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
        if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
        if ( info%stat /= 0 ) return

        ! initialize the execution information
        info%iteration = 0
        info%residual_norms = ZERO
        info%err_lambda = -ONE
        info%err_X = -ONE
        info%converged = 0
        info%left = 0
        info%next_left = huge(ONE)

        ! nothing to do
        if ( nep <= 0 ) then
          rci%job = SSMFE_DONE
          return
        end if

        ! save the input data
        keep%options%err_est = options%err_est
        keep%left = nep
        keep%right = 0
        keep%options%minAprod = options%minAprod
        keep%options%minBprod = options%minBprod

        ! compute the number of extra vectors iterated for the sake
        ! of better convergence
        if ( options%extra_left < 0 ) then
          keep%options%extra_left = max(10, nep/10) ! use default
        else
          keep%options%extra_left = options%extra_left ! take user's value
        end if

      end if

      ! if the total number of eigenvalues to be computed (wanted + extra)
      ! exceeds the block size m, they are computed in portions
      ! of sizes keep%left (those left of sigma);
      ! once a portion is computed, iterations are restarted (see below)
      ! (note that more than wanted + extra eigenpairs are needed if gaps in
      ! the spectrum separating wanted eigenvalues from the rest are sought)

      i = nep + keep%options%extra_left
      if ( block_size < i .or. gap /= ZERO ) then
        keep%options%extra_left = min(keep%options%extra_left, block_size/2)
        keep%left = block_size - keep%options%extra_left
      else
        keep%left = nep
      end if
      ! if gaps in the spectrum are sought, the number of computed
      ! eigenpairs must be larger than requested
      if ( keep%left == nep .and. gap /= ZERO &
          .and. keep%options%extra_left > 0 ) then
        keep%left = keep%left + 1
        keep%options%extra_left = keep%options%extra_left - 1
      end if

      if ( rci%job == SSMFE_START ) then

        keep%options%min_gap = 0.05 ! see ssmfe_core spec

        keep%lcon = 0 ! number of converged eigenpairs

        keep%left_converged = .false.

        if ( options%max_left < 0 ) then ! the number of eigenvalues not known
          keep%max_left = max_nep + 1     ! set to no-impact value
        else
          keep%max_left = options%max_left
        end if

      end if

      ! allocate eigenvalues array of ssmfe_core solver
      if ( allocated(keep%lmd) ) deallocate ( keep%lmd )
      allocate ( keep%lmd(block_size), stat = info%stat )
      if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
      if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
      if ( info%stat /= 0 ) return

      keep%first = 1          ! first non-converged eigenpair index
      keep%last  = block_size ! last non-converged eigenpair index

      ! print problem summary
      call ssmfe_msg( problem, options, nep, 0, block_size )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) rci%job = SSMFE_START

    else if ( rci%job == SSMFE_SAVE_CONVERGED .and. rci%i < 0 ) then

      ! some converged eigenpairs just saved to storage,
      ! decide on whether to go for more

      if ( .not. keep%left_converged ) then

        if ( info%iteration >= options%max_iterations ) then
          ! check the number of iterations performed
          rci%job = SSMFE_QUIT
          info%flag = MAX_NUM_ITERATIONS_EXCEEDED
        else if ( keep%lcon == max_nep ) then
          ! check the remaining storage
          rci%job = SSMFE_QUIT
          info%flag = OUT_OF_STORAGE
        end if
      end if

      if ( rci%job /= SSMFE_QUIT ) then
        if ( keep%left_converged ) then
          rci%job = SSMFE_DONE
        else
          if ( keep%first > keep%left ) then
            ! the current portion of eigenpairs converged,
            ! but not all wanted eigenpairs have been computed,
            ! restart to compute the next portion
            rci%job = SSMFE_RESTART
            rci%jx = keep%first
            rci%nx = keep%last - keep%first + 1
            rci%k = 0
            return
          end if
        end if
      end if

      if ( rci%job < 0 ) then ! computation terminated
        info%non_converged = max(0, nep - keep%lcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
        end if
        if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )
        return
      end if

    end if

    ! call the core solver
1   call ssmfe &
      ( rci, problem, keep%left, keep%right, &
        block_size, keep%lmd, rr_matrices, ind, &
        keep%keep, keep%options, keep%info )

    select case (rci%job)

    case (SSMFE_SAVE_CONVERGED)

      if ( rci%i < 0 ) return

      ! if the number of converged eigenpairs exceeds the storage size
      ! reduce the number of newly converged eigenpairs to be saved
      if ( keep%lcon + rci%nx > max_nep ) then
        rci%nx = max_nep - keep%lcon
      end if

      do k = 0, rci%nx - 1
        j = keep%lcon + k + 1
        info%residual_norms(j) = keep%info%residual_norms(rci%jx + k)
        info%err_lambda(j) = keep%info%err_lambda(rci%jx + k)
        info%err_X(j) = keep%info%err_X(rci%jx + k)
        l = keep%info%converged(rci%jx + k)
        if ( l > 0 ) then
          info%converged(j) = max(1, info%iteration)
        else
          info%converged(j) = -max(1, info%iteration)
        end if
        lambda(j) = keep%lmd(rci%jx + k)
      end do
      j = rci%jx + rci%nx
      info%next_left = huge(ONE)
      if ( j <= block_size ) info%next_left = keep%lmd(j)
      keep%first = j

      if ( options%print_level == 2 .and. u_diag > NONE ) then

        frmt = '(i7, i8, es23.14, a, es10.1, 2es11.1)'

        first = rci%jx
        last = rci%jx + rci%nx - 1
        if ( .not. keep%left_converged .and. last < block_size &
              .and. keep%lcon + rci%nx < keep%max_left ) last = last + 1

        do i = first, last
          if ( keep%info%converged(i) /= 0 ) then
            word = '  yes'
          else
            word = '   no'
          end if
          write( u_diag, frmt ) &
            info%iteration, keep%lcon + i - first + 1, &
            keep%lmd(i), word, &
            keep%info%residual_norms(i), &
            keep%info%err_lambda(block_size + i), &
            keep%info%err_X(block_size + i)
          if ( keep%info%converged(i) == 0 ) exit
        end do

      end if

      keep%lcon = keep%lcon + rci%nx

      if ( options%print_level > 2 .and. u_diag > NONE ) then

        write( u_diag, '(/a, i5, a, i4)' ) &
          'Iteration: ', info%iteration, ', not converged: ', &
          max(0, nep - keep%lcon)

        head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
        neck = &
 '                       |           |           |   errors   |    errors'
        line = &
 '-------------------------------------------------------------------------'

        frmt = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'

        write( u_diag, '(a/a)' ) &
          trim(line), trim(head)
        write( u_diag, '(a)' ) trim(neck)
        write( u_diag, '(a)' ) trim(line)

        do i = 1, block_size
          if ( keep%info%converged(i) > 0 ) then
            word = '   yes'
          else
            word = '    no'
          end if
          write( u_diag, frmt ) &
            keep%lmd(i), ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            keep%info%err_lambda(block_size + i), '  |', &
            keep%info%err_X(block_size + i)
        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line)

      end if ! print_level > 2

      if ( gap == ZERO ) then

        keep%left_converged = keep%lcon >= nep
        if ( keep%left_converged ) then
          info%left = keep%lcon
          info%next_left = huge(ONE)
          if ( keep%first <= block_size ) then
            info%next_left = keep%lmd(keep%first)
          end if
        end if

      else

        if ( gap > 0 ) then
          t = gap
          q = ZERO
        else
          t = abs(gap) * keep%av_dist
          q = MIN_REL_GAP
        end if

        keep%left_converged = keep%left == 0 .or. keep%lcon >= keep%max_left
        do i = nep + 1, keep%lcon
          j = maxloc(lambda(1 : i - 1), 1)
          r = lambda(i) - lambda(j)
          if ( r >= t .and. r > q * abs(lambda(i)) ) then
            info%left = i - 1
            info%next_left = lambda(i)
            keep%left_converged = .true.
            exit
          end if
        end do
        if ( .not. keep%left_converged .and. keep%first <= block_size ) then
          if ( keep%lcon >= nep ) then
            s =  keep%lmd(keep%first) - lambda(keep%lcon)
            r = keep%info%err_lambda(keep%first)
            if ( r > ZERO .and. r < s*sqrt(epsilon(ONE)) .and. s >= t ) then
              info%left = keep%lcon
              info%next_left = keep%lmd(keep%first)
              keep%left_converged = .true.
            end if
          end if
        end if

      end if

    case (SSMFE_CHECK_CONVERGENCE)

      delta = max(options%tol_X, sqrt(options%rel_tol_lambda))
      delta = sqrt(max(delta, epsilon(ONE)))

      l = keep%lcon
      if ( keep%lcon > 0 ) then
        t = lambda(keep%lcon)
      else
        t = keep%lmd(1)
      end if
      first = keep%first
      last = keep%last
      s = t
      do i = first, last
        if ( keep%info%err_X(block_size + i) > delta .and. i /= first ) exit
        l = l + 1
        s = keep%lmd(i)
      end do

      if ( keep%lcon > 0 ) then
        t = lambda(1)
      else
        t = keep%lmd(1)
      end if
      keep%av_dist = (s - t)/l

      r = ZERO
      if ( keep%lcon > 0 ) r = abs(lambda(1))
      r = r*epsilon(ONE)

      check_res = options%abs_tol_residual > 0 &
             .or. options%rel_tol_residual > 0
      check_lmd = options%abs_tol_lambda > 0 &
             .or. options%rel_tol_lambda > 0

      t = 10*epsilon(ONE)

      do i = 1, block_size

        if ( keep%info%converged(i) /= 0 ) cycle

        converged = check_res .or. check_lmd .or. options%tol_X /= 0

        if ( check_res ) then
          s = abs(keep%lmd(i)) * max(options%rel_tol_residual, t)
          s = max(s, options%abs_tol_residual)
          converged = keep%info%residual_norms(i) <= s
        end if

        if ( check_lmd ) then
          if ( keep%info%err_lambda(i) >= ZERO ) then
            s = keep%av_dist * max(options%rel_tol_lambda, t*t)
            s = max(s, options%abs_tol_lambda)
            converged = converged .and. keep%info%err_lambda(i) <= s
          else
            converged = .false.
          end if
        end if

        if ( options%tol_X > ZERO ) then
          converged = converged &
          .and. keep%info%err_X(i) >= ZERO &
          .and. keep%info%err_X(i) <= max(options%tol_X, t)
        else if ( options%tol_X < ZERO ) then
          converged = converged &
          .and. keep%info%err_X(i) >= ZERO &
          .and. keep%info%err_X(i) <= sqrt(epsilon(ONE))
        end if

        converged = converged .and. abs(keep%info%err_X(block_size + i)) < 0.05

        if ( converged ) keep%info%converged(i) = 1

      end do ! i = 1, block_size

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1

      go to 1

    end select

    if ( rci%job < 0 ) then

      info%non_converged = 0
      if ( rci%job /= SSMFE_DONE ) then
        info%non_converged = max(0, nep - keep%lcon)
        info%left = keep%lcon
      end if
      info%flag = keep%info%flag

      if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )

    end if

  end subroutine ssmfe_direct_rci_double

  subroutine ssmfe_expert_free_keep_double( keep )
    type(ssmfe_expert_keep), intent(inout) :: keep

    if ( allocated(keep%lmd) ) deallocate( keep%lmd )
    call ssmfe_free( keep%keep, keep%info )

  end subroutine ssmfe_expert_free_keep_double

  subroutine ssmfe_expert_free_double( keep, info )
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: info

    call ssmfe_free( keep )
    call ssmfe_free( info )

  end subroutine ssmfe_expert_free_double

  subroutine ssmfe_inverse_rci_double_complex &
    ( problem, sigma, left, right, &
      max_nep, lambda, block_size, rr_matrices, ind, &
      rci, keep, options, info )
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_APPLY_A           = 1
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE  = 9
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0

    real(kind = PRECISION), parameter :: MIN_REL_GAP = 1E-3

    integer, intent(in) :: problem
    real(PRECISION), intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: max_nep
    real(PRECISION), intent(inout) :: lambda(max_nep)
    integer, intent(in) :: block_size
    complex(PRECISION), intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)
    integer, intent(inout) :: ind(block_size)
    type(ssmfe_rciz       ), intent(inout) :: rci
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_options    ), intent(in   ) :: options
    type(ssmfe_inform     ), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, frmt, line
    character(7) :: word

    integer :: u_diag
    integer :: m, mep, new, ncon, first, last, step
    integer :: i, j, k, l

    real(PRECISION) :: left_gap
    real(PRECISION) :: right_gap
    real(PRECISION) :: delta
    real(PRECISION) :: q, r, s, t

    u_diag = options%unit_diagnostic

    left_gap = options%left_gap
    right_gap = options%right_gap

    info%flag = 0

    if ( left < 0 &
        .or. options%max_left >= 0 .and. options%max_left < left ) then
      info%flag = WRONG_LEFT
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( right < 0 &
        .or. options%max_right >= 0 .and. options%max_right < right ) then
      info%flag = WRONG_RIGHT
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( left + right > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( block_size < 2 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( problem == 2 .and. sigma == ZERO ) then
      info%flag = WRONG_SIGMA
      rci%job = SSMFE_QUIT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        mep = max_nep
        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda    ) ) deallocate ( info%err_lambda     )
        if ( allocated(info%err_X         ) ) deallocate ( info%err_X          )
        if ( allocated(info%converged     ) ) deallocate ( info%converged      )
        allocate &
          ( info%residual_norms(mep), info%err_lambda(mep), &
            info%err_X(mep), info%converged(mep), stat = info%stat )
        if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
        if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
        if ( info%stat /= 0 ) return

        info%iteration = 0
        info%residual_norms = ZERO
        info%err_lambda = -ONE
        info%err_X = -ONE
        info%converged = 0
        info%left = 0
        info%right = 0
        info%next_left = sigma
        info%next_right = sigma
        if ( left <= 0 .and. right <= 0 ) then
          rci%job = SSMFE_DONE
          return
        end if

        keep%problem = -abs(problem)
        keep%options%err_est = options%err_est
        keep%left = left
        keep%right = right

        ! compute the number of extra vectors iterated for the sake
        ! of better convergence
        if ( left < 1 ) then
          keep%options%extra_left = 0
        else if ( options%extra_left < 0 ) then
          keep%options%extra_left = max(10, left/10) ! use default
        else
          keep%options%extra_left = options%extra_left ! take user's value
        end if
        if ( right < 1 ) then
          keep%options%extra_right = 0
        else if ( options%extra_right < 0 ) then
          keep%options%extra_right = max(10, right/10) ! use default
        else
          keep%options%extra_right = options%extra_right ! take user's value
        end if

        if ( options%max_left >= 0 ) keep%options%extra_left = &
          min(keep%options%extra_left, options%max_left - left)
        if ( options%max_right >= 0 ) keep%options%extra_right = &
          min(keep%options%extra_right, options%max_right - right)

        keep%options%minAprod = options%minAprod
        keep%options%minBprod = options%minBprod

      end if

      m = block_size
      ! if the total number of eigenvalues to be computed (wanted + extra)
      ! exceeds the block size m, they are computed in portions
      ! of sizes keep%left (those left of sigma) and keep%right (right of it);
      ! once a portion is computed, iterations are restarted (see below)
      ! (note that more than wanted + extra eigenpairs are needed if gaps in
      ! the spectrum separating wanted eigenvalues from the rest are sought)

      if ( keep%left > 0 ) then
        ! the total numers of eigenvalues left of sigma to be computed
        i = left + keep%options%extra_left
        if ( keep%right > 0 ) then
          ! the total numers of eigenvalues right of sigma to be computed
          j = right + keep%options%extra_right
          ! workspace for eigenvectors corresponding to left eigenvalues
          l = m*i/(i + j)
          keep%left = left ! assumingly can be computed in one go
          if ( l < i .or. left_gap /= ZERO ) then ! not enough space
            keep%options%extra_left = min(keep%options%extra_left, l/2)
            keep%left = l - keep%options%extra_left ! portion size
            ! limiting extras number to l/2 ensures that the portion size
            ! is substantial thus reducing the number of restarts
          end if
          keep%right = right ! assumingly can be computed in one go
          if ( m - l < j .or. right_gap /= ZERO ) then ! not enough space
            keep%options%extra_right = min(keep%options%extra_right, (m - l)/2)
            keep%right = m - l - keep%options%extra_right ! portion size
          end if
        else
          keep%left = left
          if ( m < i .or. left_gap /= ZERO ) then
            keep%options%extra_left = min(keep%options%extra_left, m/2)
            keep%left = m - keep%options%extra_left
          end if
          keep%options%extra_right = 0
          keep%right = 0
        end if
      else
        keep%options%extra_left = 0
        if ( m < right + keep%options%extra_right .or. right_gap /= ZERO ) then
          keep%options%extra_right = min(keep%options%extra_right, m/2)
          keep%right = m - keep%options%extra_right
        end if
      end if
      ! if gaps in the spectrum are sought, the number of computed
      ! eigenpairs must be larger than requested
      if ( left > 0 .and. keep%left == left .and. left_gap /= ZERO &
          .and. keep%options%extra_left > 0 ) then
        keep%left = keep%left + 1
        keep%options%extra_left = keep%options%extra_left - 1
      end if
      if ( right > 0 .and. keep%right == right .and. right_gap /= ZERO &
          .and. keep%options%extra_right > 0 ) then
        keep%right = keep%right + 1
        keep%options%extra_right = keep%options%extra_right - 1
      end if

      if ( rci%job == SSMFE_START ) then

        keep%options%min_gap = 0.05

        keep%lcon = 0
        keep%rcon = 0

        ! initialize the convergence flags
        keep%left_converged = left == 0
        keep%right_converged = right == 0

        if ( options%max_left < 0 ) then
          keep%max_left = max_nep + 1
        else
          keep%max_left = options%max_left
        end if
        if ( options%max_right < 0 ) then
          keep%max_right = max_nep + 1
        else
          keep%max_right = options%max_right
        end if

      else

        if ( allocated(keep%lmd) ) deallocate ( keep%lmd )

      end if

      ! allocate work arrays
      ! keep%lmd: eigenvalues array of ssmfe solver
      ! keep%v: workspace for ssmfe solver matrices
      m = block_size
      allocate ( keep%lmd(m), stat = info%stat )
      if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
      if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
      if ( info%stat /= 0 ) return

      keep%first = 1
      keep%last = m

      call ssmfe_msg( problem, options, left, right, m )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) &
        rci%job = SSMFE_START

    else if ( rci%job == SSMFE_SAVE_CONVERGED .and. rci%i < 0 ) then

      mep = max_nep
      m = block_size
      ncon = keep%lcon + keep%rcon

      if ( .not. (keep%left_converged .and. keep%right_converged) ) then
        if ( info%iteration >= options%max_iterations ) then
          rci%job = SSMFE_QUIT
          info%flag = MAX_NUM_ITERATIONS_EXCEEDED
        else if ( ncon == mep &
          .or. .not. keep%left_converged &
               .and. keep%lcon >= mep - max(right, keep%rcon) &
          .or. .not. keep%right_converged &
               .and. keep%rcon >= mep - max(left, keep%lcon) &
            ) then
          rci%job = SSMFE_QUIT
          info%flag = OUT_OF_STORAGE
        end if
      end if

      if ( rci%job /= SSMFE_QUIT ) then

        if ( keep%left_converged ) keep%left = 0
        if ( keep%right_converged ) keep%right = 0

        if ( keep%left_converged .and. keep%right_converged ) then
          rci%job = SSMFE_DONE
        else
          if ( (keep%left > 0 .and. keep%first > keep%left .or. &
                keep%right > 0 .and. block_size - keep%last >= keep%right) &
              ) then
            rci%job = SSMFE_RESTART
            rci%jx = keep%first
            rci%nx = keep%last - keep%first + 1
            rci%k = 0
            return
          end if
        end if

      end if

      if ( rci%job < 0 ) then
        info%non_converged = &
          max(0, left - keep%lcon) + max(0, right - keep%rcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
          info%right = keep%rcon
        end if
        if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )
        return
      end if

    end if

    ! call the engine
1   call ssmfe &
      ( rci, keep%problem, keep%left, keep%right, block_size, &
        keep%lmd, rr_matrices, ind, keep%keep, keep%options, keep%info )

    select case (rci%job)

    case (SSMFE_APPLY_A)

      rci%job = SSMFE_DO_SHIFTED_SOLVE

    case (SSMFE_SAVE_CONVERGED)

      ! the number of newly converged eigenpairs
      new = rci%nx
      ! number of eigenpairs converged previously
      ncon = keep%lcon + keep%rcon
      mep = max_nep

      ! if the number of converged eigenpairs exceeds the storage size
      ! reduce the number of newly converged eigenpairs to be saved
      if ( rci%i > 0 ) then
        if ( left_gap == ZERO ) new = min(new, left - keep%lcon)
        k = max(right, keep%rcon) ! space reserved for the right eigenpairs
        if ( new > mep - k - keep%lcon ) new = mep - k - keep%lcon
      else
        if ( right_gap == ZERO ) new = min(new, right - keep%rcon)
        k = max(left, keep%lcon) ! space reserved for the left eigenpairs
        if ( new > mep - k - keep%rcon ) new = mep - k - keep%rcon
      end if
      rci%nx = new

      ! copy the newly converged eigenpairs and related data
      ! to lambda, X and info arrays
      do i = 1, new
        k = (i - 1)*rci%i
        if ( rci%i > 0 ) then
          j = keep%lcon + i
        else
          j = mep - keep%rcon - i + 1
        end if
        s = keep%lmd(rci%jx + k)
        t = keep%info%err_lambda(rci%jx + k)
        info%residual_norms(j) = keep%info%residual_norms(rci%jx + k)
        if ( t < 0 ) then
          info%err_lambda(j) = t
        else
          info%err_lambda(j) = si_map_err(problem, sigma, s, t)
        end if
        info%err_lambda(j) = keep%info%err_lambda(rci%jx + k)
        info%err_X(j) = keep%info%err_X(rci%jx + k)
        l = keep%info%converged(rci%jx + k)
        if ( l > 0 ) then
          info%converged(j) = max(1, info%iteration)
        else
          info%converged(j) = -max(1, info%iteration)
        end if
        lambda(j) = si_map(problem, sigma, keep%lmd(rci%jx + k))
      end do
      if ( rci%i > 0 ) then
        j = rci%jx + new
        m = block_size
        if ( .not. keep%left_converged ) then
          info%next_left = sigma
          if ( j <= m ) then
            if ( keep%lcon + new > 0 .and. keep%lmd(j)*keep%lmd(1) > 0 ) &
              info%next_left = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%first = j
      else
        j = rci%jx - new
        m = block_size
        if ( .not. keep%right_converged ) then
          info%next_right = sigma
          if ( j > 0 ) then
            if ( keep%rcon + new > 0 .and. keep%lmd(j)*keep%lmd(m) > 0 ) &
              info%next_right = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%last = j
      end if

      ! report the convergence situation
      if ( options%print_level == 2 .and. u_diag > NONE ) then

        frmt = '(i7, i8, es23.14, a, es10.1, 2es11.1)'

        ! the first and last converged eigenpair
        first = rci%jx
        last = rci%jx + (rci%nx - 1)*rci%i
        ! print one non-convergent if needed and possible
        if ( rci%i > 0 .and. (.not. keep%left_converged) &
          .and. keep%lcon + new < keep%max_left .and. last < block_size ) then
          if ( keep%lmd(last + 1 ) < ZERO ) last = last + 1
        end if
        if ( rci%i < 0 .and. (.not. keep%right_converged) &
          .and. keep%rcon + new < keep%max_right .and. last > 1 ) then
          if ( keep%lmd(last - 1) > ZERO ) last = last - 1
        end if

        step = rci%i

        do i = first, last, step

          if ( keep%info%converged(i) /= 0 ) then
            word = '  yes'
          else
            word = '   no'
          end if
          s = keep%lmd(i)
          ! t gets the eigenvalue error estimate when available,
          ! otherwise a rough guess, see the core solver keep%err_lambda meaning
          t = keep%info%err_lambda(block_size + i)
          ! recomputing the estimate for wanted eigenvalues
          ! from that for shifted-inverted
          if ( t >= ZERO ) t = si_map_err(problem, sigma, s, t)
          ! same for eigenvalues
          s = si_map(problem, sigma, keep%lmd(i))

          ! use negative indices for eigenvalues left of sigma
          ! and positive for those on the right
          if ( rci%i > 0 ) then
            j = -(keep%lcon + i - first + 1)
          else
            j = (keep%rcon + first - i + 1)
          end if
          write( u_diag, frmt ) &
            info%iteration, j, s, word, &
            keep%info%residual_norms(i), t, keep%info%err_X(block_size + i)

          ! quit printing after the first non-convergent eigenpair
          if ( keep%info%converged(i) == 0 ) exit

        end do

      end if

      ! update the numbers of converged eigenpairs
      if ( rci%i > 0 ) then
        keep%lcon = keep%lcon + new
      else
        keep%rcon = keep%rcon + new
      end if
      ncon = keep%lcon + keep%rcon

      if ( options%print_level > 2 .and. rci%i < 0 .and. u_diag > NONE ) then

        write( u_diag, '(/a, i5, a, i4, i4)' ) &
          'Iteration: ', info%iteration, ', not converged: ', &
          max(0, left - keep%lcon), max(0, right - keep%rcon)

        head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
        neck = &
 '                       |           |           |   errors   |    errors'
        line = &
 '-------------------------------------------------------------------------'

        frmt = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'

        write( u_diag, '(a/a)' ) &
          trim(line), trim(head)
        write( u_diag, '(a)' ) trim(neck)
        write( u_diag, '(a)' ) trim(line)

        do i = 1, block_size

          if ( keep%info%converged(i) /= 0 ) then
            word = '   yes'
          else
            word = '    no'
          end if

          s = keep%lmd(i)
          t = keep%info%err_lambda(block_size + i)
          if ( t >= ZERO ) t = si_map_err(problem, sigma, s, t)
          s = si_map(problem, sigma, keep%lmd(i))

          write( u_diag, frmt ) &
            s, ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            t, '  |', keep%info%err_X(block_size + i)

        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line)

      end if ! print_level > 2

      select case ( rci%i )

      case ( 1 )

        if ( left_gap == ZERO ) then

          keep%left_converged = keep%lcon >= left
          if ( keep%left_converged ) then
            info%left = keep%lcon
            info%next_left = sigma
            if ( keep%first <= block_size ) then
              if ( keep%lmd(keep%first) < ZERO ) &
                info%next_left = si_map(problem, sigma, keep%lmd(keep%first))
            end if
          end if

        else

          if ( left_gap > 0 ) then
            t = left_gap
            q = ZERO
          else
            t = abs(left_gap) * keep%av_dist
            q = MIN_REL_GAP
          end if

          keep%left_converged = keep%left == 0 .or. keep%lcon >= keep%max_left
          do i = keep%lcon, left + 1, -1
            j = minloc(lambda(1 : i - 1), 1)
            s = lambda(j) - lambda(i)
            if ( s >= t .and. s > q*(sigma - lambda(i)) ) then
              info%next_left = lambda(i)
              info%left = i - 1
              keep%left_converged = .true.
              exit
            end if
          end do
          if ( keep%first <= block_size ) then
            if ( .not. keep%left_converged .and. keep%lcon >= left &
              .and. si_map(problem, sigma, keep%lmd(keep%first)) < sigma ) then
              r = si_map_err &
                (problem, sigma, keep%lmd(keep%first), &
                  keep%info%err_lambda(keep%first))
              s = lambda(keep%lcon) &
                - si_map(problem, sigma, keep%lmd(keep%first))
              if ( r > ZERO .and. r < s*sqrt(epsilon(ONE)) .and. s >= t ) then
                info%left = keep%lcon
                i = keep%first
                info%next_left = si_map(problem, sigma, keep%lmd(keep%first))
                keep%left_converged = .true.
              end if
            end if
          end if

        end if

      case ( -1 )

        if ( right_gap == ZERO ) then

          keep%right_converged = keep%rcon >= right
          if ( keep%right_converged ) then
            info%right = keep%rcon
            info%next_right = sigma
            if ( keep%last > 0 ) then
              if ( keep%lmd(keep%last) > ZERO ) &
                info%next_right = si_map(problem, sigma, keep%lmd(keep%last))
            end if
          end if

        else

          if ( right_gap > 0 ) then
            t = right_gap
            q = ZERO
          else
            t = abs(right_gap) * keep%av_dist
            q = MIN_REL_GAP
          end if

          keep%right_converged = keep%right == 0 .or. keep%rcon >= keep%max_right
          do i = mep - keep%rcon + 1, mep - right
            j = i + maxloc(lambda(i + 1 : mep), 1)
            s = lambda(i) - lambda(j)
            if ( s >= t .and. s > q*(lambda(i) - sigma) ) then
              info%right = mep - i
              info%next_right = lambda(i)
              keep%right_converged = .true.
              exit
            end if
          end do
          if ( keep%last >= 1 ) then
            if ( .not. keep%right_converged .and. keep%rcon >= right &
              .and. si_map(problem, sigma, keep%lmd(keep%last)) > sigma ) then
              r = si_map_err &
                (problem, sigma, keep%lmd(keep%last), &
                  keep%info%err_lambda(keep%last))
              s = si_map(problem, sigma, keep%lmd(keep%last)) &
                - lambda(mep - keep%rcon + 1)
              if ( r > ZERO .and. r < s*sqrt(epsilon(ONE)) .and. s >= t ) then
                info%right = keep%rcon
                info%next_right = si_map(problem, sigma, keep%lmd(keep%last))
                keep%right_converged = .true.
              end if
            end if
          end if

        end if

      end select ! case ( rci%i )

    case (SSMFE_CHECK_CONVERGENCE)

      mep = max_nep
      m = block_size

      delta = max(options%tol_X, sqrt(options%rel_tol_lambda))
      delta = sqrt(max(delta, epsilon(ONE)))

      r = ZERO
      k = 0

      l = keep%lcon
      if ( keep%lcon > 0 ) then
        t = lambda(keep%lcon)
      else
        t = si_map(problem, sigma, keep%lmd(1))
      end if
      first = keep%first
      if ( keep%left /= 0 ) then
        last = keep%last
      else
        last = first - 1
      end if
      do step = -1, 1, 2
        s = t
        do i = first, last, -step
          if ( step*(si_map(problem, sigma, keep%lmd(i)) - sigma) < ZERO ) exit
          if ( keep%info%err_X(m + i) > delta .and. i /= first ) exit
          l = l + 1
          s = si_map(problem, sigma, keep%lmd(i))
        end do
        if ( step > 0 ) exit
        k = l
        r = s
        l = keep%rcon
        if ( keep%rcon > 0 ) then
          t = lambda(mep - keep%rcon + 1)
        else
          t = si_map(problem, sigma, keep%lmd(m))
        end if
        first = keep%last
        if ( keep%right > 0 ) then
          last = keep%first
        else
          last = first + 1
        end if
      end do

      if ( k == 0 .and. l == 0 ) then
        keep%av_dist = ZERO
      else if ( k == 0 ) then
        if ( keep%rcon > 0 ) then
          t = lambda(mep)
        else
          t = si_map(problem, sigma, keep%lmd(m))
        end if
        keep%av_dist = (s - t)/l
      else if ( l == 0 ) then
        if ( keep%lcon > 0 ) then
          t = lambda(1)
        else
          t = si_map(problem, sigma, keep%lmd(1))
        end if
        keep%av_dist = (t - r)/k
      else
        keep%av_dist = (s - r)/(k + l)
      end if

      check_res = options%abs_tol_residual > 0 &
             .or. options%rel_tol_residual > 0
      check_lmd = options%abs_tol_lambda > 0 &
             .or. options%rel_tol_lambda > 0

      t = 10*epsilon(ONE)

      do i = 1, m

        if ( keep%info%converged(i) /= 0 ) cycle

        converged = check_res .or. check_lmd .or. options%tol_X /= 0

        if ( check_res ) then
          s = abs(keep%lmd(i)) * max(options%rel_tol_residual, t)
          s = max(s, options%abs_tol_residual)
          converged = keep%info%residual_norms(i) <= s
        end if

        if ( check_lmd ) then
          if ( keep%info%err_lambda(i) >= ZERO ) then
            s = keep%av_dist * max(options%rel_tol_lambda, t*t)
            s = max(s, options%abs_tol_lambda)
            r = si_map_err(problem, sigma, keep%lmd(i), keep%info%err_lambda(i))
            converged = converged .and. r <= s
          else
            converged = .false.
          end if
        end if

        if ( options%tol_X > ZERO ) then
          converged = converged &
          .and. keep%info%err_X(i) >= ZERO &
          .and. keep%info%err_X(i) <= max(options%tol_X, t)
        else if ( options%tol_X < ZERO ) then
          converged = converged &
          .and. keep%info%err_X(i) >= ZERO &
          .and. keep%info%err_X(i) <= sqrt(epsilon(ONE))
        end if

        converged = converged .and. abs(keep%info%err_X(m + i)) < 0.05

        if ( converged ) keep%info%converged(i) = 1

      end do ! i = 1, m

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1

      go to 1

    end select

    if ( rci%job < 0 ) then

      info%non_converged = 0
      if ( rci%job /= SSMFE_DONE ) then
        info%left = keep%lcon
        info%right = keep%rcon
        info%non_converged = &
          max(0, left - keep%lcon) + max(0, right - keep%rcon)
      end if
      info%flag = keep%info%flag

      if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )

    end if

  contains

    real(PRECISION) function si_map( problem, sigma, mu )
      integer, intent(in) :: problem
      real(PRECISION), intent(in) :: sigma, mu
      if ( problem == 2 ) then
        si_map = sigma*mu/(mu - 1)
      else
        si_map = sigma + 1/mu
      end if
    end function si_map

    real(PRECISION) function si_map_err( problem, sigma, mu, delta )
      integer, intent(in) :: problem
      real(PRECISION), intent(in) :: sigma, mu, delta
      if ( problem == 2 ) then
        si_map_err = sigma*delta/(mu - 1)**2
      else
        si_map_err = delta/mu**2
      end if
    end function si_map_err

  end subroutine ssmfe_inverse_rci_double_complex

  subroutine ssmfe_direct_rci_double_complex &
    ( problem, nep, max_nep, lambda, block_size, rr_matrices, ind, &
      rci, keep, options, info )
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0

    real(kind = PRECISION), parameter :: MIN_REL_GAP = 1E-3

    type(ssmfe_rciz), intent(inout) :: rci
    integer, intent(in) :: problem
    integer, intent(in) :: nep
    integer, intent(in) :: max_nep
    integer, intent(in) :: block_size
    real(PRECISION), intent(inout) :: lambda(max_nep)
    complex(PRECISION), intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)
    integer, intent(inout) :: ind(block_size)
    type(ssmfe_options    ), intent(in) :: options
    type(ssmfe_expert_keep), intent(inout) :: keep
    type(ssmfe_inform     ), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, frmt, line
    character(7) :: word

    integer :: u_diag
    integer :: first, last
    integer :: i, j, k, l

    real(PRECISION) :: gap
    real(PRECISION) :: delta
    real(PRECISION) :: q, r, s, t

    u_diag = options%unit_diagnostic

    gap = options%left_gap

    info%flag = 0

    if ( nep < 0 &
        .or. options%max_left >= 0 .and. options%max_left < nep ) then
      info%flag = WRONG_LEFT
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( nep > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( block_size < 1 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      call ssmfe_errmsg( options, info )
      return
    end if

    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda    ) ) deallocate ( info%err_lambda     )
        if ( allocated(info%err_X         ) ) deallocate ( info%err_X          )
        if ( allocated(info%converged     ) ) deallocate ( info%converged      )
        allocate &
          ( info%residual_norms(max_nep), info%err_lambda(max_nep), &
            info%err_X(max_nep), info%converged(max_nep), stat = info%stat )
        if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
        if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
        if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
        if ( info%stat /= 0 ) return

        info%iteration = 0
        info%residual_norms = ZERO
        info%err_lambda = -ONE
        info%err_X = -ONE
        info%converged = 0
        info%left = 0
        info%next_left = huge(ONE)
        if ( nep <= 0 ) then
          rci%job = SSMFE_DONE
          return
        end if

        ! copy some interface optionss to the engine conttrols
        keep%options%err_est = options%err_est
        keep%left = nep
        keep%right = 0

        ! compute the number of extra vectors iterated for the sake
        ! of better convergence
        if ( options%extra_left < 0 ) then
          keep%options%extra_left = max(10, nep/10) ! use default
        else
          keep%options%extra_left = options%extra_left ! take user's value
        end if
        keep%options%minAprod = options%minAprod
        keep%options%minBprod = options%minBprod

      end if

      ! if the total number of eigenvalues to be computed (wanted + extra)
      ! exceeds the block size m, they are computed in portions
      ! of sizes keep%left (those left of sigma) and keep%right (right of it);
      ! once a portion is computed, iterations are restarted (see below)
      ! (note that more than wanted + extra eigenpairs are needed if gaps in
      ! the spectrum separating wanted eigenvalues from the rest are sought)

      i = nep + keep%options%extra_left
      if ( block_size < i .or. gap /= ZERO ) then
        keep%options%extra_left = min(keep%options%extra_left, block_size/2)
        keep%left = block_size - keep%options%extra_left
      else
        keep%left = nep
      end if
      ! if gaps in the spectrum are sought, the number of computed
      ! eigenpairs must be larger than requested
      if ( keep%left == nep .and. gap /= ZERO &
          .and. keep%options%extra_left > 0 ) then
        keep%left = keep%left + 1
        keep%options%extra_left = keep%options%extra_left - 1
      end if

      if ( rci%job == SSMFE_START ) then

        keep%options%min_gap = 0.05
        keep%lcon = 0
        ! initialize the convergence flags
        keep%left_converged = .false.
        if ( options%max_left < 0 ) then
          keep%max_left = max_nep + 1
        else
          keep%max_left = options%max_left
        end if

      end if

      if ( allocated(keep%lmd) ) deallocate ( keep%lmd )
      allocate ( keep%lmd(block_size), stat = info%stat )
      if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
      if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( info%stat /= 0 ) call ssmfe_errmsg( options, info )
      if ( info%stat /= 0 ) return

      keep%first = 1
      keep%last = block_size

      call ssmfe_msg( problem, options, nep, 0, block_size )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) rci%job = SSMFE_START

    else if ( rci%job == SSMFE_SAVE_CONVERGED .and. rci%i < 0 ) then

      if ( .not. keep%left_converged ) then
        if ( info%iteration >= options%max_iterations ) then
          rci%job = SSMFE_QUIT
          info%flag = MAX_NUM_ITERATIONS_EXCEEDED
        else if ( keep%lcon == max_nep ) then
          rci%job = SSMFE_QUIT
          info%flag = OUT_OF_STORAGE
        end if
      end if

      if ( rci%job /= SSMFE_QUIT ) then

        if ( keep%left_converged ) then
          rci%job = SSMFE_DONE
        else
          if ( keep%first > keep%left ) then
            rci%job = SSMFE_RESTART
            rci%jx = keep%first
            rci%nx = keep%last - keep%first + 1
            rci%k = 0
            return
          end if
        end if

      end if

      if ( rci%job < 0 ) then
        info%non_converged = max(0, nep - keep%lcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
        end if
        if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )
        return
      end if

    end if

    ! call the engine
1   call ssmfe &
      ( rci, problem, keep%left, keep%right, &
        block_size, keep%lmd, rr_matrices, ind, &
        keep%keep, keep%options, keep%info )

    select case (rci%job)

    case (SSMFE_SAVE_CONVERGED)

      if ( rci%i < 0 ) return

      ! if the number of converged eigenpairs exceeds the storage size
      ! reduce the number of newly converged eigenpairs to be saved
      if ( keep%lcon + rci%nx > max_nep ) then
        rci%nx = max_nep - keep%lcon
      end if

      do k = 0, rci%nx - 1
        j = keep%lcon + k + 1
        info%residual_norms(j) = keep%info%residual_norms(rci%jx + k)
        info%err_lambda(j) = keep%info%err_lambda(rci%jx + k)
        info%err_X(j) = keep%info%err_X(rci%jx + k)
        l = keep%info%converged(rci%jx + k)
        if ( l > 0 ) then
          info%converged(j) = max(1, info%iteration)
        else
          info%converged(j) = -max(1, info%iteration)
        end if
        lambda(j) = keep%lmd(rci%jx + k)
      end do
      j = rci%jx + rci%nx
      info%next_left = huge(ONE)
      if ( j <= block_size ) info%next_left = keep%lmd(j)
      keep%first = j

      if ( options%print_level == 2 .and. u_diag > NONE ) then

        frmt = '(i7, i8, es23.14, a, es10.1, 2es11.1)'

        first = rci%jx
        last = rci%jx + rci%nx - 1
        if ( .not. keep%left_converged .and. last < block_size &
              .and. keep%lcon + rci%nx < keep%max_left ) last = last + 1

        do i = first, last
          if ( keep%info%converged(i) /= 0 ) then
            word = '  yes'
          else
            word = '   no'
          end if
          write( u_diag, frmt ) &
            info%iteration, keep%lcon + i - first + 1, &
            keep%lmd(i), word, &
            keep%info%residual_norms(i), &
            keep%info%err_lambda(block_size + i), &
            keep%info%err_X(block_size + i)
          if ( keep%info%converged(i) == 0 ) exit
        end do

      end if

      keep%lcon = keep%lcon + rci%nx

      if ( options%print_level > 2 .and. u_diag > NONE ) then

        write( u_diag, '(/a, i5, a, i4)' ) &
          'Iteration: ', info%iteration, ', not converged: ', &
          max(0, nep - keep%lcon)

        head = &
 '      eigenvalues      |   locked  | residuals | eigenvalue | eigenvector'
        neck = &
 '                       |           |           |   errors   |    errors'
        line = &
 '-------------------------------------------------------------------------'

        frmt = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'

        write( u_diag, '(a/a)' ) &
          trim(line), trim(head)
        write( u_diag, '(a)' ) trim(neck)
        write( u_diag, '(a)' ) trim(line)

        do i = 1, block_size
          if ( keep%info%converged(i) /= 0 ) then
            word = '   yes'
          else
            word = '    no'
          end if
          write( u_diag, frmt ) &
            keep%lmd(i), ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            keep%info%err_lambda(block_size + i), '  |', &
            keep%info%err_X(block_size + i)
        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line)

      end if ! print_level > 2

      if ( gap == ZERO ) then

        keep%left_converged = keep%lcon >= nep
        if ( keep%left_converged ) then
          info%left = keep%lcon
          info%next_left = huge(ONE)
          if ( keep%first <= block_size ) then
            info%next_left = keep%lmd(keep%first)
          end if
        end if

      else

        if ( gap > 0 ) then
          t = gap
          q = ZERO
        else
          t = abs(gap) * keep%av_dist
          q = MIN_REL_GAP
        end if

        keep%left_converged = keep%left == 0 .or. keep%lcon >= keep%max_left
        do i = keep%lcon, nep + 1, -1
          j = maxloc(lambda(1 : i - 1), 1)
          r = lambda(i) - lambda(j)
          if ( r >= t .and. r > q * abs(lambda(i)) ) then
            info%left = i - 1
            info%next_left = lambda(i)
            keep%left_converged = .true.
            exit
          end if
        end do
        if ( keep%first <= block_size ) then
          if ( .not. keep%left_converged .and. keep%lcon >= nep ) then
            s =  keep%lmd(keep%first) - lambda(keep%lcon)
            r = keep%info%err_lambda(keep%first)
            if ( r > ZERO .and. r < s*sqrt(epsilon(ONE)) .and. s >= t ) then
              info%left = keep%lcon
              info%next_left = keep%lmd(keep%first)
              keep%left_converged = .true.
            end if
          end if
        end if

      end if

    case (SSMFE_CHECK_CONVERGENCE)

      delta = max(options%tol_X, sqrt(options%rel_tol_lambda))
      delta = sqrt(max(delta, epsilon(ONE)))

      l = keep%lcon
      if ( keep%lcon > 0 ) then
        t = lambda(keep%lcon)
      else
        t = keep%lmd(1)
      end if
      first = keep%first
      last = keep%last
      s = t
      do i = first, last
        if ( keep%info%err_X(block_size + i) > delta .and. i /= first ) exit
        l = l + 1
        s = keep%lmd(i)
      end do

      if ( keep%lcon > 0 ) then
        t = lambda(1)
      else
        t = keep%lmd(1)
      end if
      keep%av_dist = (s - t)/l

      r = ZERO
      if ( keep%lcon > 0 ) r = abs(lambda(1))
      r = r*epsilon(ONE)

      check_res = options%abs_tol_residual > 0 &
             .or. options%rel_tol_residual > 0
      check_lmd = options%abs_tol_lambda > 0 &
             .or. options%rel_tol_lambda > 0

      t = 10*epsilon(ONE)

      do i = 1, block_size

        if ( keep%info%converged(i) /= 0 ) cycle

        converged = check_res .or. check_lmd .or. options%tol_X /= 0

        if ( check_res ) then
          s = abs(keep%lmd(i)) * max(options%rel_tol_residual, t)
          s = max(s, options%abs_tol_residual)
          converged = keep%info%residual_norms(i) <= s
        end if

        if ( check_lmd ) then
          if ( keep%info%err_lambda(i) >= ZERO ) then
            s = keep%av_dist * max(options%rel_tol_lambda, t*t)
            s = max(s, options%abs_tol_lambda)
            converged = converged .and. keep%info%err_lambda(i) <= s
          else
            converged = .false.
          end if
        end if

        if ( options%tol_X > ZERO ) then
          converged = converged &
          .and. keep%info%err_X(i) >= ZERO &
          .and. keep%info%err_X(i) <= max(options%tol_X, t)
        else if ( options%tol_X < ZERO ) then
          converged = converged &
          .and. keep%info%err_X(i) >= ZERO &
          .and. keep%info%err_X(i) <= sqrt(epsilon(ONE))
        end if

        converged = converged .and. abs(keep%info%err_X(block_size + i)) < 0.05

        if ( converged ) keep%info%converged(i) = 1

      end do ! i = 1, block_size

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1

      go to 1

    end select

    if ( rci%job < 0 ) then

      info%non_converged = 0
      if ( rci%job /= SSMFE_DONE ) then
        info%non_converged = max(0, nep - keep%lcon)
        info%left = keep%lcon
      end if
      info%flag = keep%info%flag

      if ( info%flag /= 0 ) call ssmfe_errmsg( options, info )

    end if

  end subroutine ssmfe_direct_rci_double_complex

  subroutine ssmfe_errmsg( options, inform )
    integer, parameter :: NONE = -1

    type(ssmfe_options) :: options
    type(ssmfe_inform) :: inform

    logical :: oom

    integer :: u_errr, u_warn

    u_errr    = options%unit_error
    u_warn    = options%unit_warning

    oom = inform%flag == OUT_OF_MEMORY .and. u_errr > NONE
    if ( oom ) write( u_errr, '(/a/)' ) '??? Out of memory'

    select case ( inform%flag )

    case ( WRONG_LEFT )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong number of left eigenpairs'

    case ( WRONG_RIGHT )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong number of right eigenpairs'

    case ( WRONG_STORAGE_SIZE )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong eigenvalue storage size'

    case ( WRONG_SIGMA )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Zero sigma in buckling mode'

    case ( WRONG_BLOCK_SIZE )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong block size'

    case ( WRONG_ERR_EST )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong err_est'

    case ( WRONG_MINPROD )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error: minAprod and minBprod must be true'

    case ( B_NOT_POSITIVE_DEFINITE )

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong B or linear depended initial vectors'

    case ( NO_SEARCH_DIRECTIONS_LEFT )

      if ( u_warn > NONE ) &
        write( u_warn, '(/a,a/)' ) &
          '??? WARNING: iterations terminated because no further progress ', &
          'is possible'

    case ( MAX_NUM_ITERATIONS_EXCEEDED )

      if ( u_warn > NONE ) &
        write( u_warn, '(/a/)' ) &
          '??? WARNING: maximum number of iterations exceeded'

    case ( OUT_OF_STORAGE )

      if ( u_warn > NONE ) &
        write( u_warn, '(/a/)' ) &
          '??? WARNING: out of storage for converged eigenpairs'

    end select

  end subroutine ssmfe_errmsg

  subroutine ssmfe_msg( problem, options, left, right, m )
    type(ssmfe_options) :: options
    integer :: problem, left, right, m
    intent(in) :: options, left, right, m

    integer, parameter :: NONE = -1

    logical :: minAprod, minBprod
    integer :: print_lev, max_it
    integer :: u_diag
    real(kind = PRECISION) :: abs_tol, rel_tol, tol, abs_res, rel_res

    print_lev = options%print_level
    u_diag    = options%unit_diagnostic
    max_it    = options%max_iterations
    abs_tol   = options%abs_tol_lambda
    rel_tol   = options%rel_tol_lambda
    abs_res   = options%abs_tol_residual
    rel_res   = options%rel_tol_residual
    tol       = options%tol_X
    minAprod  = options%minAprod
    minBprod  = options%minBprod

    if ( print_lev <= NONE ) u_diag = NONE
    if ( u_diag <= NONE .and. print_lev > NONE ) print_lev = 0

    if ( print_lev > 0 ) then
      if ( problem == 0 ) then
        write( u_diag, '(/a)' ) &
          'Solving the standard eigenvalue problem A x = lambda x'
      else
        write( u_diag, '(/a)' ) &
          'Solving the generalized eigenvalue problem A x = lambda B x'
      end if
      if ( left > 0 ) &
        write ( u_diag, '(a, i4)' ) &
          'leftmost eigenpairs requested:', left
      if ( left >= 0 .and. right > 0  ) &
        write ( u_diag, '(a, i4)' ) &
          'rightmost eigenpairs requested:', right
      write( u_diag, '(a, i4 )' ) 'iterated subspace dimension:', m
      if ( abs_res >= 0 .and. rel_res >= 0 .and. abs_res + rel_res > 0 ) &
        write( u_diag, '(a, es8.0, a, es8.0 )' ) &
          'residual tolerances: absolute =', abs_res, &
          ', relative = ', rel_res
      if ( abs_tol >= 0 .and. rel_tol >= 0 .and. abs_tol + rel_tol > 0 ) &
        write( u_diag, '(a, es8.0, a, es8.0 )' ) &
          'eigenvalue error tolerances: absolute =', abs_tol, &
          ', relative = ', rel_tol
      if ( tol > 0.0 ) then
        write( u_diag, '(a, es8.0)' ) &
          'eigenvector error tolerance:', tol
      else if ( tol < 0.0 ) then
        write( u_diag, '(a, es8.0)' ) &
          'eigenvector error tolerance:', sqrt(epsilon(1.0D0))
      end if
      if ( minAprod ) &
        write ( u_diag, '(a)' ) &
          'the number of multiplications by A is minimized'
      if ( minBprod .and. problem /= 0 ) &
        write ( u_diag, '(a)' ) &
          'the number of multiplications by B is minimized'
    end if

    if ( print_lev == 2 .and. max_it > 0 ) &
      write( u_diag, '(/60x,a/a,2x,a,7x,a,6x,a,2x,a,2x,a,1x,a/)' ) &
        'Estimated Errors', &
        'Iteration', 'Index', 'Eigenvalue', 'Locked', 'Residual', &
        'Eigenvalue', 'Eigenvector'

  end subroutine ssmfe_msg

end module SPRAL_ssmfe_expert

