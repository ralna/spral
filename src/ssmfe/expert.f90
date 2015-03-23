! COPYRIGHT (c) 2014 The Science and Technology Facilities Council (STFC)
! Original date 1 December 2014, Version 1.0.0
!
! Written by: Evgueni Ovtchinnikov
!
module SPRAL_ssmfe_expert

  use SPRAL_ssmfe_core, ssmfe_work => ssmfe_keep, ssmfe_opts => ssmfe_options

  integer, parameter, private :: PRECISION = kind(1.0D0)
  integer, parameter, private :: SSMFE_ITYPE = 1
  character, parameter, private :: SSMFE_JOBZ = 'V', SSMFE_UPLO = 'U'
  
  integer, parameter, private :: WRONG_RCI_JOB      = -1
  integer, parameter, private :: WRONG_BLOCK_SIZE   = -2
  integer, parameter, private :: WRONG_ERR_EST      = -3
  integer, parameter, private :: WRONG_MINPROD      = -4
  integer, parameter, private :: WRONG_EXTRAS       = -5
  integer, parameter, private :: WRONG_MIN_GAP      = -6
  integer, parameter, private :: WRONG_CF_MAX       = -7
  integer, parameter, private :: WRONG_LEFT         = -11
  integer, parameter, private :: WRONG_RIGHT        = -12
  integer, parameter, private :: WRONG_STORAGE_SIZE = -13
  integer, parameter, private :: WRONG_SIGMA        = -14
  
  integer, parameter, private :: OUT_OF_MEMORY           = -100
  integer, parameter, private :: B_NOT_POSITIVE_DEFINITE = -200

  integer, parameter, private :: NO_SEARCH_DIRECTIONS_LEFT = 1

  integer, parameter, private :: SSMFE_DONE  = -1
  integer, parameter, private :: SSMFE_QUIT  = -2
  integer, parameter, private :: SSMFE_ABORT = -3

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

  interface ssmfe_terminate
    module procedure &
      ssmfe_terminate_double, &
      ssmfe_delete_keep_double
  end interface 
  
  ! error estimation schemes
  integer, parameter, private :: SSMFE_RESIDUAL = 1
  integer, parameter, private :: SSMFE_KINEMATIC = 2
  
  type ssmfe_options
  
    integer :: print_level = 0
    integer :: unit_error = 6
    integer :: unit_warning = 6 
    integer :: unit_diagnostic = 6
        
    integer :: max_iterations = 100
    integer :: user_X = 0
    integer :: err_est = SSMFE_KINEMATIC

    double precision :: abs_tol_lambda   =  0.0
    double precision :: rel_tol_lambda   =  0.0
    double precision :: abs_tol_residual =  0.0
    double precision :: rel_tol_residual =  0.0
    double precision :: tol_X = -1.0
    
    double precision :: left_gap = 0.0
    double precision :: right_gap = 0.0

    integer :: extra_left = -1, extra_right = -1    
    integer :: max_left = -1, max_right = -1

    logical :: minAprod = .true.
    logical :: minBprod = .true.

  end type ssmfe_options
  
  ! private subroutines and types

  type ssmfe_keep
  
    private
    
    logical :: left_converged, right_converged

    integer :: problem
    integer :: user_X = 0
    integer :: work_space = 0
    integer :: left, right
    integer :: lft, rgt
    integer :: max_left  = -1
    integer :: max_right = -1
    integer :: lcon = 0
    integer :: rcon = 0
    integer :: first, last
    integer :: step = 0
    
    double precision :: av_dist
    double precision, allocatable :: lmd(:)

    type(ssmfe_inform) :: info
    type(ssmfe_work) :: keep
    type(ssmfe_opts) :: options

  end type ssmfe_keep

contains

  subroutine ssmfe_expert_double &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep

    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    real(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_rci_double &
      ( 0, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_double

  subroutine ssmfe_expert_double_complex &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    complex(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_rci_double_complex &
      ( 0, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_double_complex

  subroutine ssmfe_expert_gen_double &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep

    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    real(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_rci_double &
      ( 1, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_double

  subroutine ssmfe_expert_gen_double_complex &
      ( rci, left, mep, lambda, m, rr, ind, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: m
    complex(PRECISION), intent(inout) :: rr(2*m, 2*m, 3)
    integer, intent(inout) :: ind(m)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_rci_double_complex &
      ( 1, left, mep, lambda, m, rr, ind, rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_double_complex

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
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_inverse_rci_double &
      ( 0, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_shift_double

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
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_inverse_rci_double_complex &
      ( 0, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_shift_double_complex

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
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_inverse_rci_double &
      ( 1, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_shift_double

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
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_inverse_rci_double_complex &
      ( 1, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_gen_shift_double_complex

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
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_inverse_rci_double &
      ( 2, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_buckling_double

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
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_inverse_rci_double_complex &
      ( 2, sigma, left, right, mep, lambda, m, rr, ind, &
        rci, keep, options, inform )

  end subroutine ssmfe_expert_buckling_double_complex

  subroutine ssmfe_inverse_rci_double &
    ( problem, sigma, left, right, &
      max_nep, lambda, block_size, rr_matrices, ind, &
      rci, keep, options, info )

    implicit none
    
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_APPLY_A           = 1
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE  = 9
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    double precision, parameter :: ZERO = 0.0D0, ONE = 1.0D0
    double precision, parameter :: NIL = ZERO
    double precision, parameter :: UNIT = ONE
    
    real(kind = PRECISION), parameter :: GAP_REL_A = 1E-3, GAP_REL_D = 1E-1

    integer, intent(in) :: problem
    double precision, intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: max_nep
    double precision, intent(inout) :: lambda(max_nep)
    integer, intent(in) :: block_size
    double precision, intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)
    integer, intent(inout) :: ind(block_size)
    type(ssmfe_rcid    ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, form, line
    character(7) :: word

    integer :: u_diag, u_errr, verb
    integer :: m, mm, mep, new, ncon, first, last, step
    integer :: i, j, k, l
    
    double precision :: left_gap
    double precision :: right_gap
    double precision :: delta
    double precision :: q, r, s, t
    
    u_diag = options%unit_diagnostic
    u_errr = options%unit_error
    verb   = options%print_level
    
    left_gap  = options%left_gap
    right_gap = options%right_gap

    info%flag = 0
    
    if ( options%max_left >= 0 .and. options%max_left < left ) then
      info%flag = WRONG_LEFT
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong number of left eigenpairs'
      return
    end if

    if ( options%max_right >= 0 .and. options%max_right < right ) then
      info%flag = WRONG_RIGHT
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong number of right eigenpairs'
      return
    end if

    if ( left + right > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong eigenvalue storage size'
      return
    end if
    
    if ( block_size < 2 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a, i5/)' ) &
          '??? Wrong block size', block_size
      return
    end if
    
    if ( problem == 2 .and. sigma == ZERO ) then
      info%flag = WRONG_SIGMA
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Zero sigma in buckling mode'
      return
    end if
    
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        info%iteration = 0
        mep = max_nep
        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda) ) deallocate ( info%err_lambda )
        if ( allocated(info%err_X) ) deallocate ( info%err_X )
        if ( allocated(info%converged) ) deallocate ( info%converged )
        allocate ( info%residual_norms(mep), info%err_lambda(mep), &
                   info%err_X(mep), info%converged(mep) )
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

        keep%options%minAprod = .true.
        keep%options%minBprod = problem /= 0
        
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
          if ( l < i .or. left_gap /= ZERO ) then ! not enough space
            keep%options%extra_left = min(keep%options%extra_left, l/2)
            keep%left = l - keep%options%extra_left ! portion size
            ! limiting extras number to l/2 ensures that the portion size
            ! is substantial thus reducing the number of restarts
          else
            keep%left = left ! can be computed in one go
          end if
          if ( m - l < j .or. right_gap /= ZERO ) then ! not enough space
            keep%options%extra_right = min(keep%options%extra_right, (m - l)/2)
            keep%right = m - l - keep%options%extra_right ! portion size
          else
            keep%right = right ! can be computed in one go
          end if
        else
          if ( m < i .or. left_gap /= ZERO ) then
            keep%options%extra_left = min(keep%options%extra_left, m/2)
            keep%left = m - keep%options%extra_left
          else
            keep%left = left
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

        keep%user_X = options%user_X
        keep%options%min_gap = 0.05

        keep%lcon = 0
        keep%rcon = 0
        
        ! initialize the convergence flags
        keep%left_converged = left == 0
        keep%right_converged = right == 0

        ! initialize the number of converged eigenpairs before gaps
        keep%lft = 0
        keep%rgt = 0
        
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

        deallocate ( keep%lmd )

      end if

      ! allocate work arrays
      ! keep%lmd: eigenvalues array of ssmfe solver
      ! keep%v: workspace for ssmfe solver matrices
      m = block_size
      mm = m + m
      allocate ( keep%lmd(m) ) !, keep%ind(m), keep%v(mm, mm, 3) )

      keep%first = 1
      keep%last = m

      call ssmfe_msg( problem, options, left, right, m, 0 )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) &
        rci%job = SSMFE_START
      
    else if ( rci%job == SSMFE_SAVE_CONVERGED ) then

      mep = max_nep
      ncon = keep%lcon + keep%rcon
      
      if ( .not. (keep%left_converged .and. keep%right_converged) ) then
        if ( info%iteration >= options%max_iterations ) then
          rci%job = SSMFE_QUIT
          info%flag = 2
        else if ( ncon == mep &
          .or. .not. keep%left_converged &
               .and. keep%lcon >= mep - max(right, keep%rcon) &
          .or. .not. keep%right_converged &
               .and. keep%rcon >= mep - max(left, keep%lcon) &
            ) then
          rci%job = SSMFE_QUIT
          info%flag = 3
        end if
      end if

      if ( rci%job /= SSMFE_QUIT ) then
        if ( keep%left_converged ) keep%left = 0
        if ( keep%right_converged ) keep%right = 0
        if ( keep%left_converged .and. keep%right_converged ) then
          rci%job = SSMFE_DONE
        else
          if ( rci%i < 0 .and. & ! restart after printout
               (keep%left > 0 .and. keep%first > keep%left .or. &
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
        info%data = max(0, left - keep%lcon) + max(0, right - keep%rcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
          info%right = keep%rcon
        end if
        if ( info%flag /= 0 ) &
          call ssmfe_msg &
            ( problem, options, keep%left, keep%right, m, info%flag )
        return
      end if

    end if
    
    ! call the engine
    call ssmfe &
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
        if ( .not. keep%left_converged ) then
          if ( j > m ) then
            info%next_left = sigma
          else
            if ( keep%lcon + new > 0 .and. keep%lmd(j)*keep%lmd(1) > 0 ) &
              info%next_left = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%first = j
      else
        j = rci%jx - new
        m = block_size
        if ( .not. keep%right_converged ) then
          if ( j < 1 ) then
            info%next_right = sigma
          else 
            if ( keep%rcon + new > 0 .and. keep%lmd(j)*keep%lmd(m) > 0 ) &
              info%next_right = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%last = j
      end if
      
      ! report the convergence situation
      if ( options%print_level == 2 .and. u_diag > NONE ) then

        if ( PRECISION == kind(1.0) ) then
          form = '(i7, i8, es18.6, 4x, a, 3es11.1)'
        else
          form = '(i7, i8, es23.14, a, es10.1, 2es11.1)'
        end if

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
          write( u_diag, form ) &
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

        if ( PRECISION == kind(1.0) ) then
          form = '(es17.6,5x,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        else
          form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        end if

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

          write( u_diag, form ) &
            s, ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            t, '  |', keep%info%err_X(block_size + i)

        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line) 

      end if ! print_level > 2

      select case ( rci%i )
      
      case ( 1 )

        if ( keep%left == 0 ) then

          keep%left_converged = .true.

        else if ( left_gap == ZERO ) then

          keep%left_converged = keep%lcon >= left
          if ( keep%left_converged ) then
            info%left = keep%lcon
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
            q = GAP_REL_A
          end if

          keep%left_converged = (keep%lcon >= keep%max_left)
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
          if ( .not. keep%left_converged .and. keep%first <= block_size ) then
            if ( keep%lcon >= left &
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

        if ( keep%right == 0 ) then

          keep%right_converged = .true.

        else if ( right_gap == ZERO ) then

          keep%right_converged = keep%rcon >= right
          if ( keep%right_converged ) then
            info%right = keep%rcon
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
            q = GAP_REL_A
          end if

          keep%right_converged = (keep%rcon >= keep%max_right)
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

          if ( .not. keep%right_converged .and. keep%last >= 1 ) then
            if ( keep%rcon >= right &
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

        converged = converged .and. abs(keep%info%err_X(m + i)) < 0.01

        if ( converged ) keep%info%converged(i) = 1
        
      end do ! i = 1, m

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1
        
    end select

    if ( rci%job < 0 ) then

      info%left = keep%lcon
      info%right = keep%rcon
      if ( rci%job /= SSMFE_DONE ) then
        info%data = max(0, left - keep%lcon) + max(0, right - keep%rcon)
      else
        info%data = 0
      end if
      info%flag = keep%info%flag

      if ( info%flag /= 0 ) &
        call ssmfe_msg &
          ( problem, options, keep%left, keep%right, m, info%flag )

    end if

  contains

    double precision function si_map( problem, sigma, mu )
      integer, intent(in) :: problem
      double precision, intent(in) :: sigma, mu
      if ( problem == 2 ) then
        si_map = sigma*mu/(mu - 1)
      else
        si_map = sigma + 1/mu
      end if
    end function si_map

    double precision function si_map_err( problem, sigma, mu, delta )
      integer, intent(in) :: problem
      double precision, intent(in) :: sigma, mu, delta
      if ( problem == 2 ) then
        si_map_err = sigma*delta/(mu - 1)**2
      else
        si_map_err = delta/mu**2
      end if
    end function si_map_err

  end subroutine ssmfe_inverse_rci_double

  subroutine ssmfe_direct_rci_double &
    ( problem, nep, max_nep, lambda, block_size, rr_matrices, ind, &
      rci, keep, options, info )

    implicit none
    
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    double precision, parameter :: ZERO = 0.0D0, ONE = 1.0D0
    double precision, parameter :: NIL = ZERO
    double precision, parameter :: UNIT = ONE
    
    real(kind = PRECISION), parameter :: GAP_REL_A = 1E-3, GAP_REL_D = 1E-1

    type(ssmfe_rcid), intent(inout) :: rci
    integer, intent(in) :: problem
    integer, intent(in) :: nep
    integer, intent(in) :: max_nep
    integer, intent(in) :: block_size
    double precision, intent(inout) :: lambda(max_nep)
    double precision, intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)
    integer, intent(inout) :: ind(block_size)
    type(ssmfe_options), intent(in) :: options
    type(ssmfe_keep), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, form, line
    character(7) :: word

    integer :: u_diag, u_errr, verb
    integer :: mm, first, last
    integer :: i, j, k, l
    
    double precision :: gap
    double precision :: delta
    double precision :: q, r, s, t
    
    u_diag = options%unit_diagnostic
    u_errr = options%unit_error
    verb   = options%print_level
    
    gap = options%left_gap

    info%flag = 0

    if ( nep > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error in ssmfe_direct_solve: wrong eigenvalue storage size'
      return
    end if
    
    if ( block_size < 1 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a, i5/)' ) &
          '??? Error in ssmfe_direct_solve: wrong block size', block_size
      return
    end if
    
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        info%iteration = 0
        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda) ) deallocate ( info%err_lambda )
        if ( allocated(info%err_X) ) deallocate ( info%err_X )
        if ( allocated(info%converged) ) deallocate ( info%converged )
        allocate ( info%residual_norms(max_nep), info%err_lambda(max_nep), &
                   info%err_X(max_nep), info%converged(max_nep) )
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

        keep%user_X = options%user_X
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
      mm = 2*block_size
      allocate ( keep%lmd(block_size) )

      keep%first = 1
      keep%last = block_size

      call ssmfe_msg( problem, options, nep, 0, block_size, 0 )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) rci%job = SSMFE_START
      
    else if ( rci%job == SSMFE_SAVE_CONVERGED ) then

      if ( .not. keep%left_converged ) then
        if ( info%iteration >= options%max_iterations ) then
          rci%job = SSMFE_QUIT
          info%flag = 2
        else if ( keep%lcon == max_nep ) then
          rci%job = SSMFE_QUIT
          info%flag = 3
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
        info%data = max(0, nep - keep%lcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
        end if
        if ( info%flag /= 0 ) &
          call ssmfe_msg &
            ( problem, options, keep%left, 0, block_size, info%flag )
        return
      end if

    end if
    
    ! call the engine
    call ssmfe &
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
      if ( j > block_size ) then
        info%next_left = huge(ONE)
      else 
        info%next_left = keep%lmd(j)
      end if
      keep%first = j

      if ( options%print_level == 2 .and. u_diag > NONE ) then

        if ( PRECISION == kind(1.0) ) then
          form = '(i7, i8, es18.6, 4x, a, 3es11.1)'
        else
          form = '(i7, i8, es23.14, a, es10.1, 2es11.1)'
        end if

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
          write( u_diag, form ) &
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

        if ( PRECISION == kind(1.0) ) then
          form = '(es17.6,5x,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        else
          form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        end if

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
          write( u_diag, form ) &
            keep%lmd(i), ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            keep%info%err_lambda(block_size + i), '  |', &
            keep%info%err_X(block_size + i)
        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line) 

      end if ! print_level > 2

      if ( keep%left == 0 ) then

        keep%left_converged = .true.

      else if ( gap == ZERO ) then

        keep%left_converged = keep%lcon >= nep
        if ( keep%left_converged ) then
          info%left = keep%lcon
          if ( keep%first <= block_size ) then
            info%next_left = keep%lmd(keep%first)
          else
            info%next_left = huge(ONE)
          end if
        end if

      else
          
        if ( gap > 0 ) then
          t = gap
          q = ZERO
        else
          t = abs(gap) * keep%av_dist
          q = GAP_REL_A
        end if

        keep%left_converged = (keep%lcon >= keep%max_left)
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
        
      if ( l == 0 ) then
        keep%av_dist = ZERO
      else
        if ( keep%lcon > 0 ) then
          t = lambda(1)
        else
          t = keep%lmd(1)
        end if
        keep%av_dist = (s - t)/l
      end if
      
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

        converged = converged .and. abs(keep%info%err_X(block_size + i)) < 0.01

        if ( converged ) keep%info%converged(i) = 1
        
      end do ! i = 1, block_size

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1
        
    end select

    if ( rci%job < 0 ) then

      if ( rci%job /= SSMFE_DONE ) then
        info%data = max(0, nep - keep%lcon)
        info%left = keep%lcon
      else
        info%data = 0
      end if
      info%flag = keep%info%flag

      if ( info%flag /= 0 ) &
        call ssmfe_msg &
          ( problem, options, keep%left, 0, block_size, info%flag )

    end if
      
  end subroutine ssmfe_direct_rci_double

  subroutine ssmfe_delete_keep_double( keep )

    implicit none
    
    type(ssmfe_keep), intent(inout) :: keep

    if ( allocated(keep%lmd) ) deallocate( keep%lmd )
    call ssmfe_delete_work_double( keep%keep )
    call ssmfe_delete_info_double( keep%info )

  end subroutine ssmfe_delete_keep_double

  subroutine ssmfe_terminate_double( keep, info )

    implicit none
    
    type(ssmfe_keep), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: info
    
    call ssmfe_delete_keep_double( keep )
    call ssmfe_delete_info_double( info )

  end subroutine ssmfe_terminate_double

  subroutine ssmfe_inverse_rci_double_complex &
    ( problem, sigma, left, right, &
      max_nep, lambda, block_size, rr_matrices, ind, &
      rci, keep, options, info )

    implicit none
    
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_APPLY_A           = 1
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE  = 9
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    double precision, parameter :: ZERO = 0.0D0, ONE = 1.0D0
    complex(PRECISION), parameter :: NIL = ZERO
    complex(PRECISION), parameter :: UNIT = ONE
    
    real(kind = PRECISION), parameter :: GAP_REL_A = 1E-3, GAP_REL_D = 1E-1

    integer, intent(in) :: problem
    double precision, intent(in) :: sigma
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: max_nep
    double precision, intent(inout) :: lambda(max_nep)
    integer, intent(in) :: block_size
    complex(PRECISION), intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)
    integer, intent(inout) :: ind(block_size)
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, form, line
    character(7) :: word

    integer :: u_diag, u_errr, verb
    integer :: m, mm, mep, new, ncon, first, last, step
    integer :: i, j, k, l
    
    double precision :: left_gap
    double precision :: right_gap
    double precision :: delta
    double precision :: q, r, s, t
    
    u_diag = options%unit_diagnostic
    u_errr = options%unit_error
    verb   = options%print_level
    
    left_gap = options%left_gap
    right_gap = options%right_gap

    info%flag = 0
    
    if ( options%max_left >= 0 .and. options%max_left < left ) then
      info%flag = WRONG_LEFT
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error in ssmfe_solve: wrong number of left eigenpairs'
      return
    end if

    if ( options%max_right >= 0 .and. options%max_right < right ) then
      info%flag = WRONG_RIGHT
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error in ssmfe_solve: wrong number of right eigenpairs'
      return
    end if

    if ( left + right > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error in ssmfe_solve: wrong eigenvalue storage size'
      return
    end if
    
    if ( block_size < 2 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error in ssmfe_solve: wrong block size'
      return
    end if
    
    if ( problem == 2 .and. sigma == ZERO ) then
      info%flag = WRONG_SIGMA
      rci%job = SSMFE_QUIT
      if ( options%unit_error > NONE .and. options%print_level > NONE ) &
        write( options%unit_error, '(/a, i9/)' ) &
          '??? Zero sigma in buckling mode'
      return
    end if
    
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        info%iteration = 0
        mep = max_nep
        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda) ) deallocate ( info%err_lambda )
        if ( allocated(info%err_X) ) deallocate ( info%err_X )
        if ( allocated(info%converged) ) deallocate ( info%converged )
        allocate ( info%residual_norms(mep), info%err_lambda(mep), &
                   info%err_X(mep), info%converged(mep) )
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

        keep%options%minAprod = .true.
        keep%options%minBprod = problem /= 0
        
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
          if ( l < i .or. left_gap /= ZERO ) then ! not enough space
            keep%options%extra_left = min(keep%options%extra_left, l/2)
            keep%left = l - keep%options%extra_left ! portion size
            ! limiting extras number to l/2 ensures that the portion size
            ! is substantial thus reducing the number of restarts
          else
            keep%left = left ! can be computed in one go
          end if
          if ( m - l < j .or. right_gap /= ZERO ) then ! not enough space
            keep%options%extra_right = min(keep%options%extra_right, (m - l)/2)
            keep%right = m - l - keep%options%extra_right ! portion size
          else
            keep%right = right ! can be computed in one go
          end if
        else
          if ( m < i .or. left_gap /= ZERO ) then
            keep%options%extra_left = min(keep%options%extra_left, m/2)
            keep%left = m - keep%options%extra_left
          else
            keep%left = left
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

        keep%user_X = options%user_X
        keep%options%min_gap = 0.05

        keep%lcon = 0
        keep%rcon = 0
        
        ! initialize the convergence flags
        keep%left_converged = left == 0
        keep%right_converged = right == 0

        ! initialize the number of converged eigenpairs before gaps
        keep%lft = 0
        keep%rgt = 0
        
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

        deallocate ( keep%lmd )

      end if

      ! allocate work arrays
      ! keep%lmd: eigenvalues array of ssmfe solver
      ! keep%v: workspace for ssmfe solver matrices
      m = block_size
      mm = m + m
      allocate ( keep%lmd(m) ) !, keep%ind(m), keep%v(mm, mm, 3) )

      keep%first = 1
      keep%last = m

      call ssmfe_msg( problem, options, left, right, m, 0 )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) &
        rci%job = SSMFE_START
      
    else if ( rci%job == SSMFE_SAVE_CONVERGED ) then

      mep = max_nep
      ncon = keep%lcon + keep%rcon
      
      if ( rci%job < 0 ) then
        info%data = max(0, left - keep%lcon) + max(0, right - keep%rcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
          info%right = keep%rcon
        end if
        info%flag = keep%info%flag
        if ( info%flag /= 0 ) &
          call ssmfe_msg &
            ( problem, options, keep%left, keep%right, m, info%flag )
        return
      end if

      if ( .not. (keep%left_converged .and. keep%right_converged) ) then
        if ( info%iteration >= options%max_iterations ) then
          rci%job = SSMFE_QUIT
          info%flag = 2
        else if ( ncon == mep &
          .or. .not. keep%left_converged &
               .and. keep%lcon >= mep - max(right, keep%rcon) &
          .or. .not. keep%right_converged &
               .and. keep%rcon >= mep - max(left, keep%lcon) &
            ) then
          rci%job = SSMFE_QUIT
          info%flag = 3
        end if
      end if

      if ( rci%job /= SSMFE_QUIT ) then

        if ( keep%left_converged ) keep%left = 0
        if ( keep%right_converged ) keep%right = 0
        
        if ( keep%left_converged .and. keep%right_converged ) then
          rci%job = SSMFE_DONE
        else
          if ( & !block_size < keep%left + keep%right .and. &
               (keep%left > 0 .and. keep%first > keep%left .or. &
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
      
      if ( rci%job < 0 ) return

    end if
    
    ! call the engine
    call ssmfe &
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
        if ( .not. keep%left_converged ) then
          if ( j > m ) then
            info%next_left = sigma
          else
            if ( keep%lcon + new > 0 .and. keep%lmd(j)*keep%lmd(1) > 0 ) &
              info%next_left = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%first = j
      else
        j = rci%jx - new
        m = block_size
        if ( .not. keep%right_converged ) then
          if ( j < 1 ) then
            info%next_right = sigma
          else 
            if ( keep%rcon + new > 0 .and. keep%lmd(j)*keep%lmd(m) > 0 ) &
              info%next_right = si_map(problem, sigma, keep%lmd(j))
          end if
        end if
        keep%last = j
      end if

      ! report the convergence situation
      if ( options%print_level == 2 .and. u_diag > NONE ) then

        if ( PRECISION == kind(1.0) ) then
          form = '(i7, i8, es18.6, 4x, a, 3es11.1)'
        else
          form = '(i7, i8, es23.14, a, es10.1, 2es11.1)'
        end if

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
          write( u_diag, form ) &
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

        if ( PRECISION == kind(1.0) ) then
          form = '(es17.6,5x,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        else
          form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        end if

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

          write( u_diag, form ) &
            s, ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            t, '  |', keep%info%err_X(block_size + i)

        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line) 

      end if ! print_level > 2

      select case ( rci%i )
      
      case ( 1 )

        if ( keep%left == 0 ) then

          keep%left_converged = .true.

        else if ( keep%lcon >= keep%max_left ) then

          info%left = keep%max_left
          info%next_left = sigma
          keep%left_converged = .true.

        else if ( left_gap == ZERO ) then

          keep%left_converged = keep%lcon >= left
          if ( keep%left_converged ) then
            info%left = keep%lcon
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
            q = GAP_REL_A
          end if

          keep%left_converged = .false.
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

        if ( keep%right == 0 ) then

          keep%right_converged = .true.

        else if ( keep%rcon >= keep%max_right ) then

          info%right = keep%max_right
          info%next_right = sigma
          keep%right_converged = .true.

        else if ( right_gap == ZERO ) then

          keep%right_converged = keep%rcon >= right
          if ( keep%right_converged ) then
            info%right = keep%rcon
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
            q = GAP_REL_A
          end if

          keep%right_converged = .false.
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

        converged = converged .and. abs(keep%info%err_X(m + i)) < 0.01

        if ( converged ) keep%info%converged(i) = 1
        
      end do ! i = 1, m

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1
        
    end select

    if ( rci%job < 0 ) then

      info%left = keep%lcon
      info%right = keep%rcon
      if ( rci%job /= SSMFE_DONE ) then
        info%data = max(0, left - keep%lcon) + max(0, right - keep%rcon)
      else
        info%data = 0
      end if
      info%flag = keep%info%flag

      if ( info%flag /= 0 ) &
        call ssmfe_msg &
          ( problem, options, keep%left, keep%right, m, info%flag )

    end if

  contains

    double precision function si_map( problem, sigma, mu )
      integer, intent(in) :: problem
      double precision, intent(in) :: sigma, mu
      if ( problem == 2 ) then
        si_map = sigma*mu/(mu - 1)
      else
        si_map = sigma + 1/mu
      end if
    end function si_map

    double precision function si_map_err( problem, sigma, mu, delta )
      integer, intent(in) :: problem
      double precision, intent(in) :: sigma, mu, delta
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

    implicit none
    
    integer, parameter :: SSMFE_START             = 0
    integer, parameter :: SSMFE_CHECK_CONVERGENCE = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED    = 5
    integer, parameter :: SSMFE_RESTART           = 999

    integer, parameter :: NONE = -1

    double precision, parameter :: ZERO = 0.0D0, ONE = 1.0D0
    complex(PRECISION), parameter :: NIL = ZERO
    complex(PRECISION), parameter :: UNIT = ONE
    
    real(kind = PRECISION), parameter :: GAP_REL_A = 1E-3, GAP_REL_D = 1E-1

    type(ssmfe_rciz), intent(inout) :: rci
    integer, intent(in) :: problem
    integer, intent(in) :: nep
    integer, intent(in) :: max_nep
    integer, intent(in) :: block_size
    double precision, intent(inout) :: lambda(max_nep)
    complex(PRECISION), intent(inout) :: &
      rr_matrices(2*block_size, 2*block_size, 3)
    integer, intent(inout) :: ind(block_size)
    type(ssmfe_options), intent(in) :: options
    type(ssmfe_keep), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: info

    logical :: converged
    logical :: check_res, check_lmd

    character(78) :: head, neck, form, line
    character(7) :: word

    integer :: u_diag, u_errr, verb
    integer :: mm, first, last
    integer :: i, j, k, l
    
    double precision :: gap
    double precision :: delta
    double precision :: q, r, s, t
    
    u_diag = options%unit_diagnostic
    u_errr = options%unit_error
    verb   = options%print_level
    
    gap = options%left_gap

    info%flag = 0
    
    if ( nep > max_nep ) then
      info%flag = WRONG_STORAGE_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error in ssmfe_inverse_solve: wrong eigenvalue storage size'
      return
    end if
    
    if ( block_size < 2 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      if ( u_errr > NONE .and. verb > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Error in ssmfe_inverse_solve: wrong block size'
      return
    end if
    
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
        .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) ) then

      if ( rci%job == SSMFE_START ) then

        info%iteration = 0
        if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
        if ( allocated(info%err_lambda) ) deallocate ( info%err_lambda )
        if ( allocated(info%err_X) ) deallocate ( info%err_X )
        if ( allocated(info%converged) ) deallocate ( info%converged )
        allocate ( info%residual_norms(max_nep), info%err_lambda(max_nep), &
                   info%err_X(max_nep), info%converged(max_nep) )
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

        keep%user_X = options%user_X
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
      mm = 2*block_size
      allocate ( keep%lmd(mm) )

      keep%first = 1
      keep%last = block_size

      call ssmfe_msg( problem, options, nep, 0, block_size, 0 )

      if ( rci%job == SSMFE_RESTART .and. rci%k == 0 ) rci%job = SSMFE_START
      
    else if ( rci%job == SSMFE_SAVE_CONVERGED ) then

      if ( .not. keep%left_converged ) then
        if ( info%iteration >= options%max_iterations ) then
          rci%job = SSMFE_QUIT
          info%flag = 2
        else if ( keep%lcon == max_nep ) then
          rci%job = SSMFE_QUIT
          info%flag = 3
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
        info%data = max(0, nep - keep%lcon)
        if ( rci%job /= SSMFE_DONE ) then
          info%left = keep%lcon
        end if
        if ( info%flag /= 0 ) &
          call ssmfe_msg &
            ( problem, options, keep%left, 0, block_size, info%flag )
        return
      end if

    end if
    
    ! call the engine
    call ssmfe &
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
      if ( j > block_size ) then
        info%next_left = huge(ONE)
      else 
        info%next_left = keep%lmd(j)
      end if
      keep%first = j

      if ( options%print_level == 2 .and. u_diag > NONE ) then

        if ( PRECISION == kind(1.0) ) then
          form = '(i7, i8, es18.6, 4x, a, 3es11.1)'
        else
          form = '(i7, i8, es23.14, a, es10.1, 2es11.1)'
        end if

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
          write( u_diag, form ) &
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

        if ( PRECISION == kind(1.0) ) then
          form = '(es17.6,5x,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        else
          form = '(es22.14,a,1x,a,2x,a,es9.1,a,es10.1,a,es10.1)'
        end if

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
          write( u_diag, form ) &
            keep%lmd(i), ' |', word, ' |', keep%info%residual_norms(i), '  |', &
            keep%info%err_lambda(block_size + i), '  |', &
            keep%info%err_X(block_size + i)
        end do ! i = 1, block_size

        write( u_diag, '(a)' ) trim(line) 

      end if ! print_level > 2

      if ( keep%left == 0 ) then

        keep%left_converged = .true.

      else if ( keep%lcon >= keep%max_left ) then

        info%left = keep%max_left
        info%next_left = huge(ONE)
        keep%left_converged = .true.

      else if ( gap == ZERO ) then

        keep%left_converged = keep%lcon >= nep
        if ( keep%left_converged ) then
          info%left = keep%lcon
          if ( keep%first <= block_size ) then
            info%next_left = keep%lmd(keep%first)
          else
            info%next_left = huge(ONE)
          end if
        end if

      else
          
        if ( gap > 0 ) then
          t = gap
          q = ZERO
        else
          t = abs(gap) * keep%av_dist
          q = GAP_REL_A
        end if
        
        keep%left_converged = .false.
        do i = keep%lcon, nep + 1, -1
!!        do i = nep + 1, keep%lcon
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
        
      if ( l == 0 ) then
        keep%av_dist = ZERO
      else
        if ( keep%lcon > 0 ) then
          t = lambda(1)
        else
          t = keep%lmd(1)
        end if
        keep%av_dist = (s - t)/l
      end if
      
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

        converged = converged .and. abs(keep%info%err_X(block_size + i)) < 0.01

        if ( converged ) keep%info%converged(i) = 1
        
      end do ! i = 1, block_size

      if ( keep%info%iteration > 0 ) info%iteration = info%iteration + 1
        
    end select

    if ( rci%job < 0 ) then

      if ( rci%job /= SSMFE_DONE ) then
        info%data = max(0, nep - keep%lcon)
        info%left = keep%lcon
      else
        info%data = 0
      end if
      info%flag = keep%info%flag

      if ( info%flag /= 0 ) &
        call ssmfe_msg &
          ( problem, options, keep%left, 0, block_size, info%flag )

    end if
      
  end subroutine ssmfe_direct_rci_double_complex

  subroutine ssmfe_msg( problem, options, left, right, m, flag )
  
    implicit none
  
    type(ssmfe_options) :: options
    integer :: problem, left, right, m, flag
    intent(in) :: options, left, right, m, flag
    
    integer, parameter :: NONE = -1

    logical :: minAprod, minBprod
    integer :: print_lev, max_it, err_est
    integer :: u_errr, u_warn, u_diag
    real(kind = PRECISION) :: abs_tol, rel_tol, tol, abs_res, rel_res

    print_lev = options%print_level
    u_errr    = options%unit_error
    u_warn    = options%unit_warning
    u_diag    = options%unit_diagnostic
    max_it    = options%max_iterations
    err_est   = options%err_est
    abs_tol   = options%abs_tol_lambda
    rel_tol   = options%rel_tol_lambda
    abs_res   = options%abs_tol_residual
    rel_res   = options%rel_tol_residual
    tol       = options%tol_X
    minAprod  = options%minAprod
    minBprod  = options%minBprod

    if ( print_lev <= NONE ) then
      u_errr = NONE
      u_warn = NONE
      u_diag = NONE
    end if
    if ( u_diag <= NONE .and. print_lev > NONE ) print_lev = 0

    select case ( flag )
    
    case (0)
    
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
        if ( left < 0 ) &
          write ( u_diag, '(a, i4)' ) &
            'extreme eigenpairs requested:', right
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

    case (WRONG_BLOCK_SIZE)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong block size'

    case (WRONG_RCI_JOB)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong rci%job'

    case (WRONG_ERR_EST)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong err_est'

    case (WRONG_EXTRAS)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong extra_left or extra_right'

    case (WRONG_MINPROD)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a,a/)' ) &
          '??? Error: minAprod and minBprod must be true'

    case (WRONG_MIN_GAP)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a,a/)' ) &
          '??? Wrong value of min_gap'

    case (WRONG_CF_MAX)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a,a/)' ) &
          '??? Wrong value of cf_max'

    case (B_NOT_POSITIVE_DEFINITE)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Wrong B or linear depended initial vectors'

    case (OUT_OF_MEMORY)

      if ( u_errr > NONE ) &
        write( u_errr, '(/a/)' ) &
          '??? Out of memory'

    case (NO_SEARCH_DIRECTIONS_LEFT)
    
      if ( u_warn > NONE ) &
        write( u_warn, '(/a,a/)' ) &
          '??? WARNING: iterations terminated because no further progress ', &
          'is possible'

    end select

  end subroutine ssmfe_msg

end module SPRAL_ssmfe_expert

