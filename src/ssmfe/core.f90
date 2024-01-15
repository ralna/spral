! COPYRIGHT (c) 2013      Evgueni Ovtchinnikov
! COPYRIGHT (c) 2014,2015 The Science and Technology Facilities Council (STFC)
! Version 1.0.0
!
! Written by: Evgueni Ovtchinnikov
!
module spral_ssmfe_core

  implicit none

  private

  ! double precision to be used
  integer, parameter :: PRECISION = kind(1.0D0)

  ! BLAS/LAPACK flags
  character, parameter :: SSMFE_JOBZ = 'V'
  character, parameter :: SSMFE_UPLO = 'U'
  integer, parameter :: SSMFE_ITYPE = 1

  ! input error flags
  integer, parameter :: WRONG_BLOCK_SIZE = -2
  integer, parameter :: WRONG_ERR_EST    = -3
  integer, parameter :: WRONG_MINPROD    = -4
  integer, parameter :: WRONG_EXTRAS     = -5
  integer, parameter :: WRONG_MIN_GAP    = -6
  integer, parameter :: WRONG_CF_MAX     = -7

  ! execution error flags
  integer, parameter :: OUT_OF_MEMORY    = -100
  integer, parameter :: INDEFINITE_B     = -200

  ! warning flag
  integer, parameter :: NO_SEARCH_DIRECTIONS_LEFT = 1

  ! error estimation schemes
  integer, parameter :: SSMFE_RESIDUAL  = 1
  integer, parameter :: SSMFE_KINEMATIC = 2

  interface ssmfe
    module procedure ssmfe_double, ssmfe_double_complex
  end interface

  interface ssmfe_largest
    module procedure ssmfe_largest_double, ssmfe_largest_double_complex
  end interface

  interface ssmfe_free
    module procedure &
      ssmfe_core_free_double, &
      ssmfe_core_free_keep_double, &
      ssmfe_free_info_double
  end interface

  type ssmfe_core_options ! see the spec

    integer :: extra_left  = 0
    integer :: extra_right = 0
    integer :: err_est = SSMFE_KINEMATIC
    logical :: minAprod = .true.
    logical :: minBprod = .true.
    real(PRECISION) :: min_gap = 0.0
    real(PRECISION) :: cf_max = 1.0

  end type ssmfe_core_options

  type ssmfe_rcid ! see the spec

    integer :: job = 0
    integer :: nx = 0
    integer :: jx = 0
    integer :: kx = 0
    integer :: ny = 0
    integer :: jy = 0
    integer :: ky = 0
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0

    real(PRECISION) :: alpha, beta

    real(PRECISION), dimension(:,:), pointer :: x => null()
    real(PRECISION), dimension(:,:), pointer :: y => null()

  end type ssmfe_rcid

  type ssmfe_rciz ! see the spec

    integer :: job = 0
    integer :: nx = 0
    integer :: jx = 0
    integer :: kx = 0
    integer :: ny = 0
    integer :: jy = 0
    integer :: ky = 0
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0

    complex(PRECISION) :: alpha = 0
    complex(PRECISION) :: beta = 0

    complex(PRECISION), dimension(:,:), pointer :: x => null()
    complex(PRECISION), dimension(:,:), pointer :: y => null()

  end type ssmfe_rciz

  type ssmfe_inform ! see the spec

    integer :: flag          = 0
    integer :: stat          = 0
    integer :: non_converged = 0
    integer :: iteration     = 0
    integer :: left = 0
    integer :: right = 0

    integer, dimension(:), allocatable :: converged

    real(PRECISION) :: next_left  = 1
    real(PRECISION) :: next_right = -1

    real(PRECISION), dimension(:), allocatable :: residual_norms
    real(PRECISION), dimension(:), allocatable :: err_lambda
    real(PRECISION), dimension(:), allocatable :: err_x

  end type ssmfe_inform

  type ssmfe_core_keep

    private

    ! problem type
    !  0: A x = lambda x
    ! >0: A x = lambda B x
    ! <0: A B x = lambda x
    integer :: problem = 0

    integer :: step = 0 ! computational step
    integer :: stage = 0 ! computational stage
    integer :: iteration = 0 ! iteration number

    integer :: sizeX = 0 ! no. currently iterated vectors
    integer :: sizeY = 0 ! no. current search directions
    integer :: sizeZ = 0 ! no. previous search directions
    integer :: sizeXn = 0 ! updated sizeX

    integer :: firstX  = 1 ! first iterated vector position in workspace block 0
    integer :: firstXn = 1 ! updated firstXn
    integer :: leftX  = 0 ! no. currently iterated left vectors
    integer :: rightX = 0 ! no. currently iterated right vectors
    integer :: leftXn  = 0 ! updated leftX
    integer :: rightXn = 0 ! updated rightX

    integer :: left_cnv  = 0 ! no. converged left eigenpairs
    integer :: right_cnv = 0 ! no. converged right eigenpairs
    integer :: new_left  = 0 ! no. newly converged left eigenpairs
    integer :: new_right = 0 ! no. newly converged right eigenpairs

    ! indices of blocks of the workspace array
    ! zero block holds eigenvector iterates X
    integer :: kY = 1 ! current search directions Y
    integer :: kZ = 2 ! previous search directions Z
    integer :: kW = 3 ! auxiliary block W
    integer :: kAX  = 4 ! A*X
    integer :: kAYZ = 5 ! A*Y or A*Z
    integer :: kBX  = 6 ! B*X
    integer :: kBYZ = 7! B*Y or B*Z

    integer :: lwork = 0 ! LAPACK dsyev/zheev workspace size

    integer :: rec = 0 ! record number in eigenvalue decrements history
    integer :: RECORDS = 30 ! max number of records

    ! estimate for the condition number of [X Y]'*B*[X, Y]
    real(PRECISION) :: cond = 1
    ! estimate for the error of multiplying by A
    real(PRECISION) :: err_A = 0
    ! estimate for the error of multiplying by B
    real(PRECISION) :: err_B = 0

    ! eigenvalues decrements history
    real(PRECISION), dimension(:,:), allocatable :: dlmd
    ! estimates for eigenvalues convergence factors
    real(PRECISION), dimension(:), allocatable :: q
    ! approximate eigenvectors shifts after an iteration
    real(PRECISION), dimension(:), allocatable :: dX
    ! eigenvalues array
    real(PRECISION), dimension(:), allocatable :: lambda
    ! real work array for LAPACK dsyev/zheev
    real(PRECISION), dimension(:), allocatable :: dwork
    ! complex work array for LAPACK zheev
    complex(PRECISION), dimension(:), allocatable :: zwork

    ! search directions reordering index
    integer, dimension(:), allocatable :: ind

    integer :: err_est = 2 ! error estimation scheme
    logical :: minAprod = .true. ! A-multiplications minimization flag
    logical :: minBprod = .true. ! B-multiplications minimization flag

  end type ssmfe_core_keep

  public :: ssmfe_core_options, ssmfe_rcid, ssmfe_rciz, ssmfe_inform, &
            ssmfe_core_keep
  public :: ssmfe, ssmfe_largest, ssmfe_free

contains

!**************************************************************************
! real RCI for computing specified numbers of eigenvalues on both sides of
! the spectrum and corresponding eigenvectors,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_double &
    ( rci, problem, left, right, m, lambda, rr_matrices, ind, &
      keep, options, info )
    integer, intent(in) :: problem
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: m
    real(PRECISION), dimension(m), intent(inout) :: lambda
    real(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rcid        ), intent(inout) :: rci
    type(ssmfe_core_keep   ), intent(inout) :: keep
    type(ssmfe_core_options), intent(in   ) :: options
    type(ssmfe_inform   ), intent(inout) :: info

    call ssmfe_engine_double &
      ( problem, left, right, m, lambda, rr_matrices, ind, rci, keep, options, &
        info )

  end subroutine ssmfe_double

!**************************************************************************
! complex RCI for computing specified numbers of eigenvalues on both sides
! of the spectrum and corresponding eigenvectors,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_double_complex &
    ( rci, problem, left, right, m, lambda, rr_matrices, ind, &
      keep, options, info )
    integer, intent(in) :: problem
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: m
    real(PRECISION), dimension(m), intent(inout) :: lambda
    complex(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rciz        ), intent(inout) :: rci
    type(ssmfe_core_keep   ), intent(inout) :: keep
    type(ssmfe_core_options), intent(in   ) :: options
    type(ssmfe_inform      ), intent(inout) :: info

    call ssmfe_engine_double_complex &
      ( problem, left, right, m, lambda, rr_matrices, ind, rci, keep, options, &
        info )

  end subroutine ssmfe_double_complex

!**************************************************************************
! real RCI for computing a specified number of largest eigenvalues and
! corresponding eigenvectors,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_largest_double &
    ( rci, problem, nep, m, lambda, rr_matrices, ind, keep, control, info )
    integer, intent(in) :: problem
    integer, intent(in) :: nep
    integer, intent(in) :: m
    real(PRECISION), dimension(m), intent(inout) :: lambda
    real(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rcid        ), intent(inout) :: rci
    type(ssmfe_core_keep   ), intent(inout) :: keep
    type(ssmfe_core_options), intent(in   ) :: control
    type(ssmfe_inform      ), intent(inout) :: info

    call ssmfe_engine_double &
      ( problem, -1, nep, m, lambda, rr_matrices, ind, rci, keep, control, &
        info )

  end subroutine ssmfe_largest_double

!**************************************************************************
! complex RCI for computing a specified number of largest eigenvalues and
! corresponding eigenvectors,
! see ssmfe spec for the arguments
!
  subroutine ssmfe_largest_double_complex &
    ( rci, problem, nep, m, lambda, rr_matrices, ind, keep, control, info )
    integer, intent(in) :: problem
    integer, intent(in) :: nep
    integer, intent(in) :: m
    real(PRECISION), dimension(m), intent(inout) :: lambda
    complex(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rciz        ), intent(inout) :: rci
    type(ssmfe_core_keep   ), intent(inout) :: keep
    type(ssmfe_core_options), intent(in   ) :: control
    type(ssmfe_inform      ), intent(inout) :: info

    call ssmfe_engine_double_complex &
      ( problem, -1, nep, m, lambda, rr_matrices, ind, rci, keep, control, &
        info )

  end subroutine ssmfe_largest_double_complex

!**************************************************************************
! real core solver RCI
!
  subroutine ssmfe_engine_double &
      ( problem, left, right, m, lambda, rr_matrices, ind, rci, keep, control, &
        info )
    use spral_blas_iface, &
      copy => dcopy, &
      norm => dnrm2, &
      dot  => ddot, &
      axpy => daxpy, &
      gemm => dgemm, &
      trsm => dtrsm
    use spral_lapack_iface, chol => dpotrf
    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0
    real(PRECISION), parameter :: NIL = ZERO
    real(PRECISION), parameter :: UNIT = ONE
    character, parameter :: TRANS = 'T'

    ! arguments

    ! problem type
    !   0: A x = lambda x
    ! > 0: A x = lambda B x
    ! < 0: A B x = lambda x
    ! A = A', B = B' > 0
    integer, intent(in) :: problem

    ! >= 0: number of wanted eigenvalues left of sigma
    !  < 0: indicates that largest eigenvalues are wanted (see <right>)
    integer, intent(in) :: left

    ! left >= 0: number of wanted eigenvalues right of sigma
    ! left  < 0: number of wanted largest eigenvalues
    integer, intent(in) :: right

    ! BJCG block size
    integer, intent(in) :: m

    ! eigenvalues array
    real(PRECISION), dimension(m), intent(inout) :: lambda

    ! matrices for the Rayleigh-Ritz procedure
    real(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, *)

    ! work matrix column reordering index
    integer, intent(inout), dimension(m) :: ind

    ! ssmfe_core types - see the spec
    type(ssmfe_core_keep   ), intent(inout) :: keep
    type(ssmfe_rcid        ), intent(inout) :: rci
    type(ssmfe_core_options), intent(in   ) :: control
    type(ssmfe_inform      ), intent(inout) :: info

    ! rci jobs
    integer, parameter :: SSMFE_START              = 0
    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_CHECK_CONVERGENCE  = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_COPY_VECTORS       = 11
    integer, parameter :: SSMFE_COMPUTE_DOTS       = 12
    integer, parameter :: SSMFE_SCALE_VECTORS      = 13
    integer, parameter :: SSMFE_COMPUTE_YMXD       = 14
    integer, parameter :: SSMFE_COMPUTE_XY         = 15
    integer, parameter :: SSMFE_COMPUTE_XQ         = 16
    integer, parameter :: SSMFE_TRANSFORM_X        = 17
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999
    integer, parameter :: SSMFE_FINISH             = -1
    integer, parameter :: SSMFE_STOP               = -2
    integer, parameter :: SSMFE_ABORT              = -3

    ! computation stages
    integer, parameter :: INITIAL = 10
    integer, parameter :: CG_LOOP = 20

    ! computation steps
    integer, parameter :: BEGIN = 0
    integer, parameter :: QUIT  = -1

    integer, parameter :: COMPUTE_BX           = 100
    integer, parameter :: ORTHOG_X_TO_Xc       = 150
    integer, parameter :: COMPUTE_AX           = 200
    integer, parameter :: COMPUTE_INIT_XAX     = 210
    integer, parameter :: COMPUTE_INIT_XBX     = 220
    integer, parameter :: RR_IN_X              = 300
    integer, parameter :: TRANSFORM_X          = 305
    integer, parameter :: COMPUTE_INIT_AX      = 310
    integer, parameter :: COMPUTE_INIT_BX      = 320
    integer, parameter :: COMPUTE_XBX          = 330
    integer, parameter :: COPY_AX_TO_W         = 350
    integer, parameter :: COMPUTE_XAX          = 360
    integer, parameter :: COMPUTE_RESIDUALS    = 400
    integer, parameter :: ORTHOG_R_TO_Xc       = 410
    integer, parameter :: COMPUTE_RR           = 420
    integer, parameter :: COPY_W_TO_Y          = 430
    integer, parameter :: APPLY_B_TO_Y         = 440
    integer, parameter :: ORTHOG_Y_TO_Xc       = 450
    integer, parameter :: ESTIMATE_ERRORS      = 500
    integer, parameter :: CHECK_CONVERGENCE    = 600
    integer, parameter :: SAVE_RIGHT_CONVERGED = 700
    integer, parameter :: APPLY_PRECONDITIONER = 800
    integer, parameter :: COMPUTE_AZ           = 900
    integer, parameter :: COMPUTE_ZAY          = 910
    integer, parameter :: COMPUTE_BY           = 920
    integer, parameter :: COMPUTE_ZBY          = 930
    integer, parameter :: COMPUTE_YBY_DIAG     = 950
    integer, parameter :: CONJUGATE_Y          = 1000
    integer, parameter :: RECOMPUTE_BY         = 1100
    integer, parameter :: UPDATE_BY            = 1110
    integer, parameter :: COPY_W_TO_BYZ        = 1120
    integer, parameter :: APPLY_CONSTRAINTS    = 1200
    integer, parameter :: COMPUTE_YBY          = 1300
    integer, parameter :: SCALE_Y              = 1320
    integer, parameter :: COMPUTE_XBY          = 1400
    integer, parameter :: CLEANUP_Y            = 1500
    integer, parameter :: COMPUTE_AY           = 1600
    integer, parameter :: RR_IN_XY             = 2000
    integer, parameter :: COMPUTE_XAY          = 2100
    integer, parameter :: PREPARE_MATRICES     = 2200
    integer, parameter :: ANALYSE_RR_RESULTS   = 2300
    integer, parameter :: PUT_AZQ_IN_Z         = 3100
    integer, parameter :: ADD_AXQ_TO_Z         = 3200
    integer, parameter :: PUT_AXQ_IN_W         = 3300
    integer, parameter :: ADD_AZQ_TO_W         = 3400
    integer, parameter :: PUT_W_IN_AX          = 3500
    integer, parameter :: PUT_Z_IN_AZ          = 3600
    integer, parameter :: PUT_BZQ_IN_Z         = 4100
    integer, parameter :: ADD_BXQ_TO_Z         = 4200
    integer, parameter :: PUT_BXQ_IN_W         = 4300
    integer, parameter :: ADD_BZQ_TO_W         = 4400
    integer, parameter :: PUT_W_IN_BX          = 4500
    integer, parameter :: PUT_Z_IN_BZ          = 4600
    integer, parameter :: PUT_YQ_IN_Z          = 5100
    integer, parameter :: ADD_XQ_TO_Z          = 5200
    integer, parameter :: PUT_XQ_IN_W          = 5300
    integer, parameter :: ADD_YQ_TO_W          = 5400
    integer, parameter :: PUT_W_IN_X           = 5500
    integer, parameter :: CHECK_THE_GAP        = 6000

    ! round-off error level
    real(PRECISION) :: TOO_SMALL

    ! heuristic constants
    real(PRECISION) :: A_SMALL_FRACTION, A_MULTIPLE
    real(PRECISION) :: REL_GAP, MAX_BETA, Q_MAX

    ! maximal admissible condition number of the B-gram matrix
    real(PRECISION) :: GRAM_RCOND_MIN

    ! flag for 'not available yet'
    real(PRECISION), parameter :: NO_VALUE = -1.0

    ! locals
    logical :: skip, doit
    integer :: kY, kZ, kW, kAX, kAYZ, kBX, kBYZ
    integer :: mm, left_cnv, right_cnv, first, last, step, go, nsmall
    integer :: iX, jX, iY, jY, sizeXY
    integer :: i, j, k, l
    real(PRECISION) :: delta, theta
    real(PRECISION) :: q, r, s, t

    ! shortcuts for the options
    logical :: minAprod, minBprod
    integer :: err_est

    ! check for the input data errors

    if ( control%extra_left < 0 .or. control%extra_right < 0 ) then
      info%flag = WRONG_EXTRAS
      rci%job = SSMFE_ABORT
      return
    end if

    if ( control%min_gap < ZERO .or. control%min_gap > ONE ) then
      info%flag = WRONG_MIN_GAP
      rci%job = SSMFE_ABORT
      return
    end if

    if ( control%cf_max < ONE/2 .or. control%cf_max > ONE ) then
      info%flag = WRONG_CF_MAX
      rci%job = SSMFE_ABORT
      return
    end if

    if ( m < 1 &
      .or. (left < 0 .or. left > 0 .and. right > 0) .and. m < 2 ) then
      info%flag = WRONG_BLOCK_SIZE
      rci%job = SSMFE_ABORT
      return
    end if

    err_est = control%err_est
    select case ( err_est )
    case ( SSMFE_RESIDUAL, SSMFE_KINEMATIC )
    case default
      info%flag = WRONG_ERR_EST
      rci%job = SSMFE_ABORT
      return
    end select

    if ( problem < 0 ) then
      if ( .not. control%minAprod ) then
        info%flag = WRONG_MINPROD
        rci%job   = SSMFE_ABORT
        return
      end if
    end if

    A_SMALL_FRACTION = 0.01D0
    A_MULTIPLE = 100
    GRAM_RCOND_MIN = 1D-4
    MAX_BETA = 100
    TOO_SMALL = A_MULTIPLE*epsilon(ONE)

    REL_GAP = control%min_gap
    q_max = control%cf_max

    mm = 2*m

if_rci: &
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
      .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) &
      ) then
      ! things to do on (re)start

      ! reallocate the inform and keep arrays
      if ( allocated(keep%ind           ) ) deallocate ( keep%ind            )
      if ( allocated(keep%lambda        ) ) deallocate ( keep%lambda         )
      if ( allocated(keep%q             ) ) deallocate ( keep%q              )
      if ( allocated(keep%dX            ) ) deallocate ( keep%dX             )
      if ( allocated(keep%dlmd          ) ) deallocate ( keep%dlmd           )
      if ( allocated(keep%dwork         ) ) deallocate ( keep%dwork          )
      if ( allocated(info%err_lambda    ) ) deallocate ( info%err_lambda     )
      if ( allocated(info%err_X         ) ) deallocate ( info%err_X          )
      if ( allocated(info%converged     ) ) deallocate ( info%converged      )
      if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )

      keep%lwork = lwork_sevp( mm )

      allocate &
        ( info%residual_norms(m), info%err_lambda(mm), info%err_X(mm), &
          info%converged(m), keep%lambda(3*m), &
          keep%ind(2*m), keep%q(mm), keep%dX(mm), &
          keep%dwork(keep%lwork), &
          keep%dlmd(m, keep%RECORDS), stat = info%stat )
      if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
      if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( info%stat /= 0 ) return

      ! initialize the execution information
      info%flag       = 0
      info%iteration  = 0
      info%residual_norms = NO_VALUE
      info%err_lambda = NO_VALUE
      info%err_X      = NO_VALUE
      info%converged  = 0

      ! save the input data in keep
      keep%problem  = problem
      keep%err_est  = control%err_est
      keep%minAprod = control%minAprod
      keep%minBprod = control%minBprod .and. keep%problem /= 0

      ! mark workspace blocks
      ! zero block stores approximate eigenvectors   (X)
      keep%kY = 1 ! new search directions block      (Y)
      keep%kZ = 2 ! previous search directions block (Z)
      keep%kW = 3 ! auxiliary block                  (W)
      if ( keep%minAprod ) then
        keep%kAX  = keep%kW + 1  ! block for A*X
        keep%kAYZ = keep%kAX + 1 ! block for A*Y and A*Z
      else
        ! no permanent storage for A*X, A*Y and A*Z, store in auxiliary block
        keep%kAX  = keep%kW
        keep%kAYZ = keep%kW
      end if
      if ( keep%minBprod .or. keep%problem < 0 ) then
        keep%kBX  = keep%kAYZ + 1 ! block for B*X
        keep%kBYZ = keep%kBX + 1  ! block for B*Y and B*Z
      else
        ! no permanent storage for A*X, A*Y and A*Z, store in other blocks
        if ( keep%problem == 0 ) then
          keep%kBX = 0
        else
          keep%kBX = keep%kY
        end if
        keep%kBYZ = keep%kY
      end if

      ! divide the workspace into parts for storing eigenvectors
      ! corresponding to left and right margin of the spectrum
      ! (referred below to as left and right eigenvectors)

      ! determine keep%leftX, the size of the part for left eigenvectors
      if ( left >= 0 ) then
        ! divide proportionally to <left> and <right>
        if ( left == 0 ) then
          keep%leftX = 0
        else if ( right == 0 ) then
          keep%leftX = m
        else
          i = left + control%extra_left
          j = right + control%extra_right
          if ( i + j <= m ) then
            keep%leftX = i + int((m - i - j)*i*ONE/(i + j))
          else
            keep%leftX = int(m*i*ONE/(i + j))
          end if
          keep%leftX = max(keep%leftX, 1)
          keep%leftX = min(keep%leftX, m - 1)
        end if
      else
        ! divide evenly
        keep%leftX = m/2
      end if
      ! the size of the part for the right eigenvectors
      keep%rightX = m - keep%leftX

      keep%firstX = 1 ! first non-converged eigenpair index
      keep%sizeX  = m ! # of non-converged eigenpairs
      keep%sizeY  = 0 ! # of current search directions
      keep%sizeZ  = 0 ! # of previous search directions

      keep%firstXn   = 1            ! keep%firstX after update
      keep%leftXn    = keep%leftX   ! keep%leftX after update
      keep%rightXn   = keep%rightX  ! keep%rightX after update
      keep%sizeXn    = m            ! keep%sizeX after update
      keep%new_left  = 0            ! keep%leftX update increment
      keep%new_right = 0            ! keep%rightX update increment
      if ( rci%job == SSMFE_START ) then
        keep%left_cnv  = 0 ! number of converged left eigenvectors
        keep%right_cnv = 0 ! number of converged right eigenvectors
      end if

      keep%cond  = ONE  ! estimate for the condition number of X'*B*X
      keep%err_A = ZERO ! estimated error of multiplication by A
      keep%err_B = ZERO ! estimated error of multiplication by B

      keep%rec  = 0    ! record number in eigenvalue decrement history
      keep%dlmd = ZERO ! eigenvalue decrements history
      keep%q    = ONE  ! eigenvalue convergence factors
      keep%dX   = ONE  ! eigenvector shifts

      ! initial stage of computation consists in the orthogonalization
      ! of the initial approximate eigenvectors X to constraints
      ! and the Rayleigh-Ritz procedure in the trial subspace spanned
      ! by the orthogonalized initial vectors
      keep%stage = INITIAL

      ! the first step of the initial stage
      keep%step = BEGIN

      keep%iteration = 0

    end if if_rci

    ! short names for elements of keep

    minAprod = keep%minAprod
    minBprod = keep%minBprod
    err_est  = keep%err_est

    ! indices of blocks of the workspace array
    kY   = keep%kY    ! current search directions Y
    kZ   = keep%kZ    ! previous search directions Z
    kW   = keep%kW    ! general purpose workspace
    kAX  = keep%kAX   ! A*X
    kAYZ = keep%kAYZ  ! A*Y or A*Z
    kBX  = keep%kBX   ! B*X
    kBYZ = keep%kBYZ  ! B*Y or B*Z

do_select_step: &
    do ! This do loop immediately surrounds the select case construct and
       ! allows any of the cases to reset keep%step and cause the corresponging
       ! case to execute next.

      if ( left <= 0 .and. right <= 0 ) keep%step = QUIT

select_step: &
      select case ( keep%step )

      case (BEGIN) select_step

        ! computation starts with the orthogonalization of the initial
        ! approximate eigenvectors X to constraints, which, if B is not
        ! identity, requires multiplying these vectors by B

        keep%step = COMPUTE_BX

      case (COMPUTE_BX) select_step

        ! select next step
        if ( keep%stage == CG_LOOP ) then
          ! compute A*X to form residual vectors
          keep%step = COMPUTE_AX
        else
          ! orthogonalize X to constraints
          keep%step = ORTHOG_X_TO_Xc
        end if

        ! check whether the multiplication by B is actually needed:
        ! it is skipped if B*X is computed implicitly or B is identity,
        ! but is enforced every 20 iterations to eliminate posiible
        ! accumulation of errors caused by implicit computation of B*X
        if ( .not. (minBprod .and. keep%stage == CG_LOOP &
                    .and. mod(keep%iteration, 20) > 0 ) &
              .and. problem /= 0 ) then
          ! multiplication needed
          rci%job = SSMFE_APPLY_B ! instruct the user to apply B
          rci%nx = keep%sizeX     ! number of vectors to be multiplied
          rci%kx = 0              ! block index
          rci%jx = keep%firstX    ! position of the first vector in the block
          rci%ky = kBX            ! product receiving block index
          rci%jy = keep%firstX    ! pos. of the first product in receiving block
          return
        end if

      case (ORTHOG_X_TO_Xc) select_step

        rci%job = SSMFE_APPLY_CONSTRAINTS
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%jy = keep%firstX
        rci%ky = kBX

        ! next step is to compute A*X or A*B*X
        keep%step = COMPUTE_AX
        return

      case (COMPUTE_AX) select_step

        ! select next step
        if ( keep%stage /= INITIAL ) then
          keep%step = COMPUTE_XBX
        else
          keep%step = COMPUTE_INIT_XAX
        end if

        ! back to the user for A X
        if ( .not. (minAprod .and. keep%stage == CG_LOOP &
                    .and. mod(keep%iteration, 20) > 0) ) then
          rci%job = SSMFE_APPLY_A
          rci%nx = keep%sizeX
          rci%jx = keep%firstX

          if ( problem < 0 ) then
            ! in the case of the problem A B x = lambda x
            ! the Rayleigh-Ritx procedure requires the computation of
            ! X'*B*A*B*X instead of X'*A*X, hence compute A*(B*X)
            rci%kx = kBX
          else
            ! compute A*X
            rci%kx = 0
          end if

          rci%jy = keep%firstX
          rci%ky = kAX
          return
        end if

      case (COMPUTE_INIT_XAX) select_step

        ! compute H = X'*(A*X) or H = (B*X)'*(A*(B*X))

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX

        if ( problem < 0 ) then
          ! multiply (B*X)' and (A*(B*X))
          rci%kx = kBX
        else
          ! multiply X' and A*X
          rci%kx = 0
        end if

        rci%jy = keep%firstX
        rci%ky = kAX
        rci%ny = keep%sizeX
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_INIT_XBX
        return

      case (COMPUTE_INIT_XBX) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%jy = keep%firstX
        if ( problem == 0 ) then
          ! compute G = X'*X
          rci%ky = 0
        else
          ! compute G = X'*(B*X)
          rci%ky = kBX
        end if
        rci%ny = keep%sizeX
        rci%i = 1
        rci%j = 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = RR_IN_X
        return

      case (RR_IN_X) select_step ! Rayleigh-Ritz procedure in span(X)

        ! solve the eigenvalue problem H q = lambda G q
        call solve_gsevp &
          ( keep%sizeX, rr_matrices(1,1,2), mm, rr_matrices, mm, &
            lambda, keep%lwork, keep%dwork, i )
        if ( i /= 0 ) then
          info%flag = INDEFINITE_B
          rci%job = SSMFE_ABORT
          return
        end if

        if ( left < 0 ) then
          ! largest eigenvalues are computed:
          ! redistribute work space for left and right eigenvectors
          ! (in this case, left eigenvalues are negative and right
          ! eigenvalues positive)
          j = 1
          do while ( j < m )
            if ( lambda(j) > 0 ) exit
            j = j + 1
          end do
          keep%leftX = j - 1
          keep%rightX = m - keep%leftX
        end if

        keep%lambda(1 : keep%sizeX) = lambda(1 : keep%sizeX)

        keep%step = TRANSFORM_X

      case (TRANSFORM_X) select_step

        ! update X := X*Q, where the columns of Q are the eigenvectors of
        ! H q = lambda G q
        rci%job = SSMFE_TRANSFORM_X
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        rci%ky = kZ
        rci%i = 1
        rci%j = 1
        rci%k = 2
        keep%step = COMPUTE_INIT_AX
        return

      case (COMPUTE_INIT_AX) select_step

        ! update A*X := (A*X)*Q as a cheaper alternative to A*X := A*(X*Q)
        ! (BLAS level 3 can be used)

        rci%job = SSMFE_TRANSFORM_X
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kAX
        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        rci%ky = kZ
        rci%i = 1
        rci%j = 1
        rci%k = 2
        keep%step = COMPUTE_INIT_BX
        return

      case (COMPUTE_INIT_BX) select_step

        ! this is the last step of the initial stage
        if ( keep%stage == INITIAL ) keep%stage = CG_LOOP
        keep%step = COMPUTE_XBX

        if ( problem /= 0 ) then
          ! update B*X := (B*X)*Q as a cheaper alternative to B*X := B*(X*Q)
          rci%job = SSMFE_TRANSFORM_X
          rci%nx = keep%sizeX
          rci%jx = keep%firstX
          rci%kx = kBX
          rci%ny = keep%sizeX
          rci%jy = keep%firstX
          rci%ky = kZ
          rci%i = 1
          rci%j = 1
          rci%k = 2
          return
        end if

      case (COMPUTE_XBX) select_step

        ! compute X'*B*X to be used later in the Rayleigh-Ritz procedure
        ! in the trial subspace spanned by X and new search directions Y

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        if ( problem == 0 ) then
          ! compute X'*X
          rci%ky = 0
        else
          ! compute X'*(B*X)
          rci%ky = kBX
        end if
        rci%i = 1
        rci%j = 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL

        ! next, start computing residual vectors
        keep%step = COPY_AX_TO_W
        return

      case (COPY_AX_TO_W) select_step

        keep%step = COMPUTE_XAX
        if ( minAprod ) then
          rci%job = SSMFE_COPY_VECTORS
          rci%nx = keep%sizeX
          rci%jx = keep%firstX
          rci%kx = kAX
          rci%jy = keep%firstX
          rci%ky = kW
          rci%i = 0
          return
        end if

      case (COMPUTE_XAX) select_step

        ! compute X'*(A*X) or (B*X)'*A*(B*X)

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0

        if ( problem < 0 ) then
          ! compute (B*X)'*A*(B*X)
          rci%kx = kBX
        else
          ! compute X'*(A*X)
          rci%kx = 0
        end if

        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_RESIDUALS
        return

      case (COMPUTE_RESIDUALS) select_step

        do i = 1, keep%sizeX
          ! X'*(B*X) stored in rr(:,:,1) must be identity,
          ! the deviation from identity is caused mainly by the
          ! error of the Rayleigh-Ritz procedure and the error
          ! in computing B*X, and hence can be used to estimate
          ! these errors
          keep%err_B = max(keep%err_B, abs(rr_matrices(i, i, 1) - ONE))
          do j = 1, i - 1
            keep%err_B = max(keep%err_B, abs(rr_matrices(j, i, 1)))
          end do
          iX = keep%firstX + i - 1
          s = rr_matrices(i, i, 1)
          t = rr_matrices(i, i, 2)

          ! the Rayleigh quotient for the approximate eigenvector X(:,iX)
          ! is to be used as the corresponding approximate eigenvalue
          t = t/s

          if ( keep%rec > 0 ) then
            ! try to compute the eigenvalue decrement
            if ( i > keep%leftXn - keep%new_left &
              .and. i <= keep%leftXn + keep%new_right &
              .or. left < 0 .and. t*lambda(iX) < ZERO &
              ) then
              ! previous approximate eigenvalue not available,
              ! decrement cannot be computed
              keep%dlmd(iX, keep%rec) = ZERO
            else
              ! directly computed decrement
              s = abs(lambda(iX) - t)
              r = ZERO
              do j = 1, keep%sizeX
                r = r + (rr_matrices(j, i, 2) - t*rr_matrices(j, i, 1))**2
              end do
              ! sqrt(r) is the Frobenius norm of the ressidual H q - t G q
              ! and it gives an idea of the error in t
              r = max(sqrt(r), abs(t)*epsilon(ONE))
              ! directly decrement must be sufficiently larger than the error
              ! in t to be relied upon
              if ( s > A_MULTIPLE*r ) keep%dlmd(iX, keep%rec) = s
              ! if it is small enough, Jacobi estimate is used instead
              ! (see below)
            end if
          end if
          lambda(iX) = t
        end do
        if ( keep%rec > 0 ) then
          ! in the presence of tight eigenvalue clusters, decrements may be
          ! mixed up, so replace decrements for eigenvalues in a cluster
          ! with their maximum over this cluster
          do i = 1, keep%sizeX
            iX = keep%firstX + i - 1
            s = keep%dlmd(iX, keep%rec)
            if ( s == ZERO ) cycle
            do l = 1, m
              if ( abs(lambda(l) - lambda(iX)) &
                < 10*epsilon(ONE)*max(abs(lambda(l)), abs(lambda(iX))) ) &
                s = max(s, keep%dlmd(l, keep%rec))
            end do
            keep%dlmd(iX, keep%rec) = s
          end do
        end if

        ! the number of search directions is normally equal to the number
        ! of currently iterated approximate eigenvectors (both may change later)
        keep%sizeY = keep%sizeX

        ! compute residuals

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          rr_matrices(iX, m + iX, 1) = lambda(iX)
        end do
        rci%job = SSMFE_COMPUTE_YMXD
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        if ( problem <= 0 ) then
          ! compute A x - lambda x or A B x - lambda x
          rci%kx = 0
        else
          ! compute A x - lambda B x
          rci%kx = kBX
        end if
        rci%jy = keep%firstX
        rci%ky = kW
        rci%i = keep%firstX
        rci%j = m + keep%firstX
        rci%k = 1

        ! next step is to orthogonalize residual vectors
        ! to constraints Xc for correct error estimation
        if ( problem < 0 ) then
          ! residuals r = A B x - lambda x are made B-orthogonal to Xc by
          ! update r := r - Xc*Xc'*B*r
          keep%step = COPY_W_TO_Y
        else
          ! residuals r = A x - lambda B x (B may be identity) are made
          ! orthogonal to Xc by update r := r - B*Xc*Xc'*r
          keep%step = ORTHOG_R_TO_Xc
        end if

        return

      case (COPY_W_TO_Y) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kW
        rci%jy = 1
        rci%ky = kY
        rci%i = 0
        keep%step = APPLY_B_TO_Y
        return

      case (APPLY_B_TO_Y) select_step

        rci%job = SSMFE_APPLY_B
        rci%nx = keep%sizeX
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kW
        keep%step = ORTHOG_Y_TO_Xc
        return

      case (ORTHOG_Y_TO_Xc) select_step

        rci%job = SSMFE_APPLY_CONSTRAINTS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 3
        ! next, compute Gram matrix R'*(B*R) for residuals A B x - lambda x
        keep%step = COMPUTE_RR
        return

      case (ORTHOG_R_TO_Xc) select_step

        rci%job = SSMFE_APPLY_ADJ_CONSTRS
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kW
        rci%i = 1
        rci%j = 1
        rci%k = 3
        ! next, compute Gram matrix R'*R for residuals A x - lambda B x
        keep%step = COMPUTE_RR
        return

      case (COMPUTE_RR) select_step

        rci%nx = keep%sizeX
        rci%ny = keep%sizeX

        if ( problem < 0 ) then
          rci%kx = kY
          rci%jx = 1
          rci%jy = 1
        else
          rci%jx = keep%firstX
          rci%kx = kW
          rci%jy = keep%firstX
        end if

        rci%ky = kW
        rci%i = keep%firstX
        rci%j = keep%firstX
        rci%k = 3
        if ( err_est == SSMFE_RESIDUAL ) then
          ! residual error estimation requires Gram matrix for residuals
          rci%alpha = UNIT
          rci%beta = NIL
          rci%job = SSMFE_COMPUTE_XY
        else
          ! otherwise, residual norms suffice
          rci%job = SSMFE_COMPUTE_DOTS
        end if
        keep%step = ESTIMATE_ERRORS
        return

      case (ESTIMATE_ERRORS) select_step

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          info%residual_norms(iX) = sqrt(abs(rr_matrices(iX, iX, 3)))
        end do

        ! next step
        keep%step = CHECK_CONVERGENCE

        ! estimate convergence factors for eigenvalues
        ! (the ratio of the eigenvalue errors before and after
        ! an iteration)

        l = keep%rec ! current record number in eigenvalue decrements history

        ! sufficient number of records needed for estimation
        ! based on the decrement history
        if ( l > 3 ) then

          do iX = keep%firstX, keep%firstX + keep%sizeX - 1

            ! decrement not available, hence neither is the error estimate
            if ( keep%dlmd(iX,l) == ZERO ) then
              keep%q(iX) = ONE
              cycle
            end if

            ! verify that the last eigenvector shift was small
            ! enough compared to the accuracy requirements
            if ( keep%dX(iX) > A_SMALL_FRACTION ) then ! not small enough
              keep%q(iX) = ONE
              cycle
            end if

            ! estimate the asymptotic convergence factor (a.c.f.) for this
            ! eigenvalue, the geometrical average of convergence factors over
            ! iterations performed so far
            s = ZERO
            k = 0
            ! trust region is the last 1/3 of the history
            do j = l, l - l/3, -1
              q = keep%dlmd(iX,j)
              if ( q == ZERO ) exit
              s = s + q
              k = k + 1
            end do
            ! history too short or decrements not available
            if ( k < 2 .or. s == ZERO ) cycle
            q = keep%dlmd(iX,l)/s
            q = q**(ONE/(k - 1))
            keep%q(m + iX) = q ! store a.c.f. approximation
            ! if q is too close to 1, skip refinement
            if ( q > ONE - TOO_SMALL ) cycle

            ! estimate error based on a.c.f.
            theta = q/(ONE - q)
            info%err_lambda(m + iX) = keep%dlmd(iX,l)*theta

            ! use the fact that the deviations of convergence factors
            ! from a.c.f are essentially random to find a better error
            ! bound

            ! find first available decrement to establish the trust region
            i = 0
            do j = l - 1, 1, -1
              if ( keep%dlmd(iX,j) == ZERO ) then
                i = j
                exit
              end if
            end do

            k = (l - i)/3
            if ( k < 1 ) cycle ! useful history too short

            ! set trust region to the middle 1/3 of history
            first = l - 2*k + 1
            last = l - k

            ! compute the average of theta_j = log(q_j/(1 - q_j))
            ! where q_j is the ratio of the error on (j+1)-th iteration
            ! to that on j-th
            theta = ZERO
            do j = last, first, -1
              s = info%err_lambda(m + iX)
              do k = j + 1, l
                s = s + keep%dlmd(iX,k)
              end do
              ! s approximates the eigenvalue error at iteration j + 1
              ! hence s/keep%dlmd(iX,j) = q_j/(1 - q_j)
              theta = theta + log(s) - log(keep%dlmd(iX,j))
            end do
            k = last - first + 1
            theta = theta/k

            ! the range of theta_j is (-Inf, Inf), and
            ! the deviations of theta_j from theta are essentially random

            ! compute the standard deviation for theta_j
            r = ZERO
            do j = last, first, -1
              s = info%err_lambda(m + iX)
              do k = j + 1, l
                s = s + keep%dlmd(iX,k)
              end do
              r = r + (theta - log(s) + log(keep%dlmd(iX,j)))**2
            end do
            r = sqrt(r/(last - first + 1))

            ! assuming theta_j to be normally distributed,
            ! apply 3-sigma rule to obtain the upper bound for theta_j,
            ! and thus for the convergence factor keep%q(iX),
            ! that discards outliers

            theta = exp(theta + 3*r)
            keep%q(iX) = theta/(ONE + theta)

          end do ! iX

        end if ! l > 3

        ! using eigenvalue convergence factors,
        ! estimate eigenvalue errors
        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == ZERO ) then
            ! zero-residual eigenpairs are exact
            info%err_lambda(iX) = ZERO
            cycle
          end if
          l = keep%rec
          q = keep%q(iX)
          info%err_lambda(iX) = NO_VALUE
          if ( q < ONE .and. l > 0 ) then
            if ( keep%dlmd(iX,l) > ZERO ) &
              info%err_lambda(iX) = keep%dlmd(iX,l)*q/(ONE - q)
          end if
        end do

        ! using eigenvalue convergence factors,
        ! estimate eigenvector errors
        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == ZERO ) then
            ! zero-residual eigenpairs are exact
            info%err_X(iX) = ZERO
            cycle
          end if
          ! the eigenvector error is roughly the square root
          ! of the eigenvalue error, hence use the same relationship
          ! for respective convergence factors
          q = sqrt(keep%q(iX))
          if ( q < ONE .and. keep%dX(iX) > ZERO ) then
            info%err_X(iX) = min(ONE, keep%dX(iX)*q/(ONE - q))
          else
            info%err_X(iX) = NO_VALUE
          end if
        end do

        if ( err_est == SSMFE_RESIDUAL ) then

          do iX = keep%firstX, keep%firstX + keep%sizeX - 1
            ! default: Krylov-Weinstein estimate
            info%err_lambda(iX) = info%residual_norms(iX)
          end do

          ! use Lehmann bounds for eigenvalue errors and
          ! Davis-Kahan bounds for eigenvector errors

          do go = 1, 2 ! treat each margin of the spectrum separately

            if ( keep%sizeX < 2 ) exit ! not enough data

            select case ( go )
            case ( 1 )
              if ( keep%leftX == 0 ) cycle
              first = keep%firstX
              last = keep%firstX + keep%leftX - 1
              step = 1
            case ( 2 )
              if ( keep%rightX == 0 ) cycle
              first = keep%firstX + keep%sizeX - 1
              last = keep%firstX + keep%leftX
              step = -1
            end select

            ! choose Lehmann pole
            k = last

            t = ZERO
            do while ( (k - first)*step >= 0 )
              if ( step*lambda(k) < ZERO ) exit
              k = k - step
            end do
            if ( (k - first)*step < 0 ) cycle

            doit = left >= 0 .or. last + step > m .or. last + step < 1
            if ( .not. doit ) doit = step*lambda(last + step) < ZERO

            if ( doit ) then

              s = ZERO
              do iX = first, last - step, step
                s = s + info%residual_norms(iX)**2
              end do

              do while ( (k - first)*step > 0 )
                if ( step > 0 ) then
                  r = lambda(maxloc(lambda(1:k-step),1))
                else
                  r = lambda(k + minloc(lambda(k-step:m),1))
                end if
                if ( abs(lambda(k) - r) > &
                      sqrt(s) + info%residual_norms(k) ) then
                  t = lambda(k) - step*info%residual_norms(k)
                  exit
                end if
                k = k - step
                s = max(ZERO, s - info%residual_norms(k)**2)
              end do
              if ( (k - first)*step <= 0 ) cycle
              k = k - step

            end if

            ! compute Lehmann bounds for eigenvalues
            l = (k - first)*step + 1
            do i = min(first, k), max(first, k)
              do j = i, max(first, k)
                rr_matrices(i,j,3) = -step * rr_matrices(i,j,3) &
                  /sqrt((t - lambda(i))*(t - lambda(j)))
              end do
              rr_matrices(i,i,3) = lambda(i) + rr_matrices(i,i,3)
            end do
            i = min(first, k)

            call solve_sevp &
              ( 'N', l, rr_matrices(i,i,3), mm, keep%lambda(mm + i), &
                keep%lwork, keep%dwork, j )

            ! estimate errors
            do i = first, k, step
              q = info%residual_norms(i)
              s = max(q, step*(t - lambda(i)))
              info%err_X(i) = min(ONE, q/s) ! Davis-Kahan estimate for the sine
                                          ! of the angle between the Ritz vector
                                          ! and the exact invariant subspace
              ! Lehmann
              info%err_lambda(i) = step*(lambda(i) - keep%lambda(mm + i))
              r = max(abs(lambda(first)), abs(lambda(k)))
              if ( info%err_lambda(i) <= (keep%err_B + TOO_SMALL)*r ) then
                ! Lehmann estimate too polluted by round-off errors to be used
                q = q*q
                ! asymptotic quadratic residual estimate for the eigenvalue err
                info%err_lambda(i) = q/s
              end if
            end do

          end do

        end if

        do i = 1, m
          ! if error estimates not yet available, use
          ! eigenvalue decrements keep%dlmd and eigenvector shifts keep%dX
          ! as 'rough' ones for printing
          info%err_lambda(m + i) = info%err_lambda(i)
          if ( info%err_lambda(i) == NO_VALUE &
              .and. keep%rec > 0  ) then
            if ( keep%dlmd(i, keep%rec) > 0 ) &
              info%err_lambda(m + i) = keep%dlmd(i, keep%rec)
          end if
          if ( problem == 0 ) then
            if ( info%err_lambda(m + i) == NO_VALUE &
              .or. info%err_lambda(m + i) > info%residual_norms(i) ) &
              info%err_lambda(m + i) = info%residual_norms(i)
          end if
          if ( info%err_X(i) == NO_VALUE .and. keep%dX(i) > 0 ) then
            info%err_X(m + i) = keep%dX(i)
          else
            info%err_X(m + i) = abs(info%err_X(i))
          end if
        end do

        rci%job = SSMFE_CHECK_CONVERGENCE
        return

      case (CHECK_CONVERGENCE) select_step

        ! this is just to avoid 'possibly uninitialized' warning
        skip = .false.
        first = 0
        last = 0
        left_cnv = 0
        right_cnv = 0
        r = 0

        ! parameters affecting the accuracy of Ritz vectors
        s = keep%cond*epsilon(ONE)
        t = max(abs(keep%lambda(1)), abs(keep%lambda(keep%sizeX + keep%sizeZ)))

        ! verify user's info%converged flags

        ! an eigenpair is only accepted as converged if all previous ones on the
        ! same margin (i.e. those corresponding to eigenvalues that are further
        ! away from the middle spectrum than the one in focus) have either
        ! converged or stagnated

        ! here and below middle spectrum refers to eigenvalues that are not
        ! computed

        ! do each margin one after another
        do go = 1, -1, -2

          ! determine the scope to check
          select case ( go )
          case ( 1 )
            ! left margin
            skip = left == 0
            first = keep%firstX
            last = first + keep%leftX - 1
          case ( -1 )
            ! right margin
            skip = right == 0
            first = keep%firstX + keep%sizeX - 1
            last = first - keep%rightX + 1
          end select

          k = 0
          do i = first, last, go

            if ( skip ) then

              ! mark as non-converged because there is a non-converged
              ! eigenvalue further away from the middle spectrum
              info%converged(i) = 0

            else

              if ( info%converged(i) == 0 .and. keep%rec > 0 &
                .and. keep%sizeZ > 0 ) then
                ! check for the convergence stagnation, which may happen
                ! if the wanted error tolerance is not achieavable

                ! find the middle spectrum margin
                select case ( go )
                case ( 1 )
                  r = keep%lambda(keep%leftX + 1)
                case ( -1 )
                  r = keep%lambda(keep%leftX + keep%sizeZ)
                end select

                ! find maximal previous eigenvector shift over the cluster
                ! to which the eigenvalue in focus belongs
                q = keep%dX(m + i)
                do l = 1, m
                  if ( abs(lambda(i) - lambda(l)) &
                    < 10*epsilon(ONE)*max(abs(lambda(i)), abs(lambda(l))) ) &
                    q = max(q, keep%dX(m + l))
                end do

                if ( &
                  ! convergence history is long enough
                  keep%rec >= 5 &
                  .and. &
                  ! the eigenvector shift stagnated (has not reduced tangibly)
                  (keep%dX(i) > q_max * q &
                  .and. &
                  ! the eigenvector shift is of the order of the estimated
                  ! Ritz vector error caused by the conditioning of the
                  ! matrix X'*B*X (s*t) and errors in computing Rayleigh-Ritz
                  ! matrices (keep%err_A and keep%err_B)
                  min(keep%dX(i), keep%dX(m + i))*abs(r - lambda(i)) < &
                    10*(s*t + keep%err_A + abs(lambda(i))*keep%err_B) &
                    ) &
                  ) then
                  ! no further improvement in accuracy is to be expected,
                  ! mark as stagnated
                  info%converged(i) = -info%iteration
                  if ( err_est == SSMFE_KINEMATIC ) then
                    ! use a.c.f. for the error estimates
                    q = keep%q(m + i)
                    if ( q < ONE ) then
                      if ( keep%dlmd(i, keep%rec) > ZERO ) &
                        info%err_lambda(i) = keep%dlmd(i, keep%rec)*q/(ONE - q)
                      if ( keep%dX(i) > ZERO ) &
                        info%err_X(i) = min(ONE, keep%dX(i)*q/(ONE - q))
                    end if
                  end if
                end if
              end if

              if ( info%converged(i) /= 0 ) then
                ! count converged eigenpairs
                k = k + 1
                if ( info%converged(i) > 0 ) &
                  info%converged(i) = max(1, info%iteration)
              else
                ! skip the rest on this margin
                skip = .true.
              end if

            end if

          end do

          select case ( go )
          case ( 1 )
            left_cnv = k ! the number of newly converged on the left
          case ( -1 )
            right_cnv = k ! the number of newly converged on the right
          end select

        end do ! go

        if ( left < 0 ) then
          ! use different convergence verification procedure in the case
          ! where the largest eigenvalues are wanted, notably:
          ! accept an eigenvalue marked as converged by the user as
          ! converged indeed only if all eigenvalues that are not
          ! less by absolute value have converged
          left_cnv = 0
          right_cnv = 0
          first = keep%firstX
          last = keep%firstX + keep%sizeX - 1
          l = keep%firstX + keep%leftX - 1
          skip = .false.
          do while ( first <= last )
            if ( first <= l &
              .and. (abs(lambda(first)) > abs(lambda(last)) .or. last <= l) &
              ) then
              i = first
              j = keep%left_cnv + i - keep%firstX + 1
              k = 1
              first = first + 1
            else
              i = last
              j = -(keep%right_cnv + keep%firstX + keep%sizeX - i)
              k = -1
              last = last - 1
            end if
            if ( skip ) then
              info%converged(i) = 0
            else if ( info%converged(i) == 0 ) then
              skip = .true.
            else
              if ( k > 0 ) left_cnv = left_cnv + 1
              if ( k < 0 ) right_cnv = right_cnv + 1
            end if
          end do
        end if

        ! detect zero residual vectors
        rci%k = 0
        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == 0 ) rci%k = -1
        end do

        ! update counters
        keep%left_cnv = keep%left_cnv + left_cnv
        keep%right_cnv = keep%right_cnv + right_cnv
        keep%leftX = keep%leftX - left_cnv
        keep%rightX = keep%rightX - right_cnv
        keep%sizeX = keep%leftX + keep%rightX

        ! update first non-converged eigenpair index
        keep%firstX = keep%firstX + left_cnv

        ! set next step
        keep%step = SAVE_RIGHT_CONVERGED

        if ( left_cnv > 0 ) then
          ! shift the Rayleigh-Ritz matrices up the main diagonal by the
          ! number of newly converged eigenvalues on the left to remove
          ! rows and columns that are no longer needed
          do j = 1, keep%sizeX
            call copy &
              ( keep%sizeX, rr_matrices(left_cnv + 1, left_cnv + j, 1), 1, &
                rr_matrices(1, j, 1), 1 )
            call copy &
              ( keep%sizeX, rr_matrices(left_cnv + 1, left_cnv + j, 2), 1, &
                rr_matrices(1, j, 2), 1 )
          end do
        end if

        ! instruct the user to save newly converged left eigenpairs
        rci%job = SSMFE_SAVE_CONVERGED
        rci%nx = left_cnv
        rci%ny = right_cnv
        rci%jx = keep%firstX - left_cnv
        rci%kx = 0
        rci%jy = rci%jx
        rci%ky = kBX
        rci%i = 1
        return

      case (SAVE_RIGHT_CONVERGED)

        ! set next step
        keep%step = APPLY_PRECONDITIONER

        ! instruct the user to save newly converged right eigenpairs
        rci%job = SSMFE_SAVE_CONVERGED
        rci%nx = rci%ny
        rci%jx = keep%firstX + keep%sizeX + rci%nx - 1
        rci%kx = 0
        rci%jy = rci%jx
        rci%ky = kBX
        rci%i = -1
        return

      case (APPLY_PRECONDITIONER) select_step

        if ( rci%k == -1 ) then
          ! zero residual vectors cannot be handled by the search directions
          ! clean-up procedure implemented here, so restart the computation
          rci%job = SSMFE_RESTART
          ! some vectors must remain in X, the rest overwritten with random
          rci%jx = keep%firstX ! first vector to be kept in X
          rci%nx = keep%sizeX  ! number of vectors to keep in X
          rci%k = 0 ! restart is mandatory
          return
        end if

        if ( left == 0 .and. keep%leftX > 0 ) then
          ! if left eigenpairs are no longer needed, stop computing them

          ! move the index of the first non-converged eigenpair
          keep%firstX = keep%firstX + keep%leftX
          ! only right eigenvectors are now wanted
          keep%sizeX = keep%rightX
          ! shift the Rayleigh-Ritz matrices up the main diagonal by the
          ! number of left eigenpairs to remove rows and columns that are
          ! no longer needed
          do j = 1, keep%sizeX
            call copy &
              ( keep%sizeX, rr_matrices(keep%leftX + 1, keep%leftX + j, 1), 1, &
                rr_matrices(1, j, 1), 1 )
            call copy &
              ( keep%sizeX, rr_matrices(keep%leftX + 1, keep%leftX + j, 2), 1, &
                rr_matrices(1, j, 2), 1 )
          end do
          keep%leftX = 0
        end if

        ! if right eigenpairs are no longer needed, stop computing them
        if ( right == 0 ) keep%rightX = 0

        keep%sizeX = keep%leftX + keep%rightX

        if ( keep%sizeX == 0 ) then
          ! all iterated eigenpairs have converged, restart in case not all
          ! wanted have converged yet
          rci%job = SSMFE_RESTART
          ! no eigenvectors to keep, fill X with random vectors
          rci%jx = m + 1
          rci%nx = 0
          rci%k = 0 ! restart is mandatory
          return
        end if

        ! set next step
        if ( keep%sizeZ > 0 ) then
          ! conjugate search directions
          keep%step = COMPUTE_AZ
        else
          if ( problem < 0 ) then
            ! B*Y is already available in block W
            keep%step = COPY_W_TO_BYZ
          else
            ! compute B*Y
            keep%step = RECOMPUTE_BY
          end if
        end if

        if ( problem >= 0 ) then
          ! instruct the user to apply preconditioner
          rci%job = SSMFE_APPLY_PREC
          rci%nx = keep%sizeY
          rci%jx = keep%firstXn
          rci%kx = kW
          rci%jy = 1
          rci%ky = kY
          return
        end if

      case (COMPUTE_AZ) select_step

        ! set next step
        keep%step = COMPUTE_ZAY

        ! if A*Z is not available implicitly
        if ( .not. minAprod ) then
          ! instruct user to multiply previous search directions Z by A
          rci%job = SSMFE_APPLY_A
          rci%nx = keep%sizeZ
          rci%jx = 1
          rci%kx = kZ
          rci%jy = 1
          rci%ky = kAYZ
          return
        end if

      case (COMPUTE_ZAY) select_step

        ! instruct user to compute (A*Z)'*Y or (A*B*Z)'*(B*Y)
        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kAYZ
        rci%ny = keep%sizeY
        rci%jy = 1

        if ( problem < 0 ) then
          ! compute (A*B*Z)'*(B*Y)
          rci%ky = kW
        else
          ! compute (A*Z)'*Y
          rci%ky = kY
        end if

        rci%i = m + 1
        rci%j = m + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_BY
        return

      case (COMPUTE_BY) select_step

        ! set next step
        keep%step = COMPUTE_ZBY

        if ( problem == 0 ) then
          kBYZ = kZ
          keep%kBYZ = kBYZ
        else if ( problem > 0 ) then
          ! instruct user to compute B*Y
          rci%job = SSMFE_APPLY_B
          rci%jx = 1
          rci%kx = kY
          rci%nx = keep%sizeY
          rci%jy = 1
          rci%ky = kW
          return
        end if
        ! if problem < 0, B*Y is already in block W (see step APPLY_B_TO_Y)

      case (COMPUTE_ZBY) select_step

        ! set next step
        if ( problem >= 0 ) then
          ! compute y'*B*y for each column y of Y
          keep%step = COMPUTE_YBY_DIAG
        else
          ! already computed on step COMPUTE_RR
          keep%step = CONJUGATE_Y
        end if

        ! instruct user to compute Z'*Y or Z'*(B*Y)
        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%ny = keep%sizeY
        rci%jy = 1
        if ( problem == 0 ) then
          rci%ky = kY
        else
          rci%ky = kW
        end if
        rci%i = m + 1
        rci%j = m + 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        return

      case (COMPUTE_YBY_DIAG)

        ! set next step
        keep%step = CONJUGATE_Y

        ! instruct user to compute y'*B*y for each column y of Y
        rci%job = SSMFE_COMPUTE_DOTS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        if ( problem == 0 ) then
          rci%ky = kY
        else
          rci%ky = kW
        end if
        rci%i = 1
        rci%j = 1
        rci%k = 3
        return

      case (CONJUGATE_Y) select_step
        ! compute conjugation matrix of Jacobi scheme

        do i = 1, keep%sizeY

          if ( i > keep%leftXn - keep%new_left &
            .and. i <= keep%leftXn + keep%new_right ) then
            ! previous search direction not available,
            ! leave this search direction unchanged
            do j = 1, keep%sizeZ
              rr_matrices(m + j, m + i, 2) = NIL
            end do
            cycle
          end if

          iX = keep%firstXn + i - 1
          call axpy &
            ( keep%sizeZ, -lambda(iX)*UNIT, rr_matrices(m + 1, m + i, 1), 1, &
              rr_matrices(m + 1, m + i, 2), 1 )

          ! find the position of y'*B*y in rr_matrices(:,:,3)
          if ( problem >= 0 ) then
            k = i ! cf. COMPUTE_YBY_DIAG
          else
            k = keep%firstXn + i - 1 ! cf. COMPUTE_RR
          end if
          r = sqrt(abs(rr_matrices(k, k, 3)))
          t = MAX_BETA*r ! conjugation matrix elements ceiling

          skip = .false.
          do j = 1, keep%sizeZ
            s = keep%lambda(keep%leftXn + j) - lambda(iX)
            ! check before dividing by s
            if ( abs(rr_matrices(m + j, m + i, 2)) < t*abs(s) ) then
              ! this conjugation matrix element is ok
              r = rr_matrices(m + j, m + i, 2)/s
              rr_matrices(m + j, m + i, 2) = r
              t = sqrt(max((t - r)*(t + r), ZERO))
            else
              ! this conjugation matrix element is too large, skip conjugation
              skip = .true.
              exit
            end if
          end do
          if ( skip ) then
            ! leave this search direction unchanged
            do j = 1, keep%sizeZ
              rr_matrices(m + j, m + i, 2) = NIL
            end do
          end if

        end do

        ! set next step
        if ( minBprod ) then
          ! update B*Y implicitly using current B*Y and B*Z
          keep%step = UPDATE_BY
        else
          ! compute B*Y explicitly
          keep%step = RECOMPUTE_BY
        end if

        ! instruct the user to update Y by conjugating with Z
        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kY
        rci%i = m + 1
        rci%j = m + 1
        rci%k = 2
        rci%alpha = -UNIT
        rci%beta = UNIT

        return

      case (UPDATE_BY) select_step

        ! set next step
        keep%step = COPY_W_TO_BYZ

        ! instruct the user to update B*Y using current B*Y and B*Z
        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kBYZ
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kW
        rci%i = m + 1
        rci%j = m + 1
        rci%k = 2
        rci%alpha = -UNIT
        rci%beta = UNIT
        return

      case (COPY_W_TO_BYZ) select_step

        ! set next step
        keep%step = APPLY_CONSTRAINTS

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kW
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 0
        return

      case (RECOMPUTE_BY) select_step

        ! set next step
        keep%step = APPLY_CONSTRAINTS

        if ( problem == 0 ) then
          kBYZ = kY
          keep%kBYZ = kBYZ
        else
          ! instruct user to compute B*Y
          if ( .not. minBprod ) keep%kBYZ = kW
          rci%job = SSMFE_APPLY_B
          rci%nx = keep%sizeY
          rci%jx = 1
          rci%kx = kY
          rci%jy = 1
          rci%ky = keep%kBYZ
          return
        end if

      case (APPLY_CONSTRAINTS) select_step

        ! set next step
        keep%step = SCALE_Y

        ! B-orthogonalize Y to constraints Xc

        ! instruct user to update Y := Y - Xc*Xc'*(B*Y)
        rci%job = SSMFE_APPLY_CONSTRAINTS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 1
        rci%j = 1
        rci%k = 3
        return

      case (SCALE_Y) select_step

        ! set next step
        keep%step = COMPUTE_YBY

        ! B-normalize Y

        ! instruct user to update y := y/(y'*(B*y)) for each column y of Y
        rci%job = SSMFE_SCALE_VECTORS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 1
        rci%j = 1
        rci%k = 3
        return

      case (COMPUTE_YBY) select_step

        ! set next step
        keep%step = COMPUTE_XBY

        ! instruct user to compute Y'*(B*Y)
        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeX + 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        return

      case (COMPUTE_XBY) select_step

        ! set next step
        keep%step = CLEANUP_Y

        ! instruct user to compute X'*(B*Y)
        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 1
        rci%j = keep%sizeX + 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        return

      case (CLEANUP_Y) select_step

        ! get rid of near-linearly-dependent search directions

        sizeXY = keep%sizeX + keep%sizeY
        rci%ky = keep%sizeY

        if ( problem /= 0 ) then
          ! matrix Y'*B*Y is symmetric in exact arithmetic,
          ! hence non-symmetry can be used as a measure of the error
          ! in computing B*Y, which affects the accuracy of the
          ! Rayleigh-Ritz procedure
          iY = keep%sizeX + 1
          jY = keep%sizeX + keep%sizeY
          r = ZERO
          s = ZERO
          do i = iY, jY
            do j = iY, i - 1
              t = abs(rr_matrices(i,j,1) - conjugate(rr_matrices(j,i,1)))
              r = max(r, t)
              t = t/sqrt(abs(rr_matrices(i,i,1)*rr_matrices(j,j,1)))
              s = max(s, t)
            end do
          end do
          keep%err_B = max(keep%err_B, r)
        end if

        delta = GRAM_RCOND_MIN

        ! check if the Gram matrix G =[X Y]'*B*[X Y] has eigenvalues that are
        ! too small (less than delta)

        ! perform Cholesky factorization of the Gram matrix
        call copy( mm*sizeXY, rr_matrices, 1, rr_matrices(1,1,3), 1 )
        call chol( 'U', sizeXY, rr_matrices(1,1,3), mm, i )
        if ( i /= 0 .and. i <= keep%sizeX ) info%flag = INDEFINITE_B
        if ( info%flag /= 0 ) rci%job = SSMFE_ABORT
        if ( rci%job < 0 ) return

        if ( i /= 0 ) then
          ! mass matrix is numrically degenerate
          nsmall = 1
        else
          ! perform 3 inverse iterations to compute approximations to
          ! smallest eigenvalues
          nsmall = 0
          rr_matrices(1 : sizeXY - 1, mm, 2) = ZERO
          rr_matrices(sizeXY, mm, 2) = ONE
          s = ONE
          do step = 1, 3
            call trsm &
              ( 'L', 'U', 'T', 'N', sizeXY, 1, UNIT, &
                rr_matrices(1,1,3), mm, rr_matrices(1,mm,2), mm )
            t = norm(sizeXY, rr_matrices(1,mm,2), 1)
            q = s**2/t**2 ! a Rayleigh quotient for the inverse of G
            if ( q <= delta ) exit ! small eigenvalue detected
            call trsm &
              ( 'L', 'U', 'N', 'N', sizeXY, 1, UNIT, &
                rr_matrices(1,1,3), mm, rr_matrices(1,mm,2), mm )
            s = norm(sizeXY, rr_matrices(1,mm,2), 1)
          end do
          ! order of magnitude of the condition number of G
          keep%cond = ONE/q
          if ( q <= delta ) nsmall = 1
        end if

        if ( nsmall > 0 ) then

          ! initial guess for the condition number of G
          keep%cond = ONE/delta

          ! discard some search directions

          call copy &
            ( mm * keep%sizeY, rr_matrices(1, keep%sizeX + 1, 1), 1, &
              rr_matrices(1, keep%sizeX + 1, 3), 1 )

          ! compute the off-diagonal keep%sizeX by keep%sizeY block of the
          ! upper triangular matrix U such that the Gram matrix is equal to
          ! U'*U (i.e. Cholesky U-factor of the Gram matrix).
          call trsm &
            ( 'L', 'U', 'N', 'N', keep%sizeX, keep%sizeY, UNIT, &
              rr_matrices(1, 1, 3), mm, rr_matrices(1, keep%sizeX + 1, 1), mm )

          rr_matrices(1 : sizeXY, mm, 2) = ZERO

          ! start selection procedure for the basis vectors of Rayleigh-Ritz
          ! procedure

          ! the initial selection S is: S = X

          ! in what follows, G denotes the matrix of size sizeXY stored in the
          ! array rr_matrices(:,:,3)

          ! initially, the upper triangular part of G contains keep%sizeX
          ! rows of U

          ! after each successful selection, one more row of U ends up in G

          ! indices of selected columns of Y
          forall ( j = 1 : keep%sizeY ) ind(j) = j
          iX = keep%sizeX ! number of selected basis vectors
          iY = keep%sizeY ! number of remaining vectors (search directions)
          do while ( iX < sizeXY )

            ! norms of the columns of the iX by iY off-diagonal block of G
            ! are norms of projections of remaining vectors
            ! onto the subspace spanned by selected vectors

            ! find the remaining vector v with smallest projection, i.e.
            ! the one that forms greatest angle with the selected vectors
            ! subspace
            l = iX + 1
            s = norm(iX, rr_matrices(1,l,3), 1)
            do j = iX + 2, sizeXY
              t = norm(iX, rr_matrices(1,j,3), 1)
              if ( t < s ) then
                l = j
                s = t
              end if
            end do
            k = iX + 1
            if ( l /= k ) then
              ! swap rows k and l and columns k and l of matrix G
              call copy( iX, rr_matrices(1,k,3), 1, rr_matrices(k,1,3), mm )
              call copy( iX, rr_matrices(1,l,3), 1, rr_matrices(1,k,3), 1 )
              call copy( iX, rr_matrices(k,1,3), mm, rr_matrices(1,l,3), 1 )
              call copy( iY, rr_matrices(k,k,3), 1, rr_matrices(k,1,3), 1 )
              call copy( iY, rr_matrices(k,l,3), 1, rr_matrices(k,k,3), 1 )
              call copy( iY, rr_matrices(k,1,3), 1, rr_matrices(k,l,3), 1 )
              call copy( iY, rr_matrices(k,k,3), mm, rr_matrices(k,1,3), 1 )
              call copy( iY, rr_matrices(l,k,3), mm, rr_matrices(k,k,3), mm )
              call copy( iY, rr_matrices(k,1,3), 1, rr_matrices(l,k,3), mm )
              ! the remaining parts of rows k and l are not used
              j = ind(k - keep%sizeX)
              ind(k - keep%sizeX) = ind(l - keep%sizeX)
              ind(l - keep%sizeX) = j
            end if
            ! compute the square of the sine of the angle between v
            ! and the subspace spanned by the selected vectors
            s = norm(iX, rr_matrices(1,k,3), 1)**2
            t = rr_matrices(k,k,3) - s
            if ( t <= TOO_SMALL + keep%err_B ) then
              ! v and all other remaining vectors are too close to the subspace
              ! of selected ones and must be discarded, selection complete
              exit
            end if
            ! compute the k-th diagonal element of the U factor
            s = sqrt(t)
            rr_matrices(k,k,3) = s
            ! which, together with the k-th column of G above the diagonal
            ! gives us Uv, the U-factor of [S v]'*B*[S v]

            ! perform 3 inverse iterations for Uv'*Uv = [S v]'*B*[S v]
            if ( iX == keep%sizeX ) rr_matrices(k, mm, 2) = ONE
            q = ONE
            do step = 1, 3
              call trsm &
                ( 'L', 'U', 'T', 'N', k, 1, UNIT, &
                  rr_matrices(1, 1, 3), mm, rr_matrices(1, mm, 2), mm )
              t = norm(k, rr_matrices(1, mm ,2), 1)
              r = q
              q = ONE/t**2
              if ( q <= delta ) exit
              call trsm &
                ( 'L', 'U', 'N', 'N', k, 1, UNIT, &
                  rr_matrices(1, 1, 3), mm, rr_matrices(1, mm, 2), mm )
              t = norm(k, rr_matrices(1, mm, 2), 1)
              call dscal(k, ONE/t, rr_matrices(1, mm, 2), 1)
              if ( q > 0.9*r ) exit ! iterations stagnated
            end do
            ! order of magnitude of the condition number of G
            keep%cond = ONE/q

            ! q approximates an upper bound for the smallest eigenvalue of
            ! [S v]'*B*[S v]

            ! if it is not greater than delta, v and all remaining vectors
            ! must be discarded
            if ( q <= delta ) exit

            ! compute the k-th row of U right of the diagonal
            if ( iY > 1 ) &
              call gemm &
                ( TRANS, 'N', 1, iY - 1, iX, &
                  -UNIT, rr_matrices(1, k, 3), mm, &
                  rr_matrices(1, iX + 2, 3), mm, &
                  UNIT, rr_matrices(k, iX + 2, 3), mm )
            forall ( j = 2 : iY ) &
              rr_matrices(k, iX + j, 3) = rr_matrices(k, iX + j, 3)/s

            ! update S := [S v]
            iX = k
            iY = iY - 1

          end do
          k = iX - keep%sizeX ! number of selected search directions

          if ( k < 1 ) then
            info%flag = NO_SEARCH_DIRECTIONS_LEFT
            keep%step = QUIT
            cycle do_select_step
          end if

          ! collect the U-factor of the Gram matrix for selected vectors
          call copy( mm*iX, rr_matrices(1,1,3), 1, rr_matrices, 1 )

          keep%sizeY = k

        else

          ! collect the U-factor of the whole Gram matrix
          call copy( mm*sizeXY, rr_matrices(1,1,3), 1, rr_matrices, 1 )

        end if

        ! set next step
        keep%step = COMPUTE_AY

        if ( nsmall > 0 ) then
          ! instruct user to rearrange search directions
          rci%job = SSMFE_COPY_VECTORS
          rci%kx = kY
          if ( problem /= 0 ) then
            rci%ky = kBYZ
          else
            rci%ky = kY
          end if
          rci%nx = keep%sizeY
          rci%i = 1
          return
        end if

      case (COMPUTE_AY) select_step

        ! set next step
        keep%step = RR_IN_XY

        ! instruct user to compute A*Y or A*(B*Y)
        rci%job = SSMFE_APPLY_A
        rci%nx = keep%sizeY
        rci%jx = 1
        if ( problem < 0 ) then
          ! compute A*(B*Y)
          rci%kx = kBYZ
        else
          ! compute A*Y
          rci%kx = kY
        end if
        rci%jy = 1
        rci%ky = kAYZ
        return

      case (RR_IN_XY) select_step

        ! set next step
        keep%step = COMPUTE_XAY

        ! instruct user to compute Y'*(A*Y) or (B*Y)'*(A*B*Y)
        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeY
        rci%jx = 1
        if ( problem < 0 ) then
          ! compute (B*Y)'*(A*(B*Y))
          rci%kx = kBYZ
        else
          ! compute Y'*(A*Y)
          rci%kx = kY
        end if
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kAYZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeX + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        return

      case (COMPUTE_XAY) select_step

        ! set next step
        keep%step = PREPARE_MATRICES

        ! instruct user to compute X'*(A*Y) or (A*B*X)'*(B*Y)
        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        if ( problem < 0 ) then
          ! compute (A*B*X)'*(B*Y)
          rci%kx = kAX
          rci%ky = kBYZ
        else
          ! compute X'*(A*Y)
          rci%kx = 0
          rci%ky = kAYZ
        end if
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%i = 1
        rci%j = keep%sizeX + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        return

      case (PREPARE_MATRICES) select_step

        iX = 1
        jX = keep%sizeX
        iY = jX + 1
        jY = jX + keep%sizeY
        sizeXY = jY

        ! estimate the error of multiplication by A
        ! by the deviation of Y'*A*Y from symmetry
        r = ZERO
        s = ZERO
        do i = iY, jY
          do j = iY, i - 1
            t = abs(rr_matrices(i,j,2) - conjugate(rr_matrices(j,i,2)))
            r = max(r, t)
            t = t/sqrt(abs(rr_matrices(i,i,2)*rr_matrices(j,j,2)))
            s = max(s, t)
          end do
        end do
        keep%err_A = max(keep%err_A, r)

        ! compute the lower triangular part of H, where
        ! H = [X Y]'*A*[X Y] if problem >= 0 and
        ! H = [X Y]'*B*A*B*[X Y] otherwise
        do i = 1, sizeXY
          do j = 1, i - 1
            rr_matrices(i,j,2) = conjugate(rr_matrices(j,i,2))
          end do
        end do

        ! transform H := (U')^{-1} H U^{-1}, where U is computed at step
        ! CLEANUP_Y
        call trsm &
          ( 'L', 'U', 'T', 'N', sizeXY, sizeXY, UNIT, &
            rr_matrices, mm, rr_matrices(1,1,2), mm )
        call trsm &
          ( 'R', 'U', 'N', 'N', sizeXY, sizeXY, UNIT, &
            rr_matrices, mm, rr_matrices(1,1,2), mm )

        ! save the YY block of H
        do j = 1, keep%sizeY
          call copy &
            ( keep%sizeY, rr_matrices(iY, jX + j, 2), 1, &
              rr_matrices(m + 1, j, 1), 1 )
        end do

        ! make space for the next record in eigenvalue decrements history
        if ( keep%rec < keep%RECORDS ) then
          keep%rec = keep%rec + 1
        else
          do j = 1, keep%rec - 1
            do i = 1, m
              keep%dlmd(i, j) = keep%dlmd(i, j + 1)
            end do
          end do
        end if
        l = keep%rec

        ! set next step
        keep%step = ANALYSE_RR_RESULTS

        ! compute eigenvalues and eigenvectors of H
        call solve_sevp &
          ( 'V', sizeXY, rr_matrices(1,1,2), mm, keep%lambda, &
            keep%lwork, keep%dwork, i )

      case (ANALYSE_RR_RESULTS) select_step

        iX = 1
        jX = keep%sizeX
        iY = jX + 1
        jY = jX + keep%sizeY
        sizeXY = jY

        ! select which Ritz vectors to use as approximations to which wanted
        ! eigenvectors

        if ( left < 0 ) then

          ! the number of eigenpairs yet to be computed
          l = right + control%extra_right - keep%left_cnv - keep%right_cnv

          ! if some eigenpairs have converged on this iteration, the number of
          ! iterated vectors keep%sizeX has reduced and become smaller than
          ! the block size m

          ! for as long as l > keep%sizeX, the algorithm tries to increase
          ! keep%sizeX to improve the convergence to extreme eigenpairs

          ! the new value, however, cannot be greater than either the block
          ! size m or the number of Ritz values sizeXY

          ! it should also not be greater than l, the number of eigenpairs yet
          ! to be computed, to reduce the cost of last iterations if the block
          ! size is large

          keep%sizeXn = max(min(m, l, sizeXY), keep%sizeX) ! new keep%sizeX

          ! select Ritz vectors corresponding to the largest keep%sizeXn
          ! Ritz values as new approximations to wanted eigenvectors,
          ! distributing the workspace between left and right eigenpairs
          ! according to the sign of Ritz values

          i = 1      ! index of a negative Ritz value
          j = sizeXY ! index of a positive Ritz value

          ! new number of approximate eigenvectors corresponding to left
          ! (negative) eigenvalues
          keep%leftXn = 0
          do k = 1, keep%sizeXn
            if ( keep%lambda(i) > ZERO ) exit ! no more negative Ritz values
            if ( abs(keep%lambda(i)) > abs(keep%lambda(j)) ) then
              keep%leftXn = keep%leftXn + 1
              i = i + 1
            else
              j = j - 1
            end if
          end do
          ! new number of approximate eigenvectors corresponding to right
          ! (positive) eigenvalues
          keep%rightXn = keep%sizeXn - keep%leftXn

          keep%firstXn = keep%firstX
          ! if keep%sizeXn > keep%sizeX, try to decrease keep%firstXn
          ! proportionally
          if ( keep%sizeX < m ) &
            keep%firstXn = keep%firstXn &
              - (keep%firstX - 1)*(keep%sizeXn - keep%sizeX) &
              /(m - keep%sizeX)

          keep%new_left = max(0, keep%leftXn - keep%leftX)
          keep%new_right = max(0, keep%rightXn - keep%rightX)

        else

          ! just as above, the algorithm tries to increase the number of
          ! iterated approximate eigenvectors whenever it falls below the
          ! block size m

          ! the maximal number of additional eigenvectors is keep%sizeY,
          ! proportionally allocate part of it for left eigenvectors
          i = int(keep%sizeY*left*ONE/(left + right))
          ! if possible, make this part non-zero
          if ( left > 0 .and. keep%sizeY > 1 ) i = max(i, 1)
          ! make sure it is less than keep%sizeY so that the number of right
          ! eigenvectors can be increased too
          if ( right > 0 .and. keep%sizeY > 1 ) i = min(i, keep%sizeY - 1)

          ! the number of iterated left eigenvectors cannot be larger than
          k = keep%firstX + keep%leftX - 1

          ! the number of iterated left eigenvectors cannot be larger than
          l = keep%leftX + i

          if ( left == 0 ) then
            keep%leftXn = 0
            keep%new_left = 0
          else
            ! the number of left eigenpairs yet to be computed
            j = control%extra_left + left - keep%left_cnv
            ! try to bring keep%leftX closer to j
            keep%leftXn = min(j, k, l)
            ! keep%leftX must not decrease
            keep%leftXn = max(keep%leftXn, keep%leftX)
            keep%new_left = keep%leftXn - keep%leftX
          end if
          ! change keep%firstX accordingly
          keep%firstXn = keep%firstX + keep%leftX - keep%leftXn

          ! the case of the right eigenpairs is symmetric
          if ( right == 0 ) then
            keep%rightXn = 0
            keep%new_right = 0
          else
            j = control%extra_right + right - keep%right_cnv
            k = m - k
            l = sizeXY - l
            keep%rightXn = min(j, k, l)
            keep%rightXn = max(keep%rightXn, keep%rightX)
            keep%new_right = keep%rightXn - keep%rightX
          end if

          keep%sizeXn = keep%leftXn + keep%rightXn

        end if

        ! if keep%firstX has moved, the information about left eigenpairs
        ! needs to be moved accordingly
        k = keep%firstX - keep%firstXn
        l = keep%firstX + keep%leftX - 1
        if ( k > 0 ) then
          do j = 1, l - k
            lambda(j) = lambda(j + k)
            info%residual_norms(j) = info%residual_norms(j + k)
            info%err_lambda(j) = info%err_lambda(j + k)
            info%err_X(j) = info%err_X(j + k)
            info%converged(j) = info%converged(j + k)
            keep%q(j) = keep%q(j + k)
            keep%q(m + j) = keep%q(m + j + k)
            do i = 1, keep%rec
              keep%dlmd(j,i) = keep%dlmd(j + k, i)
            end do
          end do
        end if
        ! the information about newly added left eigenpairs must be initialized
        if ( k >= 0 ) then
          do j = l - k + 1, keep%firstXn + keep%leftXn - 1
            info%err_lambda(j) = NO_VALUE
            info%err_X(j) = NO_VALUE
            info%converged(j) = 0
            keep%q(j) = ONE
            keep%q(m + j) = ONE
            do i = 1, keep%rec
              keep%dlmd(j,i) = 0
            end do
          end do
        end if

        ! the case of the right eigenpairs is symmetric
        k = keep%firstXn + keep%sizeXn - keep%firstX - keep%sizeX
        l = keep%firstX + keep%leftX
        if ( k > 0 ) then
          do j = m, l + k, -1
            lambda(j) = lambda(j - k)
            info%residual_norms(j) = info%residual_norms(j - k)
            info%err_lambda(j) = info%err_lambda(j - k)
            info%err_X(j) = info%err_X(j - k)
            info%converged(j) = info%converged(j - k)
            keep%q(j) = keep%q(j - k)
            keep%q(m + j) = keep%q(m + j - k)
            do i = 1, keep%rec
              keep%dlmd(j,i) = keep%dlmd(j - k, i)
            end do
          end do
        end if
        if ( k >= 0 ) then
          do j = keep%firstXn + keep%leftXn, l + k - 1
            info%err_lambda(j) = NO_VALUE
            info%err_X(j) = NO_VALUE
            info%converged(j) = 0
            keep%q(j) = ONE
            keep%q(m + j) = ONE
            do i = 1, keep%rec
              keep%dlmd(j,i) = 0
            end do
          end do
        end if

        if ( keep%rightXn > 0 ) then
          ! move the selected right Ritz vectors next to the selected left
          k = sizeXY - keep%rightXn
          l = keep%leftXn
          call copy &
            ( (k - l)*mm, rr_matrices(1, l + 1, 2), 1, rr_matrices(1,1,3), 1 )
          do j = 1, keep%rightXn
            call copy &
              ( sizeXY, rr_matrices(1, k + j, 2), 1, &
                rr_matrices(1, l + j, 2), 1 )
          end do
          call copy &
            ( (k - l)*mm, rr_matrices(1,1,3), 1, &
              rr_matrices(1, keep%sizeXn + 1, 2), 1 )
        end if

        ! record the last used eigenvector shifts as previous in proper place
        k = keep%firstX - keep%firstXn
        l = keep%firstX + keep%leftX - 1
        if ( k >= 0 ) then
          do j = 1, l - k
            keep%dX(m + j) = keep%dX(j + k)
          end do
          do j = l - k + 1, keep%firstXn + keep%leftXn - 1
            keep%dX(m + j) = ONE
          end do
        end if
        k = keep%firstXn + keep%sizeXn - keep%firstX - keep%sizeX
        l = keep%firstX + keep%leftX
        if ( k >= 0 ) then
          do j = m, l + k, -1
            keep%dX(m + j) = keep%dX(j - k)
          end do
          do j = keep%firstXn + keep%leftXn, l + k - 1
            keep%dX(m + j) = ONE
          end do
        end if

        ! compute the new eigenvector shifts, the sines of the angles
        ! between new approximate eigenvectors and the subspace spanned
        ! by old ones, and approximations to the new eigenvalue decrements,
        ! to be used if the direct computation of decrements by subtraction
        ! is too inaccurate because they are too close to the machine accuracy
        call gemm &
          ( 'N', 'N', keep%sizeY, keep%sizeXn, keep%sizeY, &
            UNIT, rr_matrices(m + 1, 1, 1), mm, &
            rr_matrices(keep%sizeX + 1, 1, 2), mm, &
            NIL, rr_matrices(1,1,3), mm )
        do j = 1, keep%sizeXn
          jX = keep%firstXn + j - 1
          ! the sine of the angle between the new jX-th eigenvector and
          ! the subspace spanned by the old ones
          s = norm(keep%sizeY, rr_matrices(keep%sizeX + 1, j, 2), 1)
          t = dot &
            (keep%sizeY, rr_matrices(keep%sizeX + 1, j, 2), 1, &
              rr_matrices(1, j, 3), 1)
          ! an approximation to the decrement in jX-th eigenvalue
          t = abs(t - lambda(jX)*s*s)
          keep%dX(jX) = s
          if ( s < 0.1 ) then ! approximation is acceptable
            keep%dlmd(jX, keep%rec) = t
          else ! too far from convergence, mark as not available
            keep%dlmd(jX, keep%rec) = ZERO
          end if
        end do
        do jX = 1, keep%firstXn - 1
          keep%dX(jX) = ZERO
        end do
        do jX = keep%firstXn + keep%sizeXn, m
          keep%dX(jX) = ZERO
        end do

        ! compute the eigenvectors of the generalized eigenvalue problem
        ! for the original H and G = [X Y]'*B*[X Y]
        call trsm &
          ( 'L', 'U', 'N', 'N', sizeXY, sizeXY, UNIT, &
            rr_matrices, mm, rr_matrices(1,1,2), mm )

        ! the Ritz vectors that have not been selected as new approximate
        ! eigenvectors will be used as previous search directions on the
        ! next iterations
        k = keep%sizeXn - keep%sizeX
        keep%sizeZ = keep%sizeY - k

        ! compute new eigenvector approximations X and previous search
        ! directions Z based on X, Y and Q, where the columns q of Q are
        ! eigenvectors of H q = lambda G q

        ! the first keep%sizeXn columns of Q yield X and the rest yield Z

        if ( minAprod ) then
          ! update A*X and compute A*Z implicitly from old A*X and A*Y
          ! (or A*B*X and A*B*Z, if problem < 0)
          keep%step = PUT_AZQ_IN_Z
        else if ( minBprod ) then ! update B X and compute B Z
          ! update B*X and compute B*Z implicitly from old B*X and B*Y
          keep%step = PUT_BZQ_IN_Z
        else
          ! update X and compute Z from old X and Y
          keep%step = PUT_YQ_IN_Z
        end if

      case (PUT_AZQ_IN_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kAYZ
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_AXQ_TO_Z
        return

      case (ADD_AXQ_TO_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kAX
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_AXQ_IN_W
        return

      case (PUT_AXQ_IN_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kAX
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_AZQ_TO_W
        return

      case (ADD_AZQ_TO_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kAYZ
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = keep%sizeX + 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_W_IN_AX
        return

      case (PUT_W_IN_AX) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeXn
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = keep%firstXn
        rci%ky = kAX
        rci%i = 0
        keep%step = PUT_Z_IN_AZ
        return

      case (PUT_Z_IN_AZ) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%jy = 1
        rci%ky = kAYZ
        rci%i = 0
        if ( minBprod ) then
          keep%step = PUT_BZQ_IN_Z
        else
          keep%step = PUT_YQ_IN_Z
        end if
        return

      case (PUT_BZQ_IN_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kBYZ
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_BXQ_TO_Z
        return

      case (ADD_BXQ_TO_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kBX
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_BXQ_IN_W
        return

      case (PUT_BXQ_IN_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kBX
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_BZQ_TO_W
        return

      case (ADD_BZQ_TO_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kBYZ
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = keep%sizeX + 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_W_IN_BX
        return

      case (PUT_W_IN_BX) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeXn
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = keep%firstXn
        rci%ky = kBX
        rci%i = 0
        keep%step = PUT_Z_IN_BZ
        return

      case (PUT_Z_IN_BZ) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 0
        keep%step = PUT_YQ_IN_Z
        return

      case (PUT_YQ_IN_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_XQ_TO_Z
        return

      case (ADD_XQ_TO_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_XQ_IN_W
        return

      case (PUT_XQ_IN_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_YQ_TO_W
        return

      case (ADD_YQ_TO_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%k = 2
        rci%i = keep%sizeX + 1
        rci%j = 1
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_W_IN_X
        return

      case (PUT_W_IN_X) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeXn
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = keep%firstXn
        rci%ky = 0
        rci%i = 0
        keep%step = CHECK_THE_GAP
        return

      case (CHECK_THE_GAP) select_step

        ! update the counters and the position of the first eigenpair in X
        keep%leftX = keep%leftXn
        keep%rightX = keep%rightXn
        keep%sizeX = keep%sizeXn
        keep%firstX = keep%firstXn
        keep%iteration = keep%iteration + 1
        info%iteration = keep%iteration

        ! check the gap between the computed eigenvalues and the rest
        ! of the spectrum, and suggest restart if it is too small, for
        ! the sake of faster convergence

        ! Ritz values corresponding to the Ritz vectors that have not been
        ! selected as new approximate eigenvectors are used as approximations
        ! for eigenvalues in the rest of the spectrum

        sizeXY = keep%sizeX + keep%sizeZ
        rci%i = 0
        rci%j = 0
        if ( keep%leftX > 0 ) then
          iX = keep%firstX
          q = abs(keep%lambda(1))
          r = keep%dlmd(iX, keep%rec)
          s = A_SMALL_FRACTION*REL_GAP*q
          if ( r > ZERO .and. r < s ) then
            t = keep%lambda(1) + q*REL_GAP
            do i = 1, keep%sizeZ
              if ( keep%lambda(keep%leftX + i) <= t ) then
                ! too close to computed eigenvalues, suggest increasing
                ! the room for left eigenpairs
                rci%i = rci%i + 1
              else
                exit
              end if
            end do
          end if
          if ( rci%i == 0 .and. keep%leftX > 1 .and. keep%sizeZ > 0 &
            .and. keep%lambda(keep%leftX) - keep%lambda(1) <= q*REL_GAP/10 &
            ) then
            ! computed eigenvalues too close to each other, suggest increasing
            ! the room for left eigenpairs
            rci%i = 1
          end if
        end if
        if ( keep%rightX > 0 ) then
          iX = keep%firstX + keep%sizeX - 1
          q = abs(keep%lambda(sizeXY))
          r = keep%dlmd(iX, keep%rec)
          s = A_SMALL_FRACTION*REL_GAP*q
          if ( r > ZERO .and. r < s ) then
            t = keep%lambda(sizeXY) - q*REL_GAP
            do i = keep%sizeZ, rci%i + 1, -1
              if ( keep%lambda(keep%leftX + i) >= t ) then
                ! too close to computed eigenvalues, suggest to increase
                ! the room for right eigenpairs
                rci%j = rci%j + 1
              else
                exit
              end if
            end do
          end if
          i = keep%leftX + keep%sizeZ
          if ( rci%j == 0 .and. keep%rightX > 1 .and. keep%sizeZ > 0 &
            .and. keep%lambda(sizeXY) - keep%lambda(i) <= q*REL_GAP/10 &
            ) then
            ! computed eigenvalues too close to each other, suggest increasing
            ! the room for right eigenpairs
            rci%j = 1
          end if
        end if

        keep%step = COMPUTE_BX

        if ( rci%i + rci%j > 0 ) then
          ! suggest increasing the block size
          rci%jx = keep%firstX
          rci%kx = 0
          rci%nx = keep%sizeX
          rci%job = SSMFE_RESTART
          rci%k = 1 ! restart is not mandatory
          return
        end if

      case (QUIT) select_step

        rci%job = SSMFE_FINISH
        if ( info%flag /= 0 ) rci%job = SSMFE_STOP

        return

      end select select_step ! case keep%step

    end do do_select_step

  contains

    real(PRECISION) function conjugate( z )
      real(PRECISION), intent(in) :: z

      conjugate = z

    end function conjugate

  end subroutine ssmfe_engine_double

!**************************************************************************
! deallocates allocated arrays in inform, and restores default values of
! its scalar elements
!
  subroutine ssmfe_free_info_double( info )
    type(ssmfe_inform), intent(inout) :: info

    if ( allocated(info%residual_norms) ) deallocate( info%residual_norms )
    if ( allocated(info%err_lambda    ) ) deallocate( info%err_lambda     )
    if ( allocated(info%err_X         ) ) deallocate( info%err_X          )
    if ( allocated(info%converged     ) ) deallocate( info%converged      )
    info%flag = 0
    info%stat = 0
    info%non_converged = 0
    info%iteration = 0
    info%left = 0
    info%right = 0
    info%next_left = 1
    info%next_right = -1

  end subroutine ssmfe_free_info_double

!**************************************************************************
! deallocates allocated arrays in keep
!
  subroutine ssmfe_core_free_keep_double( keep )
    type(ssmfe_core_keep), intent(inout) :: keep

    if ( allocated(keep%lambda) ) deallocate( keep%lambda )
    if ( allocated(keep%dlmd  ) ) deallocate( keep%dlmd   )
    if ( allocated(keep%q     ) ) deallocate( keep%q      )
    if ( allocated(keep%dX    ) ) deallocate( keep%dX     )
    if ( allocated(keep%dwork ) ) deallocate( keep%dwork   )
    if ( allocated(keep%zwork ) ) deallocate( keep%zwork  )
    if ( allocated(keep%ind   ) ) deallocate( keep%ind    )

  end subroutine ssmfe_core_free_keep_double

  subroutine ssmfe_core_free_double( keep, info )
    type(ssmfe_core_keep), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: info

    call ssmfe_core_free_keep_double( keep )
    call ssmfe_free_info_double( info )

  end subroutine ssmfe_core_free_double

!**************************************************************************
! returns the size of workspace for dsyev and dsygv
!
  integer function lwork_sevp( n )
    integer, intent(in) :: n

    integer :: nb, ilaenv

    nb = ilaenv(1, 'DSYTRD', SSMFE_UPLO, n, -1, -1, -1)
    lwork_sevp = max(3*n - 1, (nb + 2)*n)

  end function lwork_sevp

!**************************************************************************
! interface to dsyev
!
  subroutine solve_sevp( job, n, a, lda, lambda, lwork, work, info )

    ! solves real symmetric eigenvalue problem
    ! A x = lambda x in double precision
    character, intent(in) :: job
    integer, intent(in) :: n, lda, lwork
    real(PRECISION), intent(inout) :: a(lda, n), work(*)
    real(PRECISION), intent(out) :: lambda(n)
    integer, intent(out) :: info

    call dsyev( job, SSMFE_UPLO, n, a, lda, lambda, work, lwork, info )

  end subroutine solve_sevp

!**************************************************************************
! interface to dsygv
!
  subroutine solve_gsevp( n, a, lda, b, ldb, lambda, lwork, work, info )

    ! solves generalized real symmetric eigenvalue problem
    ! A x = lambda B x in double precision
    integer, intent(in) :: n, lda, ldb, lwork
    real(PRECISION), intent(inout) :: a(lda, n), b(ldb, n), work(*)
    real(PRECISION), intent(out) :: lambda(n)
    integer, intent(out) :: info

    call dsygv &
      ( SSMFE_ITYPE, SSMFE_JOBZ, SSMFE_UPLO, &
        n, a, lda, b, ldb, lambda, work, lwork, info )

  end subroutine solve_gsevp

!**************************************************************************
! complex core solver RCI
!
  subroutine ssmfe_engine_double_complex &
      ( problem, left, right, m, lambda, rr_matrices, ind, rci, keep, control, &
        info )
    use spral_blas_iface, &
      copy => zcopy, &
      norm => dznrm2, &
      dot  => zdotc, &
      axpy => zaxpy, &
      scal => zscal, &
      gemm => zgemm, &
      trsm => ztrsm
    use spral_lapack_iface, chol => zpotrf

    real(PRECISION), parameter :: ZERO = 0.0D0, ONE = 1.0D0
    complex(PRECISION), parameter :: NIL = ZERO
    complex(PRECISION), parameter :: UNIT = ONE
    character, parameter :: TRANS = 'C'

    integer, intent(in) :: problem
    integer, intent(in) :: left, right
    integer, intent(in) :: m
    real(PRECISION), dimension(m), intent(inout) :: lambda
    complex(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, *)
    integer, intent(inout), dimension(m) :: ind
    type(ssmfe_core_keep   ), intent(inout) :: keep
    type(ssmfe_rciz        ), intent(inout) :: rci
    type(ssmfe_core_options), intent(in   ) :: control
    type(ssmfe_inform      ), intent(inout) :: info

    integer, parameter :: SSMFE_START              = 0
    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_CHECK_CONVERGENCE  = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_COPY_VECTORS       = 11
    integer, parameter :: SSMFE_COMPUTE_DOTS       = 12
    integer, parameter :: SSMFE_SCALE_VECTORS      = 13
    integer, parameter :: SSMFE_COMPUTE_YMXD       = 14
    integer, parameter :: SSMFE_COMPUTE_XY         = 15
    integer, parameter :: SSMFE_COMPUTE_XQ         = 16
    integer, parameter :: SSMFE_TRANSFORM_X        = 17
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

    integer, parameter :: SSMFE_FINISH = -1
    integer, parameter :: SSMFE_STOP   = -2
    integer, parameter :: SSMFE_ABORT  = -3

    integer, parameter :: INITIAL = 10
    integer, parameter :: CG_LOOP = 20

    integer, parameter :: BEGIN = 0
    integer, parameter :: QUIT  = -1

    integer, parameter :: COMPUTE_BX           = 100
    integer, parameter :: ORTHOG_X_TO_Xc       = 150
    integer, parameter :: COMPUTE_AX           = 200
    integer, parameter :: COMPUTE_INIT_XAX     = 210
    integer, parameter :: COMPUTE_INIT_XBX     = 220
    integer, parameter :: RR_IN_X              = 300
    integer, parameter :: TRANSFORM_X          = 305
    integer, parameter :: COMPUTE_INIT_AX      = 310
    integer, parameter :: COMPUTE_INIT_BX      = 320
    integer, parameter :: COMPUTE_XBX          = 330
    integer, parameter :: COPY_AX_TO_W         = 350
    integer, parameter :: COMPUTE_XAX          = 360
    integer, parameter :: COMPUTE_RESIDUALS    = 400
    integer, parameter :: ORTHOG_R_TO_Xc       = 410
    integer, parameter :: COMPUTE_RR           = 420
    integer, parameter :: COPY_W_TO_Y          = 430
    integer, parameter :: APPLY_B_TO_Y         = 440
    integer, parameter :: ORTHOG_Y_TO_Xc       = 450
    integer, parameter :: ESTIMATE_ERRORS      = 500
    integer, parameter :: CHECK_CONVERGENCE    = 600
    integer, parameter :: SAVE_RIGHT_CONVERGED = 700
    integer, parameter :: APPLY_PRECONDITIONER = 800
    integer, parameter :: COMPUTE_AZ           = 900
    integer, parameter :: COMPUTE_ZAY          = 910
    integer, parameter :: COMPUTE_BY           = 920
    integer, parameter :: COMPUTE_ZBY          = 930
    integer, parameter :: COMPUTE_YBY_DIAG     = 950
    integer, parameter :: CONJUGATE_Y          = 1000
    integer, parameter :: RECOMPUTE_BY         = 1100
    integer, parameter :: UPDATE_BY            = 1110
    integer, parameter :: COPY_W_TO_BYZ        = 1120
    integer, parameter :: APPLY_CONSTRAINTS    = 1200
    integer, parameter :: COMPUTE_YBY          = 1300
    integer, parameter :: SCALE_Y              = 1320
    integer, parameter :: COMPUTE_XBY          = 1400
    integer, parameter :: CLEANUP_Y            = 1500
    integer, parameter :: COMPUTE_AY           = 1600
    integer, parameter :: RR_IN_XY             = 2000
    integer, parameter :: COMPUTE_XAY          = 2100
    integer, parameter :: PREPARE_MATRICES     = 2200
    integer, parameter :: ANALYSE_RR_RESULTS   = 2300
    integer, parameter :: PUT_AZQ_IN_Z         = 3100
    integer, parameter :: ADD_AXQ_TO_Z         = 3200
    integer, parameter :: PUT_AXQ_IN_W         = 3300
    integer, parameter :: ADD_AZQ_TO_W         = 3400
    integer, parameter :: PUT_W_IN_AX          = 3500
    integer, parameter :: PUT_Z_IN_AZ          = 3600
    integer, parameter :: PUT_BZQ_IN_Z         = 4100
    integer, parameter :: ADD_BXQ_TO_Z         = 4200
    integer, parameter :: PUT_BXQ_IN_W         = 4300
    integer, parameter :: ADD_BZQ_TO_W         = 4400
    integer, parameter :: PUT_W_IN_BX          = 4500
    integer, parameter :: PUT_Z_IN_BZ          = 4600
    integer, parameter :: PUT_YQ_IN_Z          = 5100
    integer, parameter :: ADD_XQ_TO_Z          = 5200
    integer, parameter :: PUT_XQ_IN_W          = 5300
    integer, parameter :: ADD_YQ_TO_W          = 5400
    integer, parameter :: PUT_W_IN_X           = 5500
    integer, parameter :: CHECK_THE_GAP        = 6000

    real(PRECISION) :: TOO_SMALL

    real(PRECISION) :: A_SMALL_FRACTION, A_MULTIPLE
    real(PRECISION) :: REL_GAP, MAX_BETA, Q_MAX
    real(PRECISION) :: GRAM_RCOND_MIN

    real(PRECISION), parameter :: NO_VALUE = -1.0

    logical :: skip, doit
    integer :: kY, kZ, kW, kAX, kAYZ, kBX, kBYZ
    integer :: mm, left_cnv, right_cnv, first, last, step, go, nsmall
    integer :: iX, jX, iY, jY, sizeXY
    integer :: i, j, k, l
    real(PRECISION) :: delta, theta
    real(PRECISION) :: q, r, s, t
    complex(PRECISION) :: z

    logical :: minAprod, minBprod
    integer :: err_est

    if ( control%extra_left < 0 .or. control%extra_right < 0 ) then
      info%flag = WRONG_EXTRAS
      rci%job = SSMFE_ABORT
      return
    end if

    if ( control%min_gap < ZERO .or. control%min_gap > ONE ) then
      info%flag = WRONG_MIN_GAP
      rci%job = SSMFE_ABORT
      return
    end if

    if ( control%cf_max < ONE/2 .or. control%cf_max > ONE ) then
      info%flag = WRONG_CF_MAX
      rci%job = SSMFE_ABORT
      return
    end if

    err_est = control%err_est
    select case ( err_est )
    case ( SSMFE_RESIDUAL, SSMFE_KINEMATIC )
    case default
      info%flag = WRONG_ERR_EST
      rci%job = SSMFE_ABORT
      return
    end select

    if ( problem < 0 ) then
      if ( .not. control%minAprod ) then
        info%flag = WRONG_MINPROD
        rci%job   = SSMFE_ABORT
        return
      end if
    end if

    A_SMALL_FRACTION = 0.01D0
    A_MULTIPLE = 100
    GRAM_RCOND_MIN = 1D-4
    MAX_BETA = 100
    TOO_SMALL = A_MULTIPLE*epsilon(ONE)

    REL_GAP = control%min_gap
    q_max = control%cf_max

    mm = 2*m

if_rci: &
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
      .and. (rci%i == 0 .and. rci%j == 0 .or. rci%k == 0) &
      ) then

      if ( m < 1 &
        .or. (left < 0 .or. left > 0 .and. right > 0) .and. m < 2 ) then
        info%flag = WRONG_BLOCK_SIZE
        rci%job = SSMFE_ABORT
        return
      end if

      if ( allocated(keep%ind           ) ) deallocate ( keep%ind            )
      if ( allocated(keep%lambda        ) ) deallocate ( keep%lambda         )
      if ( allocated(keep%q             ) ) deallocate ( keep%q              )
      if ( allocated(keep%dX            ) ) deallocate ( keep%dX             )
      if ( allocated(keep%dlmd          ) ) deallocate ( keep%dlmd           )
      if ( allocated(keep%dwork         ) ) deallocate ( keep%dwork          )
      if ( allocated(keep%zwork         ) ) deallocate ( keep%zwork          )
      if ( allocated(info%err_lambda    ) ) deallocate ( info%err_lambda     )
      if ( allocated(info%err_X         ) ) deallocate ( info%err_X          )
      if ( allocated(info%converged     ) ) deallocate ( info%converged      )
      if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )

      keep%lwork = lwork_hevp( mm )
      allocate &
        ( info%residual_norms(m), info%err_lambda(mm), info%err_X(mm), &
          info%converged(m), keep%lambda(3*m), &
          keep%ind(2*m), keep%q(mm), keep%dX(mm), &
          keep%zwork(keep%lwork), keep%dwork(3*mm - 2), &
          keep%dlmd(m, keep%RECORDS), stat = info%stat )
      if ( info%stat /= 0 ) info%flag = OUT_OF_MEMORY
      if ( info%stat /= 0 ) rci%job = SSMFE_ABORT
      if ( info%stat /= 0 ) return

      info%flag       = 0
      info%iteration  = 0
      info%residual_norms = NO_VALUE
      info%err_lambda = NO_VALUE
      info%err_X      = NO_VALUE
      info%converged  = 0

      keep%problem  = problem
      keep%err_est  = control%err_est
      keep%minAprod = control%minAprod
      keep%minBprod = control%minBprod .and. keep%problem /= 0

      keep%kY = 1
      keep%kZ = 2
      keep%kW = 3
      if ( keep%minAprod ) then
        keep%kAX  = keep%kW + 1
        keep%kAYZ = keep%kAX + 1
      else
        keep%kAX  = keep%kW
        keep%kAYZ = keep%kW
      end if
      if ( keep%minBprod .or. keep%problem < 0 ) then
        keep%kBX  = keep%kAYZ + 1
        keep%kBYZ = keep%kBX + 1
      else
        if ( keep%problem == 0 ) then
          keep%kBX = 0
        else
          keep%kBX = keep%kY
        end if
        keep%kBYZ = keep%kY
      end if

      if ( left >= 0 ) then
        if ( left == 0 ) then
          keep%leftX = 0
        else if ( right == 0 ) then
          keep%leftX = m
        else
          i = left + control%extra_left
          j = right + control%extra_right
          if ( i + j <= m ) then
            keep%leftX = i + int((m - i - j)*i*ONE/(i + j))
          else
            keep%leftX = int(m*i*ONE/(i + j))
          end if
          keep%leftX = max(keep%leftX, 1)
          keep%leftX = min(keep%leftX, m - 1)
        end if
      else
        keep%leftX = m/2
      end if
      keep%rightX = m - keep%leftX

      keep%firstX = 1
      keep%sizeX  = m
      keep%sizeY  = 0
      keep%sizeZ  = 0

      keep%firstXn   = 1
      keep%leftXn    = keep%leftX
      keep%rightXn   = keep%rightX
      keep%sizeXn    = m
      keep%new_left  = 0
      keep%new_right = 0
      if ( rci%job == SSMFE_START ) then
        keep%left_cnv  = 0
        keep%right_cnv = 0
      end if

      keep%cond  = ONE
      keep%err_A = ZERO
      keep%err_B = ZERO

      keep%rec  = 0
      keep%dlmd = ZERO
      keep%q    = ONE
      keep%dX   = ONE

      keep%stage = INITIAL
      keep%step = BEGIN
      keep%iteration = 0

    end if if_rci

    minAprod = keep%minAprod
    minBprod = keep%minBprod
    err_est  = keep%err_est

    kY   = keep%kY
    kZ   = keep%kZ
    kW   = keep%kW
    kAX  = keep%kAX
    kAYZ = keep%kAYZ
    kBX  = keep%kBX
    kBYZ = keep%kBYZ

do_select_step: &
    do

      if ( left <= 0 .and. right <= 0 ) keep%step = QUIT

select_step: &
      select case ( keep%step )

      case (BEGIN) select_step

        keep%step = COMPUTE_BX

      case (COMPUTE_BX) select_step

        if ( keep%stage == CG_LOOP ) then
          keep%step = COMPUTE_AX
        else
          keep%step = ORTHOG_X_TO_Xc
        end if

        if ( .not. (minBprod .and. keep%stage == CG_LOOP &
                    .and. mod(keep%iteration, 20) > 0) &
              .and. problem /= 0 ) then
          rci%job = SSMFE_APPLY_B
          rci%nx = keep%sizeX
          rci%jx = keep%firstX
          rci%kx = 0
          rci%jy = keep%firstX
          rci%ky = kBX
          return
        end if

      case (ORTHOG_X_TO_Xc) select_step

        rci%job = SSMFE_APPLY_CONSTRAINTS
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%jy = keep%firstX
        rci%ky = kBX
        keep%step = COMPUTE_AX
        return

      case (COMPUTE_AX) select_step

        if ( keep%stage /= INITIAL ) then
          keep%step = COMPUTE_XBX
        else
          keep%step = COMPUTE_INIT_XAX
        end if

        if ( .not. (minAprod .and. keep%stage == CG_LOOP &
                    .and. mod(keep%iteration, 20) > 0) ) then
          rci%job = SSMFE_APPLY_A
          rci%nx = keep%sizeX
          rci%jx = keep%firstX

          if ( problem < 0 ) then
            rci%kx = kBX
          else
            rci%kx = 0
          end if

          rci%jy = keep%firstX
          rci%ky = kAX
          return
        end if

      case (COMPUTE_INIT_XAX) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX

        if ( problem < 0 ) then
          rci%kx = kBX
        else
          rci%kx = 0
        end if

        rci%jy = keep%firstX
        rci%ky = kAX
        rci%ny = keep%sizeX
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_INIT_XBX
        return

      case (COMPUTE_INIT_XBX) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%jy = keep%firstX
        if ( problem == 0 ) then
          rci%ky = 0
        else
          rci%ky = kBX
        end if
        rci%ny = keep%sizeX
        rci%i = 1
        rci%j = 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = RR_IN_X
        return

      case (RR_IN_X) select_step

        call solve_ghevp &
          ( keep%sizeX, rr_matrices(1,1,2), mm, rr_matrices, mm, &
            lambda, keep%lwork, keep%zwork, keep%dwork, i )

        if ( i /= 0 ) then
          info%flag = INDEFINITE_B
          rci%job = SSMFE_ABORT
          return
        end if

        if ( left < 0 ) then
          j = 1
          do while ( j < m )
            if ( lambda(j) > 0 ) exit
            j = j + 1
          end do
          keep%leftX = j - 1
          keep%rightX = m - keep%leftX
        end if

        keep%lambda(1 : keep%sizeX) = lambda(1 : keep%sizeX)

        keep%step = TRANSFORM_X

      case (TRANSFORM_X) select_step

        rci%job = SSMFE_TRANSFORM_X
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        rci%ky = kZ
        rci%i = 1
        rci%j = 1
        rci%k = 2
        keep%step = COMPUTE_INIT_AX
        return

      case (COMPUTE_INIT_AX) select_step

        rci%job = SSMFE_TRANSFORM_X
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kAX
        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        rci%ky = kZ
        rci%i = 1
        rci%j = 1
        rci%k = 2
        keep%step = COMPUTE_INIT_BX
        return

      case (COMPUTE_INIT_BX) select_step

        if ( keep%stage == INITIAL ) keep%stage = CG_LOOP
        keep%step = COMPUTE_XBX

        if ( problem /= 0 ) then
          rci%job = SSMFE_TRANSFORM_X
          rci%nx = keep%sizeX
          rci%jx = keep%firstX
          rci%kx = kBX
          rci%ny = keep%sizeX
          rci%jy = keep%firstX
          rci%ky = kZ
          rci%i = 1
          rci%j = 1
          rci%k = 2
          return
        end if

      case (COMPUTE_XBX) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        if ( problem == 0 ) then
          rci%ky = 0
        else
          rci%ky = kBX
        end if
        rci%i = 1
        rci%j = 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COPY_AX_TO_W
        return

      case (COPY_AX_TO_W) select_step

        keep%step = COMPUTE_XAX
        if ( minAprod ) then
          rci%job = SSMFE_COPY_VECTORS
          rci%nx = keep%sizeX
          rci%jx = keep%firstX
          rci%kx = kAX
          rci%jy = keep%firstX
          rci%ky = kW
          rci%i = 0
          return
        end if

      case (COMPUTE_XAX) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0

        if ( problem < 0 ) then
          rci%kx = kBX
        else
          rci%kx = 0
        end if

        rci%ny = keep%sizeX
        rci%jy = keep%firstX
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_RESIDUALS
        return

      case (COMPUTE_RESIDUALS) select_step

        do i = 1, keep%sizeX
          keep%err_B = max(keep%err_B, abs(rr_matrices(i, i, 1) - ONE))
          do j = 1, i - 1
            keep%err_B = max(keep%err_B, abs(rr_matrices(j, i, 1)))
          end do
          iX = keep%firstX + i - 1
          s = real(rr_matrices(i, i, 1), PRECISION)
          t = real(rr_matrices(i, i, 2), PRECISION)

          t = t/s
          if ( keep%rec > 0 ) then
            if ( i > keep%leftXn - keep%new_left &
              .and. i <= keep%leftXn + keep%new_right &
              .or. left < 0 .and. t*lambda(iX) < ZERO &
              ) then
              keep%dlmd(iX, keep%rec) = ZERO
            else
              s = abs(lambda(iX) - t)
              r = ZERO
              do j = 1, keep%sizeX
                r = r + abs(rr_matrices(j, i, 2) - t*rr_matrices(j, i, 1))**2
              end do
              r = max(sqrt(r), abs(t)*epsilon(ONE))
              if ( s > A_MULTIPLE*r ) keep%dlmd(iX, keep%rec) = s
            end if
          end if
          lambda(iX) = t
        end do
        if ( keep%rec > 0 ) then
          do i = 1, keep%sizeX
            iX = keep%firstX + i - 1
            s = keep%dlmd(iX, keep%rec)
            if ( s == ZERO ) cycle
            do l = 1, m
              if ( abs(lambda(l) - lambda(iX)) &
                < 10*epsilon(ONE)*max(abs(lambda(l)), abs(lambda(iX))) ) &
                s = max(s, keep%dlmd(l, keep%rec))
            end do
            keep%dlmd(iX, keep%rec) = s
          end do
        end if

        keep%sizeY = keep%sizeX

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          rr_matrices(iX, m + iX, 1) = lambda(iX)
        end do
        rci%job = SSMFE_COMPUTE_YMXD
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        if ( problem <= 0 ) then
          rci%kx = 0
        else
          rci%kx = kBX
        end if
        rci%jy = keep%firstX
        rci%ky = kW
        rci%i = keep%firstX
        rci%j = m + keep%firstX
        rci%k = 1

        if ( problem < 0 ) then
          keep%step = COPY_W_TO_Y
        else
          keep%step = ORTHOG_R_TO_Xc
        end if

        return

      case (COPY_W_TO_Y) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kW
        rci%jy = 1
        rci%ky = kY
        rci%i = 0
        keep%step = APPLY_B_TO_Y
        return

      case (APPLY_B_TO_Y) select_step

        rci%job = SSMFE_APPLY_B
        rci%nx = keep%sizeX
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kW
        keep%step = ORTHOG_Y_TO_Xc
        return

      case (ORTHOG_Y_TO_Xc) select_step

        rci%job = SSMFE_APPLY_CONSTRAINTS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 3
        keep%step = COMPUTE_RR
        return

      case (ORTHOG_R_TO_Xc) select_step

        rci%job = SSMFE_APPLY_ADJ_CONSTRS
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kW
        rci%i = 1
        rci%j = 1
        rci%k = 3
        keep%step = COMPUTE_RR
        return

      case (COMPUTE_RR) select_step

        rci%nx = keep%sizeX
        rci%ny = keep%sizeX

        if ( problem < 0 ) then
          rci%kx = kY
          rci%jx = 1
          rci%jy = 1
        else
          rci%jx = keep%firstX
          rci%kx = kW
          rci%jy = keep%firstX
        end if

        rci%ky = kW
        rci%i = keep%firstX
        rci%j = keep%firstX
        rci%k = 3
        if ( err_est == SSMFE_RESIDUAL ) then
          rci%alpha = UNIT
          rci%beta = NIL
          rci%job = SSMFE_COMPUTE_XY
        else
          rci%job = SSMFE_COMPUTE_DOTS
        end if
        keep%step = ESTIMATE_ERRORS
        return

      case (ESTIMATE_ERRORS) select_step

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          info%residual_norms(iX) = sqrt(abs(rr_matrices(iX, iX, 3)))
        end do

        keep%step = CHECK_CONVERGENCE

        l = keep%rec

        if ( l > 3 ) then

          do iX = keep%firstX, keep%firstX + keep%sizeX - 1

            if ( keep%dlmd(iX,l) == ZERO ) then
              keep%q(iX) = ONE
              cycle
            end if

            if ( keep%dX(iX) > A_SMALL_FRACTION ) then
              keep%q(iX) = ONE
              cycle
            end if

            s = ZERO
            k = 0
            do j = l, l - l/3, -1
              q = keep%dlmd(iX,j)
              if ( q == ZERO ) exit
              s = s + q
              k = k + 1
            end do
            if ( k < 2 .or. s == ZERO ) cycle
            q = keep%dlmd(iX,l)/s
            q = q**(ONE/(k - 1))
            keep%q(m + iX) = q
            if ( q > ONE - TOO_SMALL ) cycle

            theta = q/(ONE - q)

            info%err_lambda(m + iX) = keep%dlmd(iX,l)*theta
            i = 0
            do j = l - 1, 1, -1
              if ( keep%dlmd(iX,j) == ZERO ) then
                i = j
                exit
              end if
            end do

            k = (l - i)/3
            if ( k < 1 ) cycle
            first = l - 2*k + 1
            last = l - k
            theta = ZERO
            do j = last, first, -1
              s = info%err_lambda(m + iX)
              do k = j + 1, l
                s = s + keep%dlmd(iX,k)
              end do
              theta = theta + log(s) - log(keep%dlmd(iX,j))
            end do
            k = last - first + 1
            theta = theta/k

            r = ZERO
            do j = last, first, -1
              s = info%err_lambda(m + iX)
              do k = j + 1, l
                s = s + keep%dlmd(iX,k)
              end do
              r = r + (theta - log(s) + log(keep%dlmd(iX,j)))**2
            end do
            r = sqrt(r/(last - first + 1))

            theta = exp(theta + 3*r)
            keep%q(iX) = theta/(ONE + theta)

          end do ! iX

        end if ! l > 3

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == ZERO ) then
            info%err_lambda(iX) = ZERO
            cycle
          end if
          l = keep%rec
          q = keep%q(iX)
          info%err_lambda(iX) = NO_VALUE
          if ( q < ONE .and. l > 0 ) then
            if ( keep%dlmd(iX,l) > ZERO ) &
              info%err_lambda(iX) = keep%dlmd(iX,l)*q/(ONE - q)
          end if
        end do

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == ZERO ) then
            info%err_X(iX) = ZERO
            cycle
          end if
          q = sqrt(keep%q(iX))
          if ( q < ONE .and. keep%dX(iX) > ZERO ) then
            info%err_X(iX) = min(ONE, keep%dX(iX)*q/(ONE - q))
          else
            info%err_X(iX) = NO_VALUE
          end if
        end do

        if ( err_est == SSMFE_RESIDUAL ) then

          do iX = keep%firstX, keep%firstX + keep%sizeX - 1
            info%err_lambda(iX) = info%residual_norms(iX)
          end do

          do go = 1, 2

            if ( keep%sizeX < 2 ) exit

            select case ( go )
            case (1)
              if ( keep%leftX == 0 ) cycle
              first = keep%firstX
              last = keep%firstX + keep%leftX - 1
              step = 1
            case (2)
              if ( keep%rightX == 0 ) cycle
              first = keep%firstX + keep%sizeX - 1
              last = keep%firstX + keep%leftX
              step = -1
            end select

            k = last

            t = ZERO
            do while ( (k - first)*step >= 0 )
              if ( step*lambda(k) < ZERO ) exit
              k = k - step
            end do
            if ( (k - first)*step < 0 ) cycle

            doit = left >= 0 .or. last + step > m .or. last + step < 1
            if ( .not. doit ) doit = step*lambda(last + step) < ZERO

            if ( doit ) then

              s = ZERO
              do iX = first, last - step, step
                s = s + info%residual_norms(iX)**2
              end do

              do while ( (k - first)*step > 0 )
                if ( step > 0 ) then
                  r = lambda(maxloc(lambda(1:k-step),1))
                else
                  r = lambda(k + minloc(lambda(k-step:m),1))
                end if
                if ( abs(lambda(k) - r) > &
                      sqrt(s) + info%residual_norms(k) ) then
                  t = lambda(k) - step*info%residual_norms(k)
                  exit
                end if
                k = k - step
                s = max(ZERO, s - info%residual_norms(k)**2)
              end do
              if ( (k - first)*step <= 0 ) cycle
              k = k - step

            end if

            l = (k - first)*step + 1
            do i = min(first, k), max(first, k)
              do j = i, max(first, k)
                rr_matrices(i,j,3) = -step * rr_matrices(i,j,3) &
                  /sqrt((t - lambda(i))*(t - lambda(j)))
              end do
              rr_matrices(i,i,3) = lambda(i) + rr_matrices(i,i,3)
            end do
            i = min(first, k)

            call solve_hevp &
              ( 'N', l, rr_matrices(i,i,3), mm, keep%lambda(mm + i), &
                keep%lwork, keep%zwork, keep%dwork, j )

            do i = first, k, step
              q = info%residual_norms(i)
              s = max(q, step*(t - lambda(i)))
              info%err_X(i) = min(ONE, q/s)
              info%err_lambda(i) = step*(lambda(i) - keep%lambda(mm + i))
              r = max(abs(lambda(first)), abs(lambda(k)))
              if ( info%err_lambda(i) <= (keep%err_B + TOO_SMALL)*r ) then
                q = q*q
                info%err_lambda(i) = q/s
              end if
            end do

          end do

        end if

        do i = 1, m
          info%err_lambda(m + i) = info%err_lambda(i)
          if ( info%err_lambda(i) == NO_VALUE &
              .and. keep%rec > 0  ) then
            if ( keep%dlmd(i, keep%rec) > 0 ) &
              info%err_lambda(m + i) = keep%dlmd(i, keep%rec)
          end if
          if ( problem == 0 ) then
            if ( info%err_lambda(m + i) == NO_VALUE &
              .or. info%err_lambda(m + i) > info%residual_norms(i) ) &
              info%err_lambda(m + i) = info%residual_norms(i)
          end if
          if ( info%err_X(i) == NO_VALUE .and. keep%dX(i) > 0 ) then
            info%err_X(m + i) = keep%dX(i)
          else
            info%err_X(m + i) = abs(info%err_X(i))
          end if
        end do

        rci%job = SSMFE_CHECK_CONVERGENCE
        return

      case (CHECK_CONVERGENCE) select_step

        skip = .false.
        first = 0
        last = 0
        left_cnv = 0
        right_cnv = 0
        r = 0
        s = keep%cond*epsilon(ONE)
        t = max(abs(keep%lambda(1)), abs(keep%lambda(keep%sizeX + keep%sizeZ)))

        do go = 1, -1, -2

          select case ( go )
          case ( 1 )
            skip = left == 0
            first = keep%firstX
            last = first + keep%leftX - 1
          case ( -1 )
            skip = right == 0
            first = keep%firstX + keep%sizeX - 1
            last = first - keep%rightX + 1
          end select

          k = 0
          do i = first, last, go

            if ( skip ) then

              info%converged(i) = 0

            else

              if ( info%converged(i) == 0 .and. keep%rec > 0 &
                .and. keep%sizeZ > 0 ) then

                select case ( go )
                case (1)
                  r = keep%lambda(keep%leftX + 1)
                case (-1)
                  r = keep%lambda(keep%leftX + keep%sizeZ)
                end select
                q = keep%dX(m + i)
                do l = 1, m
                  if ( abs(lambda(i) - lambda(l)) &
                    < 10*epsilon(ONE)*max(abs(lambda(i)), abs(lambda(l))) ) &
                    q = max(q, keep%dX(m + l))
                end do

                if ( keep%rec >= 5 .and. keep%dX(i) <= A_SMALL_FRACTION &
                  .and. (keep%dX(i) > q_max * q &
                  .and. min(keep%dX(i), keep%dX(m + i))*abs(r - lambda(i)) < &
                    10*(s*t + keep%err_A + abs(lambda(i))*keep%err_B) &
                    ) &
                  ) then
                  info%converged(i) = -info%iteration
                  if ( err_est == SSMFE_KINEMATIC ) then
                    q = keep%q(m + i)
                    if ( q < ONE ) then
                      if ( keep%dlmd(i, keep%rec) > ZERO ) &
                        info%err_lambda(i) = keep%dlmd(i, keep%rec)*q/(ONE - q)
                      if ( keep%dX(i) > ZERO ) &
                        info%err_X(i) = min(ONE, keep%dX(i)*q/(ONE - q))
                    end if
                  end if
                end if
              end if

              if ( info%converged(i) /= 0 ) then
                k = k + 1
                if ( info%converged(i) > 0 ) &
                  info%converged(i) = max(1, info%iteration)
              else
                skip = .true.
              end if

            end if

          end do

          select case ( go )
          case ( 1 )
            left_cnv = k
          case ( -1 )
            right_cnv = k
          end select

        end do ! go

        if ( left < 0 ) then
          left_cnv = 0
          right_cnv = 0
          first = keep%firstX
          last = keep%firstX + keep%sizeX - 1
          l = keep%firstX + keep%leftX - 1
          skip = .false.
          do while ( first <= last )
            if ( first <= l &
              .and. (abs(lambda(first)) > abs(lambda(last)) .or. last <= l) &
              ) then
              i = first
              j = keep%left_cnv + i - keep%firstX + 1
              k = 1
              first = first + 1
            else
              i = last
              j = -(keep%right_cnv + keep%firstX + keep%sizeX - i)
              k = -1
              last = last - 1
            end if
            if ( skip ) then
              info%converged(i) = 0
            else if ( info%converged(i) == 0 ) then
              skip = .true.
            else
              if ( k > 0 ) left_cnv = left_cnv + 1
              if ( k < 0 ) right_cnv = right_cnv + 1
            end if
          end do
        end if

        rci%k = 0
        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == 0 ) rci%k = -1
        end do

        keep%left_cnv = keep%left_cnv + left_cnv
        keep%right_cnv = keep%right_cnv + right_cnv
        keep%leftX = keep%leftX - left_cnv
        keep%rightX = keep%rightX - right_cnv
        keep%sizeX = keep%leftX + keep%rightX
        keep%firstX = keep%firstX + left_cnv

        keep%step = SAVE_RIGHT_CONVERGED

        if ( left_cnv > 0 ) then
          do j = 1, keep%sizeX
            call copy &
              ( keep%sizeX, rr_matrices(left_cnv + 1, left_cnv + j, 1), 1, &
                rr_matrices(1,j,1), 1 )
            call copy &
              ( keep%sizeX, rr_matrices(left_cnv + 1, left_cnv + j, 2), 1, &
                rr_matrices(1,j,2), 1 )
          end do
        end if
        rci%job = SSMFE_SAVE_CONVERGED
        rci%nx = left_cnv
        rci%ny = right_cnv
        rci%jx = keep%firstX - left_cnv
        rci%kx = 0
        rci%jy = rci%jx
        rci%ky = kBX
        rci%i = 1
        return

      case (SAVE_RIGHT_CONVERGED)

        keep%step = APPLY_PRECONDITIONER

        rci%job = SSMFE_SAVE_CONVERGED
        rci%nx = rci%ny
        rci%jx = keep%firstX + keep%sizeX + rci%nx - 1
        rci%kx = 0
        rci%jy = rci%jx
        rci%ky = kBX
        rci%i = -1
        return

      case (APPLY_PRECONDITIONER) select_step

        if ( rci%k == -1 ) then
          rci%job = SSMFE_RESTART
          rci%jx = keep%firstX
          rci%nx = keep%sizeX
          rci%k = 0
          return
        end if

        if ( left == 0 .and. keep%leftX > 0 ) then
          keep%firstX = keep%firstX + keep%leftX
          keep%sizeX = keep%rightX
          do j = 1, keep%sizeX
            call copy &
              ( keep%sizeX, rr_matrices(keep%leftX + 1, keep%leftX + j, 1), 1, &
                rr_matrices(1,j,1), 1 )
            call copy &
              ( keep%sizeX, rr_matrices(keep%leftX + 1, keep%leftX + j, 2), 1, &
                rr_matrices(1,j,2), 1 )
          end do
          keep%leftX = 0
        end if

        if ( right == 0 ) keep%rightX = 0

        keep%sizeX = keep%leftX + keep%rightX

        if ( keep%sizeX == 0 ) then
          rci%job = SSMFE_RESTART
          rci%jx = m + 1
          rci%nx = 0
          rci%k = 0
          return
        end if

        if ( keep%sizeZ > 0 ) then
          keep%step = COMPUTE_AZ
        else
          if ( problem < 0 ) then
            keep%step = COPY_W_TO_BYZ
          else
            keep%step = RECOMPUTE_BY
          end if
        end if

        if ( problem >= 0 ) then
          rci%job = SSMFE_APPLY_PREC
          rci%nx = keep%sizeY
          rci%jx = keep%firstXn
          rci%kx = kW
          rci%jy = 1
          rci%ky = kY
          return
        end if

      case (COMPUTE_AZ) select_step

        keep%step = COMPUTE_ZAY

        if ( .not. minAprod ) then
          rci%job = SSMFE_APPLY_A
          rci%nx = keep%sizeZ
          rci%jx = 1
          rci%kx = kZ
          rci%jy = 1
          rci%ky = kAYZ
          return
        end if

      case (COMPUTE_ZAY) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kAYZ
        rci%ny = keep%sizeY
        rci%jy = 1

        if ( problem < 0 ) then
          rci%ky = kW
        else
          rci%ky = kY
        end if

        rci%i = m + 1
        rci%j = m + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_BY
        return

      case (COMPUTE_BY) select_step

        keep%step = COMPUTE_ZBY

        if ( problem == 0 ) then
          kBYZ = kZ
          keep%kBYZ = kBYZ
        else if ( problem > 0 ) then
          rci%job = SSMFE_APPLY_B
          rci%jx = 1
          rci%kx = kY
          rci%nx = keep%sizeY
          rci%jy = 1
          rci%ky = kW
          return
        end if

      case (COMPUTE_ZBY) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%ny = keep%sizeY
        rci%jy = 1
        if ( problem == 0 ) then
          rci%ky = kY
        else
          rci%ky = kW
        end if
        rci%i = m + 1
        rci%j = m + 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        if ( problem >= 0 ) then
          keep%step = COMPUTE_YBY_DIAG
        else
          keep%step = CONJUGATE_Y
        end if
        return

      case (COMPUTE_YBY_DIAG)

        rci%job = SSMFE_COMPUTE_DOTS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        if ( problem == 0 ) then
          rci%ky = kY
        else
          rci%ky = kW
        end if
        rci%i = 1
        rci%j = 1
        rci%k = 3
        keep%step = CONJUGATE_Y
        return

      case (CONJUGATE_Y) select_step

        do i = 1, keep%sizeY

          if ( i > keep%leftXn - keep%new_left &
            .and. i <= keep%leftXn + keep%new_right ) then
            do j = 1, keep%sizeZ
              rr_matrices(m + j, m + i, 2) = NIL
            end do
            cycle
          end if

          iX = keep%firstXn + i - 1
          call axpy &
            ( keep%sizeZ, -lambda(iX)*UNIT, rr_matrices(m + 1, m + i, 1), 1, &
              rr_matrices(m + 1, m + i, 2), 1 )
          if ( problem >= 0 ) then
            k = i
          else
            k = keep%firstXn + i - 1
          end if
          r = sqrt(abs(rr_matrices(k, k, 3)))
          t = MAX_BETA*r
          skip = .false.
          do j = 1, keep%sizeZ
            s = keep%lambda(keep%leftXn + j) - lambda(iX)
            if ( abs(rr_matrices(m + j, m + i, 2)) < t*abs(s) ) then
              z = rr_matrices(m + j, m + i, 2)/s
              rr_matrices(m + j, m + i, 2) = z
              t = sqrt(max(t*t - abs(z)**2, ZERO))
            else
              skip = .true.
              exit
            end if
          end do
          if ( skip ) then
            do j = 1, keep%sizeZ
              rr_matrices(m + j, m + i, 2) = NIL
            end do
          end if

        end do

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kY
        rci%i = m + 1
        rci%j = m + 1
        rci%k = 2
        rci%alpha = -UNIT
        rci%beta = UNIT

        if ( minBprod ) then
          keep%step = UPDATE_BY
        else
          keep%step = RECOMPUTE_BY
        end if

        return

      case (UPDATE_BY) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kBYZ
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kW
        rci%i = m + 1
        rci%j = m + 1
        rci%k = 2
        rci%alpha = -UNIT
        rci%beta = UNIT
        keep%step = COPY_W_TO_BYZ
        return

      case (COPY_W_TO_BYZ) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kW
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 0
        keep%step = APPLY_CONSTRAINTS
        return

      case (RECOMPUTE_BY) select_step

        keep%step = APPLY_CONSTRAINTS

        if ( problem == 0 ) then
          kBYZ = kY
          keep%kBYZ = kBYZ
        else
          if ( .not. minBprod ) keep%kBYZ = kW
          rci%job = SSMFE_APPLY_B
          rci%nx = keep%sizeY
          rci%jx = 1
          rci%kx = kY
          rci%jy = 1
          rci%ky = keep%kBYZ
          return
        end if

      case (APPLY_CONSTRAINTS) select_step

        rci%job = SSMFE_APPLY_CONSTRAINTS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 1
        rci%j = 1
        rci%k = 3
        keep%step = SCALE_Y
        return

      case (SCALE_Y) select_step

        rci%job = SSMFE_SCALE_VECTORS
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 1
        rci%j = 1
        rci%k = 3
        keep%step = COMPUTE_YBY
        return

      case (COMPUTE_YBY) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeX + 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_XBY
        return

      case (COMPUTE_XBY) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 1
        rci%j = keep%sizeX + 1
        rci%k = 1
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = CLEANUP_Y
        return

      case (CLEANUP_Y) select_step

        sizeXY = keep%sizeX + keep%sizeY
        rci%ky = keep%sizeY

        if ( problem /= 0 ) then
          iY = keep%sizeX + 1
          jY = keep%sizeX + keep%sizeY
          r = ZERO
          s = ZERO
          do i = iY, jY
            do j = iY, i - 1
              t = abs(rr_matrices(i,j,1) - conjugate(rr_matrices(j,i,1)))
              r = max(r, t)
              t = t/sqrt(abs(rr_matrices(i,i,1)*rr_matrices(j,j,1)))
              s = max(s, t)
            end do
          end do
          keep%err_B = max(keep%err_B, r)
        end if

        delta = GRAM_RCOND_MIN

        call copy( mm*sizeXY, rr_matrices, 1, rr_matrices(1,1,3), 1 )
        nsmall = 0

        call chol( 'U', sizeXY, rr_matrices(1,1,3), mm, i )

        if ( i /= 0 .and. i <= keep%sizeX ) info%flag = INDEFINITE_B
        if ( info%flag /= 0 ) rci%job = SSMFE_ABORT
        if ( rci%job < 0 ) return

        if ( i /= 0 ) then
          nsmall = 1
        else
          rr_matrices(1 : sizeXY - 1, mm, 2) = ZERO
          rr_matrices(sizeXY, mm, 2) = ONE
          s = ONE
          do step = 1, 3
            call trsm &
              ( 'L', 'U', TRANS, 'N', sizeXY, 1, UNIT, &
                rr_matrices(1,1,3), mm, rr_matrices(1,mm,2), mm )
            t = norm(sizeXY, rr_matrices(1,mm,2), 1)
            q = s**2/t**2
            if ( q <= delta ) exit
            call trsm &
              ( 'L', 'U', 'N', 'N', sizeXY, 1, UNIT, &
                rr_matrices(1,1,3), mm, rr_matrices(1,mm,2), mm )
            s = norm(sizeXY, rr_matrices(1,mm,2), 1)
          end do
          keep%cond = ONE/q
          if ( q <= delta ) nsmall = 1
        end if

        if ( nsmall > 0 ) then

          keep%cond = ONE/delta

          call copy &
            ( mm * keep%sizeY, rr_matrices(1, keep%sizeX + 1, 1), 1, &
              rr_matrices(1, keep%sizeX + 1, 3), 1 )
          call trsm &
            ( 'L', 'U', 'N', 'N', keep%sizeX, keep%sizeY, UNIT, &
              rr_matrices(1, 1, 3), mm, rr_matrices(1, keep%sizeX + 1, 1), mm )

          rr_matrices(1 : sizeXY, mm, 2) = ZERO

          forall ( j = 1 : keep%sizeY ) ind(j) = j
          iX = keep%sizeX
          iY = keep%sizeY
          do while ( iX < sizeXY )

            l = iX + 1
            s = norm(iX, rr_matrices(1,l,3), 1)
            do j = iX + 2, sizeXY
              t = norm(iX, rr_matrices(1,j,3), 1)
              if ( t < s ) then
                l = j
                s = t
              end if
            end do
            k = iX + 1
            if ( l /= k ) then
              call copy( iX, rr_matrices(1,k,3), 1, rr_matrices(k,1,3), mm )
              call copy( iX, rr_matrices(1,l,3), 1, rr_matrices(1,k,3), 1 )
              call copy( iX, rr_matrices(k,1,3), mm, rr_matrices(1,l,3), 1 )
              call copy( iY, rr_matrices(k,k,3), 1, rr_matrices(k,1,3), 1 )
              call copy( iY, rr_matrices(k,l,3), 1, rr_matrices(k,k,3), 1 )
              call copy( iY, rr_matrices(k,1,3), 1, rr_matrices(k,l,3), 1 )
              call copy( iY, rr_matrices(k,k,3), mm, rr_matrices(k,1,3), 1 )
              call copy( iY, rr_matrices(l,k,3), mm, rr_matrices(k,k,3), mm )
              call copy( iY, rr_matrices(k,1,3), 1, rr_matrices(l,k,3), mm )
              j = ind(k - keep%sizeX)
              ind(k - keep%sizeX) = ind(l - keep%sizeX)
              ind(l - keep%sizeX) = j
            end if
            s = norm(iX, rr_matrices(1,k,3), 1)**2
            t = real(rr_matrices(k,k,3), PRECISION) - s
            if ( t <= TOO_SMALL + keep%err_B ) exit
            s = sqrt(t)
            rr_matrices(k,k,3) = s

            if ( iX == keep%sizeX ) rr_matrices(k, mm, 2) = ONE
            q = ONE
            do step = 1, 3
              call trsm &
                ( 'L', 'U', TRANS, 'N', k, 1, UNIT, &
                  rr_matrices(1, 1, 3), mm, rr_matrices(1, mm, 2), mm )
              t = norm(k, rr_matrices(1, mm ,2), 1)
              r = q
              q = ONE/t**2
              if ( q <= delta ) exit
              call trsm &
                ( 'L', 'U', 'N', 'N', k, 1, UNIT, &
                  rr_matrices(1, 1, 3), mm, rr_matrices(1, mm, 2), mm )
              t = norm(k, rr_matrices(1, mm, 2), 1)
              call scal(k, UNIT/t, rr_matrices(1, mm, 2), 1)
              if ( q > 0.9*r ) exit
            end do
            keep%cond = ONE/q
            if ( q <= delta ) exit

            if ( iY > 1 ) &
              call gemm &
                ( TRANS, 'N', 1, iY - 1, iX, &
                  -UNIT, rr_matrices(1, k, 3), mm, &
                  rr_matrices(1, iX + 2, 3), mm, &
                  UNIT, rr_matrices(k, iX + 2, 3), mm )
            forall ( j = 2 : iY ) &
              rr_matrices(k, iX + j, 3) = rr_matrices(k, iX + j, 3)/s
            iX = k
            iY = iY - 1

          end do
          k = iX - keep%sizeX
          iY = keep%sizeY

          if ( k < 1 ) then
            info%flag = NO_SEARCH_DIRECTIONS_LEFT
            keep%step = QUIT
            cycle do_select_step
          end if

          call copy( mm*iX, rr_matrices(1,1,3), 1, rr_matrices, 1 )

          keep%sizeY = k

        else

          call copy( mm*sizeXY, rr_matrices(1,1,3), 1, rr_matrices, 1 )

        end if

        keep%step = COMPUTE_AY
        if ( nsmall > 0 ) then
          rci%job = SSMFE_COPY_VECTORS
          rci%kx = kY
          if ( problem /= 0 ) then
            rci%ky = kBYZ
          else
            rci%ky = kY
          end if
          rci%nx = keep%sizeY
          rci%i = 1
          return
        end if

      case (COMPUTE_AY) select_step

        rci%job = SSMFE_APPLY_A
        rci%nx = keep%sizeY
        rci%jx = 1

        if ( problem < 0 ) then
          rci%kx = kBYZ
        else
          rci%kx = kY
        end if

        rci%jy = 1
        rci%ky = kAYZ
        keep%step = RR_IN_XY
        return

      case (RR_IN_XY) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeY
        rci%jx = 1

        if ( problem < 0 ) then
          rci%kx = kBYZ
        else
          rci%kx = kY
        end if

        rci%ny = keep%sizeY
        rci%jy = 1
        rci%ky = kAYZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeX + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = COMPUTE_XAY
        return

      case (COMPUTE_XAY) select_step

        rci%job = SSMFE_COMPUTE_XY
        rci%nx = keep%sizeX
        rci%jx = keep%firstX

        if ( problem < 0 ) then
          rci%kx = kAX
        else
          rci%kx = 0
        end if

        rci%ny = keep%sizeY
        rci%jy = 1
        if ( problem < 0 ) then
          rci%ky = kBYZ
        else
          rci%ky = kAYZ
        end if
        rci%i = 1
        rci%j = keep%sizeX + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = PREPARE_MATRICES
        return

      case (PREPARE_MATRICES) select_step

        iX = 1
        jX = keep%sizeX
        iY = jX + 1
        jY = jX + keep%sizeY
        sizeXY = jY

        r = ZERO
        s = ZERO
        do i = iY, jY
          do j = iY, i - 1
            t = abs(rr_matrices(i,j,2) - conjugate(rr_matrices(j,i,2)))
            r = max(r, t)
            t = t/sqrt(abs(rr_matrices(i,i,2)*rr_matrices(j,j,2)))
            s = max(s, t)
          end do
        end do
        keep%err_A = max(keep%err_A, r)

        do i = 1, sizeXY
          do j = 1, i - 1
            rr_matrices(i,j,2) = conjugate(rr_matrices(j,i,2))
          end do
        end do

        call trsm &
          ( 'L', 'U', TRANS, 'N', sizeXY, sizeXY, UNIT, &
            rr_matrices, mm, rr_matrices(1,1,2), mm )
        call trsm &
          ( 'R', 'U', 'N', 'N', sizeXY, sizeXY, UNIT, &
            rr_matrices, mm, rr_matrices(1,1,2), mm )
        do j = 1, keep%sizeY
          call copy &
            ( keep%sizeY, rr_matrices(iY, jX + j, 2), 1, &
              rr_matrices(m + 1, j, 1), 1 )
        end do

        if ( keep%rec < keep%RECORDS ) then
          keep%rec = keep%rec + 1
        else
          do j = 1, keep%rec - 1
            do i = 1, m
              keep%dlmd(i,j) = keep%dlmd(i, j + 1)
            end do
          end do
        end if
        l = keep%rec

        keep%step = ANALYSE_RR_RESULTS

        call solve_hevp &
          ( 'V', sizeXY, rr_matrices(1,1,2), mm, keep%lambda, &
            keep%lwork, keep%zwork, keep%dwork, i )

      case (ANALYSE_RR_RESULTS) select_step

        iX = 1
        jX = keep%sizeX
        iY = jX + 1
        jY = jX + keep%sizeY
        sizeXY = jY

        if ( left < 0 ) then

          l = right + control%extra_right - keep%left_cnv - keep%right_cnv
          keep%sizeXn = max(min(m, l, sizeXY), keep%sizeX)

          i = 1
          j = sizeXY
          keep%leftXn = 0
          do k = 1, keep%sizeXn
            if ( keep%lambda(i) > ZERO ) exit
            if ( abs(keep%lambda(i)) > abs(keep%lambda(j)) ) then
              keep%leftXn = keep%leftXn + 1
              i = i + 1
            else
              j = j - 1
            end if
          end do
          keep%rightXn = keep%sizeXn - keep%leftXn

          keep%firstXn = keep%firstX

          if ( keep%sizeX < m ) &
            keep%firstXn = keep%firstXn &
              - (keep%firstX - 1)*(keep%sizeXn - keep%sizeX) &
              /(m - keep%sizeX)

          keep%new_left = keep%leftXn - keep%leftX
          keep%new_right = keep%rightXn - keep%rightX

        else

          i = int(keep%sizeY*left*ONE/(left + right))
          if ( left > 0 .and. keep%sizeY > 1 ) i = max(i, 1)
          if ( right > 0 .and. keep%sizeY > 1 ) i = min(i, keep%sizeY - 1)

          k = keep%firstX + keep%leftX - 1
          l = keep%leftX + i

          if ( left == 0 ) then
            keep%leftXn = 0
            keep%new_left = 0
          else
            j = control%extra_left + left - keep%left_cnv
            keep%leftXn = min(j, k, l)
            keep%leftXn = max(keep%leftXn, keep%leftX)
            keep%new_left = keep%leftXn - keep%leftX
          end if
          keep%firstXn = keep%firstX + keep%leftX - keep%leftXn

          if ( right == 0 ) then
            keep%rightXn = 0
            keep%new_right = 0
          else
            j = control%extra_right + right - keep%right_cnv
            k = m - k
            l = sizeXY - l
            keep%rightXn = min(j, k, l)
            keep%rightXn = max(keep%rightXn, keep%rightX)
            keep%new_right = keep%rightXn - keep%rightX
          end if

          keep%sizeXn = keep%leftXn + keep%rightXn

        end if

        k = keep%firstX - keep%firstXn
        l = keep%firstX + keep%leftX - 1
        if ( k > 0 ) then
          do j = 1, l - k
            lambda(j) = lambda(j + k)
            info%residual_norms(j) = info%residual_norms(j + k)
            info%err_lambda(j) = info%err_lambda(j + k)
            info%err_X(j) = info%err_X(j + k)
            info%converged(j) = info%converged(j + k)
            keep%q(j) = keep%q(j + k)
            keep%q(m + j) = keep%q(m + j + k)
            do i = 1, keep%rec
              keep%dlmd(j,i) = keep%dlmd(j + k, i)
            end do
          end do
        end if
        if ( k >= 0 ) then
          do j = l - k + 1, keep%firstXn + keep%leftXn - 1
            info%err_lambda(j) = NO_VALUE
            info%err_X(j) = NO_VALUE
            info%converged(j) = 0
            keep%q(j) = ONE
            keep%q(m + j) = ONE
            do i = 1, keep%rec
              keep%dlmd(j,i) = 0
            end do
          end do
        end if

        k = keep%firstXn + keep%sizeXn - keep%firstX - keep%sizeX
        l = keep%firstX + keep%leftX
        if ( k > 0 ) then
          do j = m, l + k, -1
            lambda(j) = lambda(j - k)
            info%residual_norms(j) = info%residual_norms(j - k)
            info%err_lambda(j) = info%err_lambda(j - k)
            info%err_X(j) = info%err_X(j - k)
            info%converged(j) = info%converged(j - k)
            keep%q(j) = keep%q(j - k)
            keep%q(m + j) = keep%q(m + j - k)
            do i = 1, keep%rec
              keep%dlmd(j,i) = keep%dlmd(j - k, i)
            end do
          end do
        end if
        if ( k >= 0 ) then
          do j = keep%firstXn + keep%leftXn, l + k - 1
            info%err_lambda(j) = NO_VALUE
            info%err_X(j) = NO_VALUE
            info%converged(j) = 0
            keep%q(j) = ONE
            keep%q(m + j) = ONE
            do i = 1, keep%rec
              keep%dlmd(j,i) = 0
            end do
          end do
        end if

        if ( keep%rightXn > 0 ) then
          k = sizeXY - keep%rightXn
          l = keep%leftXn
          call copy &
            ( (k - l)*mm, rr_matrices(1, l + 1, 2), 1, rr_matrices(1,1,3), 1 )
          do j = 1, keep%rightXn
            call copy &
              ( sizeXY, rr_matrices(1, k + j, 2), 1, &
                rr_matrices(1, l + j, 2), 1 )
          end do
          call copy &
            ( (k - l)*mm, rr_matrices(1,1,3), 1, &
              rr_matrices(1, keep%sizeXn + 1, 2), 1 )
        end if

        k = keep%firstX - keep%firstXn
        l = keep%firstX + keep%leftX - 1
        if ( k >= 0 ) then
          do j = 1, l - k
            keep%dX(m + j) = keep%dX(j + k)
          end do
          do j = l - k + 1, keep%firstXn + keep%leftXn - 1
            keep%dX(m + j) = ONE
          end do
        end if
        k = keep%firstXn + keep%sizeXn - keep%firstX - keep%sizeX
        l = keep%firstX + keep%leftX
        if ( k >= 0 ) then
          do j = m, l + k, -1
            keep%dX(m + j) = keep%dX(j - k)
          end do
          do j = keep%firstXn + keep%leftXn, l + k - 1
            keep%dX(m + j) = ONE
          end do
        end if

        call gemm &
          ( 'N', 'N', keep%sizeY, keep%sizeXn, keep%sizeY, &
            UNIT, rr_matrices(m + 1, 1, 1), mm, &
            rr_matrices(keep%sizeX + 1, 1, 2), mm, &
            NIL, rr_matrices(1,1,3), mm )
        do j = 1, keep%sizeXn
          jX = keep%firstXn + j - 1
          s = norm(keep%sizeY, rr_matrices(keep%sizeX + 1, j, 2), 1)
          z = dot &
            (keep%sizeY, rr_matrices(keep%sizeX + 1, j, 2), 1, &
              rr_matrices(1, j, 3), 1)
          t = abs(z - lambda(jX)*s*s)
          keep%dX(jX) = s
          if ( s < 0.1 ) then
            keep%dlmd(jX,keep%rec) = t
          else
            keep%dlmd(jX,keep%rec) = ZERO
          end if
        end do
        do jX = 1, keep%firstXn - 1
          keep%dX(jX) = ZERO
        end do
        do jX = keep%firstXn + keep%sizeXn, m
          keep%dX(jX) = ZERO
        end do

        call trsm &
          ( 'L', 'U', 'N', 'N', sizeXY, sizeXY, UNIT, &
            rr_matrices, mm, rr_matrices(1,1,2), mm )

        k = keep%sizeXn - keep%sizeX
        keep%sizeZ = keep%sizeY - k

        if ( minAprod ) then
          keep%step = PUT_AZQ_IN_Z
        else if ( minBprod ) then
          keep%step = PUT_BZQ_IN_Z
        else
          keep%step = PUT_YQ_IN_Z
        end if

      case (PUT_AZQ_IN_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kAYZ
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_AXQ_TO_Z
        return

      case (ADD_AXQ_TO_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kAX
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_AXQ_IN_W
        return

      case (PUT_AXQ_IN_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kAX
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_AZQ_TO_W
        return

      case (ADD_AZQ_TO_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kAYZ
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = keep%sizeX + 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_W_IN_AX
        return

      case (PUT_W_IN_AX) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeXn
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = keep%firstXn
        rci%ky = kAX
        rci%i = 0
        keep%step = PUT_Z_IN_AZ
        return

      case (PUT_Z_IN_AZ) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%jy = 1
        rci%ky = kAYZ
        rci%i = 0
        if ( minBprod ) then
          keep%step = PUT_BZQ_IN_Z
        else
          keep%step = PUT_YQ_IN_Z
        end if
        return

      case (PUT_BZQ_IN_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kBYZ
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_BXQ_TO_Z
        return

      case (ADD_BXQ_TO_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kBX
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_BXQ_IN_W
        return

      case (PUT_BXQ_IN_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = kBX
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_BZQ_TO_W
        return

      case (ADD_BZQ_TO_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kBYZ
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = keep%sizeX + 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_W_IN_BX
        return

      case (PUT_W_IN_BX) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeXn
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = keep%firstXn
        rci%ky = kBX
        rci%i = 0
        keep%step = PUT_Z_IN_BZ
        return

      case (PUT_Z_IN_BZ) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeZ
        rci%jx = 1
        rci%kx = kZ
        rci%jy = 1
        rci%ky = kBYZ
        rci%i = 0
        keep%step = PUT_YQ_IN_Z
        return

      case (PUT_YQ_IN_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = keep%sizeX + 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_XQ_TO_Z
        return

      case (ADD_XQ_TO_Z) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeZ
        rci%jy = 1
        rci%ky = kZ
        rci%i = 1
        rci%j = keep%sizeXn + 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_XQ_IN_W
        return

      case (PUT_XQ_IN_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeX
        rci%jx = keep%firstX
        rci%kx = 0
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%i = 1
        rci%j = 1
        rci%k = 2
        rci%alpha = UNIT
        rci%beta = NIL
        keep%step = ADD_YQ_TO_W
        return

      case (ADD_YQ_TO_W) select_step

        rci%job = SSMFE_COMPUTE_XQ
        rci%nx = keep%sizeY
        rci%jx = 1
        rci%kx = kY
        rci%ny = keep%sizeXn
        rci%jy = keep%firstXn
        rci%ky = kW
        rci%k = 2
        rci%i = keep%sizeX + 1
        rci%j = 1
        rci%alpha = UNIT
        rci%beta = UNIT
        keep%step = PUT_W_IN_X
        return

      case (PUT_W_IN_X) select_step

        rci%job = SSMFE_COPY_VECTORS
        rci%nx = keep%sizeXn
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = keep%firstXn
        rci%ky = 0
        rci%i = 0
        keep%step = CHECK_THE_GAP
        return

      case (CHECK_THE_GAP) select_step

        keep%leftX = keep%leftXn
        keep%rightX = keep%rightXn
        keep%sizeX = keep%sizeXn
        keep%firstX = keep%firstXn
        keep%iteration = keep%iteration + 1
        info%iteration = keep%iteration

        sizeXY = keep%sizeX + keep%sizeZ
        rci%i = 0
        rci%j = 0
        if ( keep%leftX > 0 ) then
          iX = keep%firstX
          q = abs(keep%lambda(1))
          r = keep%dlmd(iX, keep%rec)
          s = A_SMALL_FRACTION*REL_GAP*q
          if ( r > ZERO .and. r < s ) then
            t = keep%lambda(1) + q*REL_GAP
            do i = 1, keep%sizeZ
              if ( keep%lambda(keep%leftX + i) <= t ) then
                rci%i = rci%i + 1
              else
                exit
              end if
            end do
          end if
          if ( rci%i == 0 .and. keep%leftX > 1 .and. keep%sizeZ > 0 &
            .and. keep%lambda(keep%leftX) - keep%lambda(1) <= q*REL_GAP/10 &
            ) then
            rci%i = 1
          end if
        end if
        if ( keep%rightX > 0 ) then
          iX = keep%firstX + keep%sizeX - 1
          q = abs(keep%lambda(sizeXY))
          r = keep%dlmd(iX, keep%rec)
          s = A_SMALL_FRACTION*REL_GAP*q
          if ( r > ZERO .and. r < s ) then
            t = keep%lambda(sizeXY) - q*REL_GAP
            do i = keep%sizeZ, rci%i + 1, -1
              if ( keep%lambda(keep%leftX + i) >= t ) then
                rci%j = rci%j + 1
              else
                exit
              end if
            end do
          end if
          i = keep%leftX + keep%sizeZ
          if ( rci%j == 0 .and. keep%rightX > 1 .and. keep%sizeZ > 0 &
            .and. keep%lambda(sizeXY) - keep%lambda(i) <= q*REL_GAP/10 &
            ) then
            rci%j = 1
          end if
        end if

        keep%step = COMPUTE_BX

        if ( rci%i + rci%j > 0 ) then
          rci%jx = keep%firstX
          rci%kx = 0
          rci%nx = keep%sizeX
          rci%job = SSMFE_RESTART
          rci%k = 1
          return
        end if

      case (QUIT) select_step

        rci%job = SSMFE_FINISH
        if ( info%flag /= 0 ) rci%job = SSMFE_STOP

        return

      end select select_step

    end do do_select_step

  contains

    complex(PRECISION) function conjugate( z )
      complex(PRECISION), intent(in) :: z

      conjugate = conjg(z)

    end function conjugate

  end subroutine ssmfe_engine_double_complex

  integer function lwork_hevp( n )
    integer, intent(in) :: n

    integer :: nb, ilaenv

    nb = ilaenv(1, 'ZHETRD', SSMFE_UPLO, n, -1, -1, -1)
    lwork_hevp = max(3*n - 1, (nb + 2)*n)

  end function lwork_hevp

  subroutine solve_hevp( job, n, a, lda, lambda, &
                         lwork, work, rwork, info )
    character, intent(in) :: job
    integer, intent(in) :: n, lda, lwork
    complex(PRECISION), intent(inout) :: a(lda, n), work(*)
    real(PRECISION), intent(out) :: lambda(n), rwork(*)
    integer, intent(out) :: info

    call zheev( job, SSMFE_UPLO, n, a, lda, lambda, work, lwork, rwork, info )

  end subroutine solve_hevp

  subroutine solve_ghevp &
    ( n, a, lda, b, ldb, lambda, lwork, work, rwork, info )
    integer, intent(in) :: n, lda, ldb, lwork
    complex(PRECISION), intent(inout) :: a(lda, n), b(ldb, n), work(*)
    real(PRECISION), intent(out) :: rwork(*), lambda(*)
    integer, intent(out) :: info

    call zhegv &
      ( SSMFE_ITYPE, SSMFE_JOBZ, SSMFE_UPLO, &
        n, a, lda, b, ldb, lambda, work, lwork, rwork, info )

  end subroutine solve_ghevp

end module spral_ssmfe_core

