! COPYRIGHT (c) 2014 The Science and Technology Facilities Council (STFC)
! Original date 1 December 2014, Version 1.0.0
!
! Written by: Evgueni Ovtchinnikov
!
module spral_ssmfe_core

  integer, parameter, private :: PRECISION = kind(1.0D0)
  integer, parameter, private :: SSMFE_ITYPE = 1
  character, parameter, private :: SSMFE_JOBZ = 'V'
  character, parameter, private :: SSMFE_UPLO = 'U'
  
  interface ssmfe
    module procedure ssmfe_double, ssmfe_double_complex
  end interface 

  interface ssmfe_largest
    module procedure ssmfe_largest_double, ssmfe_largest_double_complex
  end interface 

  interface ssmfe_terminate
    module procedure &
      ssmfe_terminate_core_double, &
      ssmfe_delete_work_double, &
      ssmfe_delete_info_double
  end interface 

  private :: lwork_sevp, solve_sevp, solve_gsevp
  private :: lwork_hevp, solve_hevp, solve_ghevp
  
  ! error estimation schemes
  integer, parameter, private :: SSMFE_RESIDUAL  = 1
  integer, parameter, private :: SSMFE_KINEMATIC = 2
  
  type ssmfe_options

    integer :: extra_left  = 0
    integer :: extra_right = 0
    integer :: err_est = SSMFE_KINEMATIC
    logical :: minAprod = .true.
    logical :: minBprod = .true.
    logical :: external_eigsol = .false.
    double precision :: min_gap = 0.0
    double precision :: cf_max = 1.0

  end type ssmfe_options
  
  type ssmfe_rcid
  
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
    integer :: op = 0

    double precision :: alpha, beta
    
    double precision, dimension(:,:), pointer :: x => null()
    double precision, dimension(:,:), pointer :: y => null()
  
  end type ssmfe_rcid
  
  type ssmfe_rciz
  
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
    integer :: op = 0

    complex(PRECISION) :: alpha, beta
    
    complex(PRECISION), dimension(:,:), pointer :: x => null()
    complex(PRECISION), dimension(:,:), pointer :: y => null()
  
  end type ssmfe_rciz
  
  type ssmfe_inform

    integer :: flag      = 0
    integer :: data      = 0
    integer :: iteration  = 0
    integer :: cc = 0, dc = 0

    integer :: left = 0, right = 0
    integer, dimension(:), allocatable :: converged

    double precision :: next_left, next_right

    double precision, dimension(:), allocatable :: residual_norms
    double precision, dimension(:), allocatable :: err_lambda, err_X
    
  end type ssmfe_inform

  ! private subroutines and types

  type ssmfe_keep
  
    private
    
    integer :: problem = 0

    integer :: step, stage, iteration
    integer :: left, right
    integer :: sizeX, sizeY, sizeZ, sizeXn, firstX, firstXn
    integer :: leftX, rightX, leftXn, rightXn, new_left, new_right
    integer :: left_cnv, right_cnv
    integer :: kY, kZ, kW
    integer :: kAX, kAYZ
    integer :: kBX, kBYZ
    integer :: lwork
    integer :: rec = 0, RECORDS = 30
    
    double precision :: cond, err_A, err_B

    double precision, dimension(:,:), allocatable :: dlmd
    double precision, dimension(:), allocatable :: q, dX
    double precision, dimension(:), allocatable :: lambda
    double precision, dimension(:), allocatable :: dwork
    complex(PRECISION), dimension(:), allocatable :: zwork
    
    integer, dimension(:), allocatable :: ind, mask
    
    ! copies of control parameters

    integer :: err_est
    logical :: minAprod, minBprod
    
  end type ssmfe_keep

contains

  subroutine ssmfe_double &
    ( rci, problem, left, right, m, lambda, rr_matrices, ind, &
      keep, options, info )

    implicit none
    
    integer, intent(in) :: problem
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: m
    double precision, dimension(m), intent(inout) :: lambda
    double precision, intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: info
    
    call ssmfe_engine_double &
      ( problem, left, right, m, lambda, rr_matrices, ind, rci, keep, options, &
        info )

  end subroutine ssmfe_double

  subroutine ssmfe_double_complex &
    ( rci, problem, left, right, m, lambda, rr_matrices, ind, &
      keep, options, info )

    implicit none
    
    integer, intent(in) :: problem
    integer, intent(in) :: left
    integer, intent(in) :: right
    integer, intent(in) :: m
    double precision, dimension(m), intent(inout) :: lambda
    double complex, intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: info
    
    call ssmfe_engine_double_complex &
      ( problem, left, right, m, lambda, rr_matrices, ind, rci, keep, options, &
        info )

  end subroutine ssmfe_double_complex

  subroutine ssmfe_largest_double &
    ( rci, problem, nep, m, lambda, rr_matrices, ind, keep, control, info )

    implicit none
    
    integer, intent(in) :: problem
    integer, intent(in) :: nep
    integer, intent(in) :: m
    double precision, dimension(m), intent(inout) :: lambda
    double precision, intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: control
    type(ssmfe_inform ), intent(inout) :: info
    
    call ssmfe_engine_double &
      ( problem, -1, nep, m, lambda, rr_matrices, ind, rci, keep, control, &
        info )

  end subroutine ssmfe_largest_double

  ! double precision real solver engine

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

    implicit none
    
    double precision, parameter :: ZERO = 0.0D0, ONE = 1.0D0
    double precision, parameter :: NIL = ZERO
    double precision, parameter :: UNIT = ONE
    character, parameter :: TRANS = 'T'
    
    ! arguments

    integer, intent(in) :: problem
    integer, intent(in) :: left, right
    integer, intent(in) :: m
    double precision, dimension(m), intent(inout) :: lambda
    double precision, intent(inout) :: rr_matrices(2*m, 2*m, *)
    integer, intent(inout), dimension(m) :: ind
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_options), intent(in   ) :: control
    type(ssmfe_inform ), intent(inout) :: info
    
    ! rci

    integer, parameter :: SSMFE_START = 0

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

    integer, parameter :: SSMFE_SOLVE_EIGENPROBLEM = 30
    integer, parameter :: SSMFE_COLLECT_Q_FOR_XQ   = 31

    integer, parameter :: SSMFE_RESTART = 999

    integer, parameter :: SSMFE_FINISH = -1
    integer, parameter :: SSMFE_STOP   = -2
    integer, parameter :: SSMFE_ABORT  = -3
  
    ! error/warning codes

    integer, parameter :: WRONG_M       = -1
    integer, parameter :: WRONG_IDO     = -2
    integer, parameter :: WRONG_ERR_EST = -3
    integer, parameter :: WRONG_MINPROD = -4
    integer, parameter :: WRONG_EXTRAS  = -5
    integer, parameter :: WRONG_MIN_GAP = -6
    integer, parameter :: WRONG_CF_MAX  = -7
    
    integer, parameter :: OUT_OF_MEMORY           = -100
    integer, parameter :: B_NOT_POSITIVE_DEFINITE = -200

    integer, parameter :: NO_SEARCH_DIRECTIONS_LEFT = 1

    ! steps and stages

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
    integer, parameter :: COMPUTE_BX_NORMS     = 340
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
    integer, parameter :: START_SCALING_Y      = 1310
    integer, parameter :: SCALE_Y              = 1320
    integer, parameter :: SCALE_BY             = 1330
    integer, parameter :: COPY_Y               = 1340
    integer, parameter :: COPY_BY              = 1350
    integer, parameter :: STOP_SCALING_Y       = 1360
    integer, parameter :: COMPUTE_XBY          = 1400
    integer, parameter :: CLEANUP_Y            = 1500
    integer, parameter :: COLLECT_Y            = 1510
    integer, parameter :: COLLECT_BY           = 1520
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
    double precision :: TOO_SMALL

    ! a precaution against overflow
    double precision :: TOO_LARGE

    ! heuristic constants
    double precision :: A_SMALL_FRACTION, A_FRACTION, A_MULTIPLE
    double precision :: REL_GAP, MAX_BETA, Q_MAX
    ! maximal admissible condition number of the B-gram matrix
    double precision :: GRAM_RCOND_MIN
    
    ! flags
    integer, parameter :: NONE = -1
    double precision, parameter :: NO_VALUE = -1.0

    ! locals
    integer :: kY, kZ, kW, kAX, kAYZ, kBX, kBYZ
    integer :: mm, left_cnv, right_cnv, first, last, step, go, nsmall
    integer :: iX, jX, iY, jY, sizeXY
    integer :: i, j, k, l
    double precision :: delta, theta
    double precision :: q, r, s, t
    logical :: skip, doit

    ! shortcuts for the parameters
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
      info%flag = WRONG_MIN_GAP
      rci%job = SSMFE_ABORT
      return
    end if

    A_FRACTION = 0.1D0
    A_SMALL_FRACTION = 0.01D0
    A_MULTIPLE = 100
    GRAM_RCOND_MIN = 1D-4
    MAX_BETA = 100    
    TOO_SMALL = A_MULTIPLE*epsilon(ONE)
    TOO_LARGE = sqrt(huge(ONE))
    
    REL_GAP = control%min_gap
    q_max = control%cf_max
    
    mm = 2*m

if_rci: &
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
      .and. rci%i == 0 .and. rci%j == 0 &
      ) then
    
      if ( m < 1 &
        .or. (left < 0 .or. left > 0 .and. right > 0) .and. m < 2 ) then
        info%flag = WRONG_M
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

      if ( allocated(keep%ind       ) ) deallocate ( keep%ind        )
      if ( allocated(keep%mask      ) ) deallocate ( keep%mask       )
      if ( allocated(keep%lambda    ) ) deallocate ( keep%lambda     )
      if ( allocated(keep%q         ) ) deallocate ( keep%q          )
      if ( allocated(keep%dX        ) ) deallocate ( keep%dX         )
      if ( allocated(keep%dwork     ) ) deallocate ( keep%dwork      )
      if ( allocated(keep%zwork     ) ) deallocate ( keep%zwork      )
      if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
      if ( allocated(info%err_lambda) ) deallocate ( info%err_lambda )
      if ( allocated(info%err_X     ) ) deallocate ( info%err_X      )
      if ( allocated(info%converged ) ) deallocate ( info%converged  )
  
      allocate &
        ( info%residual_norms(m), info%err_lambda(mm), info%err_X(mm), &
          info%converged(m), stat = i )
      allocate ( keep%lambda(3*m), stat = i )
      allocate &
        ( keep%ind(2*m), keep%mask(m), keep%q(mm), keep%dX(mm), stat = i )
      keep%lwork = lwork_sevp( mm )
      allocate( keep%dwork(keep%lwork), stat = i )
      if ( i /= 0 ) then
        info%flag = OUT_OF_MEMORY
        info%data = i
        rci%job = SSMFE_ABORT
        return
      end if

      if ( allocated(keep%dlmd) ) deallocate ( keep%dlmd )
      allocate( keep%dlmd(m, keep%RECORDS), stat = i )
      if ( i /= 0 ) then
        info%flag = OUT_OF_MEMORY
        info%data = i
        rci%job = SSMFE_ABORT
        return
      end if
      
      info%flag       = 0
      info%data      = 0      
      info%iteration  = 0
      info%residual_norms = NO_VALUE
      info%err_lambda = NO_VALUE
      info%err_X      = NO_VALUE
      info%converged  = 0
      info%cc         = 0
      info%dc         = 0

      keep%problem  = problem
      keep%err_est  = control%err_est
      keep%left     = left
      keep%right    = right
      if ( keep%problem < 0 ) then
        if ( .not. (control%minAprod .and. control%minBprod) ) then
          info%flag = WRONG_MINPROD
          rci%job   = SSMFE_ABORT
          return
        end if
        keep%minAprod = .true.
        keep%minBprod = .true.
      else
        keep%minAprod = control%minAprod
        keep%minBprod = control%minBprod .and. keep%problem > 0
      end if

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
      if ( keep%minBprod ) then
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
      
      ! divide the workspace into left and right parts
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

      keep%firstX = 1 ! first non-converged eigenpair index
      keep%sizeX  = m ! # of non-converged eigenpairs
      keep%sizeY  = 0 ! # of current search directions
      keep%sizeZ  = 0 ! # of previous search directions

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
      keep%dlmd = ZERO ! eigenvalue decrements history
      keep%q    = ONE  ! eigenvalue convergence factors
      keep%dX   = ONE  ! eigenvector shifts

      keep%stage = INITIAL
      keep%step = BEGIN
      keep%iteration = 0
      
    end if if_rci

    minAprod = keep%minAprod
    minBprod = keep%minBprod
    err_est  = keep%err_est

    ! indices of blocks of the work array W
    kY   = keep%kY    ! current search directions Y
    kZ   = keep%kZ    ! previous search directions Z
    kW   = keep%kW    ! general purpose workspace
    kAX  = keep%kAX   ! A X
    kAYZ = keep%kAYZ  ! A Y or A Z
    kBX  = keep%kBX   ! B X
    kBYZ = keep%kBYZ  ! B Y or B Z
    
do_select_step: &
    do ! This do loop immediately surrounds the select case construct and
       ! allows any of the cases to reset keep%step and cause the corresponging
       ! case to execute next.

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

        ! arranging the next step
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
            rci%kx = kBX
          else
            rci%kx = 0
          end if

          rci%jy = keep%firstX
          rci%ky = kAX
          return
        end if
        
      case (COMPUTE_INIT_XAX) select_step
      
        ! compute X^T A X

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

      case (RR_IN_X) select_step ! Rayleigh-Ritz procedure in span(X)

        call solve_gsevp &
          ( keep%sizeX, rr_matrices(1,1,2), mm, rr_matrices, mm, &
            lambda, keep%lwork, keep%dwork, i )
                
        if ( i /= 0 ) then
          info%flag = B_NOT_POSITIVE_DEFINITE
          info%data = i - keep%sizeX
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

        rci%job = SSMFE_COLLECT_Q_FOR_XQ
        rci%nx = keep%sizeX
        rci%ny = keep%sizeX
        rci%i = 0
        rci%k = 2
        keep%step = TRANSFORM_X
        return

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
        
        ! compute A X
        
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

        ! compute X^T B X

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

        ! compute X^T A X and new lambda

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
          s = rr_matrices(i, i, 1)
          t = rr_matrices(i, i, 2)
          
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
                r = r + (rr_matrices(j, i, 2) - t*rr_matrices(j, i, 1))**2
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
              
        ! compute residuals

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
          
            if ( info%residual_norms(iX) == ZERO ) then
              keep%q(iX) = 0
              cycle
            end if

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

            ! estimate the asymptotic convergence factor (q)
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

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == ZERO ) then
            ! zero-residual eigenpairs are exact
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
            info%err_lambda(iX) = info%residual_norms(iX) ! default: KW estimate
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

            ! choose Lehmann pole
            k = last
            
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

            else

              t = ZERO
              do while ( (k - first)*step >= 0 )
                if ( step*lambda(k) < ZERO ) exit
                k = k - step
              end do
              if ( (k - first)*step < 0 ) cycle

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
!                keep%lwork, keep%work, keep%rwork, j )

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
          ! if error estimates not yet available, use rough ones for printing
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

        rci%jx = keep%firstX
        if ( problem /= 0 ) then
          if ( keep%minBprod ) then
            rci%kx = kBX
          else
            rci%kx = kY
          end if
        else
          rci%kx = 0
        end if
        if ( problem < 0 ) then
          rci%jy = 1
          rci%ky = kY
        else
          rci%jy = keep%firstX
          rci%ky = kW
        end if
        rci%nx = keep%sizeX
        rci%ny = keep%sizeX
        rci%job = SSMFE_CHECK_CONVERGENCE
        return

      case (CHECK_CONVERGENCE) select_step

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
          skip = .false.
          do i = first, last, go

            if ( skip ) then

              info%converged(i) = 0

            else 

              if ( info%converged(i) == 0 .and. keep%rec > 0 &
                .and. keep%sizeZ > 0 ) then

                select case ( go )
                case ( 1 )
                  j = keep%left_cnv + i - first + 1
                  r = keep%lambda(keep%leftX + 1)
                case ( -1 )
                  j = -(keep%right_cnv + first - i + 1)
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
                    s*t + keep%err_A + abs(lambda(i))*keep%err_B &
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
          case (1)
            left_cnv = k
          case (-1)
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
        
        rci%job = SSMFE_APPLY_PREC
        rci%nx = keep%sizeY
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = 1
        rci%ky = kY
        
        rci%ny = keep%sizeX
        rci%i = keep%firstX
        
        if ( keep%sizeZ > 0 ) then
          keep%step = COMPUTE_AZ
        else
          if ( problem < 0 ) then
            if ( minBprod ) then
              keep%step = COPY_W_TO_BYZ
            else
              keep%step = APPLY_CONSTRAINTS
              kBYZ = kW
              keep%kBYZ = kW
            end if
          else
            keep%step = RECOMPUTE_BY
          end if
        end if
        
        if ( problem >= 0 ) return

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

      case (COMPUTE_ZAY) select_step ! compute Z^T A Y for Y-to-Z conjugation

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
        
      case (CONJUGATE_Y) select_step ! conjugate Y to Z using Jacobi scheme

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
              r = rr_matrices(m + j, m + i, 2)/s
              rr_matrices(m + j, m + i, 2) = r
              t = sqrt(max((t - r)*(t + r), ZERO))
            else
              r = rr_matrices(m + j, m + i, 2)/s
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
      
      ! B-orthogonalize Y to constraints

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

        ! normalize Y

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

        ! compute remaining blocks of the right-hand side Gram matrix of the
        ! Rayleigh-Ritz procedure and compute its spectral condition number

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
        if ( keep%sizeX > 0 ) then
          keep%step = COMPUTE_XBY
        else
          keep%step = CLEANUP_Y
        end if
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

        if ( i /= 0 .and. i <= keep%sizeX ) then
          info%flag = B_NOT_POSITIVE_DEFINITE
          info%data = i
          rci%job = SSMFE_ABORT
          return
        end if
        
        if ( i /= 0 ) then
          nsmall = 1
        else
          rr_matrices(1 : sizeXY - 1, mm, 2) = ZERO
          rr_matrices(sizeXY, mm, 2) = ONE
          s = ONE
          do step = 1, 3
            call trsm &
              ( 'L', 'U', 'T', 'N', sizeXY, 1, UNIT, &
                rr_matrices(1,1,3), mm, rr_matrices(1,mm,2), mm )
            t = norm(sizeXY, rr_matrices(1,mm,2), 1)
            q = s**2/t**2
            if ( q <= delta ) exit
            call trsm &
              ( 'L', 'U', 'N', 'N', sizeXY, 1, UNIT, &
                rr_matrices(1,1,3), mm, rr_matrices(1,mm,2), mm )
            s = norm(sizeXY, rr_matrices(1,mm,2), 1)
          end do
          if ( q <= delta ) nsmall = 1
        end if

        if ( nsmall > 0 ) then
        
          info%cc = info%cc + 1

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
            t = rr_matrices(k,k,3) - s
            if ( t <= TOO_SMALL + keep%err_B ) then
              exit
            end if
            s = sqrt(t)
            rr_matrices(k,k,3) = s
            
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
              if ( q > 0.9*r ) exit
            end do
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
            info%data = left - keep%left_cnv
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

      case (RR_IN_XY) select_step ! apply Rayleigh-Ritz procedure in span[X Y]

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
        if ( keep%sizeX > 0 ) then
          keep%step = COMPUTE_XAY
        else
          keep%step = PREPARE_MATRICES
        end if
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
          ( 'L', 'U', 'T', 'N', sizeXY, sizeXY, UNIT, &
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
        
        if ( control%externaL_eigsol ) then
          rci%job = SSMFE_SOLVE_EIGENPROBLEM
          rci%nx = keep%sizeX
          rci%ny = keep%sizeY
          return
        end if

        call solve_sevp &
          ( 'V', sizeXY, rr_matrices(1,1,2), mm, keep%lambda, &
            keep%lwork, keep%dwork, i )
!            keep%lwork, keep%work, keep%rwork, i )

      case (ANALYSE_RR_RESULTS) select_step

        iX = 1
        jX = keep%sizeX
        iY = jX + 1
        jY = jX + keep%sizeY
        sizeXY = jY

        if ( control%externaL_eigsol ) then
          do i = 1, sizeXY
            keep%lambda(i) = rr_matrices(i,i,3)
          end do
        end if

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
        
        if ( keep%sizeXn == 0 ) then
          keep%step = QUIT
          cycle
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
          t = dot &
            (keep%sizeY, rr_matrices(keep%sizeX + 1, j, 2), 1, &
              rr_matrices(1, j, 3), 1)
          t = abs(t - lambda(jX)*s*s)
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

        if ( minAprod ) then ! update A X and compute A Z
          keep%step = PUT_AZQ_IN_Z
        else if ( minBprod ) then ! update B X and compute B Z
          keep%step = PUT_BZQ_IN_Z
        else
          keep%step = PUT_YQ_IN_Z
        end if
        
        rci%job = SSMFE_COLLECT_Q_FOR_XQ
        rci%nx = keep%sizeX
        rci%ny = keep%sizeY
        rci%jx = keep%firstXn
        rci%i = 1
        rci%j = keep%sizeXn
        rci%k = 2
        return

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
        if ( minBprod ) then ! update B X and compute B Z
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

        ! update X and compute Z

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
            ) rci%j = 1
        end if
        
        keep%step = COMPUTE_BX

        if ( rci%i + rci%j > 0 ) then
          rci%jx = keep%firstX
          rci%kx = 0
          rci%nx = keep%sizeX
          rci%jy = 1
          rci%ky = keep%kZ
          rci%ny = keep%sizeZ
          rci%job = SSMFE_RESTART
          return
        end if

      case (QUIT) select_step

        if ( info%flag == 0 ) then
          rci%job = SSMFE_FINISH
        else
          rci%job = SSMFE_STOP
        end if

        return

      end select select_step ! case keep%step

    end do do_select_step
    
  contains
  
    double precision function conjugate( z )
    
      implicit none
      
      double precision, intent(in) :: z
!      double precision :: conjugate
      
      conjugate = z

    end function conjugate

  end subroutine ssmfe_engine_double

! deallocate auxiliary arrays

  subroutine ssmfe_delete_info_double( info )


    implicit none
    
    type(ssmfe_inform), intent(inout) :: info

    if ( allocated(info%residual_norms) ) deallocate( info%residual_norms )
    if ( allocated(info%err_lambda) ) deallocate( info%err_lambda )
    if ( allocated(info%err_X) ) deallocate( info%err_X )
    if ( allocated(info%converged ) ) deallocate( info%converged  )
    info%flag = 0
    info%data = 0
    info%iteration = 0
    info%left = 0
    info%right = 0
    info%next_left = 1
    info%next_right = -1

  end subroutine ssmfe_delete_info_double

  subroutine ssmfe_delete_work_double( keep )

    implicit none
    
    type(ssmfe_keep), intent(inout) :: keep

    if ( allocated(keep%lambda) ) deallocate( keep%lambda )
    if ( allocated(keep%dlmd  ) ) deallocate( keep%dlmd   )
    if ( allocated(keep%q     ) ) deallocate( keep%q      )
    if ( allocated(keep%dX    ) ) deallocate( keep%dX     )
    if ( allocated(keep%dwork ) ) deallocate( keep%dwork   )
    if ( allocated(keep%zwork ) ) deallocate( keep%zwork  )
    if ( allocated(keep%ind   ) ) deallocate( keep%ind    )
    if ( allocated(keep%mask  ) ) deallocate( keep%mask   )

  end subroutine ssmfe_delete_work_double

  subroutine ssmfe_terminate_core_double( keep, info )

    implicit none
    
    type(ssmfe_keep  ), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: info
    
    call ssmfe_delete_work_double( keep )
    call ssmfe_delete_info_double( info )

  end subroutine ssmfe_terminate_core_double

  integer function lwork_sevp( n )
  
    implicit none
    
    integer, intent(in) :: n
    
    integer :: nb, ilaenv

    nb = ilaenv(1, 'DSYTRD', SSMFE_UPLO, n, -1, -1, -1)
    lwork_sevp = max(3*n - 1, (nb + 2)*n)
    
  end function lwork_sevp

  subroutine solve_sevp( job, n, a, lda, lambda, lwork, work, info )

    ! solves real symmetric eigenvalue problem
    ! A x = lambda x in double precision

    implicit none

    character, intent(in) :: job
    integer, intent(in) :: n, lda, lwork
    double precision, intent(inout) :: a(lda, n), work(*)
    double precision, intent(out) :: lambda(n)
    integer, intent(out) :: info

    call dsyev( job, SSMFE_UPLO, n, a, lda, lambda, work, lwork, info )

  end subroutine solve_sevp

  subroutine solve_gsevp( n, a, lda, b, ldb, lambda, lwork, work, info )

    ! solves generalized real symmetric eigenvalue problem
    ! A x = lambda B x in double precision

    implicit none

    integer, intent(in) :: n, lda, ldb, lwork
    double precision, intent(inout) :: a(lda, n), b(ldb, n), work(*)
    double precision, intent(out) :: lambda(n)
    integer, intent(out) :: info

    call dsygv &
      ( SSMFE_ITYPE, SSMFE_JOBZ, SSMFE_UPLO, &
        n, a, lda, b, ldb, lambda, work, lwork, info )

  end subroutine solve_gsevp
  
  subroutine ssmfe_largest_double_complex &
    ( problem, nep, m, lambda, rr_matrices, ind, rci, keep, control, info )

    implicit none
    
    integer, intent(in) :: problem
    integer, intent(in) :: nep
    integer, intent(in) :: m
    double precision, dimension(m), intent(inout) :: lambda
    complex(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, 3)
    integer, intent(out), dimension(m) :: ind
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: control
    type(ssmfe_inform ), intent(inout) :: info
    
    call ssmfe_engine_double_complex &
      ( problem, -1, nep, m, lambda, rr_matrices, ind, rci, keep, control, &
        info )

  end subroutine ssmfe_largest_double_complex

  ! double complex real solver engine

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

    implicit none
    
    double precision, parameter :: ZERO = 0.0D0, ONE = 1.0D0
    complex(PRECISION), parameter :: NIL = ZERO
    complex(PRECISION), parameter :: UNIT = ONE
    character, parameter :: TRANS = 'C'
    
    ! arguments

    integer, intent(in) :: problem
    integer, intent(in) :: left, right
    integer, intent(in) :: m
    double precision, dimension(m), intent(inout) :: lambda
    complex(PRECISION), intent(inout) :: rr_matrices(2*m, 2*m, *)
    integer, intent(inout), dimension(m) :: ind
    type(ssmfe_keep   ), intent(inout) :: keep
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_options), intent(in   ) :: control
    type(ssmfe_inform ), intent(inout) :: info
    
    ! rci

    integer, parameter :: SSMFE_START = 0

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

    integer, parameter :: SSMFE_SOLVE_EIGENPROBLEM = 30
    integer, parameter :: SSMFE_COLLECT_Q_FOR_XQ   = 31

    integer, parameter :: SSMFE_RESTART = 999

    integer, parameter :: SSMFE_FINISH = -1
    integer, parameter :: SSMFE_STOP   = -2
    integer, parameter :: SSMFE_ABORT  = -3
  
    ! error/warning codes

    integer, parameter :: WRONG_M       = -1
    integer, parameter :: WRONG_IDO     = -2
    integer, parameter :: WRONG_ERR_EST = -3
    integer, parameter :: WRONG_MINPROD = -4
    integer, parameter :: WRONG_EXTRAS  = -5
    integer, parameter :: WRONG_MIN_GAP = -6
    integer, parameter :: WRONG_CF_MAX  = -7
    
    integer, parameter :: OUT_OF_MEMORY           = -100
    integer, parameter :: B_NOT_POSITIVE_DEFINITE = -200

    integer, parameter :: NO_SEARCH_DIRECTIONS_LEFT = 1

    ! steps and stages

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
    integer, parameter :: COMPUTE_BX_NORMS     = 340
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
    integer, parameter :: START_SCALING_Y      = 1310
    integer, parameter :: SCALE_Y              = 1320
    integer, parameter :: SCALE_BY             = 1330
    integer, parameter :: COPY_Y               = 1340
    integer, parameter :: COPY_BY              = 1350
    integer, parameter :: STOP_SCALING_Y       = 1360
    integer, parameter :: COMPUTE_XBY          = 1400
    integer, parameter :: CLEANUP_Y            = 1500
    integer, parameter :: COLLECT_Y            = 1510
    integer, parameter :: COLLECT_BY           = 1520
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
    double precision :: TOO_SMALL

    ! a precaution against overflow
    double precision :: TOO_LARGE

    ! heuristic constants
    double precision :: A_SMALL_FRACTION, A_FRACTION, A_MULTIPLE
    double precision :: REL_GAP, MAX_BETA, Q_MAX
    ! maximal admissible condition number of the B-gram matrix
    double precision :: GRAM_RCOND_MIN
    
    ! flags
    integer, parameter :: NONE = -1
    double precision, parameter :: NO_VALUE = -1.0

    ! locals
    logical :: skip, doit
    integer :: kY, kZ, kW, kAX, kAYZ, kBX, kBYZ
    integer :: mm, left_cnv, right_cnv, first, last, step, go, nsmall
    integer :: iX, jX, iY, jY, sizeXY
    integer :: i, j, k, l
    double precision :: delta, theta
    double precision :: q, r, s, t
    complex(PRECISION) :: z

    ! shortcuts for the parameters
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
      info%flag = WRONG_MIN_GAP
      rci%job = SSMFE_ABORT
      return
    end if

    A_FRACTION = 0.1D0
    A_SMALL_FRACTION = 0.01D0
    A_MULTIPLE = 100
    GRAM_RCOND_MIN = 1D-4
    MAX_BETA = 100    
    TOO_SMALL = A_MULTIPLE*epsilon(ONE)
    TOO_LARGE = sqrt(huge(ONE))
    
    REL_GAP = control%min_gap
    q_max = control%cf_max
    
    mm = 2*m

if_rci: &
    if ( rci%job == SSMFE_START .or. rci%job == SSMFE_RESTART &
      .and. rci%i == 0 .and. rci%j == 0 &
      ) then
    
      if ( m < 1 &
        .or. (left < 0 .or. left > 0 .and. right > 0) .and. m < 2 ) then
        info%flag = WRONG_M
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

      if ( allocated(keep%ind       ) ) deallocate ( keep%ind        )
      if ( allocated(keep%mask      ) ) deallocate ( keep%mask       )
      if ( allocated(keep%lambda    ) ) deallocate ( keep%lambda     )
      if ( allocated(keep%q         ) ) deallocate ( keep%q          )
      if ( allocated(keep%dX        ) ) deallocate ( keep%dX         )
      if ( allocated(keep%dwork     ) ) deallocate ( keep%dwork      )
      if ( allocated(keep%zwork     ) ) deallocate ( keep%zwork      )
      if ( allocated(info%residual_norms) ) deallocate ( info%residual_norms )
      if ( allocated(info%err_lambda) ) deallocate ( info%err_lambda )
      if ( allocated(info%err_X     ) ) deallocate ( info%err_X      )
      if ( allocated(info%converged ) ) deallocate ( info%converged  )
  
      allocate &
        ( info%residual_norms(m), info%err_lambda(mm), info%err_X(mm), &
          info%converged(m), stat = i )
      allocate ( keep%lambda(3*m), stat = i )
      allocate &
        ( keep%ind(2*m), keep%mask(m), keep%q(mm), keep%dX(mm), stat = i )
      keep%lwork = lwork_hevp( mm )
      allocate( keep%zwork(keep%lwork), keep%dwork(3*mm - 2), stat = i )
      if ( i /= 0 ) then
        info%flag = OUT_OF_MEMORY
        info%data = i
        rci%job = SSMFE_ABORT
        return
      end if

      if ( allocated(keep%dlmd) ) deallocate ( keep%dlmd )
      allocate( keep%dlmd(m, keep%RECORDS), stat = i )
      if ( i /= 0 ) then
        info%flag = OUT_OF_MEMORY
        info%data = i
        rci%job = SSMFE_ABORT
        return
      end if
      
      info%flag       = 0
      info%data       = 0      
      info%iteration  = 0
      info%residual_norms = NO_VALUE
      info%err_lambda = NO_VALUE
      info%err_X      = NO_VALUE
      info%converged  = 0
      info%cc         = 0
      info%dc         = 0

      keep%problem  = problem
      keep%err_est  = control%err_est
      keep%left     = left
      keep%right    = right
      if ( keep%problem < 0 ) then
        if ( .not. (control%minAprod .and. control%minBprod) ) then
          info%flag = WRONG_MINPROD
          rci%job   = SSMFE_ABORT
          return
        end if
        keep%minAprod = .true.
        keep%minBprod = .true.
      else
        keep%minAprod = control%minAprod
        keep%minBprod = control%minBprod .and. keep%problem > 0
      end if

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
      if ( keep%minBprod ) then
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
      
      ! divide the workspace into left and right parts
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

      keep%firstX = 1 ! first non-converged eigenpair index
      keep%sizeX  = m ! # of non-converged eigenpairs
      keep%sizeY  = 0 ! # of current search directions
      keep%sizeZ  = 0 ! # of previous search directions

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
      keep%dlmd = ZERO ! eigenvalue decrements history
      keep%q    = ONE  ! eigenvalue convergence factors
      keep%dX   = ONE  ! eigenvector shifts

      keep%stage = INITIAL
      keep%step = BEGIN
      keep%iteration = 0
      
    end if if_rci

    minAprod = keep%minAprod
    minBprod = keep%minBprod
    err_est  = keep%err_est

    ! indices of blocks of the work array W
    kY   = keep%kY    ! current search directions Y
    kZ   = keep%kZ    ! previous search directions Z
    kW   = keep%kW    ! general purpose workspace
    kAX  = keep%kAX   ! A X
    kAYZ = keep%kAYZ  ! A Y or A Z
    kBX  = keep%kBX   ! B X
    kBYZ = keep%kBYZ  ! B Y or B Z
    
do_select_step: &
    do ! This do loop immediately surrounds the select case construct and
       ! allows any of the cases to reset keep%step and cause the corresponging
       ! case to execute next.

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

        ! arranging the next step
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
            rci%kx = kBX
          else
            rci%kx = 0
          end if

          rci%jy = keep%firstX
          rci%ky = kAX
          return
        end if
        
      case (COMPUTE_INIT_XAX) select_step
      
        ! compute X^T A X

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

      case (RR_IN_X) select_step ! Rayleigh-Ritz procedure in span(X)

        call solve_ghevp &
          ( keep%sizeX, rr_matrices(1,1,2), mm, rr_matrices, mm, &
            lambda, keep%lwork, keep%zwork, keep%dwork, i )
                
        if ( i /= 0 ) then
          info%flag = B_NOT_POSITIVE_DEFINITE
          info%data = i - keep%sizeX
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

        rci%job = SSMFE_COLLECT_Q_FOR_XQ
        rci%nx = keep%sizeX
        rci%ny = keep%sizeX
        rci%i = 0
        rci%k = 2
        keep%step = TRANSFORM_X
        return

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
        
        ! compute A X
        
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

        ! compute X^T B X

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

        ! compute X^T A X and new lambda

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
              
        ! compute residuals

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
          
            if ( info%residual_norms(iX) == ZERO ) then
              keep%q(iX) = 0
              cycle
            end if

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

            ! estimate the asymptotic convergence factor (q)
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

        do iX = keep%firstX, keep%firstX + keep%sizeX - 1
          if ( info%residual_norms(iX) == ZERO ) then
            ! zero-residual eigenpairs are exact
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
            info%err_lambda(iX) = info%residual_norms(iX) ! default: KW estimate
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

            ! choose Lehmann pole
            k = last
            
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

            else

              t = ZERO
              do while ( (k - first)*step >= 0 )
                if ( step*lambda(k) < ZERO ) exit
                k = k - step
              end do
              if ( (k - first)*step < 0 ) cycle

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

            call solve_hevp &
              ( 'N', l, rr_matrices(i,i,3), mm, keep%lambda(mm + i), &
                keep%lwork, keep%zwork, keep%dwork, j )

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
          ! if error estimates not yet available, use rough ones for printing
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

        rci%jx = keep%firstX
        if ( problem /= 0 ) then
          if ( keep%minBprod ) then
            rci%kx = kBX
          else
            rci%kx = kY
          end if
        else
          rci%kx = 0
        end if
        if ( problem < 0 ) then
          rci%jy = 1
          rci%ky = kY
        else
          rci%jy = keep%firstX
          rci%ky = kW
        end if
        rci%nx = keep%sizeX
        rci%ny = keep%sizeX
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
                  j = keep%left_cnv + i - first + 1
                  r = keep%lambda(keep%leftX + 1)
                case (-1)
                  j = -(keep%right_cnv + first - i + 1)
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
                    s*t + keep%err_A + abs(lambda(i))*keep%err_B &
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
        
        rci%job = SSMFE_APPLY_PREC
        rci%nx = keep%sizeY
        rci%jx = keep%firstXn
        rci%kx = kW
        rci%jy = 1
        rci%ky = kY
        
        rci%ny = keep%sizeX
        rci%i = keep%firstX
        
        if ( keep%sizeZ > 0 ) then
          keep%step = COMPUTE_AZ
        else
          if ( problem < 0 ) then
            if ( minBprod ) then
              keep%step = COPY_W_TO_BYZ
            else
              keep%step = APPLY_CONSTRAINTS
              kBYZ = kW
              keep%kBYZ = kW
            end if
          else
            keep%step = RECOMPUTE_BY
          end if
        end if
        
        if ( problem >= 0 ) return

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

      case (COMPUTE_ZAY) select_step ! compute Z^T A Y for Y-to-Z conjugation

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
        
      case (CONJUGATE_Y) select_step ! conjugate Y to Z using Jacobi scheme

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
!              r = rr_matrices(m + j, m + i, 2)/s
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
      
      ! B-orthogonalize Y to constraints

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

        ! normalize Y

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

        ! compute remaining blocks of the right-hand side Gram matrix of the
        ! Rayleigh-Ritz procedure and compute its spectral condition number

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
        if ( keep%sizeX > 0 ) then
          keep%step = COMPUTE_XBY
        else
          keep%step = CLEANUP_Y
        end if
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

        if ( i /= 0 .and. i <= keep%sizeX ) then
          info%flag = B_NOT_POSITIVE_DEFINITE
          info%data = i
          rci%job = SSMFE_ABORT
          return
        end if
        
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
          if ( q <= delta ) nsmall = 1
        end if

        if ( nsmall > 0 ) then
        
          info%cc = info%cc + 1

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
            info%data = left - keep%left_cnv
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

      case (RR_IN_XY) select_step ! apply Rayleigh-Ritz procedure in span[X Y]

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
        if ( keep%sizeX > 0 ) then
          keep%step = COMPUTE_XAY
        else
          keep%step = PREPARE_MATRICES
        end if
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
        
        if ( control%externaL_eigsol ) then
          rci%job = SSMFE_SOLVE_EIGENPROBLEM
          rci%nx = keep%sizeX
          rci%ny = keep%sizeY
          return
        end if

        call solve_hevp &
          ( 'V', sizeXY, rr_matrices(1,1,2), mm, keep%lambda, &
            keep%lwork, keep%zwork, keep%dwork, i )

      case (ANALYSE_RR_RESULTS) select_step

        iX = 1
        jX = keep%sizeX
        iY = jX + 1
        jY = jX + keep%sizeY
        sizeXY = jY

        if ( control%externaL_eigsol ) then
          do i = 1, sizeXY
            keep%lambda(i) = real(rr_matrices(i,i,3), PRECISION)
          end do
        end if

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
        
        if ( keep%sizeXn == 0 ) then
          keep%step = QUIT
          cycle
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

        if ( minAprod ) then ! update A X and compute A Z
          keep%step = PUT_AZQ_IN_Z
        else if ( minBprod ) then ! update B X and compute B Z
          keep%step = PUT_BZQ_IN_Z
        else
          keep%step = PUT_YQ_IN_Z
        end if
        
        rci%job = SSMFE_COLLECT_Q_FOR_XQ
        rci%nx = keep%sizeX
        rci%ny = keep%sizeY
        rci%jx = keep%firstXn
        rci%i = 1
        rci%j = keep%sizeXn
        rci%k = 2
        return

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
        if ( minBprod ) then ! update B X and compute B Z
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

        ! update X and compute Z

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
            ) rci%j = 1
        end if
        
        keep%step = COMPUTE_BX

        if ( rci%i + rci%j > 0 ) then
          rci%jx = keep%firstX
          rci%kx = 0
          rci%nx = keep%sizeX
          rci%jy = 1
          rci%ky = keep%kZ
          rci%ny = keep%sizeZ
          rci%job = SSMFE_RESTART
          return
        end if

      case (QUIT) select_step

        if ( info%flag == 0 ) then
          rci%job = SSMFE_FINISH
        else
          rci%job = SSMFE_STOP
        end if

        return

      end select select_step ! case keep%step

    end do do_select_step
    
  contains
  
    complex(PRECISION) function conjugate( z )
    
      implicit none
      
      complex(PRECISION), intent(in) :: z
!      complex(PRECISION) :: conjugate
      
      conjugate = conjg(z)

    end function conjugate

  end subroutine ssmfe_engine_double_complex

! deallocate auxiliary arrays

  subroutine ssmfe_delete_info_double_complex( info )

    implicit none
    
    type(ssmfe_inform), intent(inout) :: info

    if ( allocated(info%residual_norms) ) deallocate( info%residual_norms )
    if ( allocated(info%err_lambda) ) deallocate( info%err_lambda )
    if ( allocated(info%err_X) ) deallocate( info%err_X )
    if ( allocated(info%converged ) ) deallocate( info%converged  )
    info%flag = 0
    info%data = 0
    info%iteration = 0
    info%left = 0
    info%right = 0
    info%next_left = 1
    info%next_right = -1

  end subroutine ssmfe_delete_info_double_complex

  subroutine ssmfe_delete_work_double_complex( keep )

    implicit none
    
    type(ssmfe_keep), intent(inout) :: keep

    if ( allocated(keep%lambda) ) deallocate( keep%lambda )
    if ( allocated(keep%dlmd  ) ) deallocate( keep%dlmd   )
    if ( allocated(keep%q     ) ) deallocate( keep%q      )
    if ( allocated(keep%dX    ) ) deallocate( keep%dX     )
    if ( allocated(keep%dwork ) ) deallocate( keep%dwork  )
    if ( allocated(keep%zwork ) ) deallocate( keep%zwork  )
    if ( allocated(keep%ind   ) ) deallocate( keep%ind    )
    if ( allocated(keep%mask  ) ) deallocate( keep%mask   )

  end subroutine ssmfe_delete_work_double_complex

  subroutine ssmfe_terminate_core_double_complex( keep, info )

    implicit none
    
    type(ssmfe_keep  ), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: info
    
    call ssmfe_delete_work_double_complex( keep )
    call ssmfe_delete_info_double_complex( info )

  end subroutine ssmfe_terminate_core_double_complex

  integer function lwork_hevp( n )
  
    implicit none
    
    integer, intent(in) :: n
    
    integer :: nb, ilaenv

    nb = ilaenv(1, 'ZHETRD', SSMFE_UPLO, n, -1, -1, -1)
    lwork_hevp = max(3*n - 1, (nb + 2)*n)
    
  end function lwork_hevp

  subroutine solve_hevp( job, n, a, lda, lambda, &
                         lwork, work, rwork, info )

    ! solves real symmetric eigenvalue problem
    ! A x = lambda x in double complex

    implicit none

    character, intent(in) :: job
    integer, intent(in) :: n, lda, lwork
    complex(PRECISION), intent(inout) :: a(lda, n), work(*)
    double precision, intent(out) :: lambda(n), rwork(*)
    integer, intent(out) :: info

    call zheev( job, SSMFE_UPLO, n, a, lda, lambda, work, lwork, rwork, info )

  end subroutine solve_hevp

  subroutine solve_ghevp &
    ( n, a, lda, b, ldb, lambda, lwork, work, rwork, info )

    ! solves generalized real symmetric eigenvalue problem
    ! A x = lambda B x in double complex

    implicit none

    integer, intent(in) :: n, lda, ldb, lwork
    complex(PRECISION), intent(inout) :: a(lda, n), b(ldb, n), work(*)
    double precision, intent(out) :: rwork(*), lambda(n)
    integer, intent(out) :: info

    call zhegv &
      ( SSMFE_ITYPE, SSMFE_JOBZ, SSMFE_UPLO, &
        n, a, lda, b, ldb, lambda, work, lwork, rwork, info )

  end subroutine solve_ghevp
  
end module spral_ssmfe_core

