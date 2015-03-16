! COPYRIGHT (c) 2014 The Science and Technology Facilities Council (STFC)
! Original date 1 December 2014, Version 1.0.0
!
! Written by: Evgueni Ovtchinnikov
!
module SPRAL_ssmfe

  use SPRAL_ssmfe_expert
  
  integer, parameter, private :: PRECISION = kind(1.0D0)

  type ssmfe_keepd
  
    private
    
    integer :: block_size = 0
    integer :: lcon = 0
    integer :: rcon = 0
    integer :: step = 0
    
    integer, dimension(:), allocatable :: ind

    real(PRECISION), dimension(:,:  ), allocatable :: BX
    real(PRECISION), dimension(:,:  ), allocatable :: U
    real(PRECISION), dimension(:,:,:), allocatable :: V
    real(PRECISION), dimension(:,:,:), allocatable :: W

    type(ssmfe_keep) :: keep

  end type ssmfe_keepd
  
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

    type(ssmfe_keep) :: keep

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
  
  interface ssmfe_terminate
    module procedure &
      ssmfe_terminate_simple_double, &
      ssmfe_delete_keep_simple_double, &
      ssmfe_terminate_simple_double_complex, &
      ssmfe_delete_keep_simple_double_complex
  end interface 

contains

  subroutine ssmfe_std_double &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    real(PRECISION), intent(inout) :: X(ldX, mep)
    integer, intent(in) :: ldX
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keepd  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_srci_double &
      ( 0, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_std_double

  subroutine ssmfe_gen_double &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    real(PRECISION), intent(inout) :: X(ldX, mep)
    integer, intent(in) :: ldX
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_keepd  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_srci_double &
      ( 1, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_gen_double

  subroutine ssmfe_direct_srci_double &
      ( problem, left, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, &
      copy => dcopy, &
      norm => dnrm2, &
      dot  => ddot , &
      axpy => daxpy, &
      gemm => dgemm

    use SPRAL_lapack_iface, mxcopy => dlacpy

    implicit none

    character, parameter :: TR = 'T'

    integer, parameter :: SSMFE_START              = 0
    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_CHECK_CONVERGENCE  = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

    real(PRECISION), parameter :: ZERO = 0.0D0
    real(PRECISION), parameter :: ONE = 1.0D0

    integer, intent(in) :: problem
    integer, intent(in) :: left
    integer, intent(in) :: max_nep
    real(PRECISION), intent(inout) :: lambda(max_nep)
    integer, intent(in) :: n
    integer, intent(in) :: ldX
    real(PRECISION), intent(inout) :: X(ldX, max_nep)
    type(ssmfe_rcid   ), intent(inout) :: rci
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    type(ssmfe_keepd  ), intent(inout), target :: keep

    integer :: kw
    integer :: ldBX
    integer :: extra
    integer :: i, j, k
    
    if ( rci%job == 0 ) then
      extra = options%extra_left
      if ( extra < 0 ) then
        if ( left == 1 ) then
          extra = 7
        else
          extra = max(left/10, 10)
        end if
      else if ( extra == 0 ) then
        if ( left == 1 ) extra = 1
      end if
      extra = min(extra, n/2 - left - 1)
      keep%block_size = left + extra
      keep%lcon = 0
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate ( keep%ind(keep%block_size) )
      allocate ( keep%U(max_nep, keep%block_size) )
      allocate ( keep%V(2*keep%block_size, 2*keep%block_size, 3) )
      if ( problem == 0 ) then
        kw = 5
      else
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep) )
        kw = 7
      end if
      allocate ( keep%W(n, keep%block_size, 0:kw) )
      call random_number( keep%W(:,:,0) )
      if ( options%user_X > 0 ) &
        call mxcopy( 'A', n, options%user_X, X, ldx, keep%W, n )
    end if

    ldBX = n
    
    call ssmfe_solve &
      ( problem, left, max_nep, lambda, keep%block_size, keep%V, keep%ind, &
        rci, keep%keep, options, inform )

    select case ( rci%job )

    case ( 11:19 )

      ! apply some standard vector operations to vectors in W
      call ssmfe_vec_ops_double &
        ( n, keep%block_size, keep%W(1,1,0), n, keep%W(1,1,1), n, &
          keep%V, rci, keep%ind, keep%U )

    case ( SSMFE_SAVE_CONVERGED )
    
      if ( rci%i < 0 ) return
      do i = 1, rci%nx
        k = (i - 1)*rci%i
        j = keep%lcon + i
        call copy( n, keep%W(1, rci%jx + k, 0), 1, X(1,j), 1 )
        if ( problem /= 0 ) &
          call copy( n, keep%W(1, rci%jy + k, rci%ky), 1, keep%BX(1,j), 1 )
      end do
      keep%lcon = keep%lcon + rci%nx

    case ( SSMFE_APPLY_ADJ_CONSTRS )
    
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

      rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
      rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
      if ( keep%lcon > 0 ) then
        call gemm &
          ( TR, 'N', keep%lcon, rci%nx, n, &
            ONE, X, ldX, rci%y, n, ZERO, keep%U, max_nep )
        call gemm &
          ( 'N', 'N', n, rci%nx, keep%lcon, &
            -ONE, X, ldX, keep%U, max_nep, ONE, rci%x, n )
        if ( problem /= 0 ) &
          call gemm &
            ( 'N', 'N', n, rci%nx, keep%lcon, &
              -ONE, keep%BX, ldBX, keep%U, max_nep, ONE, rci%y, n )
      end if

    case ( SSMFE_APPLY_A, SSMFE_APPLY_B )

      rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
      rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

    case ( SSMFE_APPLY_PREC )

      rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
      rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
      call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

    case ( SSMFE_RESTART )
    
      if ( rci%k == 0 ) then
        if ( rci%jx > 1 ) call random_number( keep%W(:, 1 : rci%jx - 1, 0) )
        if ( rci%jx + rci%nx - 1 < keep%block_size ) &
          call random_number( keep%W(:, rci%jx + rci%nx : keep%block_size, 0) )
      end if
        
    end select
      
  end subroutine ssmfe_direct_srci_double
  
  subroutine ssmfe_shift_double &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )

    implicit none

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

  subroutine ssmfe_gen_shift_double &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )

    implicit none

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

  subroutine ssmfe_buckling_double &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )

    implicit none

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

  subroutine ssmfe_inverse_srci_double &
      ( problem, sigma, left, right, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, &
      copy => dcopy, &
      norm => dnrm2, &
      dot  => ddot, &
      axpy => daxpy, &
      gemm => dgemm
    use SPRAL_lapack_iface, mxcopy => dlacpy

    implicit none

    integer, parameter :: SSMFE_START              = 0
    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_CHECK_CONVERGENCE  = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE   = 9
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

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
    integer :: kw
    integer :: ldV
    integer :: ldBX
    integer :: i, j, k, m
    
    real(PRECISION) :: s
    
    if ( rci%job == SSMFE_START ) then
      keep%step = 0
      keep%lcon = 0
      keep%rcon = 0
      if ( options%extra_left < 0 ) then
        extra_left = max(10, left/10)
      else
        extra_left = options%extra_left
      end if
      if ( options%max_left >= 0 ) &
        extra_left = max(0, min(extra_left, options%max_left - left))
      if ( options%extra_right < 0 ) then
        extra_right = max(10, right/10)
      else
        extra_right = options%extra_right
      end if
      if ( options%max_right >= 0 ) &
        extra_right = max(0, min(extra_right, options%max_right - right))
      keep%block_size = left + right + extra_left + extra_right
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate ( keep%ind(keep%block_size) )
      allocate ( keep%U(max_nep, 2*keep%block_size + max_nep) )
      allocate ( keep%V(2*keep%block_size, 2*keep%block_size, 3) )
      keep%U = ZERO
      do i = 1, max_nep
        keep%U(i, keep%block_size + i) = ONE
      end do
      if ( problem == 0 ) then
        kw = 5
      else
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep) )
        kw = 7
      end if
      allocate ( keep%W(n, keep%block_size, 0:kw) )
      call random_number( keep%W(:,:,0) )
      if ( options%user_X > 0 ) &
        call mxcopy( 'A', n, options%user_X, X, ldx, keep%W, n )
    end if
    
    nep = inform%left + inform%right
    ldV = 2*keep%block_size
    ldBX = n

    select case ( keep%step )
    
    case ( 0 )

      call ssmfe_solve &
        ( problem, sigma, left, right, &
          max_nep, lambda, keep%block_size, keep%V, keep%ind, &
          rci, keep%keep, options, inform )

      select case ( rci%job )

      case ( 11:19 )

        ! apply some standard vector operations to vectors in W
        call ssmfe_vec_ops_double &
          ( n, keep%block_size, keep%W(1,1,0), n, keep%W(1,1,1), n, &
            keep%V, rci, keep%ind, keep%U )

      case ( SSMFE_SAVE_CONVERGED )
    
        if ( rci%nx < 1 ) return

        do i = 1, rci%nx
          k = (i - 1)*rci%i
          if ( rci%i > 0 ) then
            j = keep%lcon + i
          else
            j = max_nep - keep%rcon - i + 1
          end if
          call copy( n, keep%W(1, rci%jx + k, 0), 1, X(1, j), 1 )
          if ( problem /= 0 ) &
            call copy( n, keep%W(1, rci%jy + k, rci%ky), 1, keep%BX(1, j), 1 )
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
        s = ZERO
        do i = 1, max_nep
          do j = 1, max_nep
            if ( i == j ) then
              s = max(s, abs(keep%U(i, m + i) - ONE))
            else
              s = max(s, abs(keep%U(i, m + j)))
            end if
          end do
        end do

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

        if ( keep%rcon > 0 ) then
          j = max_nep - keep%rcon + 1
          call gemm &
            ( TR, 'N', keep%rcon, rci%nx, n, &
              ONE, X(1, j), ldX, rci%x, n, ZERO, keep%U, max_nep )
          if ( problem == 0 ) then
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -ONE, X(1,j), ldX, keep%U, max_nep, ONE, rci%x, n )
          else
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -ONE, keep%BX(1,j), ldBX, keep%U, max_nep, ONE, rci%x, n )
          end if
        end if

      case ( SSMFE_APPLY_CONSTRAINTS )

        if ( keep%lcon == 0 .and. keep%rcon == 0 ) return

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

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

        m = keep%block_size
        k = m + max_nep + 1
        call mxcopy &
          ( 'A', max_nep, rci%nx, keep%U, max_nep, keep%U(1, k), max_nep )
        call gemm &
          ( 'N', 'N', max_nep, rci%nx, max_nep, &
            -ONE, keep%U(1, m + 1), max_nep, keep%U(1, k), max_nep, &
            2*ONE, keep%U, max_nep )

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

!      case ( SSMFE_APPLY_A )
      case ( SSMFE_DO_SHIFTED_SOLVE )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
!        rci%job = SSMFE_DO_SHIFTED_SOLVE

      case ( SSMFE_APPLY_B )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

      case ( SSMFE_APPLY_PREC )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

        call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

      case ( SSMFE_RESTART )
    
        if ( rci%k == 0 ) then
          if ( rci%jx > 1 ) call random_number( keep%W(:, 1 : rci%jx - 1, 0) )
          if ( rci%jx + rci%nx - 1 < keep%block_size ) &
            call random_number &
              ( keep%W(:, rci%jx + rci%nx : keep%block_size, 0) )
        end if
        
      case ( :-1 )

        if ( nep < 1 ) return
    
        if ( inform%left > 0 ) then
          do j = 1, inform%left/2
            s = lambda(j)
            lambda(j) = lambda(inform%left - j + 1)
            lambda(inform%left - j + 1) = s
            do i = 1, n
              s = X(i, j)
              X(i, j) = X(i, inform%left - j + 1)
              X(i, inform%left - j + 1) = s
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
              s = X(i, inform%left + j)
              X(i, inform%left + j) = X(i, nep - j + 1)
              X(i, nep - j + 1) = s
            end do
          end do
        end if

        deallocate ( keep%W )
        allocate ( keep%W(n, nep, 3) )

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

      end select
      
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
      rci%job = SSMFE_APPLY_A
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 2
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,2)
      keep%step = 3

    case ( 3 )
    
      if ( problem == 0 ) then
        call mxcopy( 'A', n, nep, keep%W, n, keep%W(1,1,3), n )
        rci%job = 100
      else
        rci%job = SSMFE_APPLY_B
        rci%nx = nep
        rci%jx = 1
        rci%kx = 0
        rci%jy = 1
        rci%ky = 4
        rci%x => keep%W(:,:,1)
        rci%y => keep%W(:,:,3)
      end if
      keep%step = 4
      
    case ( 4 )

      if ( problem == 2 ) then
        s = -ONE
      else
        s = ONE
      end if
      call gemm &
        ( TR, 'N', nep, nep, n, s, keep%W, n, keep%W(1,1,2), n, &
          ZERO, keep%V, ldV )
      call gemm &
        ( TR, 'N', nep, nep, n, ONE, keep%W, n, keep%W(1,1,3), n, &
          ZERO, keep%V(1,1,2), ldV )
      k = 4*keep%block_size**2
      call dsygv &
        ( 1, 'V', 'U', nep, keep%V, ldV, keep%V(1,1,2), ldV, lambda, &
          keep%V(1,1,3), k, i )
      call gemm &
        ( 'N', 'N', n, nep, nep, ONE, keep%W, n, keep%V, ldV, ZERO, X, ldX )

      if ( problem == 2 ) then
        do i = 1, nep
          lambda(i) = -ONE/lambda(i)
        end do
      end if

      rci%job = -1

    end select

  end subroutine ssmfe_inverse_srci_double

  subroutine ssmfe_terminate_simple_double( keep, inform )

    implicit none
    
    type(ssmfe_keepd), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: inform
    
    call ssmfe_delete_keep_simple_double( keep )
    call ssmfe_delete_info_double( inform )

  end subroutine ssmfe_terminate_simple_double

  subroutine ssmfe_delete_keep_simple_double( keep )

    implicit none
    
    type(ssmfe_keepd), intent(inout) :: keep
    
    if ( allocated(keep%ind) ) deallocate ( keep%ind )
    if ( allocated(keep%BX ) ) deallocate ( keep%BX  )
    if ( allocated(keep%U  ) ) deallocate ( keep%U   )
    if ( allocated(keep%V  ) ) deallocate ( keep%V   )
    if ( allocated(keep%W  ) ) deallocate ( keep%W   )
    call ssmfe_delete_keep_double( keep%keep )

  end subroutine ssmfe_delete_keep_simple_double

  subroutine ssmfe_std_double_complex &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    complex(PRECISION), intent(inout) :: X(ldX, mep)
    integer, intent(in) :: ldX
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keepz  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_srci_double_complex &
      ( 0, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_std_double_complex

  subroutine ssmfe_gen_double_complex &
      ( rci, left, mep, lambda, n, X, ldx, keep, options, inform )

    integer, intent(in) :: left
    integer, intent(in) :: mep
    real(PRECISION), intent(inout) :: lambda(mep)
    integer, intent(in) :: n
    complex(PRECISION), intent(inout) :: X(ldX, mep)
    integer, intent(in) :: ldX
    type(ssmfe_rciz   ), intent(inout) :: rci
    type(ssmfe_keepz  ), intent(inout) :: keep
    type(ssmfe_options), intent(in   ) :: options
    type(ssmfe_inform ), intent(inout) :: inform
    
    call ssmfe_direct_srci_double_complex &
      ( 1, left, mep, lambda, n, X, ldX, rci, keep, options, inform )

  end subroutine ssmfe_gen_double_complex

  subroutine ssmfe_direct_srci_double_complex &
      ( problem, left, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, &
      copy => zcopy, &
      norm => dznrm2, &
      dot  => zdotc, &
      axpy => zaxpy, &
      gemm => zgemm

    use SPRAL_lapack_iface, mxcopy => zlacpy

    implicit none

    character, parameter :: TR = 'C'

    integer, parameter :: SSMFE_START              = 0
    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_CHECK_CONVERGENCE  = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

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

    integer :: kw
    integer :: ldBX
    integer :: extra
    integer :: i, j, k
    
    real(PRECISION), allocatable :: dwork(:,:)
    
    if ( rci%job == 0 ) then
      extra = options%extra_left
      if ( extra < 0 ) then
        if ( left == 1 ) then
          extra = 7
        else
          extra = max(left/10, 10)
        end if
      else if ( extra == 0 ) then
        if ( left == 1 ) extra = 1
      end if
      extra = min(extra, n/2 - left - 1)
      keep%block_size = left + extra
      keep%lcon = 0
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate ( keep%ind(keep%block_size) )
      allocate ( keep%U(max_nep, keep%block_size) )
      allocate ( keep%V(2*keep%block_size, 2*keep%block_size, 3) )
      if ( problem == 0 ) then
        kw = 5
      else
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep) )
        kw = 7
      end if
      allocate ( keep%W(n, keep%block_size, 0:kw) )
      allocate ( dwork(n, keep%block_size) )
      call random_number( dwork )
      keep%W(:,:,0) = 2*dwork - ONE
      deallocate ( dwork )
      if ( options%user_X > 0 ) &
        call mxcopy( 'A', n, options%user_X, X, ldx, keep%W, n )
    end if

    ldBX = n
    
    call ssmfe_solve &
      ( problem, left, max_nep, lambda, keep%block_size, keep%V, keep%ind, &
        rci, keep%keep, options, inform )

    select case ( rci%job )

    case ( 11:19 )

      ! apply some standard vector operations to vectors in W
      call ssmfe_vec_ops_double_complex &
        ( n, keep%block_size, keep%W(1,1,0), n, keep%W(1,1,1), n, &
          keep%V, rci, keep%ind, keep%U )

    case ( SSMFE_SAVE_CONVERGED )
    
      if ( rci%i < 0 ) return
      do i = 1, rci%nx
        k = (i - 1)*rci%i
        j = keep%lcon + i
        call copy( n, keep%W(1, rci%jx + k, 0), 1, X(1,j), 1 )
        if ( problem /= 0 ) &
          call copy( n, keep%W(1, rci%jy + k, rci%ky), 1, keep%BX(1,j), 1 )
      end do
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
      rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

    case ( SSMFE_APPLY_PREC )

      rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
      rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
      call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

    case ( SSMFE_RESTART )
    
      if ( rci%k == 0 ) then
        allocate ( dwork(n, keep%block_size) )
        call random_number( dwork )
        if ( rci%jx > 1 ) &
          keep%W(:, 1 : rci%jx - 1, 0) = 2*dwork(:, 1 : rci%jx - 1) - ONE
        if ( rci%jx + rci%nx - 1 < keep%block_size ) &
          keep%W(:, rci%jx + rci%nx : keep%block_size, 0) = &
            2*dwork(:, rci%jx + rci%nx : keep%block_size) - ONE
        deallocate ( dwork )
      end if
        
    end select
      
  end subroutine ssmfe_direct_srci_double_complex
  
  subroutine ssmfe_shift_double_complex &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )

    implicit none

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

  subroutine ssmfe_gen_shift_double_complex &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )

    implicit none

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

  subroutine ssmfe_buckling_double_complex &
      ( rci, sigma, left, right, mep, lambda, n, X, ldx, keep, options, inform )

    implicit none

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

  subroutine ssmfe_inverse_srci_double_complex &
      ( problem, sigma, left, right, max_nep, lambda, n, X, ldX, &
        rci, keep, options, inform )

    use SPRAL_blas_iface, &
      copy => zcopy, &
      norm => dznrm2, &
      dot  => zdotc, &
      axpy => zaxpy, &
      gemm => zgemm
    use SPRAL_lapack_iface, mxcopy => zlacpy

    implicit none

    integer, parameter :: SSMFE_START              = 0
    integer, parameter :: SSMFE_APPLY_A            = 1
    integer, parameter :: SSMFE_APPLY_PREC         = 2
    integer, parameter :: SSMFE_APPLY_B            = 3
    integer, parameter :: SSMFE_CHECK_CONVERGENCE  = 4
    integer, parameter :: SSMFE_SAVE_CONVERGED     = 5
    integer, parameter :: SSMFE_DO_SHIFTED_SOLVE   = 9
    integer, parameter :: SSMFE_APPLY_CONSTRAINTS  = 21
    integer, parameter :: SSMFE_APPLY_ADJ_CONSTRS  = 22
    integer, parameter :: SSMFE_RESTART            = 999

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
    integer :: kw
    integer :: ldV
    integer :: ldBX
    integer :: i, j, k
    
    real(PRECISION) :: s
    
    real(PRECISION), allocatable :: dwork(:,:)
    
    complex(PRECISION) :: z
    
    if ( rci%job == SSMFE_START ) then
      keep%step = 0
      keep%lcon = 0
      keep%rcon = 0
      keep%block_size = left + right + max((left + right)/10, 10)
      if ( allocated(keep%ind) ) deallocate ( keep%ind )
      if ( allocated(keep%U  ) ) deallocate ( keep%U   )
      if ( allocated(keep%V  ) ) deallocate ( keep%V   )
      if ( allocated(keep%W  ) ) deallocate ( keep%W   )
      allocate ( keep%ind(keep%block_size) )
      allocate ( keep%U(max_nep, keep%block_size) )
      allocate ( keep%V(2*keep%block_size, 2*keep%block_size, 3) )
      if ( problem == 0 ) then
        kw = 5
      else
        if ( allocated(keep%BX) ) deallocate ( keep%BX )
        allocate ( keep%BX(n, max_nep) )
        kw = 7
      end if
      allocate ( keep%W(n, keep%block_size, 0:kw) )
      allocate ( dwork(n, keep%block_size) )
      call random_number( dwork )
      keep%W(:,:,0) = 2*dwork - ONE
      deallocate ( dwork )
      if ( options%user_X > 0 ) &
        call mxcopy( 'A', n, options%user_X, X, ldx, keep%W, n )
    end if
    
    nep = inform%left + inform%right
    ldV = 2*keep%block_size
    ldBX = n

    select case ( keep%step )
    
    case ( 0 )

      call ssmfe_solve &
        ( problem, sigma, left, right, &
          max_nep, lambda, keep%block_size, keep%V, keep%ind, &
          rci, keep%keep, options, inform )

      select case ( rci%job )

      case ( 11:19 )

        ! apply some standard vector operations to vectors in W
        call ssmfe_vec_ops_double_complex &
          ( n, keep%block_size, keep%W(1,1,0), n, keep%W(1,1,1), n, &
            keep%V, rci, keep%ind, keep%U )

      case ( SSMFE_SAVE_CONVERGED )
    
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

        if ( keep%rcon > 0 ) then
          j = max_nep - keep%rcon + 1
          call gemm &
            ( TR, 'N', keep%rcon, rci%nx, n, &
              UNIT, X(1, j), ldX, rci%x, n, NIL, keep%U, max_nep )
          if ( problem == 0 ) then
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -UNIT, X(1,j), ldX, keep%U, max_nep, UNIT, rci%x, n )
          else
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -UNIT, keep%BX(1,j), ldBX, keep%U, max_nep, UNIT, rci%x, n )
          end if
        end if

      case ( SSMFE_APPLY_CONSTRAINTS )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

        if ( keep%lcon > 0 ) then
          call gemm &
            ( TR, 'N', keep%lcon, rci%nx, n, &
              UNIT, X, ldX, rci%y, n, NIL, keep%U, max_nep )
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
            ( TR, 'N', keep%rcon, rci%nx, n, &
              UNIT, X(1, j), ldX, rci%y, n, NIL, keep%U, max_nep )
          call gemm &
            ( 'N', 'N', n, rci%nx, keep%rcon, &
              -UNIT, X(1, j), ldX, keep%U, max_nep, UNIT, rci%x, n )
          if ( problem /= 0 ) &
            call gemm &
              ( 'N', 'N', n, rci%nx, keep%rcon, &
                -UNIT, keep%BX(1, j), ldBX, keep%U, max_nep, UNIT, rci%y, n )
        end if

      case ( SSMFE_APPLY_A )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)
        rci%job = SSMFE_DO_SHIFTED_SOLVE

      case ( SSMFE_APPLY_B )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

      case ( SSMFE_APPLY_PREC )

        rci%x => keep%W(:, rci%jx : rci%jx + rci%nx - 1, rci%kx)
        rci%y => keep%W(:, rci%jy : rci%jy + rci%ny - 1, rci%ky)

        call copy( n*rci%nx, rci%x, 1, rci%y, 1 )

      case ( SSMFE_RESTART )
    
        if ( rci%k == 0 ) then
          allocate ( dwork(n, keep%block_size) )
          call random_number( dwork )
          if ( rci%jx > 1 ) &
            keep%W(:, 1 : rci%jx - 1, 0) = 2*dwork(:, 1 : rci%jx - 1) - ONE
          if ( rci%jx + rci%nx - 1 < keep%block_size ) &
            keep%W(:, rci%jx + rci%nx : keep%block_size, 0) = &
              2*dwork(:, rci%jx + rci%nx : keep%block_size) - ONE
          deallocate ( dwork )
        end if
        
      case ( :-1 )

        if ( nep < 1 ) return
    
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

        deallocate ( keep%W )
        allocate ( keep%W(n, nep, 3) )

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

      end select
      
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
      rci%job = SSMFE_APPLY_A
      rci%nx = nep
      rci%jx = 1
      rci%kx = 0
      rci%jy = 1
      rci%ky = 2
      rci%x => keep%W(:,:,1)
      rci%y => keep%W(:,:,2)
      keep%step = 3

    case ( 3 )
    
      if ( problem == 0 ) then
        call mxcopy( 'A', n, nep, keep%W, n, keep%W(1,1,3), n )
        rci%job = 100
      else
        rci%job = SSMFE_APPLY_B
        rci%nx = nep
        rci%jx = 1
        rci%kx = 0
        rci%jy = 1
        rci%ky = 4
        rci%x => keep%W(:,:,1)
        rci%y => keep%W(:,:,3)
      end if
      keep%step = 4
      
    case ( 4 )

      if ( problem == 2 ) then
        z = -UNIT
      else
        z = UNIT
      end if
      call gemm &
        ( TR, 'N', nep, nep, n, z, keep%W, n,keep% W(1,1,2), n, &
          NIL, keep%V, ldV )
      call gemm &
        ( TR, 'N', nep, nep, n, UNIT, keep%W, n, keep%W(1,1,3), n, &
          NIL, keep%V(1,1,2), ldV )
      k = 4*keep%block_size**2
      call zhegv &
        ( 1, 'V', 'U', nep, keep%V, ldV, keep%V(1,1,2), ldV, lambda, &
          keep%V(1,1,3), k, i )
      call gemm &
        ( 'N', 'N', n, nep, nep, UNIT, keep%W, n, keep%V, ldV, NIL, X, ldX )

      if ( problem == 2 ) then
        do i = 1, nep
          lambda(i) = -ONE/lambda(i)
        end do
      end if

      rci%job = -1

    end select

  end subroutine ssmfe_inverse_srci_double_complex

  subroutine ssmfe_terminate_simple_double_complex( keep, inform )

    implicit none
    
    type(ssmfe_keepz), intent(inout) :: keep
    type(ssmfe_inform), intent(inout) :: inform
    
    call ssmfe_delete_keep_simple_double_complex( keep )
    call ssmfe_delete_info_double( inform )

  end subroutine ssmfe_terminate_simple_double_complex

  subroutine ssmfe_delete_keep_simple_double_complex( keep )

    implicit none
    
    type(ssmfe_keepz), intent(inout) :: keep
    
    if ( allocated(keep%ind) ) deallocate ( keep%ind )
    if ( allocated(keep%BX ) ) deallocate ( keep%BX  )
    if ( allocated(keep%U  ) ) deallocate ( keep%U   )
    if ( allocated(keep%V  ) ) deallocate ( keep%V   )
    if ( allocated(keep%W  ) ) deallocate ( keep%W   )
    call ssmfe_delete_keep_double( keep%keep )

  end subroutine ssmfe_delete_keep_simple_double_complex

end module SPRAL_ssmfe


