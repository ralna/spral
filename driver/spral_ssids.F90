program run_prob
  use, intrinsic :: iso_c_binding
!$ use omp_lib
  use cuda_helper
  use spral_hw_topology
  use spral_rutherford_boeing
  use spral_ssids
  use spral_matrix_util, only : cscl_verify, SPRAL_MATRIX_REAL_SYM_INDEF
  use spral_scaling
  implicit none

  integer, parameter :: wp = kind(0d0)

  type(rb_read_options) :: rb_options
  integer :: rb_flag

  ! Matrix description
  integer :: m, n
  integer, dimension(:), allocatable :: ptr, row
  real(wp), dimension(:), allocatable :: val

  type(ssids_inform) :: inform
  type(ssids_akeep) :: akeep
  type(ssids_fkeep) :: fkeep
  type(ssids_options) :: options
  integer :: cuda_error
  double precision, dimension(:, :), allocatable :: rhs, soln
  double precision, dimension(:), allocatable :: res, scaling

  integer :: i, j, k, r

  integer :: start_t, stop_t, rate_t
  integer :: flag, more

  ! integer, parameter :: unit_rhs = 14

  real :: smanal, smfact, smaflop, smafact

  integer, parameter :: nfact = 1
  ! integer, parameter :: nfact = 50
  ! integer, parameter :: nfact = 100

  integer, parameter :: nslv = 1
  ! integer, parameter :: nslv = 10
  ! integer, parameter :: nslv = 100

  integer :: nrhs

  logical :: force_psdef, pos_def, time_scaling, flat_topology

  type(numa_region), dimension(:), allocatable :: topology
  character(len=:), allocatable :: filename

  integer(C_INT) :: cnt

  integer :: ngpus

  call proc_args(filename, options, force_psdef, pos_def, nrhs, time_scaling, &
       flat_topology, ngpus)
  if (nrhs .lt. 1) stop

  ! Read in a matrix
  write (*, "(3a)") "Reading '", filename, "'..."
  if (force_psdef) rb_options%values = -3 ! Force diagonal dominance
  rb_options%values = 2 ! make up values if necessary
  call rb_read(filename, m, n, ptr, row, val, rb_options, rb_flag)
  if (rb_flag .ne. 0) then
     print *, "Rutherford-Boeing read failed with error ", rb_flag
     stop
  end if
  write(*, "(a)") "ok"

  ! Make up a rhs
  allocate(rhs(n, nrhs), soln(n, nrhs))
  rhs = 0
  do r = 1, nrhs
     do i = 1, n
        do j = ptr(i), ptr(i+1)-1
           k = row(j)
           rhs(k, r) = rhs(k, r) + val(j)
           if (i .eq. k) cycle
           rhs(i, r) = rhs(i, r) + val(j)
        end do
     end do
  end do

  call cuda_init(cnt)
  if (cnt .lt. 0) then
     print *, "CUDA_INIT failed: ", cnt
     stop
  else if (cnt .gt. 0) then
     print *, "Number of CUDA devices: ", cnt
     ngpus = min(ngpus, cnt)
  end if

  if (flat_topology) then
     allocate(topology(1))
     topology(1)%nproc = 1
!$   topology(1)%nproc = min(omp_get_max_threads(),omp_get_thread_limit())
     print *, "Forcing topology to ", topology(1)%nproc
     print *, "Using", ngpus, "GPUs"
     if (ngpus .gt. 0) then
        allocate(topology(1)%gpus(ngpus))
        do i=1,ngpus
           topology(1)%gpus(i) = i-1
        end do
     else
        allocate(topology(1)%gpus(0))
     end if
  end if

  call cscl_verify(6, SPRAL_MATRIX_REAL_SYM_INDEF, n, n, &
       ptr, row, flag, more)
  if (flag .ne. 0) then
     print *, "CSCL_VERIFY failed: ", flag, more
     stop
  end if

  ! Analyse and factor
  call system_clock(start_t, rate_t)
  if (allocated(topology)) then
     call ssids_analyse(.false., n, ptr, row, akeep, &
          options, inform, val=val, topology=topology)
  else
     call ssids_analyse(.false., n, ptr, row, akeep, &
          options, inform, val=val)
  end if
  call system_clock(stop_t)
  print *, "Used order ", options%ordering
  if (inform%flag .lt. 0) then
     print *, "oops on analyse ", inform%flag
     stop
  end if
  write (*, "(a)") "ok"
  print *, "Analyse took ", (stop_t - start_t)/real(rate_t)
  !print *, "Used maximum memory of ", inform%maxmem
  smanal = (stop_t - start_t)/real(rate_t)
  print "(a,es10.2)", "Predict nfact = ", real(inform%num_factor)
  print "(a,es10.2)", "Predict nflop = ", real(inform%num_flops)
  print "(a6, i10)", "nparts", inform%nparts
  print "(a6, es10.2)", "cpu_flops", real(inform%cpu_flops)
  print "(a6, es10.2)", "gpu_flops", real(inform%gpu_flops)
  smaflop = real(inform%num_flops)
  smafact = real(inform%num_factor)

  if (time_scaling .and. (options%scaling .ne. 0)) then
     allocate(scaling(n))
     call do_timed_scaling(n, ptr, row, val, scaling)
     options%scaling = 0 ! We will user-supply
  end if

  write (*, "(a)") "Factorize..."
  call system_clock(start_t, rate_t)
  do i = 1, nfact
     if (allocated(scaling)) then
        call ssids_factor(pos_def, val, akeep, fkeep, &
             options, inform, ptr=ptr, row=row, scale=scaling)
     else
        call ssids_factor(pos_def, val, akeep, fkeep, &
             options, inform, ptr=ptr, row=row)
     end if
  end do
  call system_clock(stop_t)
  if (inform%flag .lt. 0) then
     print *, "oops on factorize ", inform%flag
     stop
  end if
  write (*, "(a)") "ok"
  print *, "Factor took ", (stop_t - start_t)/real(rate_t)
  smfact = (stop_t - start_t)/real(rate_t)

  ! Solve
  write (*, "(a)") "Solve..."
  call system_clock(start_t, rate_t)
  do i = 1, nslv
     soln = rhs
     call ssids_solve(nrhs,soln,n,akeep,fkeep,options,inform)
  end do
  call system_clock(stop_t)
  if (inform%flag .lt. 0) then
     print *, "oops on solve ", inform%flag
     stop
  end if
  write(*, "(a)") "ok"
  print *, "Solve took ", (stop_t - start_t)/real(rate_t)

  print *, "number bad cmp = ", count(abs(soln(1:n,1) - dble(1.0)) .ge. dble(1e-6))
  print *, "fwd error || ||_inf = ", maxval(abs(soln(1:n,1) - dble(1.0)))
  allocate(res(nrhs))
  call internal_calc_norm(n, ptr, row, val, soln, rhs, nrhs, res)
  print *, "bwd error scaled = ", res

  call ssids_free(akeep, fkeep, cuda_error)

  print "(a6, a10)", "cmp:","SMFCT"
  print "(a6, f10.2)", "anal:", smanal
  print "(a6, f10.2)", "fact:", smfact
  print "(a6, es10.2)", "afact:", smafact
  print "(a6, es10.2)", "aflop:", smaflop
  print "(a6, es10.2)", "nfact:", real(inform%num_factor)
  print "(a6, es10.2)", "nflop:", real(inform%num_flops)
  print "(a6, i10)", "delay:", inform%num_delay
  print "(a6, 3i10)", "inertia:", inform%num_neg, n-inform%matrix_rank,&
       inform%matrix_rank-inform%num_neg
  print "(a6, i10)", "2x2piv:", inform%num_two
  print "(a6, i10)", "maxfront:", inform%maxfront
  print "(a6, i10)", "maxsupernode:", inform%maxsupernode
  print "(a6, i10)", "not_first_pass:", inform%not_first_pass
  print "(a6, i10)", "not_second_pass:", inform%not_second_pass

  ! Free memory to ensure we pass leak-check tests
  deallocate(ptr, rhs, soln, res)
  if (allocated(topology)) deallocate(topology)

contains

  subroutine proc_args(filename, options, force_psdef, pos_def, nrhs, &
       time_scaling, flat_topology, ngpus)
    implicit none
    character(len=:), allocatable :: filename
    type(ssids_options), intent(inout) :: options
    logical, intent(out) :: force_psdef
    logical, intent(out) :: pos_def
    integer, intent(out) :: nrhs
    logical, intent(out) :: time_scaling
    logical, intent(out) :: flat_topology
    integer, intent(out) :: ngpus

    integer :: argnum, narg
    integer :: i
    character(len=200) :: argval
    logical :: seen_fname

    ! Defaults
    nrhs = 1
    force_psdef = .false.
    pos_def = .false.
    time_scaling = .false.
    flat_topology = .true.
    filename = "matrix.rb"
    seen_fname = .false.
    ngpus = 0

    ! Process args
    narg = command_argument_count()
    argnum = 1
    do while (argnum .le. narg)
       call get_command_argument(argnum, argval)
       argnum = argnum + 1
       select case (argval)
       case("--scale=none")
          options%scaling = 0 ! None
          print *, "Set scaling to None"
       case("--scale=mc64")
          options%scaling = 1 ! MC64
          print *, "Set scaling to MC64"
       case("--scale=auction")
          options%scaling = 2 ! Auction algorithm
          print *, "Set scaling to Auction"
       case("--scale=mc77")
          options%scaling = 4 ! MC77 algorithm
          print *, "Set scaling to MC77"
       case("--ordering=mc64-metis")
          options%ordering = 2 ! Matching-based ordering
          options%scaling = 3 ! Scaling from matching ordering
          print *, "Using matching-based ordering (scaling overwritten)"
       case("--force-posdef")
          force_psdef = .true.
          print *, "Forcing matrix to be positive definite"
       case("--posdef")
          pos_def = .true.
          print *, 'Matrix assumed positive definite'
       case("--time-scaling")
          time_scaling = .true.
       case("--timing")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) i
          if (i .gt. 0) options%print_level = -i
       case("--nrhs")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) nrhs
          print *, 'solving for', nrhs, 'right-hand sides'
       case("--nemin")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%nemin
          print *, 'Supernode amalgamation nemin = ', options%nemin
       case("--u")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%u
          print *, 'Pivoting threshold u = ', options%u
       case("--nstream")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%nstream
       case("--max-load-inbalance")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%max_load_inbalance
          print *, 'Max load inbalance = ', options%max_load_inbalance
       case("--pivot-method=app-aggressive")
          options%pivot_method = 1
          print *, 'Pivoting method APP_AGGRESSIVE'
       case("--pivot-method=app-block")
          options%pivot_method = 2
          print *, 'Pivoting method APP_BLOCK'
       case("--pivot-method=tpp")
          options%pivot_method = 3
          print *, 'Pivoting method TPP'
       case("--failed-pivot-method=tpp")
          options%failed_pivot_method = 1
          print *, 'Failed pivot method TPP'
       case("--failed-pivot-method=pass")
          options%failed_pivot_method = 2
          print *, 'Failed pivot method PASS'
       case("--flat-topology")
          flat_topology = .true.
          print *, 'Forcing flat topology'
       case("--no-flat-topology")
          flat_topology = .false.
          print *, 'Using machine topology'
       case("--disable-gpu")
          options%use_gpu = .false.
          print *, 'Disabling GPUs'
       case("--min-gpu-work")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%min_gpu_work
          print *, 'Min GPU work = ', &
               options%min_gpu_work
       case("--gpu-perf-coeff")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%gpu_perf_coeff
          print *, 'GPU Performance coefficient = ', &
               options%gpu_perf_coeff
       case("--small-subtree-threshold")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%small_subtree_threshold
          print *, 'Small subtree treshold = ', &
               options%small_subtree_threshold
       case("--cpu-block-size")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) options%cpu_block_size
          print *, 'CPU block size = ', options%cpu_block_size
       case("--no-ignore-numa")
          options%ignore_numa = .false.
          print *, 'Using separate NUMA regions'
       case("--ngpus")
          call get_command_argument(argnum, argval)
          argnum = argnum + 1
          read (argval, *) ngpus
          print *, 'NGPUS = ', ngpus
       case default
          if (seen_fname) then
             print *, "Unrecognised command line argument: ", argval
             stop
          else
             filename = trim(argval)
             seen_fname = .true.
          end if
       end select
    end do
  end subroutine proc_args

  subroutine do_timed_scaling(n, ptr, row, val, scaling)
    implicit none
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), dimension(n), intent(out) :: scaling

    type(auction_inform) :: ainform
    type(equilib_options) :: eoptions
    type(equilib_inform) :: einform
    type(hungarian_options) :: hoptions
    type(hungarian_inform) :: hinform

    write (*, "(a)") "Scaling..."

    call system_clock(start_t, rate_t)
    ! Note: we assume no checking required
    select case (options%scaling)
    case (1) ! Hungarian algorithm
       hoptions%scale_if_singular = .true.
       call hungarian_scale_sym(n, ptr, row, val, scaling, hoptions, hinform)
       if (hinform%flag .le. 0) then
          print *, "Error from hungarian_scale_sym()"
          stop
       end if
    case (2) ! Auction algorithm
       call auction_scale_sym(n, ptr, row, val, scaling, &
            options%auction, ainform)
       if (ainform%flag .ne. 0) then
          print *, "Error from auction_scale_sym() flag = ", ainform%flag
          stop
       end if
    case (3) ! From ordering
       print *, "--time-scaling not supported with matching-based ordering"
       stop
    case (4) ! MC77-like
       call equilib_scale_sym(n, ptr, row, val, scaling, eoptions, einform)
       if (einform%flag .ne. 0) then
          print *, "Error from equilib_scale_sym()"
          stop
       end if
    end select
    call system_clock(stop_t)
    write (*, "(a)") "ok"
    print *, "Scaling took ", (stop_t - start_t)/real(rate_t)
  end subroutine do_timed_scaling

  subroutine internal_calc_norm(n, ptr, row, val, x_vec, b_vec, nrhs, res)
    implicit none
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    integer, intent(in) :: nrhs
    real(wp), dimension(nrhs*n), intent(in) :: x_vec
    real(wp), dimension(nrhs*n), intent(in) :: b_vec
    real(wp), dimension(nrhs), intent(out) :: res

    integer :: i, j, k, r
    double precision, allocatable, dimension(:) :: x_norm
    real(wp), dimension(:), allocatable :: res_vec
    double precision :: temp
    double precision :: normA

    ! Find the residual
    allocate(res_vec(n*nrhs), x_norm(nrhs))
    res_vec = 0
    do i = 1, n
       do j = ptr(i), ptr(i+1)-1
          r = row(j)
          do k = 0, nrhs-1
             res_vec(i+k*n) = res_vec(i+k*n) + &
                  val(j) * x_vec(r+k*n)
          end do
          if (r .eq. i) cycle
          do k = 0, nrhs-1
             res_vec(r+k*n) = res_vec(r+k*n) + &
                  val(j) * x_vec(i+k*n)
          end do
       end do
    end do
    res_vec(:) = res_vec(:) - b_vec(:)

    ! Find matrix norm
    call matrix_inf_norm(n, ptr, row, val, normA)

    ! Find x norm
    do i = 1, nrhs
       x_norm(i) = 0
       do j = 1, n
          x_norm(i) = max(x_norm(i), abs(x_vec((i-1)*n+j)))
          if (x_vec((i-1)*n+j) .ne. x_vec((i-1)*n+j)) then ! Tests for NaN
             x_norm(i) = x_vec((i-1)*n+j)
             exit
          end if
       end do
    end do

    ! Scaled residual = ||r|| / ( ||A|| ||x|| + ||b|| )
    do i = 1, nrhs
       temp = normA * x_norm(i) + &
            maxval(abs(b_vec((i-1)*n+1:i*n)))
       if (temp .eq. dble(0.0)) then
          res(i) = maxval(abs(res_vec((i-1)*n+1:i*n)))
       else
          res(i) = maxval(abs(res_vec((i-1)*n+1:i*n))) / temp
       end if
    end do
  end subroutine internal_calc_norm

  subroutine matrix_inf_norm(n, ptr, row, val, norm)
    implicit none
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    real(wp), dimension(ptr(n+1)-1), intent(in) :: val
    real(wp), intent(out) :: norm

    real(wp), allocatable, dimension(:) :: row_norm
    integer :: i

    allocate(row_norm(n))

    row_norm = 0
    do i = 1, ptr(n+1)-1
       row_norm(row(i)) = row_norm(row(i)) + abs(val(i))
    end do

    norm = maxval(row_norm)
  end subroutine matrix_inf_norm

  ! subroutine perm_mat(matrix, perm)
  !   implicit none
  !   type(zd11_type), intent(inout) :: matrix
  !   integer, dimension(matrix%n), intent(in) :: perm

  !   integer :: ne
  !   integer :: i, j, k
  !   integer, dimension(:), allocatable :: ptr2, row2, invp
  !   double precision, dimension(:), allocatable :: val2

  !   ne = matrix%ptr(matrix%n+1)-1
  !   allocate(ptr2(matrix%n+1), row2(ne), val2(ne))

  !   allocate(invp(matrix%n))
  !   do i = 1, matrix%n
  !      invp( perm(i) ) = i
  !   end do

  !   ! Take a copy with permuted rows
  !   ptr2(:) = matrix%ptr(1:matrix%n+1)
  !   val2(:) = matrix%val(1:ne)
  !   do i = 1, ne
  !      row2(i) = invp( matrix%row(i) )
  !   end do

  !   k = 1
  !   matrix%ptr(1) = 1
  !   do i = 1, matrix%n
  !      j = order(i)
  !      k = k + ptr2(j+1)-ptr2(j)
  !      matrix%ptr(i+1) = k
  !      matrix%row(matrix%ptr(i):matrix%ptr(i+1)-1) = &
  !           row2(ptr2(j):ptr2(j+1)-1)
  !      matrix%val(matrix%ptr(i):matrix%ptr(i+1)-1) = &
  !           val2(ptr2(j):ptr2(j+1)-1)
  !   end do
  ! end subroutine perm_mat
end program
