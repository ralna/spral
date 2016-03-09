program run_prob
   use spral_core_analyse
   use spral_matrix_util
   use spral_metis_wrapper
   use spral_nd
   use spral_random
   use spral_rutherford_boeing
   use spral_scaling
   use spral_timer
   use spral_nd_types, only : FLAG_BIG_COL, FLAG_BIG_BOTH, FLAG_SMALL
   use hsl_ma97_double
   implicit none

   integer, parameter :: wp = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)

   ! RB Reader
   type(rb_reader_options) :: rb_options
   integer :: rb_flag
   integer :: flag, more, st
   character(len=3) :: type_code

   ! Matrix
   integer :: m, n
   integer, dimension(:), allocatable :: ptr, row, col
   real(wp), dimension(:), allocatable :: val

   ! ND and stats
   type(nd_options) :: options
   type(nd_inform) :: inform
   integer, dimension(:), allocatable :: perm2x2, perm, invp
   integer(long) :: nfact, nflops
   real(wp), dimension(:), allocatable :: scaling

   type(ma97_control) :: control97

   ! Timing
   type(timespec) :: t1, t2
   integer :: dummy

   ! Controls
   integer :: random, order_type
   logical :: with_ma87, with_ma97

   integer, parameter :: TYPE_ND                      = 0, &
                         TYPE_METIS                   = 1, &
                         TYPE_NUM_AWARE               = 2, &
                         TYPE_GUPTA                   = 3, &
                         TYPE_GUPTA_ZERO              = 4, &
                         TYPE_MATCH_ORDER_ND          = 5, &
                         TYPE_MATCH_ORDER_METIS       = 6, &
                         TYPE_MATCH_ORDER_ND_RELAX    = 7, &
                         TYPE_MATCH_ORDER_METIS_RELAX = 8

   call proc_args(options, random, with_ma87, with_ma97, order_type, control97)

   ! Read in a matrix
   write(*, "(a)", advance="no") "Reading..."
   rb_options%values = 2 ! make up values if necessary
   call rb_read("matrix.rb", m, n, ptr, row, col, val, rb_options, rb_flag, &
      type_code=type_code)
   if(rb_flag.ne.0) then
      print *, "Rutherford-Boeing read failed with error ", rb_flag
      stop
   endif
   write(*, "(a)") "ok"

   ! If unsymmetric, form A+A^T
   if(type_code(2:2).eq.'u') then
      if(m.ne.n) then
         print *, "Matrix is not square!"
         stop 1
      endif
      print *, "Symmetrizing unsymmetric problem..."
      call symmetrize_problem(n, ptr, row, val)
   endif
   print *, "Input matrix n = ", n
   print *, "Input matrix nz = ", ptr(n+1)-1

   ! Force to be pos-def
   if(with_ma87) call make_diagdom(n, ptr, row, val)

   ! Just to be safe... (if we've symmetrized we don't guaruntee ascending idx)
   if(type_code(2:2).ne.'u') then
      call cscl_verify(6, SPRAL_MATRIX_REAL_SYM_INDEF, n, n, &
         ptr, row, flag, more)
      if(flag.ne.0) then
         print *, "CSCL_VERIFY failed: ", flag, more
         stop
      endif
   endif

   ! Randomize order
   if(random.ne.-1) &
      call randomize_matrix(n, ptr, row, random)

   ! Order using spral_nd
   allocate(perm(n), invp(n))
   select case(order_type)
   case(TYPE_ND)
      write(*, "(a)", advance="no") "Ordering with ND..."
      dummy = clock_gettime(0, t1)
      call nd_order(0, n, ptr, row, perm, options, inform)
      dummy = clock_gettime(0, t2)
      if (inform%flag < 0) then
         print *, "oops on nd ", inform%flag
         stop
      endif
      write(*, "(a)") "ok"
   case(TYPE_METIS)
      write(*, "(a)", advance="no") "Ordering with Metis..."
      dummy = clock_gettime(0, t1)
      call metis_order(n, ptr, row, perm, invp, flag, st)
      dummy = clock_gettime(0, t2)
      if (flag .ne. 0) then
         print *, "oops on metis ", flag
         stop
      endif
      write(*, "(a)") "ok"
   !case(TYPE_GUPTA)
   !   allocate(scaling(n))
   !   call find_gupta_order(n, ptr, row, val, scaling, perm, options%u)
   case(TYPE_NUM_AWARE)
      write(*, "(a)", advance="no") "Ordering with ND..."
      allocate(perm2x2(n))
      dummy = clock_gettime(0, t1)
      call nd_order(0, n, ptr, row, perm2x2, options, inform, val=val)
      dummy = clock_gettime(0, t2)
      perm(:) = abs(perm2x2(:))
      if (inform%flag < 0) then
         print *, "oops on nd ", inform%flag
         stop
      endif
      write(*, "(a)") "ok"
   case(TYPE_GUPTA:)
      write(*, "(a)", advance="no") "Ordering with Match order..."
      dummy = clock_gettime(0, t1)
      call find_match_order(order_type, n, ptr, row, val, perm, options)
      dummy = clock_gettime(0, t2)
      write(*, "(a)") "ok"
   end select
   print *, "order took ", tdiff(t1, t2)

   ! Determine quality
   write(*, "(a)", advance="no") "Determing stats..."
   call calculate_stats(n, ptr, row, perm, nfact, nflops)
   write(*, "(a)") "ok"
   print "(a,es10.2)", "literal nfact = ", real(nfact)
   print "(a,es10.2)", "literal nflop = ", real(nflops)
   if(with_ma87) call run_ma87(n, ptr, row, val, perm)
   if(with_ma97) call run_ma97(n, ptr, row, val, perm, scaling, control97)

contains

   subroutine find_match_order(order_type, n, ptr, row, val, perm, options)
      use spral_match_order
      integer, intent(in) :: order_type
      integer, intent(in) :: n
      integer, dimension(:), intent(in) :: ptr
      integer, dimension(:), allocatable, intent(in) :: row
      real(wp), dimension(:), allocatable, intent(in) :: val
      integer, dimension(n), intent(out) :: perm
      type(nd_options), intent(in) :: options

      integer, dimension(:), allocatable :: ptr2, row2, work
      real(wp), dimension(:), allocatable :: val2, scaling

      type(mo_options) :: mooptions
      type(mo_inform) :: moinform

      allocate(scaling(n), work(n))
      allocate(ptr2(n+1), row2(2*(ptr(n+1)-1)), val2(2*(ptr(n+1)-1)))
      ptr2(1:n+1) = ptr(1:n+1)
      row2(1:ptr(n+1)-1) = row(1:ptr(n+1)-1)
      val2(1:ptr(n+1)-1) = val(1:ptr(n+1)-1)
      call half_to_full(n, row2, ptr2, work, a=val2)
      mooptions%nd_options = options
      mooptions%u = options%u
      select case(order_type)
      case(TYPE_GUPTA)
         mooptions%match_method = 3 ! Gupta
         mooptions%order_method = 2 ! ND
      case(TYPE_GUPTA_ZERO)
         mooptions%match_method = 4 ! Gupta-zero
         mooptions%order_method = 2 ! ND
      case(TYPE_MATCH_ORDER_ND)
         mooptions%order_method = 2 ! ND
      case(TYPE_MATCH_ORDER_METIS)
         ! Leave as defaults
      case(TYPE_MATCH_ORDER_ND_RELAX)
         mooptions%match_method = 2 ! mc64 relaxed
         mooptions%order_method = 2 ! ND
      case(TYPE_MATCH_ORDER_METIS_RELAX)
         mooptions%match_method = 2 ! mc64 relaxed
      end select
      call match_order(n, ptr2, row2, val2, perm, scaling, mooptions, &
         moinform)
      if(moinform%flag.ne.0) then
         print *, "match_order() returned error ", moinform%flag
         stop
      endif
      print *, "Number matched = ", moinform%nmatch
   end subroutine find_match_order

   subroutine proc_args(options, random, with_ma87, with_ma97, order_type, &
         control97)
      type(nd_options), intent(inout) :: options
      integer, intent(out) :: random
      logical, intent(out) :: with_ma87
      logical, intent(out) :: with_ma97
      integer, intent(out) :: order_type
      type(ma97_control) :: control97

      integer :: argnum, narg
      character(len=200) :: argval

      ! Defaults
      random = -1
      with_ma87 = .false.
      with_ma97 = .false.
      order_type = TYPE_ND

      ! Other settings (not changable by command line)
      control97%ordering = 0 ! user supplied
      
      ! Process args
      narg = command_argument_count()
      argnum = 1
      do while(argnum <= narg)
         call get_command_argument(argnum, argval)
         argnum = argnum + 1
         select case(argval)
         case("--ma87")
            with_ma87 = .true.
            print *, "Running ma87"
         case("--ma97")
            with_ma97 = .true.
            print *, "Running ma97"
         case("--cost")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%cost_function
            print *, "Set cost function to ", options%cost_function
         case("--refinement")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%refinement
            print *, "Set options%refinement = ", options%refinement
         case("--amd-call")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%amd_call
            print *, "Set options%amd_call = ", options%amd_call
         case("--amd-switch2")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%amd_switch2
            print *, "Set options%amd_switch2 = ", options%amd_switch2
         case("--stop-coarsening1")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%stop_coarsening1
            print *, "Set options%stop_coarsening1 = ", options%stop_coarsening1
         case("--stop-coarsening2")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%stop_coarsening2
            print *, "Set options%stop_coarsening2 = ", options%stop_coarsening2
         case("--min-reduction")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%min_reduction
            print *, "Set options%min_reduction = ", options%min_reduction
         case("--max-reduction")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%max_reduction
            print *, "Set options%max_reduction = ", options%max_reduction
         case("--nosv")
            options%find_supervariables = .false.
            print *, "Disabled supervariables"
         case("--reord=1")
            options%reord = 1
            print *, "Set Jonathan's preprocessing ordering"
         case("--reord=2")
            options%reord = 2
            print *, "Set Sue's preprocessing ordering"
         case("--partition-method")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%partition_method
            print *, "Set options%partition_method = ", options%partition_method
         case("--coarse-partition-method")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%coarse_partition_method
            print *, "Set options%coarse_partition_method = ", &
               options%coarse_partition_method
         case("--matching")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%matching
            print *, "Set options%matching = ", &
               options%matching
         case("--balance")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%balance
            print *, "Set options%balance = ", &
               options%matching
         case("--print-level")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%print_level
            print *, "Set options%print_level = ", &
               options%print_level
         case("--randomize")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) random
            print *, "Randomizing matrix row order, seed = ", random
         case("--type=num-aware")
            order_type = TYPE_NUM_AWARE
            options%find_supervariables = .false.
            print *, "Using numerically aware ordering"
            print *, "NOTE: Had to disable supervariables for num-aware order"
         case("--type=gupta")
            order_type = TYPE_GUPTA
            print *, "Using Gupta-type ordering"
         case("--type=gupta-zero")
            order_type = TYPE_GUPTA_ZERO
            print *, "Using Gupta-type ordering"
         case("--type=match-order-nd")
            order_type = TYPE_MATCH_ORDER_ND
            print *, "Using ND matching-based ordering"
         case("--type=match-order-nd-relax")
            order_type = TYPE_MATCH_ORDER_ND_RELAX
            print *, "Using ND matching-based relaxed ordering"
         case("--type=metis")
            order_type = TYPE_METIS
            print *, "Using normal MeTiS ordering"
         case("--type=match-order-metis")
            order_type = TYPE_MATCH_ORDER_METIS
            print *, "Using MeTiS matching-based ordering"
         case("--u_ord")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) options%u
            print *, "Set u_ord = ", options%u
         case("--u_num")
            call get_command_argument(argnum, argval)
            argnum = argnum + 1
            read( argval, * ) control97%u
            print *, "Set u_num = ", control97%u
         case default
            print *, "Unrecognised command line argument: ", argval
            stop
         end select
      end do
   end subroutine proc_args

   subroutine add_diag_entries(n, ptr, row, val)
      integer, intent(in) :: n
      integer, dimension(:), intent(inout) :: ptr
      integer, dimension(:), allocatable, intent(inout) :: row
      real(wp), dimension(:), allocatable, intent(inout) :: val

      integer :: i, j, k, insert
      integer, dimension(:), allocatable :: ptr_out, row_out
      real(wp), dimension(:), allocatable :: val_out

      ! Copy matrix to _out arrays, adding diagonal values as necessary
      allocate(ptr_out(n+1), row_out(ptr(n+1)-1+n), val_out(ptr(n+1)-1+n))
      insert = 1
      do i = 1, n
         ptr_out(i) = insert
         k = -1
         do j = ptr(i), ptr(i+1)-1
            k = row(j)
            if(i.eq.k) exit ! Found diagonal entry
         end do
         if(i.ne.k) then
            ! Add diagonal entry and initialise to 0.0
            row_out(insert) = i
            val_out(insert) = 0.0
            insert = insert + 1
         endif
         ! Literal copy of existing values
         row_out(insert:insert+ptr(i+1)-ptr(i)) = row(ptr(i):ptr(i+1)-1)
         val_out(insert:insert+ptr(i+1)-ptr(i)-1) = val(ptr(i):ptr(i+1)-1)
         insert = insert + ptr(i+1)-ptr(i)
      end do
      ptr_out(n+1) = insert
      if(insert.eq.ptr(n+1)) return ! No extra entries required

      !
      ! Copy *_out over original for output
      !
      deallocate(row, val)
      allocate(row(ptr_out(n+1)-1), val(ptr_out(n+1)-1))
      ptr(1:n+1) = ptr_out(1:n+1)
      row(1:ptr(n+1)-1) = row_out(1:ptr_out(n+1)-1)
      val(1:ptr(n+1)-1) = val_out(1:ptr_out(n+1)-1)
   end subroutine add_diag_entries

   subroutine make_diagdom(n, ptr, row, val)
      integer, intent(in) :: n
      integer, dimension(:), intent(inout) :: ptr
      integer, dimension(:), allocatable, intent(inout) :: row
      real(wp), dimension(:), allocatable, intent(inout) :: val

      integer :: i, j, k
      integer, dimension(:), allocatable :: dloc

      call add_diag_entries(n, ptr, row, val)

      ! Find location of diagonal entries
      allocate(dloc(n))
      dloc(:) = -1
      do i = 1, n
         do j = ptr(i), ptr(i+1)-1
            k = row(j)
            if(i.eq.k) then
               dloc(i) = j
               val(j) = 1.0 + abs(val(j)) ! force diagonal entry to be > 0.0
               exit
            endif
         end do
      end do

      ! Now add sum of row entries to diagonals
      do i = 1, n
         do j = ptr(i), ptr(i+1)-1
            k = row(j)
            if(i.eq.k) cycle ! diagonal entry
            val(dloc(i)) = val(dloc(i)) + abs(val(j))
            val(dloc(k)) = val(dloc(k)) + abs(val(j))
         end do
      end do
   end subroutine make_diagdom

   subroutine symmetrize_problem(n, ptr, row, val)
      integer, intent(in) :: n
      integer, dimension(:), intent(inout) :: ptr
      integer, dimension(:), allocatable, intent(inout) :: row
      real(wp), dimension(:), allocatable, intent(inout) :: val

      integer :: i, j, k, insert
      integer, dimension(:), allocatable :: ptr_trans, ptr_out
      integer, dimension(:), allocatable :: row_trans, row_out
      real(wp), dimension(:), allocatable :: val_trans, val_out
      integer, dimension(:), allocatable :: seen

      !
      ! Form A^T
      !
      ! count entries into ptr_trans at offset +2
      allocate(ptr_trans(n+3))
      ptr_trans(:) = 0
      do i = 1, n
         do j = ptr(i), ptr(i+1)-1
            k = row(j)
            ptr_trans(k+2) = ptr_trans(k+2) + 1
         end do
      end do
      ! calculate row starts at offset +1
      ptr_trans(1:2) = 1
      do i = 1, n
         ptr_trans(i+2) = ptr_trans(i+1) + ptr_trans(i+2)
      end do
      ! drop entries into place
      allocate(row_trans(ptr_trans(n+2)), val_trans(ptr_trans(n+2)))
      do i = 1, n
         do j = ptr(i), ptr(i+1)-1
            k = row(j)
            row_trans(ptr_trans(k+1)) = i
            val_trans(ptr_trans(k+1)) = val(j)
            ptr_trans(k+1) = ptr_trans(k+1) + 1
         end do
      end do

      !
      ! Form A+A^T
      !
      allocate(ptr_out(n+1), row_out(2*ptr(n+1)), val_out(2*ptr(n+1)))
      allocate(seen(n))
      insert = 1
      do i = 1, n
         ptr_out(i) = insert
         seen(:) = 0
         ! Add A entries
         do j = ptr(i), ptr(i+1)-1
            k = row(j)
            if(k.lt.i) cycle ! Skip upper triangle
            row_out(insert) = k
            val_out(insert) = val(j)
            seen(k) = insert
            insert = insert + 1
         end do
         ! Add A^T entries
         do j = ptr_trans(i), ptr_trans(i+1)-1
            k = row_trans(j)
            if(k.lt.i) cycle ! Skip upper triangle
            if(seen(k).eq.0) then
               ! New entry
               row_out(insert) = k
               val_out(insert) = val_trans(j)
               seen(k) = insert
               insert = insert + 1
            else
               ! Repeat entry, add to previous
               val_out(seen(k)) = val_out(seen(k)) + val_trans(j)
            endif
         end do
      end do
      ptr_out(n+1) = insert

      !
      ! Copy *_out over original for output
      !
      deallocate(row, val)
      allocate(row(ptr_out(n+1)-1), val(ptr_out(n+1)-1))
      ptr(1:n+1) = ptr_out(1:n+1)
      row(1:ptr_out(n+1)-1) = row_out(1:ptr_out(n+1)-1)
      val(1:ptr_out(n+1)-1) = val_out(1:ptr_out(n+1)-1)
   end subroutine symmetrize_problem

   subroutine randomize_list(n, list, state)
      integer, intent(in) :: n
      integer, dimension(n), intent(inout) :: list
      type(random_state), intent(inout) :: state

      integer :: i, idx1, idx2, temp

      ! do n random swaps
      do i = 1, n
         idx1 = random_integer(state, n)
         idx2 = random_integer(state, n)
         temp = list(idx1)
         list(idx1) = list(idx2)
         list(idx2) = temp
      end do
   end subroutine randomize_list

   subroutine randomize_matrix(n, ptr, row, seed)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(inout) :: row
      integer, intent(in) :: seed

      integer :: i
      type(random_state) :: state

      call random_set_seed(state, seed)
      do i = 1, n
         call randomize_list(ptr(i+1)-ptr(i), row(ptr(i)), state)
      end do
   end subroutine randomize_matrix

   subroutine calculate_stats(n, ptr, row, perm, nfact, nflops)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      integer, dimension(n), intent(in) :: perm
      integer(long), intent(out) :: nfact
      integer(long), intent(out) :: nflops

      integer :: i, realn, st
      integer, dimension(:), allocatable :: perm2, invp, parent, cc
      integer, dimension(:), allocatable :: sptr, ptr2, row2, iw

      ! Expand to a full matrix
      allocate(ptr2(n+1), row2(2*ptr(n+1)), iw(n), sptr(n+1))
      ptr2(1:n+1) = ptr(1:n+1)
      row2(1:ptr(n+1)-1) = row(1:ptr(n+1)-1)
      call half_to_full(n, row2, ptr2, iw)

      ! Limited analyse phase nemin=1
      allocate(parent(n), invp(n), perm2(n), cc(n+1))
      perm2(:) = perm(:)
      do i = 1, n
         invp(perm(i)) = i
         sptr(i) = i
      end do
      sptr(n+1) = n+1
      call find_etree(n, ptr2, row2, perm2, invp, parent, st)
      if(st.ne.0) goto 10
      call find_postorder(n, realn, ptr2, perm2, invp, parent, st)
      if(st.ne.0) goto 10
      call find_col_counts(n, ptr2, row2, perm2, invp, parent, cc, st)
      if(st.ne.0) goto 10
      call calc_stats(n, sptr, cc, nfact=nfact, nflops=nflops)

      return

      10 continue
      print *, "Allocation error in finding stats"
   end subroutine calculate_stats

   subroutine run_ma87(n, ptr, row, val, order)
      use hsl_ma87_double
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      real(wp), dimension(ptr(n+1)-1), intent(in) :: val
      integer, dimension(n), intent(inout) :: order

      type(ma87_keep) :: keep
      type(ma87_control) :: control
      type(ma87_info) :: info

      real(wp), dimension(:), allocatable :: rhs

      type(timespec) :: t1, t2
      integer :: dummy

      ! Analyse
      dummy = clock_gettime(0, t1)
      call ma87_analyse(n, ptr, row, order, keep, control, info)
      dummy = clock_gettime(0, t2)
      if(info%flag.ne.0) then
         print *, "ma87_analyse() failed with flag ", info%flag
         stop 1
      endif
      print *, "ma87 analyse took ", tdiff(t1, t2)
      print "(a,es10.2)", "ma87 nfact = ", real(info%num_factor)
      print "(a,es10.2)", "ma87 nflops = ", real(info%num_flops)

      ! Factor
      dummy = clock_gettime(0, t1)
      call ma87_factor(n, ptr, row, val, order, keep, control, info)
      dummy = clock_gettime(0, t2)
      if(info%flag.lt.0) then
         print *, "ma87_factor() failed with flag ", info%flag
         stop 1
      endif
      print *, "ma87 factor took ", tdiff(t1, t2)

      ! Solve
      allocate(rhs(n))
      rhs(1:n) = 1.0
      dummy = clock_gettime(0, t1)
      call ma87_solve(rhs, order, keep, control, info)
      dummy = clock_gettime(0, t2)
      if(info%flag.lt.0) then
         print *, "ma87_solve() failed with flag ", info%flag
         stop 1
      endif
      print *, "ma87 solve took ", tdiff(t1, t2)
   end subroutine run_ma87

   subroutine run_ma97(n, ptr, row, val, order, scaling, control)
      use hsl_mc69_double
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      real(wp), dimension(ptr(n+1)-1), intent(in) :: val
      integer, dimension(n), intent(inout) :: order
      real(wp), dimension(:), allocatable, intent(inout) :: scaling
      type(ma97_control), intent(in) :: control

      type(ma97_akeep) :: akeep
      type(ma97_fkeep) :: fkeep
      type(ma97_info) :: info

      type(hungarian_options) :: hoptions
      type(hungarian_inform) :: hinform

      integer :: i, j, k, irit
      real(wp), dimension(:), allocatable :: rhs, soln, dx
      real(wp) :: res

      type(timespec) :: t1, t2
      integer :: dummy

      real(wp), parameter :: conv_tol = 1d-14

      ! Make up a rhs
      allocate(rhs(n), soln(n))
      rhs(:) = 0
      do i = 1, n
         do j = ptr(i), ptr(i+1)-1
            k = row(j)
            rhs(k) = rhs(k) + val(j)
            if(i.eq.k) cycle
            rhs(i) = rhs(i) + val(j)
         end do
      end do

      if(.not.allocated(scaling)) then
         ! Generate MC64 scaling
         allocate(scaling(n))
         call hungarian_scale_sym(n, ptr, row, val, scaling, hoptions, hinform)
      endif

      ! Analyse
      dummy = clock_gettime(0, t1)
      call ma97_analyse(.false., n, ptr, row, akeep, control, info, order=order)
      dummy = clock_gettime(0, t2)
      if(info%flag.ne.0) then
         print *, "ma97_analyse() failed with flag ", info%flag
         stop 1
      endif
      print *, "ma97 analyse took ", tdiff(t1, t2)
      print "(a,es10.2)", "ma97 afact = ", real(info%num_factor)
      print "(a,es10.2)", "ma97 aflops = ", real(info%num_flops)

      ! Factor
      dummy = clock_gettime(0, t1)
      call ma97_factor(HSL_MATRIX_REAL_SYM_INDEF, val, akeep, fkeep, control, &
         info, ptr=ptr, row=row, scale=scaling)
      dummy = clock_gettime(0, t2)
      if(info%flag.lt.0) then
         print *, "ma97_factor() failed with flag ", info%flag
         stop 1
      endif
      print *, "ma97 factor took ", tdiff(t1, t2)
      print "(a,es10.2)", "ma97 ffact = ", real(info%num_factor)
      print "(a,es10.2)", "ma97 fflops = ", real(info%num_flops)
      print *, "ma97 ndelay = ", info%num_delay

      ! Solve
      soln(:) = rhs(:)
      dummy = clock_gettime(0, t1)
      call ma97_solve(soln, akeep, fkeep, control, info)
      dummy = clock_gettime(0, t2)
      if(info%flag.lt.0) then
         print *, "ma97_solve() failed with flag ", info%flag
         stop 1
      endif
      print *, "ma97 solve took ", tdiff(t1, t2)

      print *, "number bad cmp = ", count(abs(soln(1:n)-1.0).ge.1e-6)
      print *, "fwd error || ||_inf = ", maxval(abs(soln(1:n)-1.0))
      call internal_calc_norm(n, ptr, row, val, soln, rhs, res)
      print *, "bwd error scaled = ", res
      print *, "Inertia = ", info%matrix_rank - info%num_neg, &
         info%num_neg, n-info%matrix_rank

      ! Iterative Refinement
      irit = 1
      if(res.ge.conv_tol) then
         allocate(dx(n))
         do irit = 2, 10
            print *, "It ref iteration ", irit
            call spmv(n, ptr, row, val, soln, dx)
            dx = rhs(:) - dx
            call ma97_solve(dx,akeep,fkeep,control,info)
            soln(:) = soln(:) + dx(:)
            print *, "fwd error || ||_inf = ", maxval(abs(soln(1:n)-1.0))
            call internal_calc_norm(n, ptr, row, val, soln, rhs, res)
            print *, "bwd error scaled = ", res
            if(res.lt.conv_tol) exit
         end do
      endif
      if(res.lt.conv_tol .and. maxval(abs(soln(1:n)-1.0)).lt.1e10) then
         print *, "IR converged in ", irit, " iterations"
      else
         print *, "IR failed to converge"
      endif
   end subroutine run_ma97

   subroutine internal_calc_norm(n, ptr, row, val, x_vec, b_vec, res)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      real(wp), dimension(ptr(n+1)-1), intent(in) :: val
      real(wp), dimension(n), intent(in) :: x_vec
      real(wp), dimension(n), intent(in) :: b_vec
      real(wp), intent(out) :: res

      integer :: i, j, r
      real(wp), dimension(:), allocatable :: res_vec
      real(wp) :: temp, x_norm, normA

      ! Find the residual
      allocate(res_vec(n))
      res_vec = 0
      do i = 1, n
         do j = ptr(i), ptr(i+1)-1
            r = row(j)
            res_vec(i) = res_vec(i) + val(j) * x_vec(r)
            if(r.eq.i) cycle
            res_vec(r) = res_vec(r) + val(j) * x_vec(i)
         end do
      end do
      res_vec(:) = res_vec(:) - b_vec(:)

      ! Find matrix norm
      call matrix_inf_norm(n, ptr, row, val, normA)

      ! Find x norm
      x_norm = 0
      do j =1, n
         x_norm = max(x_norm, abs(x_vec(j)))
         if(x_vec(j).ne.x_vec(j)) then ! Tests for NaN
            x_norm = x_vec(j)
            exit
         endif
      end do

      print *, "||r|| = ", maxval(abs(res_vec(1:n)))
      print *, "||A|| = ", normA
      print *, "||x|| = ", x_norm
      print *, "||b|| = ", maxval(abs(b_vec(1:n)))

      ! Scaled residual = ||r|| / ( ||A|| ||x|| + ||b|| )
      temp = normA * x_norm + maxval(abs(b_vec(1:n)))
      if(temp .eq. 0) then
         res = maxval(abs(res_vec(1:n)))
      else
         res = maxval(abs(res_vec(1:n))) / temp
      endif
   end subroutine internal_calc_norm

   subroutine matrix_inf_norm(n, ptr, row, val, norm)
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

   ! Calculates y = Ax
   subroutine spmv(n, ptr, row, val, x, y)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      real(wp), dimension(ptr(n+1)-1), intent(in) :: val
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: y

      integer :: i, j, r

      y(:) = 0
      do i = 1, n
         do j = ptr(i), ptr(i+1)-1
            r = row(j)
            y(i) = y(i) + val(j) * x(r)
            if(r.eq.i) cycle
            y(r) = y(r) + val(j) * x(i)
         end do
      end do
   end subroutine spmv

!   subroutine find_gupta_order(n, ptr, row, val, scaling, order, u)
!      use spral_nd_numaware
!      integer, intent(in) :: n
!      integer, dimension(n+1), intent(in) :: ptr
!      integer, dimension(ptr(n+1)-1), intent(in) :: row
!      real(wp), dimension(ptr(n+1)-1), intent(in) :: val
!      real(wp), dimension(n), intent(out) :: scaling
!      integer, dimension(n), intent(out) :: order
!      real(wp), intent(in) :: u
!
!      integer :: i, k
!      integer :: jj
!      integer, dimension(:), allocatable :: bigflag, iw, match, ptr2, row2
!      real(wp), dimension(:), allocatable :: val2
!
!      type(hungarian_options) :: hoptions
!      type(hungarian_inform) :: hinform
!
!      ! Find and apply MC64 scaling
!      call hungarian_scale_sym(n, ptr, row, val, scaling, hoptions, hinform)
!      allocate(val2(2*ptr(n+1)-1))
!      do i = 1, n
!         do jj = ptr(i), ptr(i+1)-1
!            k = row(jj)
!            val2(jj) = scaling(i) * val(jj) * scaling(k)
!         end do
!      end do
!
!      ! Expand matrix
!      allocate(ptr2(n+1), row2(2*ptr(n+1)-1), iw(n))
!      ptr2(1:n+1) = ptr(1:n+1)
!      row2(1:ptr(n+1)-1) = row(1:ptr(n+1)-1)
!      call half_to_full(n, row2, ptr2, iw, a=val2)
!
!      ! Determine flags array
!      allocate(bigflag(ptr2(n+1)-1), match(n))
!      match(:) = -1 ! set everything unmatched so no artifically "big" entries
!      call nd_set_a_flags(u, n, ptr2, row2, val2, match, bigflag)
!
!      ! Determine a matching
!      call find_matching(n, ptr2, row2, val2, bigflag, match)
!
!      ! Obtain ordering from call of metis on compressed graph
!      call compress_metis_call(n, ptr2, row2, match, order)
!   end subroutine

   ! We expect both lwr and upr parts to be passed to the below
   !
   ! match(:) has the following values:
   ! -2 unmatched
   ! -1 matched as singleton
   !  0 not yet seen
   ! >0 matched with specified node
   subroutine find_matching(n, ptr, row, val, bigflag, match)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      real(wp), dimension(ptr(n+1)-1), intent(in) :: val
      integer, dimension(ptr(n+1)-1), intent(in) :: bigflag
      integer, dimension(n), intent(out) :: match

      integer :: i, k
      integer(long) :: jj, pp
      integer :: nz_extra, best_idx
      integer, dimension(:), allocatable :: pattern
      real(wp) :: score, best_score

      allocate(pattern(n))
      pattern(:) = 0

      match(:) = -2 ! Initially all unmatched
      do i = 1, n
         if(match(i).ne.-2) cycle ! Already matched
         ! Find diagonal and check if sufficiently large
         ! Otherwise create pattern(:) of column
         do jj = ptr(i), ptr(i+1)-1
            k = row(jj)
            if(i.eq.k .and. bigflag(jj).ne.FLAG_SMALL) then
               match(i) = -1 ! Matched as singleton on diagonal
               exit ! No need to work on rest of column
            endif
            pattern(k) = i
         end do
         if(match(i).eq.-1) cycle ! Diagonal was large enough
         ! Build a list of candidate merge columns
         ! Score them based on number of extra entries merging will create
         best_idx = -1
         best_score = -huge(best_score)
         do jj = ptr(i), ptr(i+1)-1
            if(bigflag(jj).ne.FLAG_BIG_COL .and. bigflag(jj).ne.FLAG_BIG_BOTH) cycle ! skip small entry
            k = row(jj)
            if(match(k).gt.0) cycle ! k already matched
            ! Count number of extra entries from merging columns i and k
            nz_extra = int(ptr(i+1)-ptr(i)+1) + int(ptr(k+1)-ptr(k)+1)
            do pp = ptr(k), ptr(k+1)-1
               if(pattern(row(pp)) .ge. i) & ! Overlapping entry
                  nz_extra = nz_extra - 2
            end do
            score = abs(val(jj)) / nz_extra
            if(score.lt.best_score) then
               best_score = score
               best_idx = k
            endif
         end do
         if(best_idx.ne.-1) then
            match(i) = best_idx
            match(best_idx) = i
         endif
      end do
   end subroutine find_matching

   subroutine compress_metis_call(n, ptr, row, match, order)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(ptr(n+1)-1), intent(in) :: row
      integer, dimension(n), intent(in) :: match
      integer, dimension(n), intent(out) :: order

      integer :: i, j, j1, j2, jj, k, krow, metis_flag, stat
      integer(long) :: klong
      integer :: ncomp, ncomp_matched
      integer, dimension(:), allocatable :: old_to_new, new_to_old
      integer, dimension(:), allocatable :: ptr3
      integer, dimension(:), allocatable :: iwork, row3, invp

      do i = 1, n
         if(match(i).gt.0) then
            if(match(i) .gt. n) then
               print *, "match(", i, ") = ", match(i), " > n = ", n
               stop
            endif
            if(match(match(i)).ne.i) then
               print *, "match(match(", i, ")=", match(i), ") = ", match(match(i))
               stop
            endif
         endif
      end do

      !
      ! Build maps for new numbering schemes
      !
      allocate(old_to_new(n), new_to_old(n))
      k = 1
      do i = 1, n
         j = match(i)
         if (j<i .and. j.gt.0) cycle
         old_to_new(i) = k
         new_to_old(k) = i ! note: new_to_old only maps to first of a pair
         if (j.gt.0) old_to_new(j) = k   
         k = k + 1
      end do
      ncomp_matched = k-1

      !
      ! Produce a condensed version of the matrix for ordering.
      ! Hold pattern using ptr3 and row3.
      !
      allocate(ptr3(ncomp_matched+1), row3(ptr(n+1)-1), iwork(n))
      ptr3(1) = 1
      iwork(:) = 0 ! Use to indicate if entry is in a paired column
      ncomp = 1
      jj = 1
      do i = 1, n
         j = match(i)
         if (j<i .and. j.gt.0) cycle ! already seen
         do klong = ptr(i), ptr(i+1)-1
            krow = old_to_new(row(klong))
            if (iwork(krow).eq.i) cycle ! already added to column
            if (krow>ncomp_matched) cycle ! unmatched row not participating
            row3(jj) = krow
            jj = jj + 1
            iwork(krow) = i
         end do
         if (j.gt.0) then
            ! Also check column match(i)
            do klong = ptr(j), ptr(j+1)-1
               krow = old_to_new(row(klong))
               if (iwork(krow).eq.i) cycle ! already added to column
               if (krow>ncomp_matched) cycle ! unmatched row not participating
               row3(jj) = krow
               jj = jj + 1
               iwork(krow) = i
            end do
         end if
         ptr3(ncomp+1) = jj
         ncomp = ncomp + 1
      end do
      ncomp = ncomp - 1

      ! store just lower triangular part for input to hsl_mc68
      ptr3(1) = 1
      jj = 1
      j1 = 1
      do i = 1, ncomp
         j2 = ptr3(i+1)
         do k = j1, j2-1
            krow = row3(k)
            if ( krow.lt.i ) cycle ! already added to column
            row3(jj) = krow
            jj = jj + 1
         end do
         ptr3(i+1) = jj
         j1 = j2
      end do

      allocate(invp(ncomp))

      ! reorder the compressed matrix using metis.
      ! switch off metis printing
      call metis_order(ncomp,ptr3,row3,order,invp,metis_flag,stat)
      select case(metis_flag)
      case(0)
         ! OK, do nothing
      case default
         ! Unknown error, should never happen
         print *, "metis_order() returned unknown error ", metis_flag
         stop
      end select

      do i = 1, ncomp
         j = order(i)
         iwork(j) = i
      end do

      !
      ! Translate inverse permutation in iwork back to 
      ! permutation for original variables.
      !
      k = 1
      do i = 1, ncomp
         j = new_to_old( iwork(i) )
         order(j) = k
         k = k + 1
         if (match(j).gt.0) then
            j = match(j)
            order(j) = k
            k = k + 1
         end if
      end do
      
   end subroutine compress_metis_call



end program
