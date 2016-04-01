program run_prob
   use spral_core_analyse
   use spral_matrix_util
   use spral_metis_wrapper
   use spral_nd
   use spral_random
   use spral_rutherford_boeing
   use spral_timer
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
   integer, dimension(:), allocatable :: perm, invp
   integer(long) :: nfact, nflops

   ! Timing
   type(timespec) :: t1, t2
   integer :: dummy

   ! Controls
   integer :: random, ndmethod
   logical :: with_metis, with_nd, with_ma87

   call proc_args(options, with_metis, with_nd, ndmethod, random, with_ma87)

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
   if(with_nd) then
      write(*, "(a)", advance="no") "Ordering with ND..."
      dummy = clock_gettime(0, t1)
      call nd_order(ndmethod, 0, n, ptr, row, perm, options, inform)
      dummy = clock_gettime(0, t2)
      if (inform%flag < 0) then
         print *, "oops on analyse ", inform%flag
         stop
      endif
      write(*, "(a)") "ok"
      print *, "nd_order() took ", tdiff(t1, t2)
      ! Determine quality
      write(*, "(a)", advance="no") "Determing stats..."
      call calculate_stats(n, ptr, row, perm, nfact, nflops)
      write(*, "(a)") "ok"
      print "(a,es10.2)", "nd nfact = ", real(nfact)
      print "(a,es10.2)", "nd nflop = ", real(nflops)
      print "(a,i10)", "nd ndense = ", inform%dense
      print "(a,i10)", "nd nsuper = ", inform%nsuper
      print "(a,i10)", "nd nzsuper = ", inform%nzsuper
      print "(a,i10)", "nd ncomp = ", inform%num_components
      print "(a,i10)", "nd n_max_component = ", inform%n_max_component
      print "(a,i10)", "nd nz_max_component = ", inform%nz_max_component
      print "(a,f10.2)", "nd band = ", inform%band
      if(with_ma87) call run_ma87(n, ptr, row, val, perm)
   endif

   ! Order using metis
   if(with_metis) then
      write(*, "(a)", advance="no") "Ordering with Metis..."
      dummy = clock_gettime(0, t1)
      call metis_order(n, ptr, row, perm, invp, flag, st)
      dummy = clock_gettime(0, t2)
      if (inform%flag < 0) then
         print *, "oops on analyse ", inform%flag
         stop
      endif
      write(*, "(a)") "ok"
      print *, "metis_order() took ", tdiff(t1, t2)
      ! Determine quality
      write(*, "(a)", advance="no") "Determing stats..."
      call calculate_stats(n, ptr, row, perm, nfact, nflops)
      write(*, "(a)") "ok"
      print "(a,es10.2)", "metis nfact = ", real(nfact)
      print "(a,es10.2)", "metis nflop = ", real(nflops)
      if(with_ma87) call run_ma87(n, ptr, row, val, perm)
   endif

contains

   subroutine proc_args(options, with_metis, with_nd, ndmethod, random, with_ma87)
      type(nd_options), intent(inout) :: options
      logical, intent(out) :: with_metis
      logical, intent(out) :: with_nd
      integer, intent(out) :: ndmethod
      integer, intent(out) :: random
      logical, intent(out) :: with_ma87

      integer :: argnum, narg
      character(len=200) :: argval

      ! Defaults
      with_metis = .false.
      with_nd = .true.
      random = -1
      ndmethod = 0 ! Non-multilevel
      with_ma87 = .false.
      
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
         case("--metis")
            with_metis = .true.
            print *, "MeTiS run requested"
         case("--nond")
            with_nd = .false.
            print *, "ND run disabled"
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
            read( argval, * ) ndmethod
            print *, "Set ndmethod = ", ndmethod
         case("--non-multilevel")
            ndmethod = 0
            print *, "Set ndmethod = ", ndmethod, "(non-multilevel)"
         case("--multilevel")
            ndmethod = 1
            print *, "Set ndmethod = ", ndmethod, "(multilevel)"
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

end program
