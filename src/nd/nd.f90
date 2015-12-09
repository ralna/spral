module spral_nd
   use spral_amd
   use spral_nd_preprocess
   use spral_nd_types
   implicit none

   private
   public :: nd_order, nd_refine_fm
   public :: nd_options, nd_inform

   ! ---------------------------------------------------
   ! Derived type definitions
   ! ---------------------------------------------------


   type nd_matrix
      integer :: m ! number rows
      integer :: n ! number columns
      integer :: ne ! number entries in matrix
      integer, allocatable, dimension(:) :: ptr ! pointer into col array
      integer, allocatable, dimension(:) :: col ! column indices
      integer, allocatable, dimension(:) :: val ! values
   end type nd_matrix

   ! *****************************************************************

   type nd_multigrid
      integer :: size ! size of this level (number of rows)
      type (nd_matrix), allocatable :: graph ! this level of matrix
      integer, allocatable, dimension(:) :: where ! where each row of
         ! this level of matrix will go (ie ordering for this level)
      integer, allocatable, dimension(:) :: row_wgt ! number of
         ! vertices this vertex of the coarse graph matrix represents
      integer :: level = 0 ! the level
      integer :: part_div(2) ! number of vertices in each part
      type (nd_multigrid), pointer :: coarse => null() ! child coarse grid
         ! (NB: owned by this instance)
      type (nd_multigrid), pointer :: fine => null() ! pointer to parent fine
         ! grid (NB: owns this instance)
      type (nd_matrix), allocatable :: p ! the prolongation operator
      type (nd_matrix), allocatable :: r ! the restriction operator
   end type nd_multigrid

   ! *****************************************************************

   ! This is a version of maxflow with no bandgraph and assumption that
   ! there are no duplicated edges.
   type network
      integer :: nnode ! number of nodes in the network
      integer :: narc ! number of arcs in the network
      integer :: source ! source node = 1
      integer :: sink ! sink node = nnode
      integer, allocatable :: inheads(:) ! for each node u first incoming arc
         ! for u
      integer, allocatable :: outheads(:) ! for each node u first outgoing arc
         ! for u
      integer, allocatable :: firsts(:) ! for each arc e first node for arc e
      integer, allocatable :: seconds(:) ! for each arc e second node for arc e
      integer, allocatable :: capacities(:) ! for each arc e capacity of arc e
      integer, allocatable :: flows(:) ! for each arc e flow through arc e
      integer, allocatable :: nextin(:) ! for each arc e next incoming arc
      integer, allocatable :: nextout(:) ! for each arc e next outgoing arc
   end type network

contains
   !
   ! Main user callable routine
   !
   subroutine nd_order(mtx,n,ptr,row,perm,options,info,seps)
      integer, intent(in) :: mtx ! 0 if lower triangular part matrix input,
         ! 1 if both upper and lower parts input
      integer, intent(in) :: n ! number of rows in the matrix
      integer, intent(in) :: ptr(n+1) ! column pointers
      integer, intent(in) :: row(ptr(n+1)-1) ! row indices
      integer, intent(out) :: perm(n) ! permutation: row i becomes row
         ! perm(i)
      type (nd_options), intent(in) :: options
      type (nd_inform), intent(inout) :: info
      integer, dimension(n), optional, intent(out) :: seps
         ! seps(i) is -1 if vertex is not in a separator; otherwise it
         ! is equal to lev, where lev is the nested dissection level at
         ! which it became part of the separator

      ! Throughout the code we will use a modified compressed sparse
      ! column/row format (same due to symmetry). However there is no
      ! a_ptr(n+1) entry, instead the number of non-zeroes is stored as
      ! a_ne. These variables store that information.
      integer :: a_n
      integer :: a_ne
      integer, allocatable, dimension(:) :: a_ptr
      integer, allocatable, dimension(:) :: a_row

      ! Other local variables
      integer :: i
      integer :: unit_error
      integer :: unit_diagnostics
      logical :: printe, printi, printd

      ! ---------------------------------------------
      ! Printing levels
      unit_diagnostics = options%unit_diagnostics
      unit_error = options%unit_error
      printe = (options%print_level>=0 .and. unit_error>=0)
      printi = (options%print_level==1 .and. unit_diagnostics>=0)
      printd = (options%print_level>=2 .and. unit_diagnostics>=0)
      ! ---------------------------------------------------

      if (printi .or. printd) then
         write (unit_diagnostics,'(a)') ' '
         write (unit_diagnostics,'(a)') 'nd_order:'
      end if

      ! Error checks
      if (n<1) then
         info%flag = ND_ERR_N
         if (printe) call nd_print_message(info%flag,unit_error, &
            'nd_order')
         return
      end if

      ! Convert matrix to internal format without diagonals
      if (mtx<1) then
         call construct_full_from_lower(n, ptr, row, a_n, a_ne, a_ptr, a_row, &
            options, info%stat)
      else
         call construct_full_from_full(n, ptr, row, a_n, a_ne, a_ptr, a_row, &
            options, info%stat)
      endif
      if(info%stat.ne.0) then
         info%flag = ND_ERR_MEMORY_ALLOC
         if (printe) &
            call nd_print_message(info%flag,unit_error, 'nd_order')
         return
      endif

      ! Output summary of input matrix post-conversion
      if (printd) then
         ! Print out a_ptr and a_row
         write (unit_diagnostics,'(a8)') 'a_ptr = '
         write (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n)
         write (unit_diagnostics,'(a8)') 'a_row = '
         write (unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne)
      else if (printi) then
         ! Print out first few entries of a_ptr and a_row
         write (unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n)) = '
         write (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,min(5,a_n))
         write (unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne)) = '
         write (unit_diagnostics,'(5i15)') (a_row(i),i=1,min(5,a_ne))
      end if

      ! Call main worker routine
      call nd_nested_both(a_n,a_ne,a_ptr,a_row,perm,options,info,seps)

      ! Output summary of results
      if (printd) then
         ! Print out perm
         write (unit_diagnostics,'(a8)') 'perm = '
         write (unit_diagnostics,'(5i15)') (perm(i),i=1,n)
      else if (printi) then
         ! Print out first few entries of perm
         write (unit_diagnostics,'(a21)') 'perm(1:min(5,n)) = '
         write (unit_diagnostics,'(5i15)') (perm(i),i=1,min(5,n))
      end if

      info%flag = 0
      if (printi .or. printd) &
         call nd_print_message(info%flag,unit_diagnostics,'nd_order')

   end subroutine nd_order



   !
   ! Main wrapper routine
   !
   subroutine nd_nested_both(a_n,a_ne,a_ptr,a_row,perm,options,info,seps)
      ! Expects input in internal CSC format
      integer, intent(inout) :: a_n
      integer, intent(inout) :: a_ne
      integer, dimension(a_n), intent(inout) :: a_ptr
      integer, dimension(a_ne), intent(inout) :: a_row
      integer, dimension(a_n), intent(out) :: perm ! row i becomes row perm(i)
      type (nd_options), intent(in) :: options
      type (nd_inform), intent(inout) :: info
      integer, dimension(a_n), optional, intent(out) :: seps
         ! seps(i) is -1 if vertex is not in a separator; otherwise it
         ! is equal to l, where l is the nested dissection level at
         ! which it became part of the separator

      integer :: i, j, k, l, ll, lwork, lirn
      integer :: a_n_orig, a_n_curr, a_ne_curr, num_zero_row
      integer, dimension(:), allocatable :: amd_order_row, amd_order_ptr, &
         amd_order_sep, amd_order_perm, amd_order_work, amd_order_iperm
      integer, dimension(:), allocatable :: work_iperm, work_seps
      integer :: unit_error
      integer :: unit_diagnostics
      integer :: nsvar
      integer :: sumweight
      integer, dimension(:), allocatable :: svar, sinvp ! supervariable info
      integer, allocatable, dimension(:) :: a_weight ! a_weight(i) will
         ! contain the weight of variable (column) i ends for the
         ! expanded matrix
      integer, allocatable, dimension(:) :: iperm ! inverse of perm(:)
      integer, allocatable, dimension(:) :: work ! space for doing work
      logical :: printe, printi, printd
      logical :: use_multilevel
      type (nd_multigrid) :: grid

      ! ---------------------------------------------
      ! Printing levels
      unit_diagnostics = options%unit_diagnostics
      unit_error = options%unit_error
      printe = (options%print_level>=0 .and. unit_error>=0)
      printi = (options%print_level==1 .and. unit_diagnostics>=0)
      printd = (options%print_level>=2 .and. unit_diagnostics>=0)
      ! ---------------------------------------------------

      if (printi .or. printd) then
         write (unit_diagnostics,'(a)') ' '
         write (unit_diagnostics,'(a)') 'nd_nested_both:'
      end if

      ! Record original size for alter use
      a_n_orig = a_n

      ! Allocate iperm and initialize to identity
      allocate (iperm(a_n),stat=info%stat)
      if (info%stat/=0) go to 10
      iperm(:) = (/ (i,i=1,a_n) /)

      ! Initialize all variables to not be in a seperator
      if (present(seps)) seps(:) = -1

      ! Remove any dense rows from matrix and modify iperm (if enabled)
      if (options%remove_dense_rows) &
         call remove_dense_rows(a_n, a_ne, a_ptr, a_row, iperm, options, info)

      ! Return if matrix is diagonal
      if (a_ne.eq.0) then ! (recall diagonal entries are not stored)
         if (printi .or. printd) then
            write (unit_diagnostics,'(a)') ' '
            write (unit_diagnostics,'(a)') 'Matrix is diagonal'
         end if

         info%nsuper = a_n
         info%nzsuper = 0
         ! Create perm from iperm
         do i = 1, a_n_orig
            j = iperm(i)
            perm(j) = i
         end do

         return
      end if

      allocate(a_weight(a_n), stat=info%stat)
      if (info%stat.ne.0) go to 10

      if (options%find_supervariables) then
         allocate(svar(a_n), sinvp(a_n), stat=info%stat)
         if (info%stat.ne.0) go to 10
         call compress_by_svar(a_n, a_ne, a_ptr, a_row, a_weight, a_n_curr, &
            a_ne_curr, nsvar, svar, sinvp, num_zero_row, options, info%stat)
         if (info%stat.ne.0) go to 10
      else
         ! Otherwise, set things up for consistency with lack of supervariables
         a_n_curr = a_n
         a_ne_curr = a_ne
         nsvar = a_n
         num_zero_row = 0
         a_weight(1:a_n_curr) = 1
      end if

      ! Proceed to main ordering now preprocessing is complete
      info%nzsuper = a_ne_curr
      info%nsuper = a_n_curr

      if (options%amd_switch2.le.0 .or. &
            a_n_curr.le.max(2,max(options%amd_call,options%amd_switch1)) ) then
         ! Apply AMD to matrix
         ! Allocate work to have length 5*a_n_curr+a_ne_curr
         if (printd) &
            write (unit_diagnostics,'(a)') ' Form AMD ordering'
         lirn = a_ne_curr + a_n_curr
         allocate(amd_order_ptr(a_n_curr), amd_order_row(lirn), &
            amd_order_sep(a_n_curr), amd_order_perm(a_n_curr), &
            amd_order_work(7*a_n_curr), amd_order_iperm(a_n))
         if (info%stat/=0) go to 10

         amd_order_ptr(1:a_n_curr) = a_ptr(1:a_n_curr)
         amd_order_row(1:a_ne_curr) = a_row(1:a_ne_curr)
         amd_order_row(a_ne_curr+1:lirn) = 0
         amd_order_sep(:) = ND_SEP_FLAG - 1
         call amd_order(a_n_curr, a_ne_curr, lirn, amd_order_ptr, &
            amd_order_row, amd_order_sep, amd_order_perm, amd_order_work)

         ! Extract perm from amd_order to apply to iperm
         if (nsvar+num_zero_row.eq.a_n) then
            do i = 1, a_n_curr
               j = amd_order_perm(i)
               amd_order_work(i) = iperm(j)
            end do
            iperm(1:a_n_curr) = amd_order_work(1:a_n_curr)
         else
            do i = 1, a_n_curr
               j = amd_order_perm(i)
               amd_order_work(i) = j
            end do
            ! Expand to matrix before supervariables detected
            k = 1
            do i = 1, a_n_curr
               j = amd_order_work(i)
               if (j.eq.1) then
                  ll = 1
               else
                  ll = svar(j-1) + 1
               end if
               do l = ll, svar(j)
                  amd_order_iperm(k) = iperm(sinvp(l))
                  k = k + 1
               end do
            end do

            iperm(1:a_n) = amd_order_iperm(1:a_n)
         end if

      else
         ! Apply ND to matrix

         if (printd) &
            write (unit_diagnostics,'(a)') ' Form ND ordering'

         if (nsvar+num_zero_row.ne.a_n) then
            ! Create shadow versions of iperm and seps that work on
            ! supervariables rather than variables
            allocate(work_iperm(a_n), stat=info%stat)
            if (info%stat.ne.0) go to 10
            work_iperm(1:a_n_curr) = (/ (i,i=1,a_n_curr) /)
            if (present(seps)) then
               allocate(work_seps(a_n_curr), stat=info%stat)
               if (info%stat.ne.0) go to 10
               work_seps(1:a_n_curr) = -1
            end if
         end if
         ! Allocate a workspace that can be reused at lower levels
         allocate (work(a_n+14*a_n_curr+a_ne_curr), stat=info%stat)
         if (info%stat/=0) go to 10

         use_multilevel = .true.
         if (nsvar+num_zero_row==a_n) then
            sumweight = sum(a_weight(1:a_n_curr))
            lwork = 12*a_n_curr + sumweight + a_ne_curr
            if (present(seps)) then
               call nd_nested_internal(a_n_curr, a_ne_curr, a_ptr(1:a_n_curr), &
                  a_row(1:a_ne_curr), a_weight(1:a_n_curr), sumweight,         &
                  iperm(1:a_n_curr), work(1:lwork),                            &
                  work(lwork+1:lwork+a_n_curr),                                &
                  work(lwork+a_n_curr+1:lwork+2*a_n_curr), 0, options, info,   &
                  .false., use_multilevel, grid, seps=seps(1:a_n_curr))
            else
               call nd_nested_internal(a_n_curr, a_ne_curr, a_ptr(1:a_n_curr), &
                  a_row(1:a_ne_curr), a_weight(1:a_n_curr), sumweight,         &
                  iperm(1:a_n_curr), work(1:lwork),                            &
                  work(lwork+1:lwork+a_n_curr),                                &
                  work(lwork+a_n_curr+1:lwork+2*a_n_curr), 0, options, info,   &
                  .false.,use_multilevel,grid)

            end if
         else
            sumweight = sum(a_weight(1:a_n_curr))
            lwork = 12*a_n_curr + sumweight + a_ne_curr
            if (present(seps)) then
               call nd_nested_internal(a_n_curr, a_ne_curr, a_ptr(1:a_n_curr), &
                  a_row(1:a_ne_curr), a_weight(1:a_n_curr), sumweight,         &
                  work_iperm, work(1:lwork), work(lwork+1:lwork+a_n_curr),     &
                  work(lwork+a_n_curr+1:lwork+2*a_n_curr), 0, options, info,   &
                  .false., use_multilevel, grid,                               &
                  seps=work_seps)
            else
               call nd_nested_internal(a_n_curr, a_ne_curr, a_ptr(1:a_n_curr), &
                  a_row(1:a_ne_curr), a_weight(1:a_n_curr), sumweight,         &
                  work_iperm, work(1:lwork), work(lwork+1:lwork+a_n_curr),     &
                  work(lwork+a_n_curr+1:lwork+2*a_n_curr), 0, options, info,   &
                  .false., use_multilevel, grid)
            end if
         end if

         if (grid%level==1) call mg_grid_destroy(grid,info%flag)

         if (nsvar+num_zero_row.eq.a_n) then
            if (present(seps)) then
               ! reorder seps
               do i = 1, a_n_curr
                  j = iperm(i)
                  work(j) = seps(i)
               end do
               seps(1:a_n_curr) = work(1:a_n_curr)
            end if
         else
            if (present(seps)) then
               ! Expand and reorder seps
               do i = 1, a_n_curr
                  j = work_iperm(i)
                  if (j==1) then
                     ll = 1
                  else
                     ll = svar(j-1) + 1
                  end if
                  do l = ll, svar(j)
                     seps(sinvp(l)) = work_seps(i)
                  end do
               end do
            end if

            ! Expand iperm to matrix before supervariables detected
            k = a_n
            do i = a_n_curr, 1, -1
               j = work_iperm(i)
               if (j==1) then
                  ll = 1
               else
                  ll = svar(j-1) + 1
               end if
               do l = ll, svar(j)
                  work_iperm(k) = iperm(sinvp(l))
                  k = k - 1
               end do
            end do

            iperm(1:a_n) = work_iperm(1:a_n)

         end if
      end if

      ! Create perm from iperm
      do i = 1, a_n_orig
         j = iperm(i)
         perm(j) = i
      end do

      info%flag = 0
      if (printi .or. printd) &
         call nd_print_message(info%flag, unit_diagnostics, 'nd_nested_both')
      return

      10 continue
      info%flag = ND_ERR_MEMORY_ALLOC
      if (printe) &
         call nd_print_message(info%flag, unit_error, 'nd_nested_both')
      return
   end subroutine nd_nested_both


   ! ---------------------------------------------------
   ! nd_nested_internal
   ! ---------------------------------------------------
   ! Does the recursive nested dissection

   recursive subroutine nd_nested_internal(a_n,a_ne,a_ptr,a_row, &
       a_weight,sumweight,iperm,work,work_comp_n,work_comp_nz,level, &
       options,info,use_amd,use_multilevel,grid,seps)

     integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
     integer, intent(in) :: a_ne ! no. nonzeros of subproblem
     integer, intent(inout) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start. This is then
     ! used to hold positions for submatrices after partitioning
     integer, intent(inout) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.This is then used to hold row indices for
     ! submatrices after partitioning
     integer, intent(inout) :: a_weight(a_n) ! On input a_weight(i)
     ! contains
     ! weight of column i. This is then
     ! used to hold the weights for submatrices after partitioning
     integer, intent(in) :: sumweight ! sum entries in a_weight
     ! (unchanged)
     integer, intent(inout) :: iperm(a_n) ! On input, iperm(i) contains
     ! the
     ! row in the original matrix (when nd_nested was called) that
     ! row i in this sub problem maps to. On output, this is updated to
     ! reflect the computed permutation.
     integer, intent(out) :: work_comp_n(a_n)
     integer, intent(out) :: work_comp_nz(a_n)
     integer, intent(out) :: work(12*a_n+sumweight+a_ne) ! Used during the
     ! algorithm to reduce need for allocations. The output is garbage.
     integer, intent(in) :: level ! which level of nested dissection is
     ! this
     type (nd_options), intent(in) :: options
     type (nd_inform), intent(inout) :: info
     logical, intent(in) :: use_amd
     logical, intent(inout) :: use_multilevel
     integer, intent(inout), optional :: seps(a_n)
     ! seps(i) is -1 if vertex i of permuted submatrix is not in a
     ! separator; otherwise it is equal to l, where l is the nested
     ! dissection level at which it became part of the separator
     type (nd_multigrid), intent(inout) :: grid

     ! ---------------------------------------------
     ! Local variables
     integer :: i, j, k, l, m, s
     integer :: unit_diagnostics ! unit on which to print diagnostics
     integer :: num_components ! Number of independent components found
     integer :: compwork, lwork
     integer :: sumweight_sub, maxdeg_max_component
     integer :: offset_ptr, offset_row
     integer :: a_n1, a_n2, a_ne1, a_ne2
     integer :: amd_order_irn, amd_order_ip, amd_order_sep, amd_order_perm, amd_order_work, lirn
     logical :: printi, printd, use_amdi

     ! ---------------------------------------------
     ! Printing levels
     unit_diagnostics = options%unit_diagnostics
     printi = (options%print_level==1 .and. unit_diagnostics>=0)
     printd = (options%print_level>=2 .and. unit_diagnostics>=0)
     use_amdi = .false.

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a,i6)') 'Nested dissection level ', level
     end if

     ! Check whether matrix is diagonal and act accordingly
     if (a_ne==0) then
       if (printi .or. printd) then
         write (unit_diagnostics,'(a)') ' '
         write (unit_diagnostics,'(a)') 'Submatrix is diagonal'
       end if
       return
     end if

     if (level==0) then
       maxdeg_max_component = a_ne + 1 - a_ptr(a_n)
       do i = 1, a_n - 1
         if (maxdeg_max_component<a_ptr(i+1)-a_ptr(i)) &
           maxdeg_max_component = a_ptr(i+1) - a_ptr(i)
       end do

     end if


     ! Check whether max number of levels has been reached or if matrix
     ! size is below
     ! nd_switch
     if (level>=options%amd_switch2 .or. a_n<=max(2,options%amd_switch1) &
       .or. use_amd) &
       go to 10
     lwork = 12*a_n + sumweight + a_ne
     call nd_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,level, &
       a_n1,a_n2,a_ne1,a_ne2,iperm,work(1:lwork),options,info, &
       use_multilevel,grid)


     if (a_n1==a_n) then
       go to 10
     end if


     if (a_n1/=0 .and. a_n2/=0 .and. a_n1+a_n2==a_n) then
       ! matrix is reducible
       if (printi .or. printd) then
         write (unit_diagnostics,'(a)') ' '
         write (unit_diagnostics,'(a)') 'Matrix is reducible'
       end if
       compwork = 0
       ! work array needs to be total length 5*a_n+a_ne
       call nd_find_indep_comps(a_n,a_ne,a_ptr,a_row,a_weight,iperm, &
         num_components,work_comp_n(1:a_n),work_comp_nz(1:a_n), &
         work(compwork+1:compwork+3*a_n+a_ne),options,info)
       if (num_components==1) then
         k = options%amd_switch2 ! Should never be reached - indep. comps.
         ! only found if it has been detected that they exist
       else
         k = level
         if (k==0) info%num_components = num_components
       end if

       ! Apply the ND to each component - do not test for indep comps
       offset_ptr = a_n + 1
       offset_row = a_ne + 1
       i = num_components
       ! work from last component so that space at end of work_comp_n and
       ! work_comp_nz can be reused without worrying about overwriting
       ! important data
       j = a_n
       do while (i>=1)
         if (level==0) use_multilevel = .true.
         l = work_comp_n(i)
         if (level==0 .and. l .le. options%amd_call) then
            use_amdi = .true.
         else
            use_amdi = .false.
         end if
         m = work_comp_nz(i)
         s = sum(a_weight(offset_ptr-l:offset_ptr-1))
         if (m>0) then
           ! Matrix not diagonal
           if (present(seps)) then
             call nd_nested_internal(l,m,a_ptr(offset_ptr-l:offset_ptr-1 &
               ),a_row(offset_row-m:offset_row-1), &
               a_weight(offset_ptr-l:offset_ptr-1),s, &
               iperm(offset_ptr-l:offset_ptr-1),work(compwork+1:compwork+12 &
               *l+s+m),work_comp_n(j-l+1:j),work_comp_nz(j-l+1:j),k, &
               options,info,use_amdi,use_multilevel,grid,seps(offset_ptr-l: &
               offset_ptr-1))
           else
             call nd_nested_internal(l,m,a_ptr(offset_ptr-l:offset_ptr-1 &
               ),a_row(offset_row-m:offset_row-1), &
               a_weight(offset_ptr-l:offset_ptr-1),s, &
               iperm(offset_ptr-l:offset_ptr-1),work(compwork+1:compwork+12 &
               *l+s+m),work_comp_n(j-l+1:j),work_comp_nz(j-l+1:j),k, &
               options,info,use_amdi,use_multilevel,grid)

           end if
         end if
         offset_ptr = offset_ptr - l
         offset_row = offset_row - m
         j = j - l
         i = i - 1
       end do
       return
     else
       if (level==0 .and. a_n>info%n_max_component) then
         info%n_max_component = a_n
         info%nz_max_component = a_ne
         info%maxdeg_max_component = maxdeg_max_component
       end if
     end if

     if (present(seps)) then
       seps(a_n1+a_n2+1:a_n) = level
     end if
     if (a_n1>max(2,options%amd_switch1)) then
       sumweight_sub = sum(a_weight(1:a_n1))
       if (present(seps)) then
         call nd_nested_internal(a_n1,a_ne1,a_ptr(1:a_n1), &
           a_row(1:a_ne1),a_weight(1:a_n1),sumweight_sub,iperm(1:a_n1), &
           work(1:12*a_n1+sumweight_sub+a_ne1),work_comp_n(1:a_n1), &
           work_comp_nz(1:a_n1),level+1,options,info,use_amdi,&
           use_multilevel,grid, &
           seps(1:a_n1))
       else
         call nd_nested_internal(a_n1,a_ne1,a_ptr(1:a_n1), &
           a_row(1:a_ne1),a_weight(1:a_n1),sumweight_sub,iperm(1:a_n1), &
           work(1:12*a_n1+sumweight_sub+a_ne1),work_comp_n(1:a_n1), &
           work_comp_nz(1:a_n1),level+1,options,info,use_amdi,&
           use_multilevel,grid)
       end if

     end if

     if (a_n2>max(2,options%amd_switch1)) then
       if (a_n1>max(2,options%amd_switch1)) then
         sumweight_sub = sum(a_weight(a_n1+1:a_n1+a_n2))
         if (present(seps)) then
           call nd_nested_internal(a_n2,a_ne2,a_ptr(a_n1+1:a_n1+a_n2), &
             a_row(a_ne1+1:a_ne1+a_ne2),a_weight(a_n1+1:a_n1+a_n2), &
             sumweight_sub,iperm(a_n1+1:a_n1+a_n2), &
             work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2), &
             work_comp_nz(1:a_n2),level+1,options,info,use_amdi,&
             use_multilevel,grid, &
             seps(a_n1+1:a_n1+a_n2))
         else
           call nd_nested_internal(a_n2,a_ne2,a_ptr(a_n1+1:a_n1+a_n2), &
             a_row(a_ne1+1:a_ne1+a_ne2),a_weight(a_n1+1:a_n1+a_n2), &
             sumweight_sub,iperm(a_n1+1:a_n1+a_n2), &
             work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2), &
             work_comp_nz(1:a_n2),level+1,options,info,use_amdi,&
             use_multilevel,grid)
         end if
       else
         sumweight_sub = sum(a_weight(a_n1+1:a_n1+a_n2))
         if (present(seps)) then
           call nd_nested_internal(a_n2,a_ne2,a_ptr(1:a_n2), &
             a_row(1:a_ne2),a_weight(a_n1+1:a_n1+a_n2),sumweight_sub, &
             iperm(a_n1+1:a_n1+a_n2),work(1:12*a_n2+sumweight_sub+a_ne2), &
             work_comp_n(1:a_n2),work_comp_nz(1:a_n2),level+1,options,info, &
             use_amdi,use_multilevel,grid,seps(a_n1+1:a_n1+a_n2))
         else
           call nd_nested_internal(a_n2,a_ne2,a_ptr(1:a_n2), &
             a_row(1:a_ne2),a_weight(a_n1+1:a_n1+a_n2),sumweight_sub, &
             iperm(a_n1+1:a_n1+a_n2),work(1:12*a_n2+sumweight_sub+a_ne2), &
             work_comp_n(1:a_n2),work_comp_nz(1:a_n2),level+1,options,info, &
             use_amdi,use_multilevel,grid)

         end if

       end if

     end if
     go to 20
     return

     ! No partition found or max number of levels have been reached
10      continue
     ! Apply AMD to matrix (note: a_n always greater than 1)
     lirn = a_ne + a_n
     amd_order_irn = 0
     amd_order_ip = amd_order_irn + lirn
     amd_order_sep = amd_order_ip + a_n
     amd_order_perm = amd_order_sep + a_n
     amd_order_work = amd_order_perm + a_n

     work(amd_order_irn+1:amd_order_irn+a_ne) = a_row(1:a_ne)
     work(amd_order_irn+a_ne+1:amd_order_irn+lirn) = 0
     work(amd_order_ip+1:amd_order_ip+a_n) = a_ptr(1:a_n)
     work(amd_order_sep+1:amd_order_sep+a_n) = ND_SEP_FLAG - 1
     call amd_order(a_n,a_ne,lirn,work(amd_order_ip+1:amd_order_ip+a_n), &
       work(amd_order_irn+1:amd_order_irn+lirn),work(amd_order_sep+1:amd_order_sep+a_n), &
       work(amd_order_perm+1:amd_order_perm+a_n),work(amd_order_work+1:amd_order_work+7*a_n))

     ! Extract perm from amd_order to apply to iperm
     do i = 1, a_n
       j = work(amd_order_perm+i)
       work(amd_order_work+i) = iperm(j)
     end do
     iperm(1:a_n) = work(amd_order_work+1:amd_order_work+a_n)

20      info%flag = 0
     if (printi .or. printd) then
       call nd_print_message(info%flag,unit_diagnostics, &
         'nd_nested_internal')
     end if
     return

   end subroutine nd_nested_internal

   ! ---------------------------------------------------
   ! nd_find_indep_comps
   ! ---------------------------------------------------
   ! Finds and forms independent components in a matrix
   subroutine nd_find_indep_comps(a_n,a_ne,a_ptr,a_row,a_weight,iperm, &
       comp_num,compsizes,compnzs,work,options,info)
     integer, intent(in) :: a_n ! size of matrix
     integer, intent(in) :: a_ne ! no. nonzeros in matrix
     integer, intent(inout) :: a_ptr(a_n) ! On entry, column ptrs for
     ! input
     ! matrix.
     ! On exit, a_ptr(1:compsizes(1)) contains column ptrs for compontent
     ! 1;
     ! a_ptr(compsizes(1)+1:compsizes(1)+compsizes(2)) contains column ptrs
     ! for compontent 2; etc.
     integer, intent(inout) :: a_row(a_ne) ! On entry, row indices for
     ! input
     ! matrix.
     ! On exit, a_row(1:compnzs(1)) contains row indices for compontent 1;
     ! a_ptr(compnzs(1)+1:compnzs(1)+compnzs(2)) contains row indices
     ! for compontent 2; etc.
     integer, intent(inout) :: a_weight(a_n) ! On entry, a_weight(i)
     ! contains
     ! weight of column i for input matrix.
     ! On exit, a_weight(1:compsizes(1)) contains column weights for
     ! compontent 1;
     ! a_weight(compsizes(1)+1:compsizes(1)+compsizes(2)) contains column
     ! weights or compontent 2; etc.
     integer, intent(inout) :: iperm(a_n) ! On input, iperm(i) contains
     ! the
     ! row in the original matrix (when nd_nested was called) that
     ! row i in this sub problem maps to. On output, this is updated to
     ! reflect the computed permutation.
     integer, intent(out) :: comp_num ! number independent components
     ! found
     integer, intent(out) :: compsizes(a_n) ! compsizes(i) will contain
     ! the
     ! size of compontent i
     integer, intent(out) :: compnzs(a_n) ! compnzs(i) will contain the
     ! number of nonzeros in compontent i
     integer, intent(out) :: work(3*a_n+a_ne) ! used as work arrays during
     ! computation
     type (nd_options), intent(in) :: options
     type (nd_inform), intent(inout) :: info

     ! ---------------------------------------------
     ! Local variables
     integer :: unit_diagnostics ! unit on which to print diagnostics
     integer :: mask, ptr_temp, row_temp, front ! pointers into work array
     integer :: root
     integer :: front_sta, front_sto ! start and end of front list
     integer :: num_assigned ! number of cols that have been assigned to
     ! a
     ! component
     integer :: i, j, l, u, v, p
     integer :: ptr_temp_pos, row_temp_pos, front_pos
     logical :: printi, printd

     ! ---------------------------------------------
     ! Printing levels
     unit_diagnostics = options%unit_diagnostics
     printi = (options%print_level==1 .and. unit_diagnostics>=0)
     printd = (options%print_level>=2 .and. unit_diagnostics>=0)

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a)') 'Start independent component search'
     end if

     mask = 0
     front = mask + a_n
     ptr_temp = front + a_n
     row_temp = ptr_temp + a_n ! total length of work array 3*a_n+a_ne

     ! Initialise arrays
     work(1:3*a_n) = 0
     compsizes(1:a_n) = 0
     compnzs(1:a_n) = 0

     ! Check from all possible roots
     num_assigned = 0
     comp_num = 0
     do root = 1, a_n
       if (work(mask+root)==0) then
         comp_num = comp_num + 1
         front_sta = num_assigned + 1
         front_sto = front_sta
         work(front+front_sta) = root
         compsizes(comp_num) = 1
         compnzs(comp_num) = 0
         work(mask+root) = compsizes(comp_num)
         num_assigned = num_assigned + 1

         do while (front_sto-front_sta>=0)
           do i = front_sta, front_sto
             ! pick vertex from front
             v = work(front+i)
             ! update compnzs
             if (v<a_n) then
               l = a_ptr(v+1)
             else
               l = a_ne + 1
             end if
             compnzs(comp_num) = compnzs(comp_num) + l - a_ptr(v)
             do j = a_ptr(v), l - 1
               ! pick a neighbour
               u = a_row(j)
               if (work(mask+u)/=0) cycle
               ! found unmasked vertex
               compsizes(comp_num) = compsizes(comp_num) + 1
               num_assigned = num_assigned + 1
               ! mask this vertex
               work(mask+u) = compsizes(comp_num)
               ! add vertex to component
               work(front+num_assigned) = u
             end do
           end do
           front_sta = front_sto + 1
           front_sto = num_assigned
         end do
       end if
     end do

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(i5,a)') comp_num, ' components found'
     end if

     if (comp_num>1) then
       ! Reorder matrix into block diagonal form
       ! Indexing will be local to each block
       ptr_temp_pos = 1
       row_temp_pos = 1
       front_pos = 1
       do l = 1, comp_num ! for each component
         work(ptr_temp+ptr_temp_pos) = 1
         ptr_temp_pos = ptr_temp_pos + 1
         do u = 1, compsizes(l)
           v = work(front+front_pos) ! for each column in the component
           front_pos = front_pos + 1
           if (v==a_n) then
             p = a_ne
           else
             p = a_ptr(v+1) - 1
           end if
           do i = a_ptr(v), p ! for each nonzero in the column
             work(row_temp+row_temp_pos) = work(mask+a_row(i))
             row_temp_pos = row_temp_pos + 1
           end do
           if (u<compsizes(l)) then
             work(ptr_temp+ptr_temp_pos) = work(ptr_temp+ptr_temp_pos-1) + &
               p + 1 - a_ptr(v)
             ptr_temp_pos = ptr_temp_pos + 1
           end if
         end do
       end do

       a_ptr(1:a_n) = work(ptr_temp+1:ptr_temp+a_n)
       a_row(1:a_ne) = work(row_temp+1:row_temp+a_ne)

       ! Reorder iperm and a_weight
       do front_pos = 1, a_n
         j = work(front+front_pos)
         work(front+front_pos) = iperm(j)
         work(mask+front_pos) = a_weight(j)
       end do
       iperm(1:a_n) = work(front+1:front+a_n)
       a_weight(1:a_n) = work(mask+1:mask+a_n)
     end if

     info%flag = 0
     if (printi .or. printd) then
       call nd_print_message(info%flag,unit_diagnostics, &
         'nd_find_indep_comps')
     end if
     return

   end subroutine nd_find_indep_comps

   ! ---------------------------------------------------
   ! nd_partition matrix
   ! ---------------------------------------------------
   ! Partition the matrix and if one (or more) of the generated submatrices
   ! is
   ! small enough, apply halo amd

   subroutine nd_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
       level,a_n1,a_n2,a_ne1,a_ne2,iperm,work,options,info,use_multilevel, &
       grid)

     integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
     integer, intent(in) :: a_ne ! no. nonzeros of subproblem
     integer, intent(inout) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start. This is then
     ! used to hold positions for submatrices after partitioning
     integer, intent(inout) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.This is then used to hold row indices for
     ! submatrices after partitioning
     integer, intent(inout) :: a_weight(a_n) ! On input a_weight(i)
     ! contains
     ! the weight of column i. This is then used to hold the weights for
     ! the submatrices after partitioning.
     integer, intent(in) :: sumweight ! Sum entries in a_weight.
     ! Unchanged.
     integer, intent(in) :: level ! Current nested dissection level
     integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
     integer, intent(out) :: a_ne1, a_ne2 ! no. nonzeros in two
     ! submatrices
     integer, intent(inout) :: iperm(a_n) ! On input, iperm(i) contains
     ! the
     ! row in the original matrix (when nd_nested was called) that
     ! row i in this sub problem maps to. On output, this is updated to
     ! reflect the computed permutation.
     logical, intent(inout) :: use_multilevel
     integer, intent(out) :: work(12*a_n+sumweight+a_ne)
     type (nd_options), intent(in) :: options
     type (nd_inform), intent(inout) :: info
     type (nd_multigrid), intent(inout) :: grid

     ! ---------------------------------------------
     ! Local variables
     integer :: unit_diagnostics ! unit on which to print diagnostics
     logical :: printi, printd
     integer :: partition_ptr ! pointer into work array
     integer :: work_ptr ! pointer into work array
     integer :: part_ptr ! pointer into work array
     integer :: a_ptr_sub_ptr ! pointer into work array
     integer :: a_row_sub_ptr ! pointer into work array
     integer :: partition_method
     integer :: i, j, k
     integer :: a_weight_1, a_weight_2, a_weight_sep
     integer :: ref_method, ref_options
     integer :: a_n1_new, a_n2_new, a_weight_1_new, a_weight_2_new, &
       a_weight_sep_new
     real(wp) :: ratio, tau_best, tau, band, depth
     logical :: imbal


     ! ---------------------------------------------
     ! Printing levels
     unit_diagnostics = options%unit_diagnostics
     printi = (options%print_level==1 .and. unit_diagnostics>=0)
     printd = (options%print_level>=2 .and. unit_diagnostics>=0)

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a)') 'Start finding a partition'
     end if

     ! If matrix is full, then don't partition
     if (a_ne==a_n*(a_n-1)) then
       a_n1 = a_n
       a_n2 = 0
       a_ne1 = a_ne
       a_ne2 = 0
       go to 10
     end if

     ! Find the partition
     if (options%coarse_partition_method<=1) then
       partition_method = 1
     else
       partition_method = 2
     end if


     partition_ptr = 0 ! length a_n
     work_ptr = partition_ptr + a_n ! max length 9*a_n+sumweight+a_ne/2
     ! do p = 1,a_n
     ! do q = p+1,a_n
     ! if (options%partition_method .ge. 2) use_multilevel = .true.
     if (partition_method==1) then
       ! Ashcraft method
       call nd_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,level, &
         a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
         work(partition_ptr+1:partition_ptr+a_n), &
         work(work_ptr+1:work_ptr+9*a_n+sumweight),options,info%flag,band, &
         depth,use_multilevel,grid)
     else
       call nd_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
         level,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
         work(partition_ptr+1:partition_ptr+a_n), &
         work(work_ptr+1:work_ptr+9*a_n+sumweight),options,info%flag, &
         band,depth,use_multilevel,grid)

     end if

     if (a_n1+a_n2==a_n) then
       return
     end if

     if (level==0) then
       if (a_n>info%n_max_component .or. (a_n==info%n_max_component .and. &
           band>info%band)) then
         info%band = band
         info%depth = depth
       end if

     end if

     if (a_n1/=0 .and. a_n2/=0 .and. a_n>=3) then
       if ( .not. use_multilevel) then

         if (options%refinement>6) then
           ref_options = 3
         else
           if (options%refinement<1) then
             ref_options = 1
           else
             ref_options = options%refinement
           end if
         end if

         select case (ref_options)
         case (1)
           ref_method = 1

         case (2)
           ref_method = 2

         case (3)
           if (real(max(a_weight_1,a_weight_2))/real(min(a_weight_1, &
               a_weight_2)+a_weight_sep)>=max(real(1.0, &
               wp),options%balance)) then
             ref_method = 2
           else
             ref_method = 1
           end if

         case (4)
           ref_method = 0

         case (5)
           ref_method = 2

         case (6)
           if (real(max(a_weight_1,a_weight_2))/real(min(a_weight_1, &
               a_weight_2)+a_weight_sep)>=max(real(1.0, &
               wp),options%balance)) then
             ref_method = 2
           else
             ref_method = 0
           end if

         end select
         if (printd) then
           write (unit_diagnostics,'(a)') 'Partition before refinement'
           write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', a_n1, &
             ',  a_n2=', a_n2, ',  a_n_sep=', a_n - a_n1 - a_n2
         end if

         select case (ref_method)

         case (0)
           call nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1, &
             a_n2,a_weight_1,a_weight_2,a_weight_sep, &
             work(partition_ptr+1:partition_ptr+a_n), &
             work(work_ptr+1:work_ptr+8),options)

         case (1)
           if (min(a_weight_1,a_weight_2)+a_weight_sep< &
               max(a_weight_1,a_weight_2)) then
             call nd_refine_block_trim(a_n,a_ne,a_ptr,a_row, &
               a_weight,sumweight,a_n1,a_n2,a_weight_1,a_weight_2, &
               a_weight_sep,work(partition_ptr+1:partition_ptr+a_n), &
               work(work_ptr+1:work_ptr+5*a_n),options)
           else
             call nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
               sumweight,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
               work(partition_ptr+1:partition_ptr+a_n), &
               work(work_ptr+1:work_ptr+3*a_n),options)
           end if


         case (2)
           call nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
             a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
             work(partition_ptr+1:partition_ptr+a_n), &
             work(work_ptr+1:work_ptr+3*a_n),options)


         end select


         if (printd) then
           write (unit_diagnostics,'(a)') 'Partition after refinement'
           write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', a_n1, &
             ',  a_n2=', a_n2, ',  a_n_sep=', a_n - a_n1 - a_n2
         end if

         if (options%max_improve_cycles>0) then
           ratio = max(real(1.0,wp),options%balance)
           if (ratio>real(sumweight-2)) then
             imbal = .false.
           else
             imbal = .true.
           end if
           call cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight, &
             ratio,imbal,options%cost_function,tau_best)
           a_n1_new = a_n1
           a_n2_new = a_n2
           a_weight_1_new = a_weight_1
           a_weight_2_new = a_weight_2
           a_weight_sep_new = a_weight_sep
         end if

         part_ptr = work_ptr + 5*a_n
         work(part_ptr+1:part_ptr+a_n) = work(partition_ptr+1:partition_ptr &
           +a_n)

         k = options%max_improve_cycles
         do i = 1, k

           call expand_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1_new, &
             a_n2_new,a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
             work(part_ptr+1:part_ptr+a_n),work(work_ptr+1:work_ptr+5*a_n))


           if (printd) then
             write (unit_diagnostics,'(a)') &
               'Partition sizes after expansion'
             write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', &
               a_n1_new, ',  a_n2=', a_n2_new, ',  a_n_sep=', &
               a_n - a_n1_new - a_n2_new
           end if

           ! call
           ! check_partition1(a_n,a_ne,a_ptr,a_row,a_n1_new,a_n2_new,work(p
           ! art_ptr+1:part_ptr+a_n))

           select case (ref_options)

           case (3)
             if (real(max(a_weight_1_new,a_weight_2_new))/real(min( &
                 a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>max(real( &
                 1.0,wp),options%balance)) then
               ref_method = 2
             else
               ref_method = 1
             end if

           case (6)
             if (real(max(a_weight_1_new,a_weight_2_new))/real(min( &
                 a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>max(real( &
                 1.0,wp),options%balance)) then
               ref_method = 2
             else
               ref_method = 0
             end if
           end select


           select case (ref_method)

           case (0)
             call nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight, &
               a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
               a_weight_sep_new,work(part_ptr+1:part_ptr+a_n), &
               work(work_ptr+1:work_ptr+8),options)

           case (1)
             if (min(a_weight_1,a_weight_2)+a_weight_sep< &
                 max(a_weight_1,a_weight_2)) then
               call nd_refine_block_trim(a_n,a_ne,a_ptr,a_row, &
                 a_weight,sumweight,a_n1_new,a_n2_new,a_weight_1_new, &
                 a_weight_2_new,a_weight_sep_new, &
                 work(part_ptr+1:part_ptr+a_n),work(work_ptr+1:work_ptr+5* &
                 a_n),options)
             else
               call nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
                 sumweight,a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
                 a_weight_sep_new,work(part_ptr+1:part_ptr+a_n), &
                 work(work_ptr+1:work_ptr+3*a_n),options)
             end if


           case (2)
             call nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
               a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
               a_weight_sep_new,work(part_ptr+1:part_ptr+a_n), &
               work(work_ptr+1:work_ptr+3*a_n),options)

           end select


           if (printd) then
             write (unit_diagnostics,'(a)') &
               'Partition sizes after refinement'
             write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', &
               a_n1_new, ',  a_n2=', a_n2_new, ',  a_n_sep=', &
               a_n - a_n1_new - a_n2_new
           end if


           call cost_function(a_weight_1_new,a_weight_2_new, &
             a_weight_sep_new,sumweight,ratio,imbal,options%cost_function,tau)
           if (tau<tau_best) then
             tau_best = tau
             work(partition_ptr+1:partition_ptr+a_n) &
               = work(part_ptr+1:part_ptr+a_n)
             a_n1 = a_n1_new
             a_n2 = a_n2_new
             a_weight_1 = a_weight_1_new
             a_weight_2 = a_weight_2_new
             a_weight_sep = a_weight_sep_new
           else
             exit
           end if
         end do

         call nd_refine_fm(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
           a_n2,a_weight_1,a_weight_2,a_weight_sep, &
           work(partition_ptr+1:partition_ptr+a_n), &
           work(work_ptr+1:work_ptr+8*a_n+sumweight),options)


       end if
     else
       go to 10
     end if

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') 'Partition found:'
       write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', a_n1, &
         ', a_n2=', a_n2, ', a_nsep=', a_n - a_n1 - a_n2
     end if

     if ((a_n1<=max(2,options%amd_switch1) .and. a_n2<=max(2, &
         options%amd_switch1))) then
       ! apply halo amd to submatrics
       call amd_order_both(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
         work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight, &
         work(work_ptr+1:work_ptr+12*a_n+a_ne))
       a_ne1 = 0
       a_ne2 = 0

     else if (a_n1<=max(2,options%amd_switch1)) then
       ! apply halo amd to [A1, B1'; B1, I] using two levels
       ! return A2 and apply ND to it
       call amd_order_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
         work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight,a_ne2, &
         work(work_ptr+1:work_ptr+12*a_n+a_ne))
       a_ne1 = 0


     else if (a_n2<=max(2,options%amd_switch1)) then
       ! apply halo amd to [A2, B2'; B2, I] using two levels
       ! return A1 and apply ND to it
       call amd_order_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
         work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight,a_ne1, &
         work(work_ptr+1:work_ptr+12*a_n+a_ne))
       a_ne2 = 0

     else
       ! return A1 and A2 and apply ND to them
       a_ptr_sub_ptr = work_ptr + a_n ! length a_n
       a_row_sub_ptr = a_ptr_sub_ptr + a_n ! length a_ne

       call extract_both_matrices(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
         work(partition_ptr+1:partition_ptr+a_n1+a_n2),a_ne1,a_ne2, &
         work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2),a_ne, &
         work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne), &
         work(work_ptr+1:work_ptr+a_n))

       ! Copy extracted matrices
       a_ptr(1:a_n1+a_n2) = work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2)
       a_row(1:a_ne1+a_ne2) = work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne1+ &
         a_ne2)

       ! Update iperm and a_weight
       do i = 1, a_n
         j = work(partition_ptr+i)
         work(a_ptr_sub_ptr+i) = iperm(j)
         work(work_ptr+i) = a_weight(j)
       end do
       do i = 1, a_n
         iperm(i) = work(a_ptr_sub_ptr+i)
         a_weight(i) = work(work_ptr+i)
       end do

     end if
     go to 20

10      if (printi .or. printd) then
       write (unit_diagnostics,'(a)') 'No partition found'
     end if

20      info%flag = 0
     if (printi .or. printd) then
       call nd_print_message(info%flag,unit_diagnostics, &
         'nd_partition')
     end if
     return

   end subroutine nd_partition

   ! ---------------------------------------------------
   ! nd_ashcraft
   ! ---------------------------------------------------
   ! Partition the matrix using the Ashcraft method
   recursive subroutine nd_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight, &
       sumweight,level,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
       partition,work,options,info,band,depth,use_multilevel,grid)

     integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
     integer, intent(in) :: a_ne ! no. nonzeros of subproblem
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i.
     integer, intent(in) :: sumweight ! sum of entries in a_weight
     integer, intent(in) :: level ! current level of nested dissection
     integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
     integer, intent(out) :: a_weight_1, a_weight_2, a_weight_sep ! Weight
     ! ed
     ! size of partitions and separator
     integer, intent(out) :: partition(a_n) ! First a_n1 entries will
     ! contain
     ! list of (local) indices in partition 1; next a_n2 entries will
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end
     integer, intent(out) :: work(9*a_n+sumweight) ! used as work array
     type (nd_options), intent(in) :: options
     integer, intent(inout) :: info
     real(wp), intent(out) :: band ! If level = 0, then on
     ! output band = 100*L/a_n, where L is the size of the
     ! largest levelset
     real(wp), intent(out) :: depth ! If level = 0, then
     ! on
     ! output depth = num_levels_nend
     logical, intent(inout) :: use_multilevel ! are we allowed to use a
     ! multilevel
     ! partitioning strategy
     type (nd_multigrid), intent(inout) :: grid

     ! ---------------------------------------------
     ! Local variables
     integer :: nstrt, nend
     integer :: i, j, dptr, p1sz, p2sz, sepsz, lwork, k
     integer :: unit_diagnostics ! unit on which to print diagnostics
     integer :: mask_p, level_p, level_ptr_p, level2_p, level2_ptr_p, &
       work_p
     integer :: num_levels_nend ! no. levels in structure rooted at nend
     integer :: num_levels_nstrt ! no. levels in structure rooted at nstrt
     integer :: num_entries ! no. entries in level set structure
     integer :: best_sep_start
     integer :: distance
     integer :: distance_ptr
     integer :: lwidth, mindeg, degree, max_search
     integer :: ww
     integer :: stop_coarsening2 ! max no. multigrid levels
     real(wp) :: bestval
     real(wp) :: val
     real(wp) :: ratio
     logical :: printi, printd
     logical :: imbal, use_multilevel_copy

     ! ---------------------------------------------
     ! Printing levels
     unit_diagnostics = options%unit_diagnostics
     printi = (options%print_level==1 .and. unit_diagnostics>=0)
     printd = (options%print_level>=2 .and. unit_diagnostics>=0)

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a)') 'Use two-sided level set method'
     end if
     ratio = max(real(1.0,wp),options%balance)
     if (ratio>real(sumweight-2)) then
       imbal = .false.
     else
       imbal = .true.
     end if
     p2sz = 0
     sepsz = 0
     use_multilevel_copy = use_multilevel


     band = -1
     depth = -1

     if (options%partition_method==1 .and. use_multilevel) go to 10

     if (options%partition_method>1 .and. level>0 .and. use_multilevel) &
       go to 10

     ! Find pseudoperipheral nodes nstart and nend, and the level structure
     ! rooted at nend
     mask_p = 0 ! size a_n
     level_ptr_p = mask_p + a_n ! size a_n
     level_p = level_ptr_p + a_n ! size a_n
     level2_ptr_p = level_p + a_n ! size a_n
     level2_p = level2_ptr_p + a_n ! size a_n
     work_p = level2_p + a_n ! size 2*a_n

     nend = -1

     ! Choose nstrt
     ! node with minimum degree
     mindeg = sumweight + 1
     do i = 1, a_n
       if (i<a_n) then
         k = a_ptr(i+1) - 1
       else
         k = a_ne
       end if
       degree = 0
       do j = a_ptr(i), k
         degree = degree + a_weight(a_row(j))
       end do
       if (degree<mindeg) then
         mindeg = degree
         nstrt = i
       end if
     end do

     max_search = 5

     call nd_find_pseudo(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
       work(level_ptr_p+1:level_ptr_p+a_n),work(level_p+1:level_p+a_n), &
       nstrt,nend,max_search,work(work_p+1:work_p+2*a_n),num_levels_nend, &
       num_entries,lwidth)

     if (num_entries<a_n) then
       ! matrix is separable
       a_n1 = num_entries
       a_n2 = a_n - a_n1
       a_weight_sep = 0
       a_weight_1 = 0
       work(work_p+1:work_p+a_n) = 0
       do i = 1, num_entries
         j = work(level_p+i)
         partition(i) = j
         work(work_p+j) = 1
         a_weight_1 = a_weight_1 + a_weight(j)
       end do
       a_weight_2 = sumweight - a_weight_1
       j = num_entries + 1
       do i = 1, a_n
         if (work(work_p+i)==0) then
           partition(j) = i
           j = j + 1
         end if

       end do
       if (level==0) then
         band = -real(lwidth,wp)
       end if

       return
     end if

     ! ********************************************************************
     ! ***
     if (level==0) then
       band = 100.0*real(lwidth,wp)/real(a_n,wp)
       ! band = max(band,real(lwidth,wp))
       ! write(*,*) sqrt(real(lwidth))/real(num_levels_nend), lwidth, &
       ! real(sumweight)/real(num_levels_nend), &
       ! real(num_levels_nend)*(real(lwidth)/real(sumweight))
       depth = 100.0*real(num_levels_nend,wp)/ &
         real(a_n,wp)
     end if
     if (options%stop_coarsening2<=0 .or. options%partition_method<1) then
       use_multilevel = .false.
     end if
     if (options%partition_method>=2 .and. use_multilevel) then
       ! if (real(lwidth,wp) .le. 2.0* &
       ! real(sumweight,wp)/real(num_levels_nend,wp) )
       ! then
       if (100.0*real(lwidth,wp)/real(sumweight,wp)<=3.0 .or. &
          a_n.lt. options%ml_call) &
           then
         use_multilevel = .false.
       else
         use_multilevel = .true.
       end if
     end if

     if (use_multilevel) go to 10


     ! Find level structure rooted at nstrt
     work(mask_p+1:mask_p+a_n) = 1
     call nd_level_struct(nstrt,a_n,a_ne,a_ptr,a_row, &
       work(mask_p+1:mask_p+a_n),work(level2_ptr_p+1:level2_ptr_p+a_n), &
       work(level2_p+1:level2_p+a_n),num_levels_nstrt,lwidth,num_entries)

     ! Calculate difference in distances from nstart and nend, and set up
     ! lists D_i of nodes with same distance
     distance_ptr = work_p
     distance = mask_p
     call nd_distance(a_n,num_levels_nend,work(level_ptr_p+1:level_ptr_p &
       +a_n),work(level_p+1:level_p+a_n),num_levels_nstrt, &
       work(level2_ptr_p+1:level2_ptr_p+a_n),work(level2_p+1:level2_p+a_n), &
       work(distance_ptr+1:distance_ptr+2*a_n-1), &
       work(distance+1:distance+a_n))

     ! Do not need the information in work(level_ptr_p+1:level2_ptr_p)
     ! Calculate total weight in each distance level
     dptr = level_ptr_p + a_n
     work(dptr+1-a_n:dptr+a_n-1) = 0
     do i = 1 - num_levels_nstrt, num_levels_nend - 2
       do j = work(distance_ptr+a_n+i), work(distance_ptr+a_n+i+1) - 1
         work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
       end do
     end do
     i = num_levels_nend - 1
     do j = work(distance_ptr+a_n+i), a_n
       work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
     end do

     ! Find first possible separator
     ww = 2
     p1sz = 0
     do i = 1 - num_levels_nstrt, num_levels_nend - ww - 2
       p1sz = p1sz + work(dptr+i)
       sepsz = work(dptr+i+1)
       do j = 1, ww - 1
         sepsz = sepsz + work(dptr+i+1+j)
       end do
       p2sz = sumweight - p1sz - sepsz
       if (p1sz>0 .and. sepsz>0) then
         exit
       end if
     end do

     if (i+ww>=num_levels_nend-1 .or. p2sz==0) then
       ! Not possible to find separator
       ! This can only be invoked for a fully connected graph. The
       ! partition
       ! subroutine checks for this case so it should never be called
       a_n1 = a_n
       a_n2 = 0
       partition(1:a_n) = (/ (i,i=1,a_n) /)
       return
     end if

     ! Entries in level i will form first possible partition
     best_sep_start = i + 1
     a_n1 = work(distance_ptr+a_n+i+1) - 1
     a_n2 = a_n - work(distance_ptr+a_n+i+1+ww) + 1
     a_weight_1 = p1sz
     a_weight_2 = p2sz
     a_weight_sep = sepsz
     call cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
         options%cost_function,bestval)

     ! Search for best separator using tau

     do j = i + 1, num_levels_nend - ww - 2
       p1sz = p1sz + work(dptr+j)
       sepsz = work(dptr+j+1)
       do k = 1, ww - 1
         sepsz = sepsz + work(dptr+j+1+k)
       end do
       p2sz = sumweight - sepsz - p1sz
       if (p2sz==0) exit
       call cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
          options%cost_function,val)
       if (val<bestval) then
         bestval = val
         best_sep_start = j + 1
         a_n1 = work(distance_ptr+a_n+j+1) - 1
         a_n2 = a_n - work(distance_ptr+a_n+j+1+ww) + 1
         a_weight_1 = p1sz
         a_weight_2 = p2sz
         a_weight_sep = sepsz
       end if
     end do

     if (imbal .and. use_multilevel_copy .and. options%partition_method>=2) &
         then
       if (real(max(a_weight_1,a_weight_2))/real(min(a_weight_1, &
           a_weight_2))>ratio) then
         use_multilevel = .true.
         go to 10

       end if
     end if

     ! Rearrange partition
     ! Entries in partition 1
     j = 1
     do i = 1, work(distance_ptr+a_n+best_sep_start) - 1
       partition(j) = work(distance+i)
       j = j + 1
     end do

     ! Entries in partition 2
     do i = work(distance_ptr+a_n+best_sep_start+ww), a_n
       partition(j) = work(distance+i)
       j = j + 1
     end do

     ! Entries in separator
     do i = work(distance_ptr+a_n+best_sep_start), &
         work(distance_ptr+a_n+best_sep_start+ww) - 1
       partition(j) = work(distance+i)
       j = j + 1
     end do
     go to 20

10      continue

     if (use_multilevel) then
       stop_coarsening2 = options%stop_coarsening2
       lwork = 9*a_n + sumweight
       call multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
         partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,options, &
         info,lwork,work(1:lwork),stop_coarsening2,grid)
       return
     end if


20      info = 0
     if (printi .or. printd) then
       call nd_print_message(info,unit_diagnostics,'nd_ashcraft')
     end if
     return

   end subroutine nd_ashcraft


   ! ---------------------------------------------------
   ! nd_level_set
   ! ---------------------------------------------------
   ! Partition the matrix using the level set method
   recursive subroutine nd_level_set(a_n,a_ne,a_ptr,a_row,a_weight, &
       sumweight,level,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
       partition,work,options,info,band,depth,use_multilevel,grid)

     integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
     integer, intent(in) :: a_ne ! no. nonzeros of subproblem
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i.
     integer, intent(in) :: sumweight ! sum of entries in a_weight
     integer, intent(in) :: level ! current nested dissection level
     integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
     integer, intent(out) :: a_weight_1, a_weight_2, a_weight_sep ! Weight
     ! ed
     ! size of partitions and separator
     integer, intent(out) :: partition(a_n) ! First a_n1 entries will
     ! contain
     ! list of (local) indices in partition 1; next a_n2 entries will
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end
     integer, intent(out) :: work(9*a_n+sumweight)
     type (nd_options), intent(in) :: options
     integer, intent(inout) :: info
     real(wp), intent(out) :: band ! If level = 0, then on
     ! output band = 100*L/a_n, where L is the size of the
     ! largest levelset
     real(wp), intent(out) :: depth ! If level = 0, then
     ! on
     ! output band = num_levels_nend
     logical, intent(inout) :: use_multilevel ! are we allowed to use a
     ! multilevel
     ! partitioning strategy
     type (nd_multigrid), intent(inout) :: grid

     ! ---------------------------------------------
     ! Local variables
     integer :: unit_diagnostics ! unit on which to print diagnostics
     integer :: nstrt, nend
     logical :: printi, printd
     integer :: level_p, level_ptr_p, work_p
     integer :: num_levels_nend ! no. levels in structure rooted at nend
     integer :: num_entries ! no. entries in level structure rooted at
     ! nend
     integer :: best_sep_start
     integer :: i, j, k, p1sz, p2sz, sepsz, lwidth
     integer :: stop_coarsening2, lwork
     integer :: mindeg, degree, max_search
     real(wp) :: bestval
     real(wp) :: val
     real(wp) :: ratio
     logical :: imbal

     ! ---------------------------------------------
     ! Printing levels
     unit_diagnostics = options%unit_diagnostics
     printi = (options%print_level==1 .and. unit_diagnostics>=0)
     printd = (options%print_level>=2 .and. unit_diagnostics>=0)

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a)') 'Use one-sided level set method'
     end if

     ratio = max(real(1.0,wp),options%balance)
     if (ratio>real(sumweight-2)) then
       imbal = .false.
     else
       imbal = .true.
     end if

     band = -1
     depth = -1

     if (options%partition_method==1 .and. use_multilevel) then
       use_multilevel = .true.
       go to 10
     end if

     if (options%partition_method>1 .and. level>0 .and. use_multilevel) &
       go to 10

     ! Find pseudoperipheral nodes nstart and nend, and the level structure
     ! rooted at nend
     level_ptr_p = 0 ! size a_n
     level_p = level_ptr_p + a_n ! size a_n
     work_p = level_p + a_n ! size 2*a_n
     mindeg = sumweight + 1
     do i = 1, a_n
       if (i<a_n) then
         k = a_ptr(i+1) - 1
       else
         k = a_ne
       end if
       degree = 0
       do j = a_ptr(i), k
         degree = degree + a_weight(a_row(j))
       end do
       if (degree<mindeg) then
         mindeg = degree
         nstrt = i
       end if
     end do
     max_search = 5

     call nd_find_pseudo(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
       work(level_ptr_p+1:level_ptr_p+a_n),work(level_p+1:level_p+a_n), &
       nstrt,nend,max_search,work(work_p+1:work_p+2*a_n),num_levels_nend, &
       num_entries,lwidth)

     if (num_entries<a_n) then
       ! matrix is separable
       a_n1 = num_entries
       a_n2 = a_n - a_n1
       work(work_p+1:work_p+a_n) = 0
       a_weight_1 = 0
       do i = 1, num_entries
         j = work(level_p+i)
         partition(i) = j
         work(work_p+j) = 1
         a_weight_1 = a_weight_1 + a_weight(j)
       end do
       j = num_entries + 1
       a_weight_2 = 0
       do i = 1, a_n
         if (work(work_p+i)==0) then
           partition(j) = i
           a_weight_2 = a_weight_2 + a_weight(j)
           j = j + 1
         end if
       end do
       a_weight_sep = 0
       if (level==0) then
         band = -real(lwidth,wp)
       end if
       return
     end if

     if (level==0) then
       band = 100.0*real(lwidth,wp)/real(sumweight,wp)
       depth = 100.0*real(num_levels_nend,wp)/ &
         real(sumweight,wp)
       ! band = max(band,real(lwidth,wp))
     end if

     if ((options%partition_method<=0) .or. (use_multilevel .and. options% &
         stop_coarsening2<=0)) then
       use_multilevel = .false.
     end if
     if (options%partition_method>=2 .and. use_multilevel) then
       if (100.0*real(lwidth,wp)/real(sumweight,wp)<=3.0 .or. &
          a_n.lt. options%ml_call) &
           then
         use_multilevel = .false.
       else
         use_multilevel = .true.
       end if
     end if

10      continue

     if (use_multilevel) then
       stop_coarsening2 = options%stop_coarsening2
       lwork = 9*a_n + sumweight

       call multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
         partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,options, &
         info,lwork,work(1:lwork),stop_coarsening2,grid)
       return
     end if

     if (num_levels_nend<=2) then
       ! Not possible to find separator
       ! This can only be invoked for a full connected graph. The partition
       ! subroutine checks for this case so it should never be called
       a_n1 = a_n
       a_n2 = 0
       partition(1:a_n) = (/ (i,i=1,a_n) /)
       return
     end if

     ! Calculate total weight in each level set
     work(work_p+1:work_p+num_levels_nend) = 0
     do i = 1, num_levels_nend - 1
       do j = work(level_ptr_p+i), work(level_ptr_p+i+1) - 1
         work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
       end do
     end do
     i = num_levels_nend
     do j = work(level_ptr_p+i), a_n
       work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
     end do


     ! First possible separator contains all of the nodes in level 2
     p1sz = work(work_p+1)
     sepsz = work(work_p+2)
     p2sz = sum(work(work_p+3:work_p+num_levels_nend))
     a_weight_1 = p1sz
     a_weight_2 = p2sz
     a_weight_sep = sepsz
     best_sep_start = 2
     a_n1 = work(level_ptr_p+2) - 1
     a_n2 = a_n - work(level_ptr_p+3) + 1
     call cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
      options%cost_function,bestval)

     ! Search for best separator using tau
     do j = 2, num_levels_nend - 4
       p1sz = p1sz + work(work_p+j)
       sepsz = work(work_p+j+1)
       p2sz = p2sz - work(work_p+j+1)
       call cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
          options%cost_function,val)
       if (val<bestval) then
         bestval = val
         best_sep_start = j + 1
         a_n1 = work(level_ptr_p+j+1) - 1
         a_n2 = a_n - work(level_ptr_p+j+2) + 1
         a_weight_1 = p1sz
         a_weight_2 = p2sz
         a_weight_sep = sepsz
       end if
     end do


     ! Rearrange partition
     ! Entries in partition 1
     j = 1
     do i = 1, work(level_ptr_p+best_sep_start) - 1
       partition(j) = work(level_p+i)
       j = j + 1
     end do

     ! Entries in partition 2
     do i = work(level_ptr_p+best_sep_start+1), a_n
       partition(j) = work(level_p+i)
       j = j + 1
     end do

     ! Entries in separator
     do i = work(level_ptr_p+best_sep_start), work(level_ptr_p+ &
         best_sep_start+1) - 1
       partition(j) = work(level_p+i)
       j = j + 1
     end do

     info = 0
     if (printi .or. printd) then
       call nd_print_message(info,unit_diagnostics,'nd_level_set')
     end if
     return

   end subroutine nd_level_set

   ! ---------------------------------------------------
   ! nd_find_pseudo
   ! ---------------------------------------------------
   ! Find pseudoperipheral pairs of nodes for irreducible graph
   subroutine nd_find_pseudo(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
       level_ptr,level,nstrt,nend,max_search,work,num_levels,num_entries, &
       lwidth)
     integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
     integer, intent(in) :: a_ne ! no. nonzeros of subproblem
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! weight of vertex i
     integer, intent(in) :: sumweight ! sum of entries in a_weight
     integer, intent(out) :: level_ptr(a_n) ! On output level_ptr(i)
     ! contains
     ! position in level that entries for level i start.
     integer, intent(out) :: level(a_n) ! On output level contains lists
     ! of
     ! rows according to the level set that they are in
     integer, intent(inout) :: nstrt ! Starting pseudoperipheral node
     integer, intent(out) :: nend ! End pseudoperipheral node
     integer, intent(in) :: max_search
     integer, intent(out) :: work(2*a_n)
     integer, intent(out) :: num_levels
     integer, intent(out) :: num_entries ! number of entries in level
     ! structure
     integer, intent(out) :: lwidth
     ! Based on MC60HD

     ! ---------------------------------------------
     ! Local variables
     integer :: i, j, k, l, ll
     integer :: mindeg, maxdep, main, lwidth1
     integer :: mask, list
     integer :: nstop ! ending pseudoperipheral node
     integer :: node ! Index of graph node
     integer :: lsize ! size of final levelset
     integer :: nlsize ! no. nodes of differing degrees in final
     ! level set
     integer :: minwid ! minimum levelset width
     integer :: nlvl ! number of levels in level structure
     integer :: minwid1

     j = 0
     mask = 0
     list = mask + a_n
     ! Initialise work(mask+1:mask+a_n) and work(list+1:list+a_n)
     work(mask+1:mask+a_n) = 1
     work(list+1:list+a_n) = 0
     level_ptr(:) = 0


     ! First guess for starting node is input nstrt

     ! Generate level structure for node nstrt
     call nd_level_struct(nstrt,a_n,a_ne,a_ptr,a_row, &
       work(mask+1:mask+a_n),level_ptr,level,maxdep,lwidth,num_entries)
     if (num_entries<a_n) then
       ! matrix is separable
       num_levels = maxdep
       return

     end if

     nstop = level(a_n)
     do main = 1, min(a_n,10)
       ! Store nodes in the final level set and their (weighted) degrees
       lsize = 0
       j = level_ptr(maxdep)
       do i = j, a_n
         node = level(i)
         lsize = lsize + 1
         work(list+lsize) = node
         level_ptr(node) = 0
         if (node==a_n) then
           k = a_ne
         else
           k = a_ptr(node+1) - 1
         end if
         do l = a_ptr(node), k
           level_ptr(node) = level_ptr(node) + a_weight(a_row(l))
         end do
       end do

       ! Choose at most max_search nodes
       do nlsize = 1, max_search
         ! Look for candiate with least degree
         mindeg = sumweight + 1
         ! mindeg = -1
         do i = nlsize, lsize
           if (level_ptr(work(list+i))<mindeg) then
             j = i
             mindeg = level_ptr(work(list+i))
           end if
         end do
         ! Jump out of loop if no candidates left
         if (mindeg==sumweight+1) go to 10
         ! if (mindeg .eq. -1) go to 55
         ! Swap chose candidate to next position
         node = work(list+j)
         work(list+j) = work(list+nlsize)
         work(list+nlsize) = node
         ! Rule out the neighbours of the chosen node
         if (node==a_n) then
           k = a_ne
         else
           k = a_ptr(node+1) - 1
         end if
         do i = a_ptr(node), k
           level_ptr(a_row(i)) = sumweight + 1
         end do
       end do
10        nlsize = nlsize - 1

       ! Loop over nodes in list
       minwid = huge(a_n)
       minwid1 = huge(a_n)

       do i = 1, nlsize
         node = work(list+i)

         ! Form rooted level structures for node
         call nd_level_struct(node,a_n,a_ne,a_ptr,a_row, &
           work(mask+1:mask+a_n),level_ptr,level,nlvl,lwidth,num_entries)
         ! if (lwidth .le. minwid) then
         lwidth1 = 0
         do k = 1, nlvl - 1
           do l = level_ptr(k), level_ptr(k+1) - 1
             ll = level(l)
             lwidth1 = lwidth1 + k*a_weight(ll)
           end do
         end do
         do l = level_ptr(nlvl), a_n
           ll = level(l)
           lwidth1 = lwidth1 + nlvl*a_weight(ll)
         end do

         if (nlvl>maxdep) then
           ! Level structure of greater depth. Begin a new iteration.
           nstrt = node
           maxdep = nlvl

           go to 20
         else
           if (lwidth1<minwid1) then

             nstop = node
             minwid = lwidth
             minwid1 = lwidth1
           end if
         end if
         ! end if
       end do
       go to 30
20        continue
     end do
30      if (nstop/=node) then
       call nd_level_struct(node,a_n,a_ne,a_ptr,a_row, &
         work(mask+1:mask+a_n),level_ptr,level,nlvl,lwidth,num_entries)
     end if
     num_levels = maxdep
     nend = nstop

     return
   end subroutine nd_find_pseudo


   ! ---------------------------------------------------
   ! nd_level_struct
   ! ---------------------------------------------------
   ! Given a root, calculate the level structure of a given graph from that
   ! root

   subroutine nd_level_struct(root,a_n,a_ne,a_ptr,a_row,mask,level_ptr, &
       level,num_levels,lwidth,num_entries)

     integer, intent(in) :: root ! Root node for level structure
     integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
     integer, intent(in) :: a_ne ! no. nonzeros of subproblem
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(inout) :: mask(a_n) ! Always restored to input value
     ! at
     ! end of the call. mask(node) > 0 for all visible nodes
     integer, intent(out) :: level_ptr(a_n) ! On output level_ptr(i)
     ! contains
     ! position in level that entries for level i start.
     integer, intent(out) :: level(a_n) ! On output level contains lists
     ! of
     ! rows according to the level set that they are in
     integer, intent(out) :: num_levels ! On output num_levels contains
     ! the
     ! number of levels
     integer, intent(out) :: lwidth ! On output, contains the width of the
     ! structure
     integer, intent(out) :: num_entries ! On output, contains number of
     ! entries in the tree structure containing root

     ! ---------------------------------------------
     ! Local variables
     integer :: lvlend ! End of previous level in level
     integer :: lnbr ! Next position in level
     integer :: lbegin ! Beginning of present level in level
     integer :: lw ! Level width
     integer :: nbr ! neighbour node
     integer :: node ! current node
     integer :: nvars
     integer :: i, j, k

     ! Based on mc60ld
     ! Establish level 1
     level_ptr(:) = 0
     mask(root) = -mask(root)
     level(1) = root
     lvlend = 0
     nvars = 0
     lnbr = 1
     lwidth = 1
     do num_levels = 1, a_n
       ! Generate next level by finding all unmasked neighbours of nodes in
       ! present level
       lbegin = lvlend + 1
       lvlend = lnbr
       level_ptr(num_levels) = lbegin
       lw = 0
       do i = lbegin, lvlend
         node = level(i)
         if (node==a_n) then
           k = a_ne
         else
           k = a_ptr(node+1) - 1
         end if
         do j = a_ptr(node), k
           nbr = a_row(j)
           if (mask(nbr)>0) then
             lnbr = lnbr + 1
             level(lnbr) = nbr
             mask(nbr) = -mask(nbr)
             ! lw = lw + a_weight(nbr)
             lw = lw + 1
           end if
         end do
       end do
       lwidth = max(lwidth,lw)
       nvars = nvars + lw
       ! If no neighbours found, we are done
       if (lnbr==lvlend) exit
       ! Abort construnction if level structure too wide
       ! if (lwidth .gt. maxwid) exit
     end do
     ! Reset mask
     do i = 1, lnbr
       mask(level(i)) = abs(mask(level(i)))
     end do
     num_entries = lnbr

   end subroutine nd_level_struct


   ! ---------------------------------------------------
   ! nd_distance
   ! ---------------------------------------------------
   ! Given two level structures, calculate the difference in each nodes
   ! distance
   ! from the start and end node
   subroutine nd_distance(a_n,num_levels_nend,level_ptr_nend,level_nend, &
       num_levels_nstrt,level_ptr_nstrt,level_nstrt,distance_ptr,distance)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: num_levels_nend ! number of levels with root
     ! nend
     integer, intent(in) :: level_ptr_nend(a_n) ! level_ptr(i) contains
     ! position in level that entries for level i start (root = nend)
     integer, intent(in) :: level_nend(a_n) ! Contains lists of rows
     ! according to the level set that they are in (root = nend)
     integer, intent(in) :: num_levels_nstrt ! no. of levels with root
     ! nstrt
     integer, intent(inout) :: level_ptr_nstrt(a_n) ! level_ptr(i)
     ! contains
     ! position in level that entries for level i start (root = nstrt)
     ! Reused during subroutine
     integer, intent(in) :: level_nstrt(a_n) ! Contains lists of rows
     ! according to the level set that they are in (root = nstrt)
     integer, intent(out) :: distance_ptr(2*a_n-1) ! distance(i) contains
     ! position in distance where entries with distance i-a_n
     integer, intent(out) :: distance(a_n) ! Contains lists of rows
     ! ordered
     ! according to their distance

     ! ---------------------------------------------
     ! Local variables
     integer :: j, k
     integer :: dptr ! used as offset into distance_ptr
     integer :: lev ! stores current level

     dptr = a_n
     distance_ptr(dptr+1-a_n:dptr+a_n-1) = 0
     distance(1:a_n) = 0

     ! set distance(i) to hold the level that i belongs to in the structure
     ! rooted at nend (= distance from end node + 1)
     do lev = 1, num_levels_nend - 1
       do j = level_ptr_nend(lev), level_ptr_nend(lev+1) - 1
         distance(level_nend(j)) = lev
       end do
     end do
     do j = level_ptr_nend(num_levels_nend), a_n
       distance(level_nend(j)) = num_levels_nend
     end do

     ! now consider level structure rooted at start node
     do lev = 1, num_levels_nstrt - 1
       do j = level_ptr_nstrt(lev), level_ptr_nstrt(lev+1) - 1
         distance(level_nstrt(j)) = distance(level_nstrt(j)) - lev
         k = distance(level_nstrt(j))
         distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
       end do
     end do
     do j = level_ptr_nstrt(num_levels_nstrt), a_n
       distance(level_nstrt(j)) = distance(level_nstrt(j)) - lev
       k = distance(level_nstrt(j))
       distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
     end do

     ! Copy distance into level_ptr_nstrt to save memory
     level_ptr_nstrt(1:a_n) = distance(1:a_n)

     ! Set distance_ptr(i) to point one place to the right of where entries
     ! with distance i will be
     do j = 1 - a_n, a_n - 1
       if (distance_ptr(dptr+j)>0) then
         exit
       else
         distance_ptr(dptr+j) = 1
       end if
     end do
     distance_ptr(dptr+j) = distance_ptr(dptr+j) + 1
     do k = j + 1, a_n - 1
       distance_ptr(dptr+k) = distance_ptr(dptr+k) + distance_ptr(dptr+k-1)
     end do

     ! Set up lists of rows with same value of distance
     do j = 1, a_n
       k = level_ptr_nstrt(j)
       distance_ptr(dptr+k) = distance_ptr(dptr+k) - 1
       distance(distance_ptr(dptr+k)) = j
     end do


   end subroutine nd_distance

   ! ---------------------------------------------------
   ! nd_convert_partition_flags
   ! ---------------------------------------------------
   ! Given a partition array, convert the partition into a flag array
   subroutine nd_convert_partition_flags(a_n,a_n1,a_n2,partition,flag_1, &
       flag_2,flag_sep,flags)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_n1 ! Size of partition 1
     integer, intent(in) :: a_n2 ! Size of partition 2
     integer, intent(in) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end.
     integer, intent(in) :: flag_1 ! flag for rows in partition 1
     integer, intent(in) :: flag_2 ! flag for rows in partition 2
     integer, intent(in) :: flag_sep ! flag for rows in separator
     integer, intent(out) :: flags(a_n) ! flags(i) contains flag for row i
     ! and indicates which partition it is in

     integer :: j, k

     do j = 1, a_n1
       k = partition(j)
       flags(k) = flag_1
     end do
     do j = a_n1 + 1, a_n1 + a_n2
       k = partition(j)
       flags(k) = flag_2
     end do
     do j = a_n1 + a_n2 + 1, a_n
       k = partition(j)
       flags(k) = flag_sep
     end do

   end subroutine nd_convert_partition_flags


   ! ---------------------------------------------------
   ! nd_flags_partition
   ! ---------------------------------------------------
   ! Given a partition array, convert the partition into a flag array
   subroutine nd_convert_flags_partition(a_n,a_n1,a_n2,flags,flag_1, &
       flag_2,partition)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_n1 ! Size of partition 1
     integer, intent(in) :: a_n2 ! Size of partition 2
     integer, intent(in) :: flags(a_n) ! flags(i) contains flag for row i
     ! and indicates which partition it is in
     integer, intent(in) :: flag_1 ! flag for rows in partition 1
     integer, intent(in) :: flag_2 ! flag for rows in partition 2
     integer, intent(out) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end.

     integer :: i, j, k, l

     i = 1
     j = a_n1 + 1
     k = a_n1 + a_n2 + 1
     do l = 1, a_n
       if (flags(l)==flag_1) then
         partition(i) = l
         i = i + 1
       else
         if (flags(l)==flag_2) then
           partition(j) = l
           j = j + 1
         else
           partition(k) = l
           k = k + 1
         end if
       end if
     end do

   end subroutine nd_convert_flags_partition

   ! ---------------------------------------------------
   ! nd_move_partition
   ! ---------------------------------------------------
   ! Given a flag array, move the separator by forming an edge separator
   ! between the input separator and the larger of P1 and P2
   subroutine nd_move_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
       a_weight_1,a_weight_2,a_weight_sep,flag_1,flag_2,flag_sep,flags)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i
     integer, intent(inout) :: a_n1 ! Size of partition 1
     integer, intent(inout) :: a_n2 ! Size of partition 2
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(in) :: flag_1 ! flag for rows in partition 1
     integer, intent(in) :: flag_2 ! flag for rows in partition 2
     integer, intent(in) :: flag_sep ! flag for rows in separator
     integer, intent(inout) :: flags(a_n) ! flags(i) contains flag for row
     ! i
     ! and indicates which partition it is in. This is updated

     integer :: part_no, a_n1_temp, a_n2_temp
     integer :: a_weight_1_temp, a_weight_2_temp, a_weight_sep_temp
     integer :: j, k, l
     logical :: stay_sep

     ! Decide whether new (large) separator is formed from
     ! intersection of the separator with partition 1 or the intersection
     ! of
     ! the separator with partition 2
     if (a_weight_1>a_weight_2) then
       ! Apply to partition 1
       part_no = flag_1
     else
       ! Apply to partition 2
       part_no = flag_2
     end if

     a_n1_temp = a_n1
     a_n2_temp = a_n2
     a_weight_1_temp = a_weight_1
     a_weight_2_temp = a_weight_2
     a_weight_sep_temp = a_weight_sep

     do j = 1, a_n
       if (flags(j)/=flag_sep) cycle
       ! j is in initial separator
       if (j==a_n) then
         k = a_ne
       else
         k = a_ptr(j+1) - 1
       end if
       stay_sep = .false.
       do l = a_ptr(j), k
         if (flags(a_row(l))==part_no) then
           stay_sep = .true.
           if (part_no==flag_1) then
             if (a_n1_temp>1) then
               a_n1_temp = a_n1_temp - 1
               flags(a_row(l)) = -1
               a_weight_sep_temp = a_weight_sep_temp + a_weight(a_row(l))
               a_weight_1_temp = a_weight_1_temp - a_weight(a_row(l))
             end if
           else
             if (a_n2_temp>1) then
               a_n2_temp = a_n2_temp - 1
               flags(a_row(l)) = -1
               a_weight_sep_temp = a_weight_sep_temp + a_weight(a_row(l))
               a_weight_2_temp = a_weight_2_temp - a_weight(a_row(l))
             end if
           end if
         else if (flags(a_row(l))==-1) then
           stay_sep = .true.
         end if
       end do
       if ( .not. stay_sep) then
         if (part_no==flag_1) then
           flags(j) = flag_2
         else
           flags(j) = flag_1
         end if
         a_weight_sep_temp = a_weight_sep_temp - a_weight(j)
         if (part_no==flag_1) then
           a_n2_temp = a_n2_temp + 1
           a_weight_2_temp = a_weight_2_temp + a_weight(j)
         else
           a_n1_temp = a_n1_temp + 1
           a_weight_1_temp = a_weight_1_temp + a_weight(j)
         end if
       end if
     end do

     do j = 1, a_n
       if (flags(j)==-1) then
         flags(j) = flag_sep
       end if
     end do

     a_n1 = a_n1_temp
     a_n2 = a_n2_temp
     a_weight_1 = a_weight_1_temp
     a_weight_2 = a_weight_2_temp
     a_weight_sep = a_weight_sep_temp

   end subroutine nd_move_partition


   ! ---------------------------------------------------
   ! nd_refine_edge
   ! ---------------------------------------------------
   ! Given a partition, refine the partition to improve the (weighted) value
   ! of the cost function. An edge separator is formed between the input
   ! separator and the larger partition, and this is then minimal using
   ! trimming or max flow
   subroutine nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
       a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,options)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i
     integer, intent(in) :: sumweight ! Sum of weights in a_weight
     integer, intent(inout) :: a_n1 ! Size of partition 1
     integer, intent(inout) :: a_n2 ! Size of partition 2
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end. This is updated to the new
     ! partition
     integer, intent(out) :: work(3*a_n) ! Work array
     type (nd_options), intent(in) :: options

     ! ---------------------------------------------
     ! Local variables
     integer :: fm_flags ! pointer into work array for start of
     ! flags from FM

     fm_flags = 0 ! length a_n

     ! Initialise work(fm_flags+1:fm_flags+a_n)
     call nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
       ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
       work(fm_flags+1:fm_flags+a_n))

     ! Create new separator by forming edge separator between input
     ! separator and largest of P1 and P2

     call nd_move_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
       a_weight_1,a_weight_2,a_weight_sep,ND_PART1_FLAG, &
       ND_PART2_FLAG,ND_SEP_FLAG,work(fm_flags+1:fm_flags+a_n))

     ! Update partition
     call nd_convert_flags_partition(a_n,a_n1,a_n2, &
       work(fm_flags+1:fm_flags+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
       partition(1:a_n))

     if (options%refinement>3) then
       call nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
         a_weight_1,a_weight_2,a_weight_sep,partition,work(1:8),options)
     else
       call nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
         a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work(1:3*a_n), &
         options)
     end if



   end subroutine nd_refine_edge

   ! ---------------------------------------------------
   ! nd_refine_fm
   ! ---------------------------------------------------
   ! Given a partition, refine the partition using FM refinement. Wrapper
   ! for code
   subroutine nd_refine_fm(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
       a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,options)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i
     integer, intent(in) :: sumweight ! Sum of weights in a_weight
     integer, intent(inout) :: a_n1 ! Size of partition 1
     integer, intent(inout) :: a_n2 ! Size of partition 2
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end. This is updated to the new
     ! partition
     integer, intent(out) :: work(8*a_n+sumweight) ! Work array
     type (nd_options), intent(in) :: options

     ! ---------------------------------------------
     ! Local variables
     integer :: fm_flags ! pointer into work array for start of
     ! flags from FM
     integer :: fm_ipart ! pointer into work array for start of
     ! ipart from FM
     integer :: fm_next ! pointer into work array for start of next
     ! from FM
     integer :: fm_last ! pointer into work array for start of last
     ! from FM
     integer :: fm_gain1 ! pointer into work array for start of
     ! gain1 from FM
     integer :: fm_gain2 ! pointer into work array for start of
     ! gain2 from FM
     integer :: fm_done ! pointer into work array for start of done
     ! from FM
     integer :: fm_head ! pointer into work array for start of head
     ! from FM
     integer :: fm_distance ! pointer into work array for start of head
     ! from FM
     integer :: icut, mult ! Used within FM refinement
     integer :: band
     real(wp) :: ratio
     logical :: imbal ! Should we check for imbalance?

     if (options%refinement_band .lt. 1) return

     ratio = max(real(1.0,wp),options%balance)
     if (ratio>real(sumweight-2)) then
       imbal = .false.
     else
       imbal = .true.
     end if
     fm_flags = 0 ! length a_n

     ! Initialise work(fm_flags+1:fm_flags+a_n)
     call nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
       ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
       work(fm_flags+1:fm_flags+a_n))

     fm_ipart = fm_flags + a_n ! length a_n
     fm_next = fm_ipart + a_n ! length a_n
     fm_last = fm_next + a_n ! length a_n
     fm_gain1 = fm_last + a_n ! length a_n
     fm_gain2 = fm_gain1 + a_n ! length a_n
     fm_done = fm_gain2 + a_n ! length a_n
     fm_head = fm_done + a_n ! length icut+mult+1
     icut = min(sumweight-1,3*(sumweight/a_n))
     icut = min(icut,5*maxVAL(a_weight))
     ! icut = sumweight/2
     ! mult = min(sumweight/20,10*sumweight/a_n) - 1
     mult = sumweight - icut - 1
     mult = min(mult,icut)
     ! mult = sumweight/2-1
     fm_distance = fm_head + icut + mult + 1

     band = min(options%refinement_band,a_n)

     call nd_fm_refinement(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,icut, &
       mult,options%cost_function,a_n1,a_n2,a_weight_1,a_weight_2,&
       a_weight_sep,band,ratio, &
       work(fm_flags+1:fm_flags+a_n),work(fm_ipart+1:fm_ipart+a_n), &
       work(fm_next+1:fm_next+a_n),work(fm_last+1:fm_last+a_n), &
       work(fm_gain1+1:fm_gain1+a_n),work(fm_gain2+1:fm_gain2+a_n), &
       work(fm_done+1:fm_done+a_n),work(fm_head+1:fm_head+icut+mult+1), &
       work(fm_distance+1:fm_distance+a_n))

     ! Update partition
     call nd_convert_flags_partition(a_n,a_n1,a_n2, &
       work(fm_flags+1:fm_flags+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
       partition(1:a_n))

   end subroutine nd_refine_fm




   ! The subroutine nd_fm_refinement uses a version of the Fiduccia-
   ! Mattheyses refinement algorithm on a tripartite partitioing of the
   ! nodes of a
   ! graph where a node in the first partition is not connected to any node
   ! in the second partition and any path between nodes in partition 1 and
   ! partition 2 must go through a node in the cutset (partition 3).
   ! The intention of the algorithm is to reduce f(P1,P2,P3), where
   ! f(P1,P2,P3) =   |P3|/(|P1||P2|) if min(|P1|,|P2|)/max(|P1|,|P2|) >=
   ! ratio
   ! =  sumweight - 2 + max(|P1|,|P2|)/min(|P1|,|P2|), otherwise
   ! This is a banded version so only nodes a distance of at most band from
   ! the
   ! input separator can be moved into the new separator

   subroutine nd_fm_refinement(n,a_ne,ptr,col,weight,sumweight,icut, &
       mult,costf,a_n1,a_n2,wnv1,wnv2,wns,band,ratio,flags,ipart,next,last,gain1, &
       gain2,done,head,dist)

     ! Matrix is held in matrix using compressed column scheme
     integer, intent(in) :: n ! size of matrix
     integer, intent(in) :: a_ne ! no. nonzeros in matrix
     integer, intent(in) :: ptr(n) ! row pointers
     integer, intent(in) :: col(a_ne) ! column indices
     ! type (nd_matrix), intent(inout) :: matrix
     ! The array weight is used to hold a weight on the vertices indicating
     ! how many vertices from the finer graphs have been combined into the
     ! current coarse graph vertex.
     integer, intent(in) :: weight(n)
     integer, intent(in) :: sumweight
     integer, intent(in) :: icut ! Used to limit search
     integer, intent(in) :: mult ! Used to bound search
     integer, intent(in) :: costf ! Determines which cost function used
     integer, intent(inout) :: a_n1 ! No. vertices partition 1
     integer, intent(inout) :: a_n2 ! No. vertices partition 2
     integer, intent(inout) :: wnv1 ! Weighted sum of vertices partition 1
     integer, intent(inout) :: wnv2 ! Weighted sum of vertices partition 2
     integer, intent(inout) :: wns ! Weighted sum of vertices separator
     integer, intent(in) :: band ! width of band around initial separator
     ! that the separator can lie in
     real(wp), intent(in) :: ratio ! ratio to determine
     ! whether
     ! partition is balanced

     ! flags holds a list of nodes stating which partition node i is in.
     ! The whole point of this routine is to return a revised partition
     ! with better properties.  Normally less nodes in the cutset while
     ! maintaining
     ! a balance between the number of nodes in the two components.
     ! flags(i) == ND_PART1_FLAG : i is in partition 1
     ! flags(i) == ND_PART2_FLAG : i is in partition 2
     ! flags(i) == ND_SEP_FLAG   : i is in separator/cutset
     integer, intent(inout) :: flags(n)
     ! info holds parameters giving information about the performance of
     ! the
     ! subroutine
     integer, intent(out) :: ipart(n), next(n), last(n)
     integer, intent(out) :: gain1(n), gain2(n), done(n), head(-mult:icut)
     integer, intent(out) :: dist(n)

     ! Number nodes in each partition
     integer :: nv1, ns, nv2, inv1, inv2, ins
     ! Weighted nodes in each partition
     integer :: winv1, winv2, wins
     integer :: i, j, ii, jj, eye, k, l
     integer :: inn, outer
     integer :: move, ming, gain, old_gain, inext, idummy
     integer :: first, tail
     real(wp) :: eval, evalc, evalo, eval1, eval2
     logical :: imbal


     if (ratio>real(sumweight-2)) then
       imbal = .false.
     else
       imbal = .true.
     end if

     ! Set-up distance to hold (min) distance of node from separator. Only
     ! find
     ! distance for nodes within band distance from separator
     first = 0
     tail = 0
     do i = 1, n
       if (flags(i)==ND_SEP_FLAG) then
         dist(i) = 0
         if (first==0) then
           first = i
           tail = i
         else
           next(tail) = i
           tail = i
         end if
       else
         dist(i) = -2
       end if
     end do

     do while (first/=0)
       j = dist(first)
       if (j==band-1) then
         if (first==n) then
           l = a_ne
         else
           l = ptr(first+1) - 1
         end if
         do i = ptr(first), l
           if (dist(col(i))==-2) then
             k = col(i)
             dist(k) = j + 1
           end if
         end do

       else
         if (first==n) then
           l = a_ne
         else
           l = ptr(first+1) - 1
         end if
         do i = ptr(first), l
           if (dist(col(i))==-2) then
             k = col(i)
             dist(k) = j + 1
             next(tail) = k
             tail = k
           end if
         end do

       end if
       if (first==tail) then
         first = 0
       else
         k = next(first)
         first = k
       end if
     end do
     next(1:n) = 0

     ! nv1,nv2,ns are the number of nodes in partitions 1, 2 and the cutset
     ! in
     ! the current partition
     ! inv1,inv2,ins,ipart are the equivalent quantities within the inner
     ! loop
     ! The same identifiers prefixed by w refer to weighted counts
     ! inner and outer are the two main loop indices

     ! Initialize nv1,nv2,ns
     nv1 = a_n1
     nv2 = a_n2
     ns = n - (a_n1+a_n2)
     ii = 1
     jj = ii + nv1

     ! Initialize ipart
     ipart(1:n) = flags(1:n)

     ! Initialize array done that flags that a node has been considered in
     ! an inner loop pass
     done = 0

     ! Compute evaluation function for current partitioning

     call cost_function(wnv1+1,wnv2+1,wns,sumweight,ratio,imbal,costf,evalc)

     ! icut is set to limit search in inner loop .. may later be a
     ! parameter
     ! we allow gains of up to max(weight)*5

     head(-mult:icut) = 0

     ! Set up doubly linked list linking nodes with same gain and headers
     ! to starts (up to cut off value icut)

     ! Compute gains for nodes in cutset
     ming = sumweight
     do i = 1, n

       if (flags(i)==ND_SEP_FLAG) then
         ! Node i is in cutset
         ! gain1(i) is change to cutset size if node i is moved to
         ! partition 1.
         ! gain2(i) is change to cutset size if node i is moved to
         ! partition 2.
         ! Run through all neighbours of node i to see what loss/gain is if
         ! node
         ! i is moved

         call compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,flags)
         gain = max(-mult,min(gain1(i),gain2(i)))

         if (gain<ming) ming = gain
         if (gain>icut) cycle
         ! New node is put at head of list
         call add_to_list(n,mult,icut,next,last,head,i,gain)
       end if
     end do

     ! Initilialization finished.  Now perform F-M algorithm in two loops.
     ! In each inner loop we choose the best obtained and if the evaluation
     ! function for this is better than previous best we perform another
     ! inner
     ! loop; otherwise we terminate.
     evalo = evalc
     do outer = 1, n
       ! Set partition that will be altered in inner loop
       inv1 = nv1
       inv2 = nv2
       ins = ns
       winv1 = wnv1
       winv2 = wnv2
       wins = wns
       ipart(1:n) = flags(1:n)
inNER:    do inn = 1, n

         ! Choose best eligible move
         do idummy = 1, n

           do gain = ming, icut
             if (head(gain)/=0) exit
           end do
           if (gain>icut) exit inNER

           ! Now cycle through nodes of least gain
           ! Currently inefficient because of re-searching linked list
           inext = head(gain)
           k = 0
10            i = inext
           if (i==0) cycle
           ! Use node if it has not been considered already
           if (done(i)<outer) go to 20
           inext = next(i)
           ! !! Extra statements to trap infinite loop
           k = k + 1
           if (k>ins) then
             ! write (*,*) 'Bug in code because of infinite loop'
             ! !! You may wish to change this to a stop
             exit
           end if
           go to 10
         end do
         exit inNER
         ! Node i has been selected as the best eligible node
         ! Set flag so only considered once in this pass
20          done(i) = outer
         ! As i will not be chosen again in this pass, remove from list
         call remove_from_list(n,mult,icut,next,last,head,i,gain)
         ! Move the node to the appropriate partition and reset partition
         ! information
         ! We will try both weighted and unweighted

         if (wnv1==0 .and. wnv2>0) then

           ! Move node i to partition 1
           move = ND_PART1_FLAG
           inv1 = inv1 + 1
           winv1 = winv1 + weight(i)
           ins = ins - 1
           wins = wins - weight(i)
         else if (wnv2==0 .and. wnv1>0) then
           ! Move node i to partition 2
           move = ND_PART2_FLAG
           inv2 = inv2 + 1
           winv2 = winv2 + weight(i)
           ins = ins - 1
           wins = wins - weight(i)

         else
           call cost_function(winv1+weight(i),winv2+1-gain1(i)-weight(i), &
             wins+gain1(i)-1,sumweight,ratio,imbal,costf,eval1)

           call cost_function(winv1+1-gain2(i)-weight(i),winv2+weight(i), &
             wins+gain2(i)-1,sumweight,ratio,imbal,costf,eval2)
           if ((eval1<eval2) .or. ((eval1==eval2) .and. (wnv1<wnv2))) then
             ! Move node i to partition 1
             move = ND_PART1_FLAG
             inv1 = inv1 + 1
             winv1 = winv1 + weight(i)
             ins = ins - 1
             wins = wins - weight(i)
           else
             ! Move node i to partition 2
             move = ND_PART2_FLAG
             inv2 = inv2 + 1
             winv2 = winv2 + weight(i)
             ins = ins - 1
             wins = wins - weight(i)
           end if
         end if
         ! Set new partition for node i
         ipart(i) = move
         ! Run through neigbours of node i to update data
         if (i==n) then
           l = a_ne
         else
           l = ptr(i+1) - 1
         end if
         do jj = ptr(i), l
           j = col(jj)
           ! Check which partition node j is in and take appropriate action
           if (ipart(j)==move) cycle
           ! If node j is in cutset, update its gain value
           if (ipart(j)==ND_SEP_FLAG) then
             ! If it has already been chosen in this pass just skip it
             if (done(j)==outer .or. dist(j)==-2) cycle
             ! old_gain is present gain

             old_gain = max(-mult,min(gain1(j),gain2(j)))
             ! old_gain = min(gain1(j),gain2(j))

             if (move==ND_PART1_FLAG) gain2(j) = gain2(j) + weight(i)
             if (move==ND_PART2_FLAG) gain1(j) = gain1(j) + weight(i)
             gain = max(-mult,min(gain1(j),gain2(j)))
             ! gain = min(gain1(j),gain2(j))

             if (old_gain==gain) cycle
             ! Remove from old list
             if (old_gain<=icut) then
               call remove_from_list(n,mult,icut,next,last,head,j,old_gain)
             end if
             ! gain has changed so move to new linked list if less than
             ! icut
             if (gain<=icut) then
               ! Reset ming if necessary
               if (gain<ming) ming = gain
               call add_to_list(n,mult,icut,next,last,head,j,gain)
             end if
           end if
           if (ipart(j)==2-move) then
             ! We have a new node in the cutset.
             ipart(j) = ND_SEP_FLAG
             ! Compute gains for this new node in the cutset and place in
             ! linked list
             ! We intentionally did not do this earlier but we do now
             ! [maybe not since won't access this node again in this pass]
             ! We use done array to record this but probably not necessary
             ! as not put
             ! in head linked list so won't be accessed
             ! First check that it was not earlier moved from cutset
             if (done(j)/=outer .and. dist(j)/=-2) then
               ! Compute gain
               call compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,j,ipart)
               gain = max(-mult,min(gain1(j),gain2(j)))
               ! gain = min(gain1(j),gain2(j))
               ! !! Just added this
               if (gain<ming) ming = gain
               ! Add to  list
               if (gain<=icut) then
                 call add_to_list(n,mult,icut,next,last,head,j,gain)
               end if
             end if
             ! Update partition and gain of any nodes in cutset connected
             ! to node j
             ins = ins + 1
             wins = wins + weight(j)
             if (move==ND_PART1_FLAG) then
               inv2 = inv2 - 1
               winv2 = winv2 - weight(j)
             end if
             if (move==ND_PART2_FLAG) then
               inv1 = inv1 - 1
               winv1 = winv1 - weight(j)
             end if
             ! Check neighbours of j since any in cut set will have gain
             ! changed
             if (j==n) then
               l = a_ne
             else
               l = ptr(j+1) - 1
             end if
             do ii = ptr(j), l
               eye = col(ii)
               if (ipart(eye)/=ND_SEP_FLAG) cycle
               if (dist(eye)==-2) cycle
               ! Neighbour is in cutset. Recompute gain and insert in
               ! linked list.
               if (done(eye)==outer) cycle
               ! old_gain is present gain
               old_gain = max(-mult,min(gain1(eye),gain2(eye)))
               ! old_gain = min(gain1(eye),gain2(eye))


               if (move==ND_PART1_FLAG) then
                 gain1(eye) = gain1(eye) - weight(j)
               end if
               if (move==ND_PART2_FLAG) then
                 gain2(eye) = gain2(eye) - weight(j)
               end if
               ! gain is new gain
               gain = max(-mult,min(gain1(eye),gain2(eye)))
               ! gain = min(gain1(eye),gain2(eye))
               if (old_gain==gain) cycle
               ! Remove from old list
               if (old_gain<=icut) then
                 call remove_from_list(n,mult,icut,next,last,head,eye,old_gain)
               end if
               ! gain has changed so move to new linked list if less than
               ! icut
               if (gain<=icut) then
                 ! Reset ming if necessary
                 if (gain<ming) ming = gain
                 call add_to_list(n,mult,icut,next,last,head,eye,gain)
               end if
             end do
           end if
           ! end of neighbours loop
         end do

         ! ii = 0
         ! do i = 1,n
         ! if (ipart(i) == 2) ii = ii + 1
         ! enddo
         ! if (ii .ne. inv2) write(6,*) 'problem in partition',ii,inv2

         ! Evaluate new partition
         call cost_function(winv1+1,winv2+1,wins,sumweight,ratio,imbal, &
           costf,eval)
         ! Compare this with best so far in inner loop and store partition
         ! information if it is the best
         if (inv1*inv2>0 .and. nv1*nv2==0) then
           ! Might have to store gains and who is in the cutset
           evalc = eval
           nv1 = inv1
           nv2 = inv2
           ns = ins
           wnv1 = winv1
           wnv2 = winv2
           wns = wins
           flags = ipart

         else if (eval<evalc .and. (inv1*inv2>0)) then
           ! Might have to store gains and who is in the cutset
           evalc = eval
           nv1 = inv1
           nv2 = inv2
           ns = ins
           wnv1 = winv1
           wnv2 = winv2
           wns = wins
           flags(1:n) = ipart(1:n)
         end if
         ! End inner loop
       end do inNER
       ! Leave loop if inner loop has not found better partition
       if (evalc>=(1.0-1.0/(LOG(real(sumweight))**2.3))*evalo) exit
       ! Otherwise we reset evalo and go back to inner loop
       evalo = evalc
       ! Recompute gains for this new partition
       ! Compute gains for nodes in cutset
       ! This is very inefficient but is in now to test functionality
       head(-mult:icut) = 0
       ming = icut + 1
       do i = 1, n
         if (flags(i)/=ND_SEP_FLAG) cycle
         if (dist(i)==-2) cycle
         ! Node i is in cutset
         ! gain1(i) is change to cutset size if node i is moved to
         ! partition 1.
         ! gain2(i) is change to cutset size if node i is moved to
         ! partition 2.
         ! Run through all neighbours of node i to see what loss/gain is if
         ! node
         ! i is moved
         call compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,flags)
         ! Recalculate doubly linked list linking nodes with same gain and
         ! headers
         ! to starts (up to cut off value icut)
         ! Initialize array done that flags that a node has been considered
         ! in
         ! an inner loop pass
         gain = max(-mult,min(gain1(i),gain2(i)))
         ! gain = min(gain1(i),gain2(i))
         if (gain>icut) cycle
         if (gain<ming) ming = gain
         ! New node is put at head of list
         call add_to_list(n,mult,icut,next,last,head,i,gain)
       end do
       ! End of outer loop
     end do
     a_n1 = nv1
     a_n2 = nv2
     return

   end subroutine nd_fm_refinement


     subroutine remove_from_list(n,mult,icut,next,last,head,irm,ig)
       integer, intent(in) :: n,mult,icut ! order matrix
       integer, intent(inout) :: next(n),last(n),head(-mult:icut)
       integer, intent(in) :: irm, ig
       integer :: inext, ilast

       inext = next(irm)
       ilast = last(irm)
       if (ilast==0) then
         head(ig) = inext
         if (inext/=0) last(inext) = 0
       else
         next(ilast) = inext
         if (inext/=0) last(inext) = ilast
       end if
     end subroutine remove_from_list


     subroutine add_to_list(n,mult,icut,next,last,head,irm,ig)
       integer, intent(in) :: n,mult,icut ! order matrix
       integer, intent(inout) :: next(n),last(n),head(-mult:icut)
       integer, intent(in) :: irm, ig
       integer :: inext

       inext = head(ig)
       head(ig) = irm
       next(irm) = inext
       if (inext/=0) last(inext) = irm
       last(irm) = 0
     end subroutine add_to_list


     subroutine compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,partit)
       integer, intent(in) :: n,a_ne,ptr(n),col(a_ne),weight(n)
       integer, intent(inout) :: gain1(n), gain2(n)
       integer, intent(in) :: i, partit(:)
       integer :: j, jj, l
       ! Initialize gain ... knowing node i will be removed from cutset
       ! The +1 is to give identical result to previous code when unit
       ! weights
       gain1(i) = -weight(i) + 1
       gain2(i) = -weight(i) + 1
       if (i==n) then
         l = a_ne
       else
         l = ptr(i+1) - 1
       end if
       do jj = ptr(i), l
         j = col(jj)
         ! Check which partition node j is in and adjust gain array
         ! appropriately
         if (partit(j)==ND_PART1_FLAG) then
           gain2(i) = gain2(i) + weight(j)
         end if
         if (partit(j)==ND_PART2_FLAG) then
           gain1(i) = gain1(i) + weight(j)
         end if
       end do
     end subroutine compute_gain

   subroutine amd_order_both(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
       iperm,a_weight,work)
     integer, intent(in) :: a_n ! order of matrix being partitioned
     integer, intent(in) :: a_ne ! no. entries in matrix being partitioned
     integer, intent(in) :: a_ptr(a_n) ! col ptrs for matrix being
     ! partitioned
     integer, intent(in) :: a_row(a_ne) ! row indices for matrix
     ! being partitioned.
     integer, intent(in) :: a_n1 ! no. rows in partition 1
     integer, intent(in) :: a_n2 ! no. rows in partition 2
     integer, intent(in) :: partition(a_n) ! the partitions
     integer, intent(inout) :: iperm(a_n) ! maps current permutation to
     ! the
     ! column indices of the matrix whose ordering is being computed
     integer, intent(inout) :: a_weight(a_n) ! weights of vertices
     integer, intent(out) :: work(12*a_n+a_ne)


     ! Local variables
     integer :: i, j
     integer :: extract_work ! pointers into work array for mask arrays
     integer :: amd_order_perm ! pointer into work array for perm array
     integer :: amd_order_work ! pointer into work array for amd_order work
     ! array
     integer :: a_ptr_sub ! pointer into work for col ptrs of
     ! submatrix
     integer :: a_irn_sub ! pointer into work for irn array of
     ! submatrix
     integer :: rows_sub ! pointer into work for rows_sub array
     integer :: a_lirn_sub ! length of irn array of submatrix
     integer :: a_n_1 ! order of submatrix 1
     integer :: a_n_2 ! order of submatrix 2
     integer :: a_n_sep ! number entries in separator
     integer :: a_ne_sub ! number entries in submatrix
     integer :: len_a_row_sub ! used when extracting submatrices


     ! Set orders of submatrices
     a_n_1 = a_n1
     a_n_2 = a_n2
     a_n_sep = 0


     ! Set pointers into work array
     amd_order_perm = 0 ! length a_n
     a_ptr_sub = amd_order_perm + a_n ! length a_n
     a_lirn_sub = a_ne + max(a_n_1,a_n_2) + 1
     ! max(a_n_1,a_n_2) + 1 .le. a_n
     a_irn_sub = a_ptr_sub + a_n ! length a_lirn_sub
     rows_sub = a_irn_sub + a_lirn_sub ! length a_n
     extract_work = rows_sub + a_n ! length a_n
     amd_order_work = rows_sub + a_n ! length 7*a_n


     ! Form submatrix 1
     ! if (a_n1 .ne. 1) then
     work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
     len_a_row_sub = a_ne
     call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,a_n_sep, &
       work(rows_sub+1:rows_sub+a_n_1),a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+ &
       a_n_1),len_a_row_sub,work(a_irn_sub+1:a_irn_sub+a_ne), &
       work(extract_work+1:extract_work+a_n))

     ! Apply amd_order
     work(rows_sub+1:rows_sub+a_n1) = ND_PART1_FLAG
     call amd_order(a_n_1,a_ne_sub,a_lirn_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_1), &
       work(a_irn_sub+1:a_irn_sub+a_lirn_sub), &
       work(rows_sub+1:rows_sub+a_n_1),work(amd_order_perm+1:amd_order_perm+a_n_1), &
       work(amd_order_work+1:amd_order_work+7*a_n_1))

     ! Overwrite first a_n1 entries of amd_order_perm with first a_n1 entries
     ! that will form new iperm. Similarly, overwrite first a_n1 entries of
     ! rows_sub with first a_n1 entries that will form new a_weight
     ! no longer need info in a_ptr
     do i = 1, a_n1
       j = work(amd_order_perm+i)
       work(a_ptr_sub+i) = iperm(partition(j))
       work(rows_sub+i) = a_weight(partition(j))
     end do
     work(amd_order_perm+1:amd_order_perm+a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)

     ! Build second submatrix
     ! if (a_n2 .ne. 1) then
     work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2) = partition(a_n1+1:a_n1+a_n2 &
       )
     len_a_row_sub = a_ne
     call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,a_n_sep, &
       work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2),a_ne_sub, &
       work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub, &
       work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
       a_n))
     ! Apply amd_order
     work(rows_sub+a_n1+1:rows_sub+a_n1+a_n2) = ND_PART1_FLAG
     call amd_order(a_n_2,a_ne_sub,a_lirn_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_2), &
       work(a_irn_sub+1:a_irn_sub+a_lirn_sub), &
       work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2), &
       work(amd_order_perm+a_n1+1:amd_order_perm+a_n),work(amd_order_work+1:amd_order_work+7* &
       a_n_2))
     ! else
     ! work(amd_order_perm+a_n1+1) = 1
     ! end if
     do i = 1, a_n_2
       j = work(amd_order_perm+a_n1+i)
       work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
       work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))
     end do
     do i = a_n_2 + 1, a_n_2 + (a_n-a_n1-a_n2)
       j = i
       work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
       work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))

     end do
     iperm(a_n1+1:a_n) = work(a_ptr_sub+1+a_n1:a_ptr_sub+a_n)
     iperm(1:a_n1) = work(amd_order_perm+1:amd_order_perm+a_n1)
     a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

   end subroutine amd_order_both


   subroutine extract_matrix(a_n,a_ne,a_ptr,a_row,a_n_part,a_n_sep, &
       rows_sub,a_ne_sub,a_ptr_sub,len_a_row_sub,a_row_sub,work)
     integer, intent(in) :: a_n ! order of matrix being partitioned
     integer, intent(in) :: a_ne ! no. entries in matrix being partitioned
     integer, intent(in) :: a_ptr(a_n) ! col ptrs for matrix being
     ! partitioned
     integer, intent(in) :: a_row(a_ne) ! row indices for matrix
     ! being partitioned.
     integer, intent(in) :: a_n_part ! no. rows in partition
     integer, intent(in) :: a_n_sep ! no. rows in partition
     integer, intent(in) :: rows_sub(a_n_part+a_n_sep) ! rows/cols of
     ! matrix
     ! to be extracted. Intersecting rows/cols of separator will be
     ! replaced
     ! by matrix of all zeros
     integer, intent(out) :: a_ne_sub ! no. entries stored in extracted
     ! matrix
     integer, intent(out) :: a_ptr_sub(a_n_part+a_n_sep) ! col ptrs for
     ! extracted matrix
     integer, intent(in) :: len_a_row_sub ! length of a_row_sub
     integer, intent(out) :: a_row_sub(len_a_row_sub) ! row indices for
     ! extracted matrix
     integer, intent(out) :: work(a_n)

     ! Local variables
     integer :: i, j, k, l, m, p
     integer :: a_n_sub ! Order of extracted matrix
     integer :: mask ! pointer into work array for mask arrays

     ! Set pointers into work array
     mask = 0 ! length a_n

     ! Set mask
     a_n_sub = a_n_part + a_n_sep
     work(mask+1:mask+a_n) = 0
     do i = 1, a_n_part
       j = rows_sub(i)
       work(mask+j) = i
     end do
     do i = a_n_part + 1, a_n_sub
       j = rows_sub(i)
       work(mask+j) = -i
     end do
     a_row_sub(:) = 0

     ! Count number of entries in  submatrix and set-up column ptrs
     a_ptr_sub(1:a_n_sub) = 0
     do j = 1, a_n_part
       a_ptr_sub(j) = 0
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)/=0) then
           a_ptr_sub(j) = a_ptr_sub(j) + 1

         end if
       end do
     end do
     do j = a_n_part + 1, a_n_sub
       a_ptr_sub(j) = 0
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)>0) a_ptr_sub(j) = a_ptr_sub(j) + 1
       end do
     end do
     a_ptr_sub(1) = a_ptr_sub(1) + 1
     do j = 2, a_n_sub
       a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
     end do
     a_ne_sub = a_ptr_sub(a_n_sub) - 1

     ! Form a_row_sub
     do j = 1, a_n_part
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)/=0) then
           a_ptr_sub(j) = a_ptr_sub(j) - 1
           p = a_ptr_sub(j)
           a_row_sub(p) = abs(work(mask+m))
         end if
       end do
     end do
     do j = a_n_part + 1, a_n_sub
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)>0) then
           a_ptr_sub(j) = a_ptr_sub(j) - 1
           p = a_ptr_sub(j)
           a_row_sub(p) = work(mask+m)
         end if
       end do
     end do
   end subroutine extract_matrix

   subroutine extract_both_matrices(a_n,a_ne,a_ptr,a_row,a_n_part1, &
       a_n_part2,rows_sub,a_ne_sub1,a_ne_sub2,a_ptr_sub,len_a_row_sub, &
       a_row_sub,work)
     integer, intent(in) :: a_n ! order of matrix being partitioned
     integer, intent(in) :: a_ne ! no. entries in matrix being partitioned
     integer, intent(in) :: a_ptr(a_n) ! col ptrs for matrix being
     ! partitioned
     integer, intent(in) :: a_row(a_ne) ! row indices for matrix
     ! being partitioned.
     integer, intent(in) :: a_n_part1 ! no. rows in partition 1
     integer, intent(in) :: a_n_part2 ! no. rows in partition 2
     integer, intent(in) :: rows_sub(a_n_part1+a_n_part2) ! rows/cols of
     ! matrices to be extracted. First a_n_part1 entries contain the
     ! rows/cols forming the first matrix to be extracted.
     integer, intent(out) :: a_ne_sub1 ! no. entries in extracted matrix 1
     integer, intent(out) :: a_ne_sub2 ! no. entries in extracted matrix 2
     integer, intent(out) :: a_ptr_sub(a_n_part1+a_n_part2) ! col ptrs for
     ! extracted matrices. First a_n_part1 are for first matrix, etc
     integer, intent(in) :: len_a_row_sub ! length of a_row_sub
     integer, intent(out) :: a_row_sub(len_a_row_sub) ! row indices for
     ! extracted matrices. First a_ne_part1 entries are for first matrix;
     ! the immediately following a_ne_part2 entries are for second matrix
     integer, intent(out) :: work(a_n)

     ! Local variables
     integer :: i, j, k, l, m, p
     integer :: mask ! pointer into work array for mask arrays


     ! Set pointers into work array
     mask = 0 ! length a_n

     ! Set mask
     work(mask+1:mask+a_n) = 0
     do i = 1, a_n_part1
       j = rows_sub(i)
       work(mask+j) = i
     end do
     do i = a_n_part1 + 1, a_n_part1 + a_n_part2
       j = rows_sub(i)
       work(mask+j) = -i + a_n_part1
     end do
     a_row_sub(:) = 0

     ! Count number of entries in each submatrix and set-up column ptrs
     a_ptr_sub(1:a_n_part1+a_n_part2) = 0
     do j = 1, a_n_part1
       a_ptr_sub(j) = 0
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)>0) then
           a_ptr_sub(j) = a_ptr_sub(j) + 1

         end if
       end do
     end do

     do j = a_n_part1 + 1, a_n_part1 + a_n_part2
       a_ptr_sub(j) = 0
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)<0) a_ptr_sub(j) = a_ptr_sub(j) + 1
       end do
     end do

     a_ptr_sub(1) = a_ptr_sub(1) + 1
     do j = 2, a_n_part1
       a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
     end do
     a_ne_sub1 = a_ptr_sub(a_n_part1) - 1

     a_ptr_sub(a_n_part1+1) = a_ptr_sub(a_n_part1+1) + 1
     do j = a_n_part1 + 2, a_n_part1 + a_n_part2
       a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
     end do
     a_ne_sub2 = a_ptr_sub(a_n_part1+a_n_part2) - 1

     ! Form a_row_sub
     do j = 1, a_n_part1
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)>0) then
           a_ptr_sub(j) = a_ptr_sub(j) - 1
           p = a_ptr_sub(j)
           a_row_sub(p) = abs(work(mask+m))
         end if
       end do
     end do

     do j = a_n_part1 + 1, a_n_part1 + a_n_part2
       i = rows_sub(j)
       if (i==a_n) then
         l = a_ne
       else
         l = a_ptr(i+1) - 1
       end if
       do k = a_ptr(i), l
         m = a_row(k)
         if (work(mask+m)<0) then
           a_ptr_sub(j) = a_ptr_sub(j) - 1
           p = a_ptr_sub(j)
           a_row_sub(p+a_ne_sub1) = -work(mask+m)
         end if
       end do
     end do
   end subroutine extract_both_matrices

   subroutine amd_order_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,iperm, &
       a_weight,a_ne_sub,work)
     ! Apply amd_order to the smaller partition and extract the other matrix
     ! into
     ! a_ptr and a_row
     integer, intent(in) :: a_n ! order of matrix being partitioned
     integer, intent(in) :: a_ne ! no. entries in matrix being partitioned
     integer, intent(inout) :: a_ptr(a_n) ! col ptrs for matrix being
     ! partitioned and, on return, for the extracted submatrix
     integer, intent(inout) :: a_row(a_ne) ! row indices for matrix
     ! being partitioned and, on return, for the extracted submatrix
     integer, intent(in) :: a_n1 ! no. rows in partition 1
     integer, intent(in) :: a_n2 ! no. rows in partition 2
     integer, intent(in) :: partition(a_n) ! the partitions
     integer, intent(inout) :: iperm(a_n) ! maps current permuation to the
     ! column indices of the matrix whose ordering is being computed
     integer, intent(inout) :: a_weight(a_n) ! weights of vertices
     integer, intent(out) :: a_ne_sub ! number entries in returned
     ! submatrix
     integer, intent(out) :: work(12*a_n+a_ne)

     ! Local variables
     integer :: i, j
     integer :: extract_work ! pointers into work array for mask arrays
     integer :: amd_order_perm ! pointer into work array for perm array
     integer :: amd_order_work ! pointer into work array for amd_order work
     ! array
     integer :: a_ptr_sub ! pointer into work for col ptrs of
     ! submatrix
     integer :: a_irn_sub ! pointer into work for irn array of
     ! submatrix
     integer :: rows_sub ! pointer into work for rows_sub array
     integer :: a_lirn_sub ! length of irn array of submatrix
     integer :: a_n_1 ! order of submatrix 1
     integer :: a_n_2 ! order of submatrix 2
     integer :: a_n_sep ! number entries in separator
     integer :: len_a_row_sub ! used when extracting submatrices

     ! Set orders of submatrices
     if (a_n1<a_n2) then
       ! Applying amd_order to first partition
       a_n_1 = a_n - a_n2
       a_n_2 = a_n2
       a_n_sep = a_n - a_n1 - a_n2
     else
       ! Applying amd_order to second partition
       a_n_2 = a_n - a_n1
       a_n_1 = a_n1
       a_n_sep = a_n - a_n1 - a_n2

     end if


     ! Set pointers into work array
     amd_order_perm = 0 ! length a_n
     a_ptr_sub = amd_order_perm + a_n ! length a_n
     if (a_n1<a_n2) then
       a_lirn_sub = a_ne + a_n_1 + 1
     else
       a_lirn_sub = a_ne + a_n_2 + 1
     end if
     ! max(a_n_1,a_n_2) + 1 .le. a_n
     a_irn_sub = a_ptr_sub + a_n ! length a_lirn_sub
     rows_sub = a_irn_sub + a_lirn_sub ! length a_n
     extract_work = rows_sub + a_n ! length a_n
     amd_order_work = rows_sub + a_n ! length 7*a_n

     ! Extract matrix that amd_order is applied to
     if (a_n1<a_n2) then

       ! Start by updating iperm and a_weight
       do i = 1, a_n
         work(amd_order_perm+i) = iperm(partition(i))
         work(rows_sub+i) = a_weight(partition(i))
       end do
       iperm(1:a_n) = work(amd_order_perm+1:amd_order_perm+a_n)
       a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

       ! Form submatrix 1
       work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
       work(rows_sub+a_n1+1:rows_sub+a_n_1) = partition(a_n1+a_n2+1:a_n)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,a_n_sep, &
         work(rows_sub+1:rows_sub+a_n_1),a_ne_sub, &
         work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub, &
         work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
         a_n))


       ! Apply amd_order
       work(rows_sub+1:rows_sub+a_n1) = ND_PART1_FLAG
       work(rows_sub+a_n1+1:rows_sub+a_n_1) = ND_SEP_FLAG
       call amd_order(a_n_1,a_ne_sub,a_lirn_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_1), &
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub), &
         work(rows_sub+1:rows_sub+a_n_1),work(amd_order_perm+1:amd_order_perm+a_n_1), &
         work(amd_order_work+1:amd_order_work+7*a_n_1))

       ! Overwrite first a_n1 entries of amd_order_perm with first a_n1 entries
       ! that will form new iperm. Similarly, overwrite first a_n1 entries
       ! of
       ! rows_sub with first a_n1 entries that will form new a_weight
       ! no longer need info in a_ptr
       ! no longer need info in a_ptr
       do i = 1, a_n1
         j = work(amd_order_perm+i)
         work(a_ptr_sub+i) = iperm(j)
         work(rows_sub+i) = a_weight(j)
       end do
       ! work(amd_order_perm+1:amd_order_perm+a_n1) =
       ! work(a_ptr_sub+1:a_ptr_sub+a_n1)


       iperm(1:a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)
       a_weight(1:a_n1) = work(rows_sub+1:rows_sub+a_n1)

       ! Build second submatrix
       work(rows_sub+1:rows_sub+a_n_2) = partition(a_n1+1:a_n1+a_n2)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,0, &
         work(rows_sub+1:rows_sub+a_n_2),a_ne_sub, &
         work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub, &
         work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
         a_n))

       ! Copy matrix to a_ptr and a_row
       a_ptr(1:a_n_2) = work(a_ptr_sub+1:a_ptr_sub+a_n_2)
       a_row(1:a_ne_sub) = work(a_irn_sub+1:a_irn_sub+a_ne_sub)

     else


       ! Form submatrix 2
       work(rows_sub+1:rows_sub+a_n_2) = partition(a_n1+1:a_n)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,a_n_sep, &
         work(rows_sub+1:rows_sub+a_n_2),a_ne_sub, &
         work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub, &
         work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
         a_n))

       ! Apply amd_order
       work(rows_sub+1:rows_sub+a_n2) = ND_PART2_FLAG
       work(rows_sub+a_n2+1:rows_sub+a_n_2) = ND_SEP_FLAG
       call amd_order(a_n_2,a_ne_sub,a_lirn_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_2), &
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub), &
         work(rows_sub+1:rows_sub+a_n_2),work(amd_order_perm+a_n1+1:amd_order_perm+ &
         a_n),work(amd_order_work+1:amd_order_work+7*a_n_2))

       do i = 1, a_n_2
         j = work(amd_order_perm+a_n1+i)
         work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
         work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))
       end do

       do i = 1, a_n_1
         work(a_ptr_sub+i) = iperm(partition(i))
         work(rows_sub+i) = a_weight(partition(i))
       end do
       iperm(1:a_n) = work(a_ptr_sub+1:a_ptr_sub+a_n)
       a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

       ! Form submatrix 1
       work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,0, &
         work(rows_sub+1:rows_sub+a_n_1),a_ne_sub, &
         work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub, &
         work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
         a_n))

       ! Copy matrix to a_ptr and a_row
       a_ptr(1:a_n_1) = work(a_ptr_sub+1:a_ptr_sub+a_n_1)
       a_row(1:a_ne_sub) = work(a_irn_sub+1:a_irn_sub+a_ne_sub)

     end if

   end subroutine amd_order_one

   subroutine multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
       partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,options, &
       info1,lwork,work,stop_coarsening2,grid)

     integer, intent(in) :: a_n ! order of matrix being partitioned
     integer, intent(in) :: a_ne ! no. entries in matrix being partitioned
     integer, intent(in) :: a_ptr(a_n) ! col ptrs for matrix being
     ! partitioned and, on return, for the extracted submatrix
     integer, intent(in) :: a_row(a_ne) ! row indices for matrix
     integer, intent(in) :: a_weight(a_n) ! weights associated with rows
     ! of matrix (useful if matrix has already been compressed)
     integer, intent(in) :: sumweight ! sum of entries in a_weight
     integer, intent(out) :: partition(a_n) ! computed partition
     integer, intent(out) :: a_n1 ! number of entries in partition 1
     integer, intent(out) :: a_n2 ! number of entries in partition 2
     integer, intent(out) :: a_weight_1, a_weight_2, a_weight_sep ! Weight
     ! ed
     ! size of partitions and separator
     type (nd_options), intent(in) :: options
     integer, intent(in) :: lwork ! length of work array: must be atleast
     ! 9a_n + sumweight
     integer, intent(out) :: work(lwork) ! work array
     integer, intent(inout) :: info1
     integer, intent(in) :: stop_coarsening2 ! no. levels in the
     ! multilevel grid (default
     ! 10)

     type (nd_multigrid), intent(inout) :: grid ! the multilevel of
     ! graphs (matrices)

     integer :: i, j, k, inv1, inv2, ins
     integer :: mp
     integer :: mglevel_cur ! current level
     integer :: err, print_level ! printing
     logical :: lerr

     info1 = 0
     ! Set up printing
     if (options%print_level<0) print_level = 0
     ! The default is options%print_level = 0
     if (options%print_level==0) print_level = 1
     if (options%print_level==1) print_level = 2
     if (options%print_level>1) print_level = 3
     mp = options%unit_diagnostics
     if (mp<0) print_level = 0
     ! Set error optionss
     lerr = options%unit_error >= 0 .and. print_level > 0
     err = options%unit_error


     if (print_level>1) then
       write (mp,'(a)') 'Start multilevel_partition:'
     end if

     ! construct the multigrid at this level

     if ( .not. allocated(grid%graph)) allocate (grid%graph)

     call nd_matrix_construct(grid%graph,a_n,a_n,a_ne,info1)
     if (info1<0) then
       if (lerr) call nd_print_message(info1,err, &
         ' multilevel_partition')
       return
     end if

     grid%graph%ptr(1:a_n) = a_ptr(1:a_n)
     grid%graph%ptr(a_n+1) = a_ne + 1
     grid%graph%col(1:a_ne) = a_row(1:a_ne)

     do i = 1, a_n - 1
       do j = a_ptr(i), a_ptr(i+1) - 1
         k = a_row(j)
         grid%graph%val(j) = a_weight(i)*a_weight(k)
       end do
     end do
     do j = a_ptr(a_n), a_ne
       k = a_row(j)
       grid%graph%val(j) = a_weight(a_n)*a_weight(k)
     end do

     grid%size = a_n
     grid%level = 1
     ! nullify (grid%p)

     call nd_assoc(grid%where,a_n,info1)
     if (info1<0) then
       if (lerr) call nd_print_message(info1,err, &
         ' multilevel_partition')
       return
     end if

     call nd_assoc(grid%row_wgt,a_n,info1)
     if (info1<0) then
       if (lerr) call nd_print_message(info1,err, &
         ' multilevel_partition')
       return
     end if

     ! Initialise row weights
     grid%row_wgt(1:a_n) = a_weight(1:a_n)

     ! initialise mglevel_cur to the maximum number of levels
     ! allowed for this bisection
     mglevel_cur = stop_coarsening2
     call multilevel(grid,options,sumweight,mglevel_cur,mp,print_level, &
       lwork,work,info1)

     if (info1/=0) then
       if (lerr) call nd_print_message(info1,err, &
         ' multilevel_partition')
       return
     end if

     inv1 = 1
     inv2 = grid%part_div(1) + 1
     ins = grid%part_div(1) + grid%part_div(2) + 1

     a_weight_1 = 0
     a_weight_2 = 0
     a_weight_sep = 0
     do i = 1, a_n
       select case (grid%where(i))
       case (ND_PART1_FLAG)
         partition(inv1) = i
         inv1 = inv1 + 1
         a_weight_1 = a_weight_1 + a_weight(i)
       case (ND_PART2_FLAG)
         partition(inv2) = i
         inv2 = inv2 + 1
         a_weight_2 = a_weight_2 + a_weight(i)
       case default
         partition(ins) = i
         ins = ins + 1
         a_weight_sep = a_weight_sep + a_weight(i)
       end select
     end do

     a_n1 = grid%part_div(1)
     a_n2 = grid%part_div(2)

     if (.false.) then
       write (*,'(a)') ' '
       write (*,'(a)') 'Multilevel partition found'
       write (*,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, ', a_n2=', a_n2, &
         ', a_n_sep=', a_n - a_n1 - a_n2
       write (*,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', a_weight_1, &
         ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
         sumweight - a_weight_1 - a_weight_2
     end if

     ! deallocate the finest level
     ! call multigrid_deallocate_first(a_n,a_n,grid,info1)
     if (info1/=0) then
       if (lerr) call nd_print_message(info1,err, &
         ' multilevel_partition')
       return
     end if

     ! deallocate (matrix%ptr,matrix%col,matrix%val,stat=st)
     ! if (st/=0) info1 = ND_ERR_MEMORY_DEALLOC
     if (info1<0) then
       if (lerr) call nd_print_message(info1,err, &
         ' multilevel_partition')
       return
     end if

     if (print_level>2) then
       write (mp,'(a)') 'multilevel_partition: successful completion'
     end if

   end subroutine multilevel_partition

   ! ********************************************************

   ! main subroutine for computing multilevel structure.
   ! Offers heavy-edge collapsing and maximal independent vertex
   ! set for coarsening. We will need to test out to see
   ! which is better.

   recursive subroutine multilevel(grid,options,sumweight,mglevel_cur,mp, &
       print_level,lwork,work,info)

     real(wp), parameter :: half = 0.5_wp
     real(wp), parameter :: one = 1.0_wp

     ! Arguments
     type (nd_multigrid), intent(inout), TARGET :: grid ! this level
     ! of matrix (grid)
     type (nd_options), intent(in) :: options
     integer, intent(in) :: sumweight ! sum of weights (unchanged between
     ! coarse and fine grid
     integer, intent(inout) :: mglevel_cur ! current grid level
     integer, intent(in) :: mp, print_level ! diagnostic printing
     integer, intent(in) :: lwork ! length of work array
     ! (.ge.9*grid%graph%n +sumweight)
     integer, intent(out) :: work(lwork) ! work array
     integer, intent(inout) :: info ! Error flag

     ! Local variables
     type (nd_multigrid), pointer :: cgrid ! the coarse level grid
     integer :: cnvtx ! number of vertices (rows) in the coarse
     ! matrix
     type (nd_matrix), pointer :: p ! the coarse grid prolongator
     type (nd_matrix), pointer :: r ! the coarse grid restrictor (= p')

     integer, dimension(:), pointer :: fwhere ! partition on fine grid
     integer, dimension(:), pointer :: cwhere ! partition on coarse grid
     type (nd_matrix), pointer :: cgraph ! the coarse graph
     type (nd_matrix), pointer :: graph ! the fine graph
     integer, dimension(:), pointer :: row_wgt ! fine
     ! graph vertex weights
     integer, dimension(:), pointer :: crow_wgt ! coarse
     ! graph vertex weights
     real(wp) :: grid_rdc_fac_min ! min grid reduction
     ! factor
     real(wp) :: grid_rdc_fac_max ! max grid reduction
     ! factor
     real(wp) :: one1
     integer :: stop_coarsening1 ! optionss when to stop coarsening
     integer :: partition_ptr, part_ptr, work_ptr, a_ne, ref_options, &
       clwork
     integer :: i, j, k, l, a_weight_1, a_weight_2, a_weight_sep, &
       ref_method, lwk
     integer :: a_n1_new, a_n2_new, a_weight_1_new, a_weight_2_new, &
       a_weight_sep_new
     logical :: imbal
     real(wp) :: tau, ratio, tau_best
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!
     info = 0
     one1 = 1.0

     stop_coarsening1 = max(2,options%stop_coarsening1)
     if (print_level.ge.2) call level_print(mp,'size of grid on level ', &
       grid%level,' is ',real(grid%size,wp))

     grid_rdc_fac_min = max(0.01_wp,options%min_reduction)
     ! max grid reduction factor must be at least half and at most one
     grid_rdc_fac_max = max(half,options%max_reduction)
     grid_rdc_fac_max = min(one,grid_rdc_fac_max)

     ! Test to see if this is either the last level or
     ! if the matrix size too small
     if (grid%level>=mglevel_cur .or. grid%size<=stop_coarsening1) then
       if (print_level.ge.2) call level_print(mp,'end of level ',grid%level)

       ! coarsest level in multilevel so compute separator
       a_ne = grid%graph%ptr(grid%graph%n+1) - 1
       call nd_coarse_partition(grid%graph%n,a_ne,grid%graph%ptr, &
         grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
         grid%part_div(2),grid%where,lwork,work,options,info)
       return
     end if

     ! Coarsest level not yet reached so carry on coarsening
     if (options%matching==1) then
       lwk = grid%size
       call coarsen_hec(grid,lwk,work(1:lwk),info)
     else
       if (options%matching>1) then
         lwk = 3*grid%size
         call coarsen_best(grid,lwk,work(1:lwk),info)
       else
         lwk = 2*grid%size
         call coarsen_cn(grid,lwk,work(1:lwk))
       end if
     end if
     if (info<0) return

     cgrid => grid%coarse
     cnvtx = cgrid%size
     ! allocate coarse grid quantities
     call nd_assoc(cgrid%where,cnvtx,info)
     if (info/=0) then
       return
     end if

     call nd_assoc(cgrid%row_wgt,cnvtx,info)
     if (info/=0) then
       return
     end if

     ! see if the grid reduction is achieved, if not, set the allowed
     ! maximum level to current level and partition this level
     ! deallocate the coarse grid quantities that haves been allocated so
     ! far
     if (real(cgrid%size)/real(grid%size)>grid_rdc_fac_max .or. &
         real(cgrid%size)/real(grid%size)<grid_rdc_fac_min .or. &
         cgrid%size<4) then

       if (print_level.ge.2) then
         ! if (.true.) then
         write (mp,'(a,i10,a,f12.4,i4)') 'at level ', grid%level, &
           ' further coarsening gives reduction factor', &
           cgrid%size/real(grid%size)
         write (mp,'(a,i10)') 'current size = ', grid%size
       end if

       ! set current grid level and recurse
       mglevel_cur = grid%level

       call multilevel(grid,options,sumweight,mglevel_cur,mp,print_level, &
         lwork,work,info)
       if (info<0) return

       return
     end if

     ! restriction ================

     ! form the coarse grid graph and matrix
     ! cmatrix = P^T*matrix = R*matrix
     p => cgrid%p
     r => cgrid%r
     graph => grid%graph
     cgraph => cgrid%graph

     ! get the coarse matrix
     lwk = 3*grid%size
     call galerkin_graph(graph,p,r,cgraph,info,lwk,work(1:lwk))
     if (info<0) return

     ! check if matrix is full
     if (real(cgrid%graph%ptr(cgrid%graph%n+1)-1)/real(cgrid%graph%n)>=real &
         (cgrid%graph%n-1)) then
       if (print_level.ge.2) then
         write (mp,'(a,i10,a)') 'at level ', grid%level, &
           ' further coarsening gives full matrix'
       end if

       ! set current grid level and recurse
       mglevel_cur = grid%level - 1
       call multilevel(grid,options,sumweight,mglevel_cur,mp,print_level, &
         lwork,work,info)
       if (info<0) return

       return
     end if

     ! row weight cw = R*w
     row_wgt => grid%row_wgt(1:grid%size)
     crow_wgt => cgrid%row_wgt(1:cgrid%size)
     call nd_matrix_multiply_vec(r,row_wgt,crow_wgt)
     clwork = 9*cgrid%graph%n + sumweight
     call multilevel(cgrid,options,sumweight,mglevel_cur,mp,print_level, &
       clwork,work(1:clwork),info)


     ! check if partition is returned
     if (cgrid%part_div(1)==0 .or. cgrid%part_div(2)==0) then
       ! Unlikely to be called because 99.999% of cases caught in full
       ! matrix check above. Follows same procedure as when full matrix found
       if (print_level.ge.2) then
         write (mp,'(a,i10,a)') 'at level ', grid%level, &
           ' no partition found'
       end if

       ! set current grid level and recurse
       mglevel_cur = grid%level - 1
       call multilevel(grid,options,sumweight,mglevel_cur,mp,print_level, &
         lwork,work,info)
       if (info<0) return

       return
     end if

     ! prolongation ================

     ! injection of the order from coarse grid to the
     ! fine grid, since cwhere(i) is the index of the
     ! i-th vertex in the new ordering, the order
     ! of this vertex should be where(i)
     ! grid%where = P*order_on_coarse_grid
     ! here P is a special matrix with only one non-zero entry per row
     fwhere => grid%where(1:grid%size)
     cwhere => cgrid%where(1:cgrid%size)
     grid%part_div(1:2) = 0
     call nd_matrix_multiply_vec(p,cwhere,fwhere)

     do i = 1, grid%size
       if (fwhere(i)==ND_PART1_FLAG) then
         grid%part_div(1) = grid%part_div(1) + 1
       else
         if (fwhere(i)==ND_PART2_FLAG) then
           grid%part_div(2) = grid%part_div(2) + 1
         end if
       end if
     end do
     a_weight_1 = 0
     a_weight_2 = 0
     a_weight_sep = 0

     ! Set partition
     partition_ptr = 0
     work_ptr = partition_ptr + grid%graph%n
     i = 1
     j = grid%part_div(1) + 1
     k = grid%part_div(1) + grid%part_div(2) + 1
     do l = 1, grid%size
       select case (grid%where(l))
       case (ND_PART1_FLAG)
         work(partition_ptr+i) = l
         a_weight_1 = a_weight_1 + grid%row_wgt(l)
         i = i + 1
       case (ND_PART2_FLAG)
         work(partition_ptr+j) = l
         a_weight_2 = a_weight_2 + grid%row_wgt(l)
         j = j + 1
       case (ND_SEP_FLAG)
         work(partition_ptr+k) = l
         a_weight_sep = a_weight_sep + grid%row_wgt(l)
         k = k + 1
       end select
     end do
     a_ne = grid%graph%ptr(grid%graph%n+1) - 1

     if (a_weight_sep>0) then
       ! Do not refine if separable graph

       if (options%refinement>6) then
         ref_options = 3
       else
         if (options%refinement<1) then
           ref_options = 1
         else
           ref_options = options%refinement
         end if
       end if

       select case (ref_options)
       case (1)
         ref_method = 1

       case (2)
         ref_method = 2

       case (3)
         if (real(max(a_weight_1,a_weight_2))/real(min(a_weight_1, &
             a_weight_2)+a_weight_sep)>max(real(1.0, &
             wp),options%balance)) then
           ref_method = 2
         else
           ref_method = 1
         end if

       case (4)
         ref_method = 0

       case (5)
         ref_method = 2

       case (6)
         if (real(max(a_weight_1,a_weight_2))/real(min(a_weight_1, &
             a_weight_2)+a_weight_sep)>max(real(1.0, &
             wp),options%balance)) then
           ref_method = 2
         else
           ref_method = 0
         end if
       end select

       select case (ref_method)
       case (0)
         call nd_refine_max_flow(grid%graph%n,a_ne,grid%graph%ptr, &
           grid%graph%col,grid%row_wgt,grid%part_div(1),grid%part_div(2), &
           a_weight_1,a_weight_2,a_weight_sep,work(partition_ptr+1: &
           partition_ptr+grid%graph%n),work(work_ptr+1:work_ptr+8),options)
       case (1)
         if (min(a_weight_1,a_weight_2)+a_weight_sep< &
             max(a_weight_1,a_weight_2)) then
           call nd_refine_block_trim(grid%graph%n,a_ne, &
             grid%graph%ptr,grid%graph%col,grid%row_wgt,sumweight, &
             grid%part_div(1),grid%part_div(2),a_weight_1,a_weight_2, &
             a_weight_sep,work(partition_ptr+1:partition_ptr+grid%graph%n), &
             work(work_ptr+1:work_ptr+5*grid%graph%n),options)
         else
           call nd_refine_trim(grid%graph%n,a_ne,grid%graph%ptr, &
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep, &
             work(partition_ptr+1:partition_ptr+grid%graph%n), &
             work(work_ptr+1:work_ptr+3*grid%graph%n),options)

         end if
       case (2)
         call nd_refine_edge(grid%graph%n,a_ne,grid%graph%ptr, &
           grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
           grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep, &
           work(partition_ptr+1:partition_ptr+grid%graph%n), &
           work(work_ptr+1:work_ptr+3*grid%graph%n),options)
       end select

       if (options%max_improve_cycles>0) then
         ratio = max(real(1.0,wp),options%balance)
         if (ratio>real(sumweight-2)) then
           imbal = .false.
         else
           imbal = .true.
         end if
         call cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight, &
           ratio,imbal,options%cost_function,tau_best)
         a_n1_new = grid%part_div(1)
         a_n2_new = grid%part_div(2)
         a_weight_1_new = a_weight_1
         a_weight_2_new = a_weight_2
         a_weight_sep_new = a_weight_sep
       end if

       part_ptr = work_ptr + 5*grid%graph%n
       work(part_ptr+1:part_ptr+grid%graph%n) = work(partition_ptr+1: &
         partition_ptr+grid%graph%n)

       k = options%max_improve_cycles
       do i = 1, k

         call expand_partition(grid%graph%n,a_ne,grid%graph%ptr, &
           grid%graph%col,grid%row_wgt,a_n1_new,a_n2_new,a_weight_1_new, &
           a_weight_2_new,a_weight_sep_new,work(part_ptr+1:part_ptr+grid% &
           graph%n),work(work_ptr+1:work_ptr+5*grid%graph%n))


         ! call
         ! check_partition1(a_n,a_ne,a_ptr,a_row,a_n1_new,a_n2_new,work(par
         ! t_ptr+1:part_ptr+a_n))

         select case (ref_options)

         case (3)
           if (real(max(a_weight_1_new,a_weight_2_new))/real(min( &
               a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>max(real( &
               1.0,wp),options%balance)) then
             ref_method = 2
           else
             ref_method = 1
           end if

         case (6)
           if (real(max(a_weight_1_new,a_weight_2_new))/real(min( &
               a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>max(real( &
               1.0,wp),options%balance)) then
             ref_method = 2
           else
             ref_method = 0
           end if
         end select


         select case (ref_method)

         case (0)
           call nd_refine_max_flow(grid%graph%n,a_ne,grid%graph%ptr, &
             grid%graph%col,grid%row_wgt,a_n1_new,a_n2_new,a_weight_1_new, &
             a_weight_2_new,a_weight_sep_new,work(part_ptr+1:part_ptr+grid% &
             graph%n),work(work_ptr+1:work_ptr+8),options)

         case (1)
           if (min(a_weight_1,a_weight_2)+a_weight_sep< &
               max(a_weight_1,a_weight_2)) then
             call nd_refine_block_trim(grid%graph%n,a_ne, &
               grid%graph%ptr,grid%graph%col,grid%row_wgt,sumweight, &
               a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
               a_weight_sep_new,work(part_ptr+1:part_ptr+grid%graph%n), &
               work(work_ptr+1:work_ptr+5*grid%graph%n),options)
           else
             call nd_refine_trim(grid%graph%n,a_ne,grid%graph%ptr, &
               grid%graph%col,grid%row_wgt,sumweight,a_n1_new,a_n2_new, &
               a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
               work(part_ptr+1:part_ptr+grid%graph%n), &
               work(work_ptr+1:work_ptr+3*grid%graph%n),options)
           end if


         case (2)
           call nd_refine_edge(grid%graph%n,a_ne,grid%graph%ptr, &
             grid%graph%col,grid%row_wgt,sumweight,a_n1_new,a_n2_new, &
             a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
             work(part_ptr+1:part_ptr+grid%graph%n), &
             work(work_ptr+1:work_ptr+3*grid%graph%n),options)

         end select


         call cost_function(a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
           sumweight,ratio,imbal,options%cost_function,tau)
         if (tau<tau_best) then
           tau_best = tau
           work(partition_ptr+1:partition_ptr+grid%graph%n) &
             = work(part_ptr+1:part_ptr+grid%graph%n)
           grid%part_div(1) = a_n1_new
           grid%part_div(2) = a_n2_new
           a_weight_1 = a_weight_1_new
           a_weight_2 = a_weight_2_new
           a_weight_sep = a_weight_sep_new
         else
           exit
         end if
       end do




       ! if (grid%level .le.2) then
       call nd_refine_fm(grid%graph%n,a_ne,grid%graph%ptr, &
         grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
         grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep, &
         work(partition_ptr+1:partition_ptr+grid%graph%n), &
         work(work_ptr+1:work_ptr+8*grid%graph%n+sumweight),options)
       ! end if

     end if

     do i = 1, grid%part_div(1)
       j = work(partition_ptr+i)
       grid%where(j) = ND_PART1_FLAG
     end do
     do i = grid%part_div(1) + 1, grid%part_div(1) + grid%part_div(2)
       j = work(partition_ptr+i)
       grid%where(j) = ND_PART2_FLAG
     end do
     do i = grid%part_div(1) + grid%part_div(2) + 1, grid%graph%n
       j = work(partition_ptr+i)
       grid%where(j) = ND_SEP_FLAG
     end do

     if (info<0) return

     if (print_level==3) call level_print(mp,' after post smoothing ', &
       grid%level)

     ! deallocate the previous level
     ! call multigrid_deallocate(cgrid,info)

   end subroutine multilevel

   ! ***************************************************************
   ! ---------------------------------------------------
   ! nd_partition matrix
   ! ---------------------------------------------------
   ! Partition the matrix and if one (or more) of the generated submatrices
   ! is
   ! small enough, apply halo amd

   subroutine nd_coarse_partition(a_n,a_ne,a_ptr,a_row,a_weight, &
       sumweight,a_n1,a_n2,where1,lwork,work,options,info)

     integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
     integer, intent(in) :: a_ne ! no. nonzeros of subproblem
     integer, intent(inout) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start. This is then
     ! used to hold positions for submatrices after partitioning
     integer, intent(inout) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.This is then used to hold row indices for
     ! submatrices after partitioning
     integer, intent(inout) :: a_weight(a_n) ! On input a_weight(i)
     ! contains
     ! the weight of column i. This is then used to hold the weights for
     ! the submatrices after partitioning.
     integer, intent(in) :: sumweight ! Sum entries in a_weight.
     ! Unchanged.
     integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
     integer, intent(out) :: where1(a_n) ! Computed partition
     integer, intent(in) :: lwork ! .ge. 9*a_n+sumweight
     integer, intent(out) :: work(lwork)
     type (nd_options), intent(in) :: options
     integer, intent(inout) :: info
     ! real(wp), optional, intent(out) :: real_work(a_n)

     ! ---------------------------------------------
     ! Local variables
     integer :: unit_diagnostics ! unit on which to print diagnostics
     logical :: printi, printd, use_multilevel
     integer :: partition_ptr ! pointer into work array
     integer :: work_ptr ! pointer into work array
     integer :: partition_method
     integer :: st
     integer :: a_weight_1, a_weight_2, a_weight_sep, ref_method, &
       ref_options
     integer, allocatable :: work1(:)
     real(wp) :: dummy, dummy1
     type (nd_multigrid) :: gridtemp

     ! ---------------------------------------------
     ! Printing levels
     unit_diagnostics = options%unit_diagnostics
     printi = (options%print_level==1 .and. unit_diagnostics>=0)
     printd = (options%print_level>=2 .and. unit_diagnostics>=0)

     if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a)') 'Start finding a coarse partition'
       write (unit_diagnostics,'(a,i10,a,i10)') 'a_n=', a_n, ', a_ne=', &
         a_ne
     end if

     ! Find the partition
     if (options%coarse_partition_method<=1) then
       partition_method = 1
     else
       partition_method = 2
     end if

     partition_ptr = 0 ! length a_n
     work_ptr = partition_ptr + a_n ! max length needed 9*a_n+a_ne

     allocate (work1(a_n),stat=st)
     if (st/=0) then
       info = ND_ERR_MEMORY_ALLOC
       return
     end if


     select case (partition_method)
     case (1)
       ! Ashcraft method
       use_multilevel = .false.
       call nd_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
         a_n2,a_weight_1,a_weight_2,a_weight_sep, &
         work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
         options,info,dummy,dummy1,use_multilevel,gridtemp)

       if (printi .or. printd) then
         write (unit_diagnostics,'(a)') ' '
         write (unit_diagnostics,'(a)') 'Initial half-level set partition'
         write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
           ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
         write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
           a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
           sumweight - a_weight_1 - a_weight_2
       end if

     case (2)
       ! Level set method
       use_multilevel = .false.
       call nd_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
         a_n2,a_weight_1,a_weight_2,a_weight_sep, &
         work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
         options,info,dummy,dummy1,use_multilevel,gridtemp)


       if (printi .or. printd) then
         write (unit_diagnostics,'(a)') ' '
         write (unit_diagnostics,'(a)') 'Initial level set partition'
         write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
           ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
         write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
           a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
           sumweight - a_weight_1 - a_weight_2
       end if

     end select



     if (a_n1/=0 .and. a_n2/=0 .and. a_n>=3) then
       if (a_n1+a_n2<a_n) then
         ! Refine the partition
         if (options%refinement>6) then
           ref_options = 3
         else
           if (options%refinement<1) then
             ref_options = 1
           else
             ref_options = options%refinement
           end if
         end if

         select case (ref_options)
         case (1)
           ref_method = 1

         case (2)
           ref_method = 2

         case (3)
           if (real(max(a_weight_1,a_weight_2))/real(min(a_weight_1, &
               a_weight_2)+a_weight_sep)>max(real(1.0, &
               wp),options%balance)) then
             ref_method = 2
           else
             ref_method = 1
           end if

         case (4)
           ref_method = 0

         case (5)
           ref_method = 2

         case (6)
           if (real(max(a_weight_1,a_weight_2))/real(min(a_weight_1, &
               a_weight_2)+a_weight_sep)>max(real(1.0, &
               wp),options%balance)) then
             ref_method = 2
           else
             ref_method = 0
           end if
         end select

         select case (ref_method)

         case (0)
           call nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1, &
             a_n2,a_weight_1,a_weight_2,a_weight_sep, &
             work1(partition_ptr+1:partition_ptr+a_n),work(1:8),options)

         case (1)
           if (min(a_weight_1,a_weight_2)+a_weight_sep< &
               max(a_weight_1,a_weight_2)) then
             call nd_refine_block_trim(a_n,a_ne,a_ptr,a_row, &
               a_weight,sumweight,a_n1,a_n2,a_weight_1,a_weight_2, &
               a_weight_sep,work1(partition_ptr+1:partition_ptr+a_n), &
               work(1:5*a_n),options)
           else
             call nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
               sumweight,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
               work1(partition_ptr+1:partition_ptr+a_n),work(1:3*a_n), &
               options)
           end if

         case (2)
           call nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
             a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
             work1(partition_ptr+1:partition_ptr+a_n),work(1:3*a_n), &
             options)


         end select

         if (printi .or. printd) then
           write (unit_diagnostics,'(a)') ' '
           write (unit_diagnostics,'(a)') 'Trimmed partition found'
           write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
             ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
           write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
             a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
             sumweight - a_weight_1 - a_weight_2
         end if

         call nd_refine_fm(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
           a_n2,a_weight_1,a_weight_2,a_weight_sep, &
           work1(partition_ptr+1:partition_ptr+a_n), &
           work(1:8*a_n+sumweight),options)

       end if
     else
       go to 10
     end if


     call nd_convert_partition_flags(a_n,a_n1,a_n2, &
       work1(partition_ptr+1:partition_ptr+a_n),ND_PART1_FLAG, &
       ND_PART2_FLAG,ND_SEP_FLAG,where1(1:a_n))

     deallocate (work1,stat=st)
     if (st/=0) then
       info = ND_ERR_MEMORY_ALLOC
       return
     end if



     if (printi .or. printd .or. .false.) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a)') 'Initial coarse partition found'
       write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
         ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
       write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
         a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
         sumweight - a_weight_1 - a_weight_2
     end if
     go to 20

10      if (printi .or. printd) then
       write (unit_diagnostics,'(a)') ' '
       write (unit_diagnostics,'(a)') 'No partition found'
     end if

20      info = 0
     if (printi .or. printd) then
       call nd_print_message(info,unit_diagnostics, &
         'nd_coarse_partition')
     end if
     return

   end subroutine nd_coarse_partition

   ! *****************************************************************

   recursive subroutine mg_grid_destroy(grid,info)
     ! deallocate a grid structure
     type (nd_multigrid) :: grid
     integer :: info

     if (associated(grid%coarse)) then

       call mg_grid_destroy(grid%coarse,info)

       if (grid%level/=1) then

         call multigrid_deallocate(grid,info)

       else

         call multigrid_deallocate_first(grid,info)

       end if

     else

       if (grid%level/=1) then

         call multigrid_deallocate_last(grid,info)

       else

         call multigrid_deallocate_first(grid,info)

       end if

     end if

   end subroutine mg_grid_destroy


   ! *****************************************************************
   subroutine multigrid_deallocate(grid,info)
     ! deallocate a grid (at given level between last and first)
     type (nd_multigrid) :: grid
     integer :: info

     call nd_matrix_destruct(grid%graph,info)
     if (info/=0) then
       return
     end if

     call nd_matrix_destruct(grid%p,info)
     if (info/=0) then
       return
     end if



     call nd_matrix_destruct(grid%r,info)
     if (info/=0) then
       return
     end if

     if(associated(grid%coarse)) deallocate(grid%coarse)
     deallocate (grid%graph,grid%p,grid%r,grid%where,grid%row_wgt)
     nullify (grid%coarse)

   end subroutine multigrid_deallocate

   ! *****************************************************************
   subroutine multigrid_deallocate_last(grid,info)

     ! deallocate a grid (at the last level). In this case the matrix
     ! grid%graph
     ! has not been formed yet
     type (nd_multigrid) :: grid
     integer, intent(inout) :: info


     integer :: ierr

     call nd_matrix_destruct(grid%p,ierr)
     if (ierr/=0) then
       info = ND_ERR_MEMORY_DEALLOC
       return
     end if

     call nd_matrix_destruct(grid%r,ierr)
     if (ierr/=0) then
       info = ND_ERR_MEMORY_DEALLOC
       return
     end if
     if(associated(grid%coarse)) deallocate(grid%coarse)
     deallocate (grid%graph,grid%p,grid%r,grid%where,grid%row_wgt)
     nullify (grid%coarse)

   end subroutine multigrid_deallocate_last
   ! *****************************************************************
   subroutine multigrid_deallocate_first(grid,info)
     ! deallocate a grid (at the first level). In this case the matrix
     ! grid%p
     ! does not exist
     type (nd_multigrid) :: grid
     integer, intent(inout) :: info
     integer :: ierr

     if (allocated(grid%graph)) then
       call nd_matrix_destruct(grid%graph,ierr)
       if (ierr/=0) then
         info = ND_ERR_MEMORY_DEALLOC
         return
       end if
     end if

     deallocate (grid%where,grid%row_wgt,stat=ierr)
     if (ierr/=0) info = ND_ERR_MEMORY_DEALLOC

   end subroutine multigrid_deallocate_first

   ! ***************************************************************
   subroutine coarsen_hec(grid,lwork,work,info)
     ! coarsen the grid using heavy-edge collapsing and set up the
     ! coarse grid equation, the prolongator and restrictor

     type (nd_multigrid), intent(inout), TARGET :: grid
     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)
     integer, intent(inout) :: info


     if ( .not. associated(grid%coarse)) allocate (grid%coarse)

     grid%coarse%fine => grid

     ! find the prolongator
     call prolng_heavy_edge(grid,lwork,work,info)


     grid%coarse%level = grid%level + 1

   end subroutine coarsen_hec


   ! ***************************************************************
   subroutine coarsen_cn(grid,lwork,work)
     ! coarsen the grid using common neighbours collapsing and set up the
     ! coarse grid equation, the prolongator and restrictor

     type (nd_multigrid), intent(inout), TARGET :: grid

     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)

     if ( .not. associated(grid%coarse)) allocate (grid%coarse)

     grid%coarse%fine => grid

     ! find the prolongator

     call prolng_common_neigh(grid,lwork,work)

     grid%coarse%level = grid%level + 1

   end subroutine coarsen_cn

   ! ***************************************************************
   subroutine coarsen_best(grid,lwork,work,info)
     ! coarsen the grid using common neighbours collapsing and set up the
     ! coarse grid equation, the prolongator and restrictor

     integer, intent(inout) :: info
     type (nd_multigrid), intent(inout), TARGET :: grid

     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)

     if ( .not. associated(grid%coarse)) allocate (grid%coarse)

     grid%coarse%fine => grid

     ! find the prolongator

     call prolng_best(grid,lwork,work,info)

     grid%coarse%level = grid%level + 1

   end subroutine coarsen_best

   ! ***************************************************************
   subroutine nd_matrix_multiply_vec(matrix,x,y)
     ! subroutine nd_matrix_multiply_vec(matrix,x,y)

     ! y = matrix*x where x and y are integer vectors. Entries of
     ! matrix is assumed to be one. Dimension of y
     ! is checked and returned if it is smaller than the row dimension
     ! of x    !

     ! matrix: of the derived type nd_matrix, intent(in),
     ! the sparse matrix in compressed sparse row format
     type (nd_matrix), intent(in) :: matrix

     ! x: integer array of intent(in), a vector to be
     ! multiplied with the matrix
     integer, intent(in), dimension(*) :: x

     ! y: integer array of intent(out), the result of
     ! matrix*x or matrix^T*x
     integer, intent(out), dimension(*) :: y

     ! local ==========
     integer :: m, n, i, l1, l2

     m = matrix%m
     n = matrix%n

     do i = 1, m
       l1 = matrix%ptr(i)
       l2 = matrix%ptr(i+1) - 1
       y(i) = sum(x(matrix%col(l1:l2)))
     end do
   end subroutine nd_matrix_multiply_vec


   ! ***************************************************************
   subroutine nd_matrix_destruct(matrix,info,stat)
     ! subroutine nd_matrix_destruct(matrix,info):

     ! destruct the matrix object by deallocating all
     ! space occupied by
     ! matrix. including matrix%ptr, matrix%col and matrix%val.

     ! matrix: is of the derived type nd_matrix,
     ! with intent(inout). It
     ! the sparse matrix object to be destroyed.
     type (nd_matrix), intent(inout) :: matrix

     ! info: is an integer scaler of intent(out).
     ! = 0 if successful
     ! = nd_ERR_MEMORY_DEALLOC if memory deallocation failed
     integer, intent(out) :: info

     ! stat: is an integer scaler of intent(out). If supplied,
     ! on exit it holds the error tag for memory allocation
     integer, optional, intent(out) :: stat

     ! ===================== local variables =============
     ! ierr: error tag for deallocation
     integer :: ierr

     info = 0
     if (present(stat)) stat = 0
     deallocate (matrix%col,matrix%ptr,stat=ierr)
     if (present(stat)) stat = ierr
     if (ierr/=0) then
       info = ND_ERR_MEMORY_DEALLOC
       return
     end if

     deallocate (matrix%val,stat=ierr)
     if (present(stat)) stat = ierr
     if (ierr/=0) then
       info = ND_ERR_MEMORY_DEALLOC
       return
     end if

   end subroutine nd_matrix_destruct


   ! ***************************************************************

   subroutine nd_matrix_construct(p,m,n,ne,info)
     ! Construct data structure for storing sparse matrix
     ! Arrays in nd_matrix will only be (re)allocated if they are not
     ! long
     ! enough. On exit,
     ! size(p%val) <-  max(ne, size(p%val)
     ! size(p%col) <-  max(ne, size(p%col)
     ! size(p%ptr) <-  max(m+1, size(p%ptr)
     type (nd_matrix), intent(inout) :: p ! matrix being formed using
     ! CSR
     integer, intent(in) :: m ! number of rows
     integer, intent(in) :: n ! number of columns
     integer, intent(in) :: ne ! number entries
     integer, intent(out) :: info

     info = 0

     p%m = m
     p%n = n
     p%ne = ne

     call nd_alloc(p%ptr,m+1,info)
     if (info<0) then
       return
     end if

     call nd_alloc(p%col,ne,info)
     if (info<0) then
       return
     end if

     call nd_alloc(p%val,ne,info)
     if (info<0) then
       return
     end if

   end subroutine nd_matrix_construct

   ! ***************************************************************

   subroutine nd_alloc(v,n,info)
     integer, intent(inout), allocatable :: v(:)
     integer, intent(in) :: n
     integer, intent(out) :: info

     integer :: st

     info = 0

     if (allocateD(v)) then
       if (SIZE(v)<n) then
         deallocate (v,stat=st)
         if (st<0) then
           info = ND_ERR_MEMORY_ALLOC
           return
         end if
       else
         return
       end if
     end if

     allocate (v(n),stat=st)
     if (st<0) then
       info = ND_ERR_MEMORY_DEALLOC
     end if


   end subroutine nd_alloc


   ! ********************************************************

   !
   ! If array has size at least sz, do nothing. Otherwise, create/resize array
   ! arr of size sz.
   !
   subroutine nd_assoc(array,sz,info)
      integer, allocatable, dimension(:), intent(inout) :: array
      integer, intent(in) :: sz
      integer, intent(out) :: info

      integer :: st

      info = 0

      if (allocated(array)) then
        if(size(array).ge.sz) return ! All is well, immediate return
        ! Otherwise deallocate
        deallocate (array)
      endif

      ! If we reach thsi point, arr is now deallocated: allocate to correct size
      allocate (array(sz),stat=st)
      if (st/=0) info = ND_ERR_MEMORY_ALLOC

   end subroutine nd_assoc

   ! ********************************************
   subroutine prolng_heavy_edge(grid,lwork,work,info)

     ! calculate the prolongator for heavy-edge collapsing:
     ! match the vertices of the heaviest edges

     integer, intent(inout) :: info
     ! input fine grid
     type (nd_multigrid), target, intent(inout) :: grid
     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)

     ! coarse grid based on the fine grid
     type (nd_multigrid), pointer :: cgrid

     ! the fine grid row connectivity graph
     type (nd_matrix), pointer :: graph

     ! the coarse grid prolongator
     type (nd_matrix), pointer :: p

     ! the coarse grid restrictor
     type (nd_matrix), pointer :: r

     ! the number of fine and coarse grid vertices
     integer :: nvtx, cnvtx

     ! working variables
     integer :: v, u, j, i, k
     integer :: nz

     ! whether a vertex is matched already
     integer, parameter :: unmatched = -1

     ! matching status of each vertex
     integer :: ptr_match

     ! maximum weight and index of edges connected to the current vertex
     integer :: maxwgt
     integer :: maxind

     ! allocate the prolongation matrix pointers
     cgrid => grid%coarse
     graph => grid%graph

     ! allocate the graph and matrix pointer and the mincut pointer
     ! so that everything is defined
     if ( .not. allocated(cgrid%graph)) allocate (cgrid%graph)

     nvtx = graph%n

     ! prolongator start here ================================

     ! initialise the matching status and randomly permute the vertex order
     ptr_match = 0

     work(ptr_match+1:ptr_match+nvtx) = unmatched

     ! loop over each vertex and match along the heaviest edge
     cnvtx = 0
     nz = 0
     do i = 1, nvtx
       v = i
       ! If already matched, next vertex please
       if (work(ptr_match+v)/=unmatched) cycle
       maxwgt = -huge(0)
       ! in the case no match is found then match itself
       maxind = v
       ! Loop over entries in row v
       do j = graph%ptr(v), graph%ptr(v+1) - 1
         ! u is col index of en entry in row v (so u is neighbor of v)
         u = graph%col(j)
         ! heavy edge matching
         ! if u is unmatched and value of the entry in col. u is greater
         ! than maxwgt, select u as the matching.
         if (work(ptr_match+u)==unmatched .and. maxwgt<abs(graph%val(j))) &
             then
           maxwgt = abs(graph%val(j))
           maxind = u
         end if
       end do
       ! notE: maxind .ge. v
       ! the neighbor with heaviest weight
       work(ptr_match+v) = maxind
       ! mark maxind as having been matched
       work(ptr_match+maxind) = v
       ! increase number of vertices in coarse graph by 1
       cnvtx = cnvtx + 1
       ! construct the prolongation matrix: find vertex v and maxind is
       ! linked
       ! with the coarse grid vertex cnvtx
       nz = nz + 1
       if (maxind/=v) then
         nz = nz + 1
       end if
     end do

     ! storage allocation for col. indices and values of prolongation
     ! matrix P (nvtx * cnvtx)
     if ( .not. allocated(cgrid%p)) then
       allocate (cgrid%p)
       p => cgrid%p
       call nd_matrix_construct(p,nvtx,cnvtx,nz,info)
     else
       p => cgrid%p
       call nd_matrix_construct(p,nvtx,cnvtx,nz,info)
     end if


     ! storage allocation for col. indices and values of restiction
     ! matrix R (cnvtx * nvtx)
     if ( .not. allocated(cgrid%r)) then
       allocate (cgrid%r)
       r => cgrid%r
       call nd_matrix_construct(r,cnvtx,nvtx,nz,info)
     else
       r => cgrid%r
       call nd_matrix_construct(r,cnvtx,nvtx,nz,info)
     end if

     r%val(1:nz) = 1

     ! store restriction matrix
     r%ptr(cnvtx+1) = nz + 1

     j = 1
     k = 1
     do i = 1, nvtx
       if (work(ptr_match+i)==i) then
         r%ptr(k) = j
         r%col(j) = i
         j = j + 1
         k = k + 1
       else
         if (work(ptr_match+i)>i) then
           r%ptr(k) = j
           r%col(j) = i
           r%col(j+1) = work(ptr_match+i)
           j = j + 2
           k = k + 1
         end if
       end if
     end do

     ! store prolongation matrix

     p%ptr(1) = 1
     do i = 1, nvtx
       p%ptr(i+1) = p%ptr(i) + 1
     end do

     p%val(1:nz) = 1

     j = 1
     do i = 1, nvtx
       k = work(ptr_match+i)
       if (k==i) then
         p%col(p%ptr(i)) = j
         j = j + 1
       else
         if (k>i) then
           p%col(p%ptr(i)) = j
           p%col(p%ptr(k)) = j
           j = j + 1
         end if
       end if
     end do

     ! size of coarse grid
     cgrid%size = cnvtx

   end subroutine prolng_heavy_edge

   ! *******************************************************************

   subroutine prolng_common_neigh(grid,lwork,work)

     ! calculate the prolongator:
     ! match the vertices of with most neighbours in common

     ! input fine grid
     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)
     type (nd_multigrid), target, intent(inout) :: grid

     ! coarse grid based on the fine grid
     type (nd_multigrid), pointer :: cgrid

     ! the fine grid row connectivity graph
     type (nd_matrix), pointer :: graph

     ! the coarse grid prolongator
     type (nd_matrix), pointer :: p

     ! the fine grid restrictor
     type (nd_matrix), pointer :: r

     ! the number of fine and coarse grid vertices
     integer :: nvtx, cnvtx

     ! working variables
     integer :: v, u, w, j, i, k

     ! whether a vertex is matched already
     integer, parameter :: unmatched = -1

     ! matching status of each vertex
     integer :: ptr_match
     ! flag array to flag up neighbours of a  node
     integer :: ptr_flag

     ! maximum no. neighbours and index of edges connected to the current
     ! vertex
     integer :: max_neigh, maxind, num

     integer :: info
     integer :: nz

     ! allocate the prolongation matrix pointers
     cgrid => grid%coarse
     graph => grid%graph

     ! allocate the graph and matrix pointer and the mincut pointer
     ! so that everything is defined
     if ( .not. allocated(cgrid%graph)) allocate (cgrid%graph)

     nvtx = graph%n

     ! prolongator start here ================================

     ! initialise the matching status
     ptr_match = 0
     ptr_flag = ptr_match + nvtx

     work(ptr_match+1:ptr_match+nvtx) = unmatched

     work(ptr_flag+1:ptr_flag+nvtx) = 0

     ! loop over each vertex and match based on number of neighbours in
     ! common
     cnvtx = 0
     nz = 0
     do i = 1, nvtx
       v = i
       ! If already matched, next vertex please
       if (work(ptr_match+v)/=unmatched) cycle
       ! access the col. indices of row v

       ! in the case no match is found then match itself
       maxind = v
       ! Loop over entries in row v and set flag for each entry
       work(ptr_flag+v) = i
       do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
         ! u is col index of en entry in row v (so u is neighbor of v)
         u = grid%graph%col(j)
         work(ptr_flag+u) = i
       end do
       ! For each unmatched neighbour of v, count the number of
       ! neighbours it has in common with v
       max_neigh = 0
       do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
         u = grid%graph%col(j)
         ! cycle is u is already matched
         if (work(ptr_match+u)/=unmatched) cycle
         num = 0
         do k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
           w = grid%graph%col(k)
           if (work(ptr_flag+w)==i) num = num + 1
         end do
         if (num>max_neigh) then
           max_neigh = num
           maxind = u
         end if
       end do

       ! the neighbor with largest number of neighbours in common with v
       work(ptr_match+v) = maxind
       ! mark maxind as having been matched
       work(ptr_match+maxind) = v
       ! increase number of vertices in coarse graph by 1
       cnvtx = cnvtx + 1
       ! construct the prolongation matrix: find vertex v and maxind is
       ! linked
       ! with the coarse grid vertex cnvtx
       nz = nz + 1
       if (maxind/=v) then
         nz = nz + 1
       end if
     end do

     ! storage allocation for col. indices and values of prolongation
     ! matrix P (order nvtx * cnvtx)
     if ( .not. allocated(cgrid%p)) then
       allocate (cgrid%p)
       p => cgrid%p
       call nd_matrix_construct(p,nvtx,cnvtx,nz,info)
     else
       p => cgrid%p
       call nd_matrix_construct(p,nvtx,cnvtx,nz,info)
     end if
     p%val(1:nz) = 0

     ! storage allocation for col. indices and values of restiction
     ! matrix R (cnvtx * nvtx)
     if ( .not. allocated(cgrid%r)) then
       allocate (cgrid%r)
       r => cgrid%r
       call nd_matrix_construct(r,cnvtx,nvtx,nz,info)
     else
       r => cgrid%r
       call nd_matrix_construct(r,cnvtx,nvtx,nz,info)
     end if

     r%val(1:nz) = 1

     ! store restriction matrix
     r%ptr(cnvtx+1) = nz + 1

     j = 1
     k = 1
     do i = 1, nvtx
       if (work(ptr_match+i)==i) then
         r%ptr(k) = j
         r%col(j) = i
         j = j + 1
         k = k + 1
       else
         if (work(ptr_match+i)>i) then
           r%ptr(k) = j
           r%col(j) = i
           r%col(j+1) = work(ptr_match+i)
           j = j + 2
           k = k + 1
         end if
       end if
     end do


     ! store prolongation matrix

     p%ptr(1) = 1
     do i = 1, nvtx
       p%ptr(i+1) = p%ptr(i) + 1
     end do

     p%val(1:nz) = 1

     j = 1
     do i = 1, nvtx
       k = work(ptr_match+i)
       if (k==i) then
         p%col(p%ptr(i)) = j
         j = j + 1
       else
         if (k>i) then
           p%col(p%ptr(i)) = j
           p%col(p%ptr(k)) = j
           j = j + 1
         end if
       end if
     end do


     ! size of coarse grid
     cgrid%size = cnvtx

   end subroutine prolng_common_neigh

   ! ********************************************
   subroutine prolng_best(grid,lwork,work,info)

     ! calculate the prolongator for heavy-edge collapsing:
     ! match the vertices of the heaviest edges

     integer, intent(inout) :: info
     ! input fine grid
     type (nd_multigrid), target, intent(inout) :: grid
     integer, intent(in) :: lwork
     integer, intent(out), TARGET :: work(lwork)

     ! coarse grid based on the fine grid
     type (nd_multigrid), pointer :: cgrid

     ! the fine grid row connectivity graph
     type (nd_matrix), pointer :: graph

     ! the coarse grid prolongator
     type (nd_matrix), pointer :: p

     ! the coarse grid restrictor
     type (nd_matrix), pointer :: r

     ! the number of fine and coarse grid vertices
     integer :: nvtx, cnvtx, cnvtx1

     ! working variables
     integer :: v, u, j, i, k
     integer :: nz

     ! whether a vertex is matched already
     integer, parameter :: unmatched = -1

     ! matching status of each vertex
     integer :: ptr_match, ptr_match1, ptr_flag, max_neigh, num, w

     ! maximum weight and index of edges connected to the current vertex
     integer :: maxwgt
     integer :: maxind

     integer, pointer, dimension(:) :: matching

     ! allocate the prolongation matrix pointers
     cgrid => grid%coarse
     graph => grid%graph

     ! allocate the graph and matrix pointer and the mincut pointer
     ! so that everything is defined
     if ( .not. allocated(cgrid%graph)) allocate (cgrid%graph)

     nvtx = graph%n

     ptr_match = 0
     ptr_match1 = ptr_match + nvtx

     ! -----------------------------------------------------------------
     ! Find heavy-edge matching

     ! initialise the matching status and randomly permute the vertex order

     work(ptr_match+1:ptr_match+nvtx) = unmatched

     ! loop over each vertex and match along the heaviest edge
     cnvtx = 0
     do i = 1, nvtx
       v = i
       ! If already matched, next vertex please
       if (work(ptr_match+v)/=unmatched) cycle
       maxwgt = -huge(0)
       ! in the case no match is found then match itself
       maxind = v
       ! Loop over entries in row v
       do j = graph%ptr(v), graph%ptr(v+1) - 1
         ! u is col index of en entry in row v (so u is neighbor of v)
         u = graph%col(j)
         ! heavy edge matching
         ! if u is unmatched and value of the entry in col. u is greater
         ! than maxwgt, select u as the matching.
         if (work(ptr_match+u)==unmatched .and. maxwgt<abs(graph%val(j))) &
             then
           maxwgt = abs(graph%val(j))
           maxind = u
         end if
       end do
       ! notE: maxind .ge. v
       ! the neighbor with heaviest weight
       work(ptr_match+v) = maxind
       ! mark maxind as having been matched
       work(ptr_match+maxind) = v
       ! increase number of vertices in coarse graph by 1
       cnvtx = cnvtx + 1
       ! construct the prolongation matrix: find vertex v and maxind is
       ! linked
       ! with the coarse grid vertex cnvtx
     end do
     nz = nvtx

     ! -----------------------------------------------------------------
     ! Find common neighbours matching

     ! initialise the matching status
     ptr_match1 = ptr_match + nvtx
     ptr_flag = ptr_match1 + nvtx

     work(ptr_match1+1:ptr_match1+nvtx) = unmatched

     work(ptr_flag+1:ptr_flag+nvtx) = 0

     ! loop over each vertex and match based on number of neighbours in
     ! common
     cnvtx1 = 0
     nz = 0
     do i = 1, nvtx
       v = i
       ! If already matched, next vertex please
       if (work(ptr_match1+v)/=unmatched) cycle
       ! access the col. indices of row v

       ! in the case no match is found then match itself
       maxind = v
       ! Loop over entries in row v and set flag for each entry
       work(ptr_flag+v) = i
       do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
         ! u is col index of en entry in row v (so u is neighbor of v)
         u = grid%graph%col(j)
         work(ptr_flag+u) = i
       end do
       ! For each unmatched neighbour of v, count the number of
       ! neighbours it has in common with v
       max_neigh = 0
       do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
         u = grid%graph%col(j)
         ! cycle is u is already matched
         if (work(ptr_match1+u)/=unmatched) cycle
         num = 0
         do k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
           w = grid%graph%col(k)
           if (work(ptr_flag+w)==i) num = num + 1
         end do
         if (num>max_neigh) then
           max_neigh = num
           maxind = u
         end if
       end do

       ! the neighbor with largest number of neighbours in common with v
       work(ptr_match1+v) = maxind
       ! mark maxind as having been matched
       work(ptr_match1+maxind) = v
       ! increase number of vertices in coarse graph by 1
       cnvtx1 = cnvtx1 + 1
       ! construct the prolongation matrix: find vertex v and maxind is
       ! linked
       ! with the coarse grid vertex cnvtx1
       nz = nz + 1
       if (maxind/=v) then
         nz = nz + 1
       end if
     end do

     ! --------------------------------------------------------------------
     ! -
     if (cnvtx<=cnvtx1) then
       ! use heavy-edge matching
       matching => work(ptr_match+1:ptr_match+nvtx)
     else
       ! use common neighbours matching
       matching => work(ptr_match1+1:ptr_match1+nvtx)
       cnvtx = cnvtx1
     end if


     ! storage allocation for col. indices and values of prolongation
     ! matrix P (nvtx * cnvtx)
     if ( .not. allocated(cgrid%p)) then
       allocate (cgrid%p)
       p => cgrid%p
       call nd_matrix_construct(p,nvtx,cnvtx,nz,info)
     else
       p => cgrid%p
       call nd_matrix_construct(p,nvtx,cnvtx,nz,info)
     end if


     ! storage allocation for col. indices and values of restiction
     ! matrix R (cnvtx * nvtx)
     if ( .not. allocated(cgrid%r)) then
       allocate (cgrid%r)
       r => cgrid%r
       call nd_matrix_construct(r,cnvtx,nvtx,nz,info)
     else
       r => cgrid%r
       call nd_matrix_construct(r,cnvtx,nvtx,nz,info)
     end if

     r%val(1:nz) = 1

     ! store restriction matrix
     r%ptr(cnvtx+1) = nz + 1

     j = 1
     k = 1
     do i = 1, nvtx
       if (matching(i)==i) then
         r%ptr(k) = j
         r%col(j) = i
         j = j + 1
         k = k + 1
       else
         if (matching(i)>i) then
           r%ptr(k) = j
           r%col(j) = i
           r%col(j+1) = matching(i)
           j = j + 2
           k = k + 1
         end if
       end if
     end do

     ! store prolongation matrix

     p%ptr(1) = 1
     do i = 1, nvtx
       p%ptr(i+1) = p%ptr(i) + 1
     end do

     p%val(1:nz) = 1

     j = 1
     do i = 1, nvtx
       k = matching(i)
       if (k==i) then
         p%col(p%ptr(i)) = j
         j = j + 1
       else
         if (k>i) then
           p%col(p%ptr(i)) = j
           p%col(p%ptr(k)) = j
           j = j + 1
         end if
       end if
     end do

     ! size of coarse grid
     cgrid%size = cnvtx
     nullify (matching)

   end subroutine prolng_best


   ! *******************************************************************
   subroutine level_print(mp,title1,level,title2,res)

     character (len=*), intent(in) :: title1
     integer, intent(in) :: mp, level
     real(wp), optional, intent(in) :: res
     character (len=*), optional, intent(in) :: title2
     integer :: char_len1, char_len2

     char_len1 = len_trim(title1)

     if (present(res) .and. present(title2)) then
       char_len2 = len_trim(title2)
       write (mp,'(a,i4,a,g14.3)') title1, level, title2, res
     else
       write (mp,'(a,i4)') title1, level
     end if

   end subroutine level_print



   ! *************************************************
   subroutine galerkin_graph(matrix,p,r,cmatrix,info,lwork,work)

     ! Given matrix on fine grid and a prolongation operator p,
     ! find the coarse matrix R*A*P

     ! matrix: fine grid matrix
     type (nd_matrix), intent(in) :: matrix
     ! p: prolongation operator
     type (nd_matrix), intent(in) :: p
     ! r: restriction operator
     type (nd_matrix), intent(in) :: r
     ! cmatrix: coarse grid matrix
     type (nd_matrix), intent(inout) :: cmatrix
     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)

     ! nvtx,cnvtx: size of fine and coarse grid
     integer :: nvtx, cnvtx
     integer :: nz

     integer, intent(inout) :: info

     ! call mc65_matrix_transpose(p,r,info65)
     ! if (info65<0) then
     ! info = info65
     ! return
     ! end if
     nvtx = matrix%n
     cnvtx = p%n

     ! get the size of the coarse matrix first
     call galerkin_graph_rap_size(nvtx,cnvtx,nz,p%ptr(nvtx+1)-1,p%col, &
       p%ptr,matrix%ptr(nvtx+1)-1,matrix%col,matrix%ptr,r%ptr(cnvtx+1)-1, &
       r%col,r%ptr,lwork,work(1:lwork))
     if (info<0) return

     call nd_matrix_construct(cmatrix,cnvtx,cnvtx,nz,info)
     if (info<0) then
       return
     end if

     call galerkin_graph_rap(nvtx,cnvtx,p%ptr(nvtx+1)-1,p%val,p%col,p%ptr, &
       matrix%ptr(nvtx+1)-1,matrix%val,matrix%col,matrix%ptr, &
       r%ptr(cnvtx+1)-1,r%val,r%col,r%ptr,nz,cmatrix%val,cmatrix%col, &
       cmatrix%ptr,lwork,work(1:lwork))
     if (info<0) return

   end subroutine galerkin_graph

   ! *************************************************

   subroutine galerkin_graph_rap_size(nvtx,cnvtx,nz,nzp,pcol,pptr,nzaa, &
       acol,aptr,nzr,rcol,rptr,lwork,work)
     ! get the number of nonzeros in R*A*P
     ! nvtx: size of aa matrix
     ! cnvtx: size of ca matrix
     integer, intent(in) :: nvtx, cnvtx
     ! nz: number of nonzeros in R*A*P
     integer, intent(out) :: nz

     ! P: matrix
     integer, intent(in) :: nzp
     integer, intent(in), dimension(nzp) :: pcol
     integer, intent(in), dimension(nvtx+1) :: pptr
     ! aa: matrix
     integer, intent(in) :: nzaa
     integer, intent(in), dimension(nzaa) :: acol
     integer, intent(in), dimension(nvtx+1) :: aptr
     ! R: matrix
     integer, intent(in) :: nzr
     integer, intent(in), dimension(nzr) :: rcol
     integer, intent(in), dimension(cnvtx+1) :: rptr

     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)

     ! mask: masking array to see if an entry has been seen before
     integer :: ptr_mask
     ! i,j,k: loop index
     integer :: i, j, k
     ! nz: number of nonzeros so far in ca
     integer :: nz1
     ! various neighbors
     integer :: neigh, neighneigh

     ! col: column index of a row of r*matrix
     integer :: ptr_col

     ptr_mask = 0
     ptr_col = ptr_mask + nvtx
     work(ptr_mask+1:ptr_mask+nvtx) = 0
     nz = 0
     ! loop over coarse grid points
     do i = 1, cnvtx
       ! first form row i of (r*matrix)
       nz1 = 0
       ! for each vertex D that restricts to C (including itself).
       do j = rptr(i), rptr(i+1) - 1
         neigh = rcol(j)
         ! find D's neighbor
         do k = aptr(neigh), aptr(neigh+1) - 1
           neighneigh = acol(k)
           if (work(ptr_mask+neighneigh)/=i) then
             nz1 = nz1 + 1
             work(ptr_col+nz1) = neighneigh
             work(ptr_mask+neighneigh) = i
           end if
         end do
       end do
       ! form row i of (r*matrix)*p
       do j = 1, nz1
         neigh = work(ptr_col+j)
         do k = pptr(neigh), pptr(neigh+1) - 1
           neighneigh = pcol(k)
           if (work(ptr_mask+neighneigh)/=-i .and. neighneigh/=i) then
             nz = nz + 1
             work(ptr_mask+neighneigh) = -i
           end if
         end do
       end do
     end do

   end subroutine galerkin_graph_rap_size
   ! ******************************************************
   subroutine galerkin_graph_rap(nvtx,cnvtx,nzp,pa,pcol,pptr,nzaa,aa,acol, &
       aptr,nzr,ra,rcol,rptr,nzca,ca,ccol,cptr,lwork,work)
     ! multiply R*A*P to get CA
     ! nvtx: size of aa matrix
     ! cnvtx: size of ca matrix
     integer, intent(in) :: nvtx, cnvtx
     ! p: matrix
     integer, intent(in) :: nzp
     integer, intent(in), dimension(nzp) :: pa
     integer, intent(in), dimension(nzp) :: pcol
     integer, intent(in), dimension(nvtx+1) :: pptr
     ! aa: matrix
     integer, intent(in) :: nzaa
     integer, intent(in), dimension(:) :: aa
     integer, intent(in), dimension(nzaa) :: acol
     integer, intent(in), dimension(:) :: aptr
     ! r: matrix
     integer, intent(in) :: nzr
     integer, intent(in), dimension(nzr) :: ra
     integer, intent(in), dimension(nzr) :: rcol
     integer, intent(in), dimension(cnvtx+1) :: rptr
     ! ca: matrix
     integer, intent(in) :: nzca
     integer, intent(inout), dimension(nzca) :: ca
     integer, intent(inout), dimension(nzca) :: ccol
     integer, intent(inout), dimension(cnvtx+1) :: cptr

     integer, intent(in) :: lwork
     integer, intent(out) :: work(lwork)


     ! mask: masking array to see if an entry has been seen before
     integer :: ptr_mask
     ! i,j,k,l: loop index
     integer :: i, j, k
     ! nz: number of nonzeros so far in ca
     integer :: nz, nzz, nz1
     ! various neighbors
     integer :: neigh, neighneigh
     ! r_ij: (i,j) element of r
     integer :: r_ij
     ! col: column index of a row of r*matrix
     ! a: values of a row of r*matrix
     integer :: ptr_col, ptr_a

     ptr_mask = 0
     ptr_col = ptr_mask + nvtx
     ptr_a = ptr_col + nvtx
     ! now get the entries of the coarse matrix
     cptr(1) = 1
     work(ptr_mask+1:ptr_mask+nvtx) = 0
     nz = 0
     ! loop over every coarse grid point
     do i = 1, cnvtx
       ! first form row i of (r*matrix)
       nz1 = 0
       ! foreach each vertex D that restricts to C (including itself).
       do j = rptr(i), rptr(i+1) - 1
         neigh = rcol(j)
         r_ij = ra(j)
         ! find D's neighbor
         do k = aptr(neigh), aptr(neigh+1) - 1
           neighneigh = acol(k)
           nzz = work(ptr_mask+neighneigh)
           if (nzz==0) then
             nz1 = nz1 + 1
             work(ptr_col+nz1) = neighneigh
             work(ptr_a+nz1) = r_ij*aa(k)
             work(ptr_mask+neighneigh) = nz1
           else
             work(ptr_a+nzz) = work(ptr_a+nzz) + r_ij*aa(k)
           end if
         end do
       end do
       do j = 1, nz1
         work(ptr_mask+work(ptr_col+j)) = 0
       end do

       ! form row i of (r*matrix)*p
       do j = 1, nz1
         neigh = work(ptr_col+j)
         r_ij = work(ptr_a+j)
         do k = pptr(neigh), pptr(neigh+1) - 1
           neighneigh = pcol(k)
           if (neighneigh==i) cycle
           nzz = work(ptr_mask+neighneigh)
           if (nzz==0) then
             nz = nz + 1
             work(ptr_mask+neighneigh) = nz
             ca(nz) = r_ij*pa(k)
             ccol(nz) = neighneigh
           else
             ca(nzz) = ca(nzz) + r_ij*pa(k)
           end if
         end do
       end do

       do j = cptr(i), nz
         work(ptr_mask+ccol(j)) = 0
       end do
       cptr(i+1) = nz + 1
     end do


   end subroutine galerkin_graph_rap

   ! ---------------------------------------------------
   ! nd_refine_trim
   ! ---------------------------------------------------
   ! Given a partition, trim the partition to make it minimal
   subroutine nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
       a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,options)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i
     integer, intent(in) :: sumweight ! Sum of weights in a_weight
     integer, intent(inout) :: a_n1 ! Size of partition 1
     integer, intent(inout) :: a_n2 ! Size of partition 2
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end. This is updated to the new
     ! partition
     integer, intent(out) :: work(3*a_n) ! Work array
     type (nd_options), intent(in) :: options

     ! ---------------------------------------------
     ! Local variables
     ! ---------------------------------------------
     integer :: work_part, work_next, work_prev, a_n1_orig, a_n2_orig
     integer :: head1, head2, tail1, tail2
     integer, parameter :: sep1 = -1
     integer, parameter :: sep2 = -2
     integer, parameter :: sep3 = -3
     integer :: i, j, k, l, m, p, q, w1, w2
     logical :: next1, next2, imbal
     real(wp) :: t1, t2
     real(wp) :: ratio

     ratio = max(real(1.0,wp),options%balance)
     if (ratio>real(sumweight-2)) then
       imbal = .false.
     else
       imbal = .true.
     end if

     ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what
     ! part of the partition the nodes are in
     work_part = 0
     call nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
       ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
       work(work_part+1:work_part+a_n))

     a_n1_orig = a_n1
     a_n2_orig = a_n2


     ! Create two lists
     work_next = work_part + a_n
     work_prev = work_next + a_n
     head1 = 0
     head2 = 0
     tail1 = 0
     tail2 = 0
     do i = a_n1 + a_n2 + 1, a_n
       next1 = .false.
       next2 = .false.
       j = partition(i)
       if (j<a_n) then
         k = a_ptr(j+1) - 1
       else
         k = a_ne
       end if
       do l = a_ptr(j), k
         m = a_row(l)
         if (work(work_part+m)==ND_PART1_FLAG) then
           next1 = .true.
         else if (work(work_part+m)==ND_PART2_FLAG) then
           next2 = .true.
         end if
       end do
       if ((next1 .and. .not. next2) .or. ( .not. next1 .and. .not. next2 &
           .and. a_n1==0)) then
         ! Add to list 1
         if (head1==0) then
           head1 = j
           work(work_next+j) = 0
           work(work_prev+j) = 0
           tail1 = j
         else
           work(work_next+tail1) = j
           work(work_prev+j) = tail1
           work(work_next+j) = 0
           tail1 = j
         end if
         work(work_part+j) = sep1
       else if ((next2 .and. .not. next1) .or. ( .not. next1 .and. .not. &
           next2 .and. a_n2==0)) then
         ! Add to list 2
         if (head2==0) then
           head2 = j
           work(work_next+j) = 0
           work(work_prev+j) = 0
           tail2 = j
         else
           work(work_next+tail2) = j
           work(work_prev+j) = tail2
           work(work_next+j) = 0
           tail2 = j
         end if
         work(work_part+j) = sep2
       else if (next1 .and. next2) then
         work(work_part+j) = sep3
       else
         continue
       end if
     end do

     do while (head1>0 .or. head2>0)
       if (head1>0 .and. head2>0) then
         w1 = a_weight(head1)
         w2 = a_weight(head2)
         call cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1, &
           sumweight,ratio,imbal,options%cost_function,t1)
         call cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2, &
           sumweight,ratio,imbal,options%cost_function,t2)

         if (t1<t2) then
           go to 10
         else
           go to 20
         end if

       else if (head1>0) then
         go to 10
       else
         go to 20
       end if


       ! move entry from separator to partition1
10        i = head1
       work(work_part+i) = ND_PART1_FLAG
       head1 = work(work_next+i)
       work(work_next+i) = 0
       a_n1 = a_n1 + 1
       a_weight_1 = a_weight_1 + a_weight(i)
       a_weight_sep = a_weight_sep - a_weight(i)
       ! update list
       if (i<a_n) then
         k = a_ptr(i+1) - 1
       else
         k = a_ne
       end if
       do l = a_ptr(i), k
         j = a_row(l)
         m = work(work_part+j)
         select case (m)
         case (ND_SEP_FLAG)
           ! Add to list 1
           work(work_next+tail1) = j
           work(work_prev+j) = tail1
           work(work_next+j) = 0
           tail1 = j
           work(work_part+j) = sep1
           if (head1==0) then
             head1 = j
           end if

         case (sep2)
           ! Remove from list 2
           p = work(work_prev+j)
           q = work(work_next+j)

           if (j/=head2 .and. j/=tail2) then
             work(work_prev+q) = p
             work(work_next+p) = q
             work(work_prev+j) = 0
             work(work_next+j) = 0
           else if (j/=head2 .and. j==tail2) then
             work(work_next+p) = 0
             work(work_prev+j) = 0
             tail2 = p
           else if (j/=tail2 .and. j==head2) then
             work(work_prev+q) = p
             work(work_next+j) = 0
             head2 = q
           else
             head2 = 0
             tail2 = 0

           end if
           work(work_part+j) = sep3

         end select
       end do
       go to 30

       ! move entry from separator to partition 2
20        i = head2
       work(work_part+i) = ND_PART2_FLAG
       head2 = work(work_next+i)
       work(work_next+i) = 0
       a_n2 = a_n2 + 1
       a_weight_2 = a_weight_2 + a_weight(i)
       a_weight_sep = a_weight_sep - a_weight(i)
       ! update list
       if (i<a_n) then
         k = a_ptr(i+1) - 1
       else
         k = a_ne
       end if
       do l = a_ptr(i), k
         j = a_row(l)
         m = work(work_part+j)
         select case (m)
         case (ND_SEP_FLAG)
           ! Add to list 2
           work(work_next+tail2) = j
           work(work_prev+j) = tail2
           work(work_next+j) = 0
           tail2 = j
           work(work_part+j) = sep2
           if (head2==0) then
             head2 = j
           end if

         case (sep1)
           ! Remove from list 1
           p = work(work_prev+j)
           q = work(work_next+j)

           if (j/=head1 .and. j/=tail1) then
             work(work_prev+q) = p
             work(work_next+p) = q
             work(work_prev+j) = 0
             work(work_next+j) = 0
           else if (j/=head1 .and. j==tail1) then
             work(work_next+p) = 0
             work(work_prev+j) = 0
             tail1 = p
           else if (j/=tail1 .and. j==head1) then
             work(work_prev+q) = p
             work(work_next+j) = 0
             head1 = q
           else
             head1 = 0
             tail1 = 0
           end if
           work(work_part+j) = sep3
         end select
       end do

30        continue

     end do

     ! Check for any entries in separator that are still inside boundary
     ! and
     ! move into a partition
     work(work_next+a_n1_orig+a_n2_orig+1:work_next+a_n) = 0
     do i = a_n1_orig + a_n2_orig + 1, a_n
       j = partition(i)
       if (work(work_part+j)==ND_SEP_FLAG) then
         ! j is not on the boundary
         if (a_weight_1<a_weight_2) then
           ! Move j into partition 1
           work(work_part+j) = ND_PART1_FLAG
           a_n1 = a_n1 + 1
           a_weight_1 = a_weight_1 + a_weight(j)
           a_weight_sep = a_weight_sep - a_weight(j)

           head1 = j
           tail1 = j
           do while (head1>0)
             q = head1
             if (q<a_n) then
               k = a_ptr(q+1) - 1
             else
               k = a_ne
             end if
             do l = a_ptr(q), k
               p = a_row(l)
               if (work(work_part+p)==ND_SEP_FLAG) then
                 work(work_part+p) = ND_PART1_FLAG
                 a_n1 = a_n1 + 1
                 a_weight_1 = a_weight_1 + a_weight(p)
                 a_weight_sep = a_weight_sep - a_weight(p)
                 work(work_next+tail1) = p
                 tail1 = p
               end if
             end do
             if (head1==tail1) then
               head1 = 0
               tail1 = 0
             else
               head1 = work(work_next+q)
               work(work_next+q) = 0
             end if

           end do

         else
           ! Move j into partition 2
           work(work_part+j) = ND_PART2_FLAG
           a_n2 = a_n2 + 1
           a_weight_2 = a_weight_2 + a_weight(j)
           a_weight_sep = a_weight_sep - a_weight(j)
           head2 = j
           tail2 = j

           do while (head2>0)
             q = head2
             if (q<a_n) then
               k = a_ptr(q+1) - 1
             else
               k = a_ne
             end if
             do l = a_ptr(q), k
               p = a_row(l)
               if (work(work_part+p)==ND_SEP_FLAG) then
                 work(work_part+p) = ND_PART2_FLAG
                 a_n2 = a_n2 + 1
                 a_weight_2 = a_weight_2 + a_weight(p)
                 a_weight_sep = a_weight_sep - a_weight(p)
                 work(work_next+tail2) = p
                 tail2 = p
               end if
             end do
             if (head2==tail2) then
               head2 = 0
               tail2 = 0
             else
               head2 = work(work_next+q)
               work(work_next+q) = 0
             end if
           end do
         end if
       end if
     end do

     a_weight_sep = sumweight - a_weight_1 - a_weight_2

     ! Reset partition matrix
     call nd_convert_flags_partition(a_n,a_n1,a_n2, &
       work(work_part+1:work_part+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
       partition(1:a_n))

     ! call
     ! check_partition1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,work(work_part+1:wor
     ! k_part+a_n),a_weight_1,a_weight_2,a_weight)

   end subroutine nd_refine_trim


   ! ---------------------------------------------------
   ! nd_refine_block_trim
   ! ---------------------------------------------------
   ! Given a partition, trim the partition using blocks to make it minimal
   subroutine nd_refine_block_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
       sumweight,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,partition, &
       work,options)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i
     integer, intent(in) :: sumweight ! Sum of weights in a_weight
     integer, intent(inout) :: a_n1 ! Size of partition 1
     integer, intent(inout) :: a_n2 ! Size of partition 2
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end. This is updated to the new
     ! partition
     integer, intent(out) :: work(5*a_n) ! Work array
     type (nd_options), intent(in) :: options

     ! ---------------------------------------------
     ! Local variables
     integer :: work_part, work_next1, work_next2, work_level1, work_level2
     integer :: a_n1_orig, a_n2_orig
     integer :: head1, head2, tail1, tail2, maxlevel1, maxlevel2
     integer :: currlevel1, currlevel2
     integer :: i, j, k, l, m, w1, w2, l1, l2
     logical :: next1, next2, imbal
     real(wp) :: t1, t2
     real(wp) :: ratio

     ratio = max(real(1.0,wp),options%balance)
     if (ratio>real(sumweight-2)) then
       imbal = .false.
     else
       imbal = .true.
     end if

     ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what
     ! part of the partition the nodes are in
     work_part = 0
     call nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
       ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
       work(work_part+1:work_part+a_n))
     a_n1_orig = a_n1
     a_n2_orig = a_n2
     a_weight_sep = sumweight - a_weight_1 - a_weight_2


     ! Set work(work_next1+1:work_next1+a_n) to hold list pointers
     ! Set work(work_next2+1:work_next2+a_n) to hold list pointers
     ! Set work(work_level1+1:work_level1+a_n) to hold distance from orig
     ! partition 1
     ! Set work(work_level2+1:work_level2+a_n) to hold distance from orig
     ! partition 2
     work_next1 = work_part + a_n
     work_next2 = work_next1 + a_n
     work_level1 = work_next2 + a_n
     work_level2 = work_level1 + a_n
     work(work_next1+1:work_next1+a_n) = 0
     work(work_next2+1:work_next2+a_n) = 0
     work(work_level1+1:work_level1+a_n) = 0
     work(work_level2+1:work_level2+a_n) = 0

     ! Create two lists
     head1 = 0
     head2 = 0
     do i = a_n1 + a_n2 + 1, a_n
       next1 = .false.
       next2 = .false.
       j = partition(i)
       if (j<a_n) then
         k = a_ptr(j+1) - 1
       else
         k = a_ne
       end if
       do l = a_ptr(j), k
         m = a_row(l)
         if (work(work_part+m)==ND_PART1_FLAG) then
           next1 = .true.
         else if (work(work_part+m)==ND_PART2_FLAG) then
           next2 = .true.
         end if
       end do
       if (next1) then
         ! Add to list 1
         if (head1==0) then
           head1 = j
         else
           work(work_next1+tail1) = j
         end if
         tail1 = j
         work(work_level1+j) = 1
       end if
       if (next2) then
         ! Add to list 2
         if (head2==0) then
           head2 = j
         else
           work(work_next2+tail2) = j
         end if
         tail2 = j
         work(work_level2+j) = 1
       end if
     end do

     ! Breadth first search of separator from entries adjacent to partition
     ! 1
     l1 = head1
     do while (l1>0)
       if (l1<a_n) then
         k = a_ptr(l1+1) - 1
       else
         k = a_ne
       end if
       do l = a_ptr(l1), k
         m = a_row(l)
         if (work(work_part+m)==ND_SEP_FLAG .and. &
             work(work_level1+m)==0) then
           ! Add to list (note list is non-empty)
           work(work_next1+tail1) = m
           tail1 = m
           work(work_level1+m) = work(work_level1+l1) + 1
         end if
       end do
       l1 = work(work_next1+l1)
     end do
     maxlevel1 = work(work_level1+tail1)

     ! Breadth first search of separator from entries adjacent to partition
     ! 2
     l1 = head2
     do while (l1>0)
       if (l1<a_n) then
         k = a_ptr(l1+1) - 1
       else
         k = a_ne
       end if
       do l = a_ptr(l1), k
         m = a_row(l)
         if (work(work_part+m)==ND_SEP_FLAG .and. &
             work(work_level2+m)==0) then
           ! Add to list (note list is non-empty)
           work(work_next2+tail2) = m
           tail2 = m
           work(work_level2+m) = work(work_level2+l1) + 1
         end if
       end do
       l1 = work(work_next2+l1)
     end do
     maxlevel2 = work(work_level2+tail2)

     ! Check for any entries in separator only reachable from one partition
     do i = a_n1 + a_n2 + 1, a_n
       j = partition(i)
       if (work(work_level2+j)==0) then
         work(work_level2+j) = maxlevel2 + 1
       else
         if (work(work_level1+j)==0) then
           work(work_level1+j) = maxlevel1 + 1
         end if
       end if
     end do

     ! Trim the separator
     currlevel1 = 1
     currlevel2 = 1
     l1 = head1
     l2 = head2
     do while (currlevel1<=maxlevel1 .or. currlevel2<=maxlevel2)
       if (currlevel1>maxlevel1) then
         t1 = huge(1.0_wp)
       else
         w1 = 0
         j = l1
         do while (work(work_level1+j)==currlevel1)
           if (work(work_level2+j)>currlevel2) then
             w1 = w1 + a_weight(j)
           end if
           j = work(work_next1+j)
           if (j==0) exit
         end do
         if (w1==0) then
           currlevel1 = currlevel1 + 1
           l1 = j
           cycle
         else
           call cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1, &
             sumweight,ratio,imbal,options%cost_function,t1)
         end if
       end if

       if (currlevel2>maxlevel2) then
         t2 = huge(1.0_wp)
       else
         w2 = 0
         j = l2
         do while (work(work_level2+j)==currlevel2)
           if (work(work_level1+j)>currlevel1) then
             w2 = w2 + a_weight(j)
           end if
           j = work(work_next2+j)
           if (j==0) exit
         end do
         if (w2==0) then
           currlevel2 = currlevel2 + 1
           l2 = j
           cycle
         else
           call cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2, &
             sumweight,ratio,imbal,options%cost_function,t2)
         end if
       end if

       ! Add entries to relevant partition and update a_n1, a_n2 etc
       if (t1<t2) then
         j = l1
         do while (work(work_level1+j)==currlevel1)
           if (work(work_level2+j)>currlevel2) then
             work(work_part+j) = ND_PART1_FLAG
             a_n1 = a_n1 + 1
           end if
           j = work(work_next1+j)
           if (j==0) exit
         end do
         a_weight_1 = a_weight_1 + w1
         a_weight_sep = a_weight_sep - w1
         l1 = j
         if (j==0) then
           currlevel1 = maxlevel1 + 1
         else
           currlevel1 = (work(work_level1+l1))
         end if

       else
         j = l2
         do while (work(work_level2+j)==currlevel2)
           if (work(work_level1+j)>currlevel1) then
             work(work_part+j) = ND_PART2_FLAG
             a_n2 = a_n2 + 1
           end if
           j = work(work_next2+j)
           if (j==0) exit
         end do
         a_weight_2 = a_weight_2 + w2
         a_weight_sep = a_weight_sep - w2
         l2 = j
         if (j==0) then
           currlevel2 = maxlevel2 + 1
         else
           currlevel2 = (work(work_level2+l2))
         end if
       end if
     end do

     ! Reset partition matrix
     call nd_convert_flags_partition(a_n,a_n1,a_n2, &
       work(work_part+1:work_part+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
       partition(1:a_n))
     a_weight_sep = sumweight - a_weight_1 - a_weight_2
     ! call check_partition1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition)

   end subroutine nd_refine_block_trim


   ! ---------------------------------------------------
   ! nd_refine_max_flow
   ! ---------------------------------------------------
   ! Given a partition, trim the partition using blocks to make it minimal
   subroutine nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1, &
       a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,options)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i
     integer, intent(inout) :: a_n1 ! Size of partition 1
     integer, intent(inout) :: a_n2 ! Size of partition 2
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end. This is updated to the new
     ! partition
     integer, intent(out) :: work(8) ! Work array
     type (nd_options), intent(in) :: options

     ! ---------------------------------------------
     ! Local variables
     integer :: msglvl
     real(wp) :: cost, ratio

     msglvl = 0
     if (options%print_level==1 .and. options%unit_diagnostics>=0) &
       msglvl = 1
     if (options%print_level>=2 .and. options%unit_diagnostics>=0) &
       msglvl = 3


     if (a_n-a_n1-a_n2>1) then
       ratio = max(real(1.0,wp),options%balance)
       call nd_maxflow(a_n,a_ne,a_ptr,a_row,a_weight,options%cost_function,&
         a_n1,a_n2, &
         a_weight_1,a_weight_2,a_weight_sep,partition,ratio,msglvl, &
         work(1:8),cost)
     end if


   end subroutine nd_refine_max_flow

   subroutine expand_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
       a_weight_1,a_weight_2,a_weight_sep,partition,work)
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix
     integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
     ! indices of the non-zero rows. Diagonal entries have been removed
     ! and the matrix expanded.
     integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
     ! the weight of column i
     integer, intent(inout) :: a_n1 ! Size of partition 1
     integer, intent(inout) :: a_n2 ! Size of partition 2
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end. This is updated to the new
     ! partition
     integer, intent(out) :: work(a_n) ! Work array

     ! Local variables
     integer :: i, j, k, l, m, w
     integer :: work_part, a_weight_sep_orig
     ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what
     ! part of the partition the nodes are in
     work_part = 0
     do i = 1, a_n1
       j = partition(i)
       work(work_part+j) = ND_PART1_FLAG
     end do
     do i = a_n1 + 1, a_n1 + a_n2
       j = partition(i)
       work(work_part+j) = ND_PART2_FLAG
     end do
     do i = a_n1 + a_n2 + 1, a_n
       j = partition(i)
       work(work_part+j) = ND_SEP_FLAG
     end do

     ! if (a_weight_1 .lt. a_weight_2) then
     ! side = ND_PART2_FLAG
     ! else if (a_weight_1 .gt. a_weight_2) then
     ! side = ND_PART1_FLAG
     ! else
     ! side = ND_SEP_FLAG
     ! end if
     a_weight_sep_orig = a_weight_sep

     do i = a_n1 + a_n2 + 1, a_n
       j = partition(i)
       ! search neighbours of j and add to separator
       if (j==a_n) then
         k = a_ne
       else
         k = a_ptr(j+1) - 1
       end if
       do l = a_ptr(j), k
         m = a_row(l)
         if (work(work_part+m)==ND_PART1_FLAG .and. a_n1>1) then
           ! if (side .eq. ND_PART1_FLAG .or. side .eq. ND_SEP_FLAG)
           ! then
           work(work_part+m) = ND_SEP_FLAG
           a_n1 = a_n1 - 1
           w = a_weight(m)
           a_weight_1 = a_weight_1 - w
           a_weight_sep = a_weight_sep + w
           ! end if
         else if (work(work_part+m)==ND_PART2_FLAG .and. a_n2>1) then
           ! if (side .eq. ND_PART2_FLAG .or. side .eq. ND_SEP_FLAG)
           ! then
           work(work_part+m) = ND_SEP_FLAG
           a_n2 = a_n2 - 1
           w = a_weight(m)
           a_weight_2 = a_weight_2 - w
           a_weight_sep = a_weight_sep + w
           ! end if
         end if
       end do
     end do
     j = 1
     k = j + a_n1
     l = k + a_n2
     do i = 1, a_n
       m = work(work_part+i)
       select case (m)
       case (ND_PART1_FLAG)
         partition(j) = i
         j = j + 1
       case (ND_PART2_FLAG)
         partition(k) = i
         k = k + 1
       case default
         partition(l) = i
         l = l + 1
       end select
     end do

   end subroutine expand_partition


   subroutine cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight, &
       ratio,imbal,costf,tau)

     integer, intent(in) :: a_weight_1, a_weight_2, a_weight_sep ! Weighte
     ! d
     ! size of partitions and separator
     integer, intent(in) :: sumweight
     real(wp), intent(in) :: ratio
     logical, intent(in) :: imbal ! Use penalty function?
     integer, intent(in) :: costf
     real(wp), intent(out) :: tau
     real(wp) :: beta, a_wgt1, a_wgt2

     beta = 0.5
     a_wgt1 = max(1,a_weight_1)
     a_wgt2 = max(1,a_weight_2)

     if (costf.le.1) then
       tau = ((real(a_weight_sep)**1.0)/real(a_wgt1))/real(a_wgt2)
       if (imbal .and. real(max(a_wgt1,a_wgt2))/real(min(a_wgt1, &
           a_wgt2))>=ratio) then
         tau = real(sumweight-2) + tau
       end if
     else
       if (imbal .and. real(max(a_wgt1,a_wgt2))/real(min(a_wgt1, &
           a_wgt2))>=ratio) then
         tau = real(sumweight)*(1.0+beta) + real(a_weight_sep)*( &
           1.0_wp+beta*real(abs(a_wgt1-a_wgt2))/real(sumweight))
       else
         tau = real(a_weight_sep)*(1.0_wp+beta*real(abs(a_wgt1- &
           a_wgt2))/real(sumweight))
       end if

     end if

   end subroutine cost_function

   ! ---------------------------------------------------
   ! nd_maxflow
   ! ---------------------------------------------------
   ! Given a partition, get better partition using maxflow algorithm
   subroutine nd_maxflow(a_n,a_ne,a_ptr,a_row,a_weight,costf,a_n1,a_n2, &
       a_weight_1,a_weight_2,a_weight_sep,partition,alpha,msglvl,stats, &
       cost)

     ! Input matrix: a_n, a_ne, a_ptr, a_row
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix (lower and
     ! upper triangle)
     integer, intent(in) :: a_ptr(:) ! On input, a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(:) ! On input, a_row contains row
     ! indices of the nonzero entries. Diagonal entries have been
     ! removed and the matrix expanded.
     ! At the moment weights are not used at all
     integer, intent(in) :: a_weight(a_n) ! On input, a_weight(i) contains
     ! the weight of column i
     integer, intent(in) :: costf ! Determines which cost function is used
     ! Data on partition a_n1, a_n2, partition ... will be updated
     integer, intent(inout) :: a_n1 ! Size of partition 1 (ie B)
     integer, intent(inout) :: a_n2 ! Size of partition 2 (ie W)
     integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
     ! hted
     ! size of partitions and separator
     integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end. This is updated to the new
     ! partition.

     ! Parameters alpha (for balance) for cost function
     real(wp), intent(in) :: alpha
     integer, intent(in) :: msglvl
     ! output --
     ! stats[1] -- weight of vertices in S
     ! stats[2] -- weight of vertices in B
     ! stats[3] -- weight of vertices in W
     ! stats[4] -- weight of edges in A_{S,S}
     ! stats[5] -- weight of edges in A_{S,B}
     ! stats[6] -- weight of edges in A_{S,W}
     ! stats[7] -- weight of edges in A_{B,B}
     ! stats[8] -- weight of edges in A_{B,W}
     ! cost     -- cost of new partition
     integer, intent(out) :: stats(8)
     real(wp), intent(out) :: cost

     type (network) :: netw

     ! Work arrays
     ! map,mapL,mapR of length a_n
     ! dmapL,dmapR,vwts of length a_ns
     ! mark1,mark2,pred,list of length number of nodes in network that
     ! is bounded by 2*a_ns+2

     ! Eventually these arrays will be allocated higher up the chain and
     ! will be passed as parameters
     integer, allocatable :: map(:), mapl(:), mapr(:)
     integer, allocatable :: dmapl(:), dmapr(:)
     integer, allocatable :: vwts(:)
     integer, allocatable :: sedge(:,:)
     integer, allocatable :: mark1(:), mark2(:), pred(:), list(:)
     real(wp), allocatable :: imb(:)

     ! Local variables
     integer :: a_ns, i, istart_s, j1, k, lp, wtw, wtb, statsr(9), &
       statsl(9)
     integer nedge, matsiz
     real(wp) :: costr, costl

     lp = 6

     ! write(9,*) 'Entering maxflow'
     ! write(0,*) 'Entering maxflow'
     allocate (map(a_n),mapl(a_n),mapr(a_n))


     ! Number vertices in separator
     a_ns = a_n - a_n1 - a_n2

     allocate (dmapl(2*a_ns),dmapr(2*a_ns),vwts(a_ns))

     ! Allocate network work arrays.  Length is upper bound
     ! allocate (mark1(2*a_ns+2),mark2(2*a_ns+2),pred(2*a_ns+2), &
     ! list(2*a_ns+2))
     matsiz = max(a_n,2*a_ns+2)
     allocate (mark1(matsiz),mark2(matsiz),pred(2*a_ns+2),list(matsiz))
     allocate (imb(2*a_ns))

     ! Set up map array to define in what partition each vertex lies
     ! At same time set weights for partition (can check with Sue's input)

     call nd_convert_partition_flags(a_n,a_n1,a_n2,partition,1,2,0, &
       map(1:a_n))

     wtb = a_weight_1
     wtw = a_weight_2

     do i = a_n1 + a_n2 + 1, a_n
       k = partition(i)
       vwts(i-a_n1-a_n2) = a_weight(k)
     end do


     ! Count edges to get upper bound on size of sedge array
     nedge = 0
     do k = 1, a_ns
       i = partition(a_n1+a_n2+k)
       j1 = a_ptr(i)
       if (i==a_n) then
         nedge = nedge + a_ne + 1 - a_ptr(i)
       else
         nedge = nedge + a_ptr(i+1) - a_ptr(i)
       end if
     end do
     allocate (sedge(nedge,2))


     ! Generate network graph.  The structure for our maxflow algorithm

     ! Work arrays dmapL and dmapR used to hold isadjsource and isadjsink
     ! Work array  mapL used to hold sep_map
     ! Source is associated with partition B (size a_n1)
     ! Sink is associated with partition W   (size a_n2)
     ! write(9,*) 'Calling mk_network'
     ! write(0,*) 'Calling mk_network'
     call mk_network(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,map,a_ns, &
       msglvl,netw,vwts,wtb,wtw,sedge,mapl,dmapl,dmapr,list,mark1, &
       mark2,imb)
     ! write(9,*) 'Leaving mk_network'
     ! write(0,*) 'Leaving mk_network'


     ! solve a max flow problem to find the two new maps dmapL and dmapR

     call solvemaxflow(netw,dmapl,dmapr,mark1,mark2,pred,list)


     if (msglvl>2) then
       write (lp,*) 'dmapL ...'
       write (lp,'(10I4)') dmapl
       write (lp,*) 'dmapR ...'
       write (lp,'(10I4)') dmapr
     end if

     mapl = map
     mapr = map
     istart_s = a_n1 + a_n2
     do i = 1, a_ns
       mapl(partition(istart_s+i)) = dmapl(i)
       mapr(partition(istart_s+i)) = dmapr(i)
     end do

     if (msglvl>2) then
       write (lp,*) 'mapL ...'
       write (lp,'(10I4)') mapl
       write (lp,*) 'mapR ...'
       write (lp,'(10I4)') mapr
     end if

     ! Use evaluation function to choose best partition from among these
     ! two
     ! Use Sue's weighted code
     call evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,mapl,alpha,costf,statsl,costl)
     call evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,mapr,alpha,costf,statsr,costr)


     ! Find the better of the two partitions


     if (statsl(9)==1 .and. statsr(9)==1) then
       if (msglvl>0) write (lp,'(A)') 'both maps are acceptable'
       if (costl<=costr) then
         map = mapl
         stats = statsl(1:8)
         cost = costl
         if (msglvl>0) write (lp,'(A)') 'left map accepted'
       else
         map = mapr
         stats = statsr(1:8)
         cost = costr
         if (msglvl>0) write (lp,'(A)') 'right map accepted'
       end if
     else if (statsl(9)==1) then
       map = mapl
       stats = statsl(1:8)
       cost = costl
       if (msglvl>0) write (lp,'(A)') &
         'right map not acceptable, left map accepted'
     else if (statsr(9)==1) then
       map = mapr
       stats = statsr(1:8)
       cost = costr
       if (msglvl>0) write (lp,'(A)') &
         'left map not acceptable, right map accepted'
     else
       if (msglvl>0) write (lp,'(A)') 'NEITHER map acceptable'
       if (costl<=costr) then
         map = mapl
         stats = statsl(1:8)
         cost = costl
         if (msglvl>0) write (lp,'(A)') 'left map accepted'
       else
         map = mapr
         stats = statsr(1:8)
         cost = costr
         if (msglvl>0) write (lp,'(A)') 'right map accepted'
       end if
     end if
     a_weight_1 = stats(2)
     a_weight_2 = stats(3)
     a_weight_sep = stats(1)

     ! Now update partition
     ! First count number of vertices in each part
     a_n1 = 0
     a_n2 = 0
     a_ns = 0
     do i = 1, a_n
       if (map(i)==1) a_n1 = a_n1 + 1
       if (map(i)==2) a_n2 = a_n2 + 1
       if (map(i)==0) a_ns = a_ns + 1
     end do


     call nd_convert_flags_partition(a_n,a_n1,a_n2,map(1:a_n),1,2, &
       partition(1:a_n))

     deallocate (map)
     deallocate (dmapl,dmapr)
     deallocate (mapl,mapr)
     deallocate (mark1,mark2,pred,list,vwts)
     deallocate (sedge)


   end subroutine nd_maxflow

   subroutine solvemaxflow(netw,maps1,maps2,mark1,mark2,pred,list)
     ! Find two partitions of a wide separator graph by solving a max flow
     ! problem




     ! output --

     ! mapS1[n_S] -- first map from wide separator to {0,1,2} = {S,B,W}
     ! mapS2[n_S] -- second map from wide separator to {0,1,2} = {S,B,W}


     ! input/output --

     ! network graph

     ! work arrays

     ! mark1, mark2, pred, list of length nnode


     type (network), intent(inout) :: netw
     integer, intent(out) :: maps1(:), maps2(:)

     integer :: mark1(:), mark2(:), pred(:), list(:)

     ! Local variables
     integer ii, lp, narc, nnode, u

     lp = 6


     nnode = netw%nnode
     narc = netw%narc


     ! Find maxflow through the network using the Ford-Fulkerson algorithm

     call findmaxflow(netw,pred,list,mark1,mark2)



     ! Find the two mincuts

     call findmincut(netw,mark1,mark2,list)


     ! Use mark1 and mark2 to generate maps1 and maps2


     maps1 = -1
     maps2 = -1

     do ii = 2, nnode - 1, 2
       u = ii/2
       if (mark1(ii)==1) then
         if (mark1(ii+1)==1) then
           maps1(u) = 1
         else
           maps1(u) = 0
         end if
       end if
       if (mark1(ii)==2) then
         if (mark1(ii+1)==2) then
           maps1(u) = 2
         else
           maps1(u) = 0
         end if
       end if
     end do

     do ii = 2, nnode - 1, 2
       u = ii/2
       if (mark2(ii)==1) then
         if (mark2(ii+1)==1) then
           maps2(u) = 1
         else
           maps2(u) = 0
         end if
       end if
       if (mark2(ii)==2) then
         if (mark2(ii+1)==2) then
           maps2(u) = 2
         else
           maps2(u) = 0
         end if
       end if
     end do

   end subroutine solvemaxflow


   subroutine mk_network(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,map,nvtx, &
       msglvl,netw,vwts,wtb,wtw,sedge,sep_map,isadjtosource,isadjtosink, &
       list,mark1,mark2,imb)
     ! Create and return a network structure


     ! Input matrix: a_n, a_ne, a_ptr, a_row
     integer, intent(in) :: a_n ! order of matrix
     integer, intent(in) :: a_ne ! number of entries in matrix (lower and
     ! upper triangle)
     integer, intent(in) :: a_ptr(:) ! On input, a_ptr(i) contains
     ! position in a_row that entries for column i start.
     integer, intent(in) :: a_row(:) ! On input, a_row contains row
     ! indices of the nonzero entries. Diagonal entries have been
     ! removed and the matrix expanded.

     integer, intent(in) :: partition(a_n) ! First a_n1 entries contain
     ! list of (local) indices in partition 1; next a_n2 entries
     ! contain list of (local) entries in partition 2; entries in
     ! separator are listed at the end.
     integer, intent(in) :: map(a_n) ! First a_n1 entries contain
     integer, intent(in) :: msglvl, nvtx
     integer, intent(in) :: vwts(:)
     integer, intent(inout) :: wtb, wtw
     integer, intent(inout) :: a_n1 ! Size of partition 1 (ie B)
     integer, intent(inout) :: a_n2 ! Size of partition 2 (ie W)

     ! Note that we still do the allocations here.  Doing it further up the
     ! call tree would need to allocate much more storage.
     type (network), intent(out) :: netw

     integer i, iarc, ii, jj, lp, narc1, narc2, narc3, narc4, narc, nnode, &
       nedge, u, v, wts
     integer j1, j2, j, k

     ! Work arrays sedge(nedge,2),sep_map(a_n),isAdjToSource(nvtx),
     ! isAdjToSink(nvtx)
     integer :: sedge(:,:)
     integer :: sep_map(:)
     integer :: isadjtosource(:), isadjtosink(:)
     real(wp) :: imb(:)
     integer :: mark1(:), mark2(:), list(:)
     logical :: augcap

     augcap = .true.

     lp = 0

     if (msglvl>0) then
       write (lp,'(/A)') '### inside mknetwork()'
       write (lp,'(A,I8,A,I8)') 'nvtx', nvtx
     end if


     ! write(0,*) 'a_n,a_n1,a_n2,a_ne',   &
     ! a_n,a_n1,a_n2,a_ne
     ! write(0,*) 'a_ptr',a_ptr(1:a_n)
     ! write(0,*) 'a_row',a_row(1:a_ne)
     ! write(0,*) 'partition',partition(1:a_n)
     ! write(0,*) 'map',map(1:a_n)
     ! write(0,*) 'vwts',vwts(1:nvtx)
     ! write(0,*) 'wtb',wtb
     ! write(0,*) 'wtw',wtw

     isadjtosource = 0
     isadjtosink = 0

     ! Determine mapping of global variables of matrix to separator set
     do k = 1, nvtx
       i = partition(a_n1+a_n2+k)
       sep_map(i) = k
     end do

     ! We could use a single array although the logic for generating it is
     ! marginally more complicated.
     ! For nodes in separator, set isadj as
     ! 1 if only connected to source
     ! 2 if only connected to sink
     ! 3 if connected to source and sink
     ! isadj = 0


     ! Run through nodes in separator S and generate edges
     nedge = 0
     do k = 1, nvtx
       i = partition(a_n1+a_n2+k)
       j1 = a_ptr(i)
       if (i==a_n) then
         j2 = a_ne
       else
         j2 = a_ptr(i+1) - 1
       end if

       ! Run through vertices connected to vertex i
       do jj = j1, j2
         j = a_row(jj)
         ! Find out in which partition node j lies using map array
         if (map(j)==1) then
           ! If in partition B add vertex k to AdjToSource
           isadjtosource(k) = 1
         end if
         if (map(j)==2) then
           ! If in partition W add vertex k to AdjToSink
           isadjtosink(k) = 1
         end if
         if (map(j)==0) then
           ! If in separator add edge accumulating number of edges
           nedge = nedge + 1
           ! Edge has this orientation to emulate matlab code
           sedge(nedge,2) = sep_map(i)
           sedge(nedge,1) = sep_map(j)
         end if
       end do

     end do


     ! narc1 is number of vertices in separator and is equal to the number
     ! of
     ! added edges u- to u+
     narc1 = nvtx
     ! narc2 is number of vertices in separator connected to source
     narc2 = 0
     ! narc3 is number of vertices in separator connected to sink
     narc3 = 0


     ! count the number of arcs
     ! Can't do easily in above loop because of multiple edges from
     ! source/sink
     ! to vertex in the separator set

     do k = 1, nvtx
       if (isadjtosource(k)==1) narc2 = narc2 + 1
       if (isadjtosink(k)==1) narc3 = narc3 + 1
     end do

     narc4 = 0
     do ii = 1, nedge
       u = sedge(ii,1)
       v = sedge(ii,2)
       ! --- ignore self edges ---
       if (u==v) cycle
       ! --- omit edges with essential vertices ---
       if (isadjtosink(u)==1 .and. isadjtosource(u)==1) cycle
       if (isadjtosink(v)==1 .and. isadjtosource(v)==1) cycle
       ! --- omit pairs both adjacent to source ---
       if ((isadjtosource(u)==1) .and. (isadjtosource(v)==1)) cycle
       ! --- omit pairs both adjacent to sink ---
       if ((isadjtosink(u)==1) .and. (isadjtosink(v)==1)) cycle
       ! Qualifying arc found
       narc4 = narc4 + 1
     end do

     nnode = 2*nvtx + 2
     netw%nnode = nnode
     narc = narc1 + narc2 + narc3 + narc4
     netw%narc = narc
     netw%source = 1
     netw%sink = 2*nvtx + 2

     if (msglvl>0) then
       write (lp,'(I8,A)') narc1, ' internal arcs'
       write (lp,'(I8,A)') narc2, ' arcs from source'
       write (lp,'(I8,A)') narc3, ' arcs from sink'
       write (lp,'(I8,A)') narc4, ' edge arcs'
       write (lp,'(I8,A)') narc, ' total arcs'
     end if


     ! create the arc arrays



     ! Allocations done here but could be easily moved up the path although
     ! values very dependent on separator size.
     allocate (netw%firsts(narc),netw%seconds(narc),netw%capacities(narc))
     allocate (netw%flows(narc))
     allocate (netw%inheads(nnode),netw%outheads(nnode),netw%nextin(narc), &
       netw%nextout(narc))

     netw%firsts = -1
     netw%seconds = -1


     ! (u-,u+) arcs first

     iarc = 0
     do u = 1, nvtx
       iarc = iarc + 1
       netw%firsts(iarc) = 2*u
       netw%seconds(iarc) = 2*u + 1
       ! We set capacities after computing imbalance penalty
       ! netw%capacities(iarc) = vwts(u)
     end do

     if (msglvl>0) write (lp,'(A,I8)') 'after (u-,u+) arcs, iarc = ', iarc


     ! (source,u) arcs

     do u = 1, nvtx
       if (isadjtosource(u)==1) then
         iarc = iarc + 1
         netw%firsts(iarc) = netw%source
         netw%seconds(iarc) = 2*u
         netw%capacities(iarc) = huge(1)/2
       end if
     end do

     if (msglvl>0) write (lp,'(A,I8)') 'after (source,u-) arcs, iarc = ', &
       iarc


     ! (u,sink) arcs

     do u = 1, nvtx
       if (msglvl>5) write (lp,'(A,I4,A,I8)') 'isAdjToSink(', u, ')= ', &
         isadjtosink(u)
       if (isadjtosink(u)==1) then
         iarc = iarc + 1
         netw%firsts(iarc) = 2*u + 1
         netw%seconds(iarc) = netw%sink
         netw%capacities(iarc) = huge(1)/2
       end if
     end do

     if (msglvl>0) write (lp,'(A,I8)') 'after (u,sink) arcs, iarc = ', iarc


     ! (u+,v-) arcs

     do ii = 1, nedge
       u = sedge(ii,1)
       v = sedge(ii,2)
       if ((u/=v) .and. (isadjtosource(u)/=1 .or. isadjtosource( &
           v)/=1) .and. (isadjtosink(u)/=1 .or. isadjtosink( &
           v)/=1) .and. (isadjtosource(u)/=1 .or. isadjtosink( &
           u)/=1) .and. (isadjtosource(v)/=1 .or. isadjtosink(v)/=1)) then
         iarc = iarc + 1
         netw%firsts(iarc) = 2*u + 1
         netw%seconds(iarc) = 2*v
         netw%capacities(iarc) = huge(1)/2
       end if
     end do

     if (msglvl>0) write (lp,'(A,I8)') 'after (u+,v-) arcs, iarc = ', iarc


     ! Generate the head vectors for in/out edges
     ! and the in/out link vectors for the arcs


     netw%inheads = -1
     netw%outheads = -1
     netw%nextin = -1
     netw%nextout = -1
     do ii = narc, 1, -1
       u = netw%firsts(ii)
       v = netw%seconds(ii)
       if (msglvl>1) write (lp,'(A,I8,A,I8,A,I8,A)') 'ii', ii, 'arc (', u, &
         ',', v, ')'
       netw%nextin(ii) = netw%inheads(v)
       netw%inheads(v) = ii
       netw%nextout(ii) = netw%outheads(u)
       netw%outheads(u) = ii
     end do

     if (augcap) then
       ! Generate wtS
       wts = 0
       do i = 1, nvtx
         wts = wts + vwts(i)
       end do
       if (msglvl>0) write (lp,*) 'Calling findpenalty'
       ! Compute network capacities for separator arcs.
       call findpenalty(msglvl,a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
         map,vwts,wtb,wtw,wts,isadjtosource,isadjtosink,mark1,mark2,list, &
         imb)
       if (msglvl>0) write (lp,*) 'Exiting findpenalty'
       do i = 1, nvtx
         netw%capacities(i) = mark1(i)
       end do

     else
       do i = 1, nvtx
         netw%capacities(i) = vwts(i)
       end do
     end if

     ! ISD Would this not be initialized in routine that computes flows?
     netw%flows = 0

     if (msglvl>0) write (lp,'(A/)') '### leaving mknetwork()'

   end subroutine mk_network


   subroutine findmaxflow(netw,pred,list,tags,deltas)
     ! Find a maximum flow through the network


     type (network), intent(inout) :: netw

     integer avail, iarc, lp, nnode, sink, source, stats(2), tag

     ! Work arrays ... all of length nnode

     integer :: pred(:), list(:), tags(:), deltas(:)

     lp = 6

     nnode = netw%nnode
     source = netw%source
     sink = netw%sink

     ! Network flows initialized to zero
     netw%flows = 0

     ! tag is just used to count which arc from source is being used to
     ! start augmenting path
     tag = 0
     iarc = netw%outheads(source)
     ! Run through all nodes connected to source
     do
       if (iarc==-1) exit

       do
         avail = netw%capacities(iarc) - netw%flows(iarc)
         if (avail>0) then
           ! Set in findaugpath

           ! use BFS to find path

           ! avail is an output from findaugpath giving available flow on
           ! augmenting
           ! path.  Not to be confused with dummy avail above. The
           ! augmenting
           ! path is given through the array pred.

           call findaugpath(netw,iarc,avail,pred,stats,list,tags,deltas)


           ! Go to next arc from source node if no augmenting path has been
           ! found
           if ((avail==0) .or. (pred(sink)==0)) exit

           ! Update flows
           call augmentpath(netw,avail,pred)

         else
           exit
         end if
       end do

       iarc = netw%nextout(iarc)
       tag = tag + 1

     end do

   end subroutine findmaxflow

   subroutine findmincut(netw,mark1,mark2,list)

     ! Finds one or two mincuts, one nearest the source, one nearest the
     ! sink

     ! Input parameters

     ! netw    -- network object
     ! msglvl  -- message level

     ! Output parameters
     ! mark1 --- to identify cut set nearest source
     ! mark2 --- to identify cut set nearest sink


     ! Workspace
     ! list --- to hold list of nodes being searched

     integer, intent(out) :: mark1(:), mark2(:), list(:)
     type (network), intent(inout) :: netw

     ! Local variables
     integer iarc, last, lp, nnode, now, sink, source, x, z

     lp = 6

     nnode = netw%nnode
     source = netw%source
     sink = netw%sink


     ! breadth first traversal from source


     mark1 = 2
     mark1(source) = 1

     list = 0
     now = 1
     last = 1
     list(now) = source
     ! while now <= last
     do
       if (now>last) exit
       x = list(now)
       now = now + 1
       iarc = netw%outheads(x)
       ! while iarc ~= -1
       ! Run through all arcs starting at node x and putting node at end of
       ! arc
       ! on list if there is spare capacity on arc
       do
         if (iarc==-1) exit
         z = netw%seconds(iarc)
         if (mark1(z)==1) then
         else
           if (netw%flows(iarc)<netw%capacities(iarc)) then
             last = last + 1
             list(last) = z
             mark1(z) = 1
           end if
         end if
         iarc = netw%nextout(iarc)
       end do
       iarc = netw%inheads(x)
       ! while iarc ~= -1
       ! Run through all arcs terminating at node x and putting node at
       ! start of arc
       ! on list if there is spare capacity on arc
       do
         if (iarc==-1) exit
         z = netw%firsts(iarc)
         if (mark1(z)==1) then
         else
           if (netw%flows(iarc)>0) then
             last = last + 1
             list(last) = z
             mark1(z) = 1
           end if
         end if
         iarc = netw%nextin(iarc)
       end do
     end do

     ! breadth first traversal from sink


     mark2 = 1
     mark2(sink) = 2

     list = 0
     now = 1
     last = 1
     list(now) = sink
     ! while now <= last
     do
       if (now>last) exit
       x = list(now)
       now = now + 1
       iarc = netw%outheads(x)
       ! while iarc ~= -1
       ! Run through all arcs starting at node x and putting node at end of
       ! arc
       ! on list if there is spare capacity on arc
       do
         if (iarc==-1) exit
         z = netw%seconds(iarc)
         if (mark2(z)==2) then
         else
           if (netw%flows(iarc)>0) then
             last = last + 1
             list(last) = z
             mark2(z) = 2
           end if
         end if
         iarc = netw%nextout(iarc)
       end do
       iarc = netw%inheads(x)
       ! while iarc ~= -1
       ! Run through all arcs terminating at node x and putting node at
       ! start of arc
       ! on list if there is spare capacity on arc
       do
         if (iarc==-1) exit
         z = netw%firsts(iarc)
         if (mark2(z)==2) then
         else
           if (netw%flows(iarc)<netw%capacities(iarc)) then
             last = last + 1
             list(last) = z
             mark2(z) = 2
           end if
         end if
         iarc = netw%nextin(iarc)
       end do
     end do


   end subroutine findmincut

   subroutine findpenalty(msglvl,a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
       map,vwts,wtb,wtw,wts,count,head,mark,mark1,list,imb)

     ! count is isAdjToSource in call
     ! head is isAdjToSink in call

     ! Computes augmented capacities using weighting based on balance

     ! Input parameters

     ! msglvl  -- message level

     ! Output parameters
     ! mark --- records level set of vertex in cutset and two-level set
     ! also
     ! mark1 --- records level set of vertex in cutset starting from sink


     integer, intent(in) :: msglvl, a_n, a_ne, a_n1, a_n2
     integer, intent(in) :: a_ptr(:), a_row(:)
     integer, intent(in) :: vwts(:), partition(:), map(:)
     integer, intent(inout) :: count(:), head(:)
     integer, intent(inout) :: wtb, wtw, wts
     integer, intent(out) :: mark(:), mark1(:)

     ! Workspace
     ! list --- to hold list of nodes being searched
     ! length bounded by a_n
     integer :: list(:)
     real(wp) :: imb(:)

     ! Local variables
     integer inode, last, lp, maxl, minl, x, z
     integer i, j1, j2, jj, k, nvtx
     integer begin_lev, end_lev, ilev, penp

     lp = 0

     ! Set imbalance weighting
     penp = 100

     ! nvtx is number of vertices in separator
     nvtx = a_n - a_n1 - a_n2

     if (msglvl>0) then
       write (lp,'(A)') ''
       write (lp,'(A)') '### inside findpenalty()'
     end if


     ! Breadth first traversal from source
     ! Source is defined as all vertices in black partition

     if (msglvl>0) write (lp,'(A)') 'breadth first traversal from source'

     mark = 0
     ! Black vertices at level 0
     ! Run through separator vertices
     last = 0
     do k = a_n1 + a_n2 + 1, a_n
       z = partition(k)
       ! Check if adjacent to source
       if (count(k-a_n1-a_n2)==1) then
         last = last + 1
         list(last) = z
         mark(z) = 1
       end if
     end do
     end_lev = 0

     ! Each pass through this loop determines all nodes in level set k
     ! a_n is just a dummy end
     do k = 2, a_n
       ! Run through all nodes in the previous level set
       begin_lev = end_lev + 1
       end_lev = last
       do ilev = begin_lev, end_lev
         x = list(ilev)
         if (msglvl>1) write (lp,'(A,I10)') 'Processing vertex', x
         ! Run through vertices connected to vertex x
         j1 = a_ptr(x)
         if (x==a_n) then
           j2 = a_ne
         else
           j2 = a_ptr(x+1) - 1
         end if
         do jj = j1, j2
           z = a_row(jj)
           ! Jump if vertex not in separator
           if (map(z)/=0) cycle
           ! Jump if vertex visited already
           if (mark(z)/=0) cycle
           mark(z) = k
           ! write(0,*) 'z,mark(z)',z,mark(z)
           ! Add node z to list
           last = last + 1
           list(last) = z
         end do
       end do
       if (last==end_lev) exit
     end do ! end of processing of nodes on level k


     ! breadth first traversal from sink
     ! Sink is defined as all vertices in white partition

     if (msglvl>0) write (lp,'(A)') 'breadth first traversal from the sink'

     mark1 = 0
     ! Keep sink at level 0
     maxl = -a_n
     minl = a_n

     ! White nodes at level 0
     ! Put all separator vertices connected to sink in list
     last = 0
     do k = a_n1 + a_n2 + 1, a_n
       z = partition(k)
       ! Check if adjacent to source
       if (head(k-a_n1-a_n2)==1) then
         last = last + 1
         list(last) = z
         mark1(z) = 1
         ! k = 1
         mark(z) = mark(z) - 1
         minl = min(minl,mark(z))
         maxl = max(maxl,mark(z))
       end if
     end do
     end_lev = 0

     ! Each pass through this loop determines all nodes in level set k
     ! a_n is just a dummy end
     do k = 2, a_n
       ! Run through all nodes in the previous level set
       begin_lev = end_lev + 1
       end_lev = last
       do ilev = begin_lev, end_lev
         x = list(ilev)
         if (msglvl>1) write (lp,'(A,I10)') 'Processing vertex', x
         ! Run through vertices connected to vertex x
         j1 = a_ptr(x)
         if (x==a_n) then
           j2 = a_ne
         else
           j2 = a_ptr(x+1) - 1
         end if
         do jj = j1, j2
           z = a_row(jj)
           ! Jump if vertex not in separator
           if (map(z)/=0) cycle
           ! Jump if vertex visited already
           if (mark1(z)/=0) cycle
           mark1(z) = k
           ! write(0,*) 'z,mark1(z)',z,mark1(z)
           ! write(0,*) 'mark(z)',mark(z)
           mark(z) = mark(z) - k
           ! write(0,*) 'z,mark(z)',z,mark(z)
           minl = min(minl,mark(z))
           maxl = max(maxl,mark(z))
           ! Add node z to list
           last = last + 1
           list(last) = z
         end do
       end do
       if (last==end_lev) exit
     end do ! end of processing of nodes on level k

     ! Compute half-level sets
     if (msglvl>1) write (lp,'(A,2I4)') 'minl, maxl ', minl, maxl

     ! We will number levels from 1 to maxl-minl+1
     ! count will hold total weight of all vertices in half-level set
     ! Nodes in level set are accessed by linked list in list with headers
     ! in head
     count(1:maxl-minl+1) = 0
     head(1:maxl-minl+1) = -1

     ! Map the mark values into the local coordinates where separator
     ! vertices numbered from 1 to nvtx
     do i = 1, nvtx
       k = partition(a_n1+a_n2+i)
       mark1(i) = mark(k)
     end do

     ! Run through all notes in cutset, resetting numbering of level set
     ! putting
     ! them in level set for level and accumulating weight of half-level
     ! set
     do i = 1, nvtx
       mark1(i) = mark1(i) - minl + 1
       list(i) = head(mark1(i))
       head(mark1(i)) = i
       count(mark1(i)) = count(mark1(i)) + vwts(i)
     end do

     if (msglvl>1) then
       write (lp,'(A)') 'Number of vertices in each half-level set'
       do i = 1, maxl - minl + 1
         write (lp,'(2I10)') i, count(i)
       end do
     end if

     ! Run through half-level sets computing imbalances
     ! wtB is weight of B
     ! wtW is set to total weight of rest of network
     wtw = wtw + wts
     if (msglvl>1) write (lp,('(A,3I10)')) 'wtB,wtW,wtS', wtb, wtw, wts
     if (maxl-minl==0) then
       ! Only one level set
       imb(1) = max(real(wtb)/real(wtw),real(wtw)/real(wtb))
     else
       wtw = wtw - count(1)
      ! if (msglvl>0) write (14,'(A)') &
      !   'Half-level set   width   |B|,|W|, imbalance'
       do k = 1, maxl - minl
         wtw = wtw - count(k+1)
         imb(k) = max(real(wtb)/real(wtw),real(wtw)/real(wtb))
      !   if (msglvl>0) write (14,'(I10,4G12.2)') k, count(k), wtb, wtw, &
      !     imb(k)
         wtb = wtb + count(k)
       end do
     end if

     if (msglvl>1) then
       write (lp,'(A)') 'Imbalances'
       do i = 1, maxl - minl
         write (lp,'(I10,G12.2)') i, imb(i)
       end do
     end if

     ! Run through nodes in level set assigning penalty to them
     do k = 1, maxl - minl + 1
       inode = head(k)
       do
         if (inode==-1) exit
         if (msglvl>1) write (lp,('(A,2I4)')) 'level and node', k, inode
         if (k==1) then
           mark(inode) = floor(penp*imb(1)*vwts(inode))
         else
           if (k==maxl-minl+1) mark(inode) = floor(penp*imb(maxl-minl)*vwts &
             (inode))
           if (k>1 .and. k<maxl-minl+1) mark(inode) = floor(penp*min(imb( &
             k),imb(k-1))*vwts(inode))
         end if
         inode = list(inode)
       end do
     end do

     if (msglvl>1) then
       write (lp,'(A)') 'Computed penalties'
       do i = 1, nvtx
         write (lp,'(2I10)') i, mark(i)
       end do
     end if

     if (msglvl>0) write (lp,'(A/)') '### leaving findpenalty()'

   end subroutine findpenalty

   subroutine augmentpath(netw,delta,pred)

     ! Reset flows on augmenting path


     ! Input
     ! delta   -- increment flow
     ! pred    -- tree predecessor vector, size nnode
     ! msglvl  -- message level

     ! Input/output
     ! network -- network object

     integer, intent(in) :: delta
     type (network), intent(inout) :: netw

     integer iarc, lp, sink, source, v, w

     integer, intent(in) :: pred(:)


     lp = 6

     source = netw%source
     sink = netw%sink

     ! Should set an error flag
     if (delta<=0 .or. pred(sink)<=0) then
       write (lp,'(A,I4,A,I4)') 'ERROR : delta', delta, ', pred(sink) = ', &
         pred(sink)
       return
     end if


     ! work back from the sink resetting network flows

     w = sink
     ! while w ~= source
     do
       if (w==source) exit
       iarc = pred(w)
       if (netw%firsts(iarc)==w) then
         v = netw%seconds(iarc)
         netw%flows(iarc) = netw%flows(iarc) - delta
       else if (netw%seconds(iarc)==w) then
         v = netw%firsts(iarc)
         netw%flows(iarc) = netw%flows(iarc) + delta
       end if
       w = v
     end do

   end subroutine augmentpath

   subroutine findaugpath(netw,iarc_m,avail,pred,stats,list,tags,deltas)

     ! Find an augmenting path starting from arc iarc_m
     ! here we use a breadth first search (BFS),
     ! in findaugpath2() we use a depth first search (DFS)
     ! in findaugpath3() we will use a max distance
     ! from source to grow the tree to find a path


     ! input --
     ! netw -- network object
     ! iarc_m -- label for starting arc (u,v)


     ! output --
     ! avail -- if nonzero, available flow on augmenting path
     ! pred -- tree predecessor vector, size nnode
     ! source <-- pred^m(sink) <-- pred(pred(sink)) ...
     ! <-- pred(sink) <-- sink
     ! stats -- statistics
     ! stats(1) = # nodes visited in search
     ! stats(2) = # arcs visited


     ! working --
     ! list -- stack vector used for depth first search, size nnode
     ! tags -- mark vector, size nnode
     ! deltas -- increment flow vector, size nnode

     integer, intent(in) :: iarc_m
     integer, intent(out) :: avail, stats(2)
     type (network), intent(in) :: netw
     integer, intent(out) :: pred(:)

     integer iarc, last, lp, nnode, now, root, sink, source, v, w

     integer :: list(:), tags(:), deltas(:)

     integer narc, u, n_nodevisit, n_arcvisit

     lp = 6

     ! As input variable is intent(in) we set a local value that we will
     ! update
     iarc = iarc_m

     nnode = netw%nnode
     narc = netw%narc
     source = netw%source
     sink = netw%sink

     ! write(0,*) 'sink',sink

     stats(1) = 0
     stats(2) = 0

     list = 0

     ! Initial working storage and array pred

     tags = 0
     deltas = 0
     pred = 0


     ! check that (source,u) is an edge
     ! that is, that iarc as input is an edge from the source node

     u = netw%seconds(iarc)

     ! This will never be the case because all arcs iarc_m come from source
     ! node
     if (netw%firsts(iarc)/=source) then
       write (lp,'(A,I4,A)') 'u', u, 'is not adjacent to source'
       return
     end if


     ! check for available capacity

     avail = netw%capacities(iarc) - netw%flows(iarc)
     if (avail==0) return


     ! Find augmenting path using an alternating tree

     root = u
     now = 1
     last = 1
     list(1) = root
     ! tags is used to tell whether node has been visited on this attempt
     ! to find an augmenting path
     tags(root) = root
     tags(source) = root
     tags(sink) = -1
     pred(sink) = -1
     deltas(root) = avail
     pred(root) = iarc
     n_nodevisit = 0
     n_arcvisit = 0

     ! while now <= last
     do
       if (now>last) exit
       v = list(now)
       now = now + 1
       n_nodevisit = n_nodevisit + 1

       iarc = netw%outheads(v)
       ! while iarc ~= -1
       do
         ! Run through all edges emanating from v
         ! First is v^- to v^+
         if (iarc==-1) exit
         w = netw%seconds(iarc)
         n_arcvisit = n_arcvisit + 1

         if (tags(w)/=root) then
           ! Node w has not yet been visited

           if (netw%capacities(iarc)>netw%flows(iarc)) then
             avail = netw%capacities(iarc) - netw%flows(iarc)

             if (avail>deltas(v)) then
               avail = deltas(v)
             end if
             deltas(w) = avail
             pred(w) = iarc
             ! Flag w as being visited
             tags(w) = root

             if (w==sink) exit

             last = last + 1
             list(last) = w

           end if
         end if
         ! Go to next arc from v
         iarc = netw%nextout(iarc)
       end do
       if (w==sink) exit


       iarc = netw%inheads(v)
       ! while iarc ~= -1
       do
         ! Run through all edges coming in to v
         if (iarc==-1) exit
         w = netw%firsts(iarc)
         n_arcvisit = n_arcvisit + 1
         if (tags(w)/=root) then
           if (netw%flows(iarc)>0) then
             if (avail>netw%flows(iarc)) then
               avail = netw%flows(iarc)
             end if
             deltas(w) = avail
             pred(w) = iarc
             tags(w) = root
             last = last + 1
             list(last) = w
           end if
         end if
         iarc = netw%nextin(iarc)
       end do
       ! Don't think you can reach this statement
       if (w==sink) exit
     end do

     ! Flag to show augmenting path not found
     if (w/=sink) avail = 0

     stats(1) = n_nodevisit
     stats(2) = n_arcvisit


   end subroutine findaugpath

   subroutine evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,map,alpha,costf,stats, &
       stats10)
     ! Matlab call
     ! function stats = evalBSW ( A, map, alpha, beta, msglvl )

     ! stats = EVALBSW ( A, map, alpha, beta )

     ! input ---
     ! map[nvtx] -- map from vertices to region
     ! map[u] == 0 --> u in S
     ! map[u] == 1 --> u in B
     ! map[u] == 2 --> u in W
     ! alpha --- acceptability parameter
     ! beta  --- imbalance penalty parameter

     ! output --
     ! stats[1] -- weight of vertices in S
     ! stats[2] -- weight of vertices in B
     ! stats[3] -- weight of vertices in W
     ! stats[4] -- weight of edges in A_{S,S}
     ! stats[5] -- weight of edges in A_{S,B}
     ! stats[6] -- weight of edges in A_{S,W}
     ! stats[7] -- weight of edges in A_{B,B}
     ! stats[8] -- weight of edges in A_{B,W}
     ! stats[9] -- 1 if acceptable, 0 if not
     ! acceptable --> alpha*min(|B|,|W|) >= max(|B|,|W|)
     ! stats10 -- cost of partition
     ! cost = |S|*(1 + (beta*| |B| - |W| |)/(|B|+|S|+|W|)) ;

     ! created -- 12jan12, cca

     integer, intent(in) :: a_n
     integer, intent(in) :: a_ne
     integer, intent(in) :: map(:), a_ptr(:), a_row(:), a_weight(:)
     real(wp), intent(in) :: alpha
     integer, intent(in) :: costf
     integer, intent(out) :: stats(9)
     real(wp), intent(out) :: stats10
     integer minbw, maxbw, nss, nsb, nsw, nbb, nww, nvtx, ns, nb, nw
     integer j, j1, j2, jj, u, v
     real(wp) diffbw, beta
     logical :: imbal

     beta = 0.5_wp
     nvtx = a_n
     stats(1:9) = -1
     ns = 0
     nb = 0
     nw = 0
     do u = 1, nvtx
       if (map(u)==0) then
         ns = ns + a_weight(u)
       else if (map(u)==1) then
         nb = nb + a_weight(u)
       else if (map(u)==2) then
         nw = nw + a_weight(u)
       end if
     end do
     stats(1) = ns
     stats(2) = nb
     stats(3) = nw
     minbw = min(nb,nw)
     maxbw = max(nb,nw)
     diffbw = real(abs(nb-nw))/real(ns+nb+nw)
     if (.false.) then
       nss = 0
       nsb = 0
       nsw = 0
       nbb = 0
       nww = 0
       ! [rows, cols, ents] = find(A) ;
       ! nzA = length(rows) ;
       do j = 1, a_n
         j1 = a_ptr(j)
         if (j==a_n) then
           j2 = a_ne
         else
           j2 = a_ptr(j+1) - 1
         end if
         v = j
         do jj = j1, j2
           u = a_row(jj)
           ! v = cols(ii) ;
           if (map(u)==0) then
             if (map(v)==0) then
               nss = nss + 1
             else if (map(v)==1) then
               nsb = nsb + 1
             else if (map(v)==2) then
               nsw = nsw + 1
             end if
           else if (map(u)==1) then
             if (map(v)==1) then
               nbb = nbb + 1
             end if
           else if (map(u)==2) then
             if (map(v)==2) then
               nww = nww + 1
             end if
           end if
         end do
       end do
       stats(4) = nss
       stats(5) = nsb
       stats(6) = nsw
       stats(7) = nbb
       stats(8) = nww
     end if
     ! stats[9] -- 1 if acceptable, 0 if not
     ! acceptable --> alpha*min(|B|,|W|) >= max(|B|,|W|)
     ! stats10 -- cost of partition
     if (alpha*minbw>=maxbw) then
       stats(9) = 1
     else
       stats(9) = 0
     end if

     if (alpha > nb+nw+ns-2) then
         imbal = .false.
     else
         imbal = .true.
     end if

     call cost_function(nb,nw,ns,nb+nw+ns, &
       alpha,imbal,costf,stats10)

   end subroutine evalbsw

end module spral_nd
