module spral_nd
   use spral_nd_hamd
   use spral_nd_maxflow
   use spral_nd_partition
   use spral_nd_preprocess
   use spral_nd_refine
   use spral_nd_types
   use spral_nd_util
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
   integer :: i, j, k, l, ll, lwork, lirn
   integer :: a_n_curr, a_ne_curr, num_zero_row
   integer, dimension(:), allocatable :: amd_order_row, amd_order_ptr, &
      amd_order_sep, amd_order_perm, amd_order_work, amd_order_iperm
   integer, dimension(:), allocatable :: work_iperm, work_seps
   integer :: nsvar
   integer :: sumweight
   integer, dimension(:), allocatable :: svar, sinvp ! supervariable info
   integer, allocatable, dimension(:) :: a_weight ! a_weight(i) will
      ! contain the weight of variable (column) i ends for the
      ! expanded matrix
   integer, allocatable, dimension(:) :: iperm ! inverse of perm(:)
   integer, allocatable, dimension(:) :: work ! space for doing work
   logical :: use_multilevel
   type (nd_multigrid) :: grid

   call nd_print_diagnostic(1, options, ' ')
   call nd_print_diagnostic(1, options, 'nd_order:')

   !!!!!!!!!!!!!!!!!!!!!!
   ! Preprocessing
   !!!!!!!!!!!!!!!!!!!!!!

   ! Error checks
   if (n.lt.1) then
      info%flag = ND_ERR_N
      call nd_print_error(info%flag, options, 'nd_order')
      return
   end if

   ! Convert matrix to internal format without diagonals
   if (mtx.lt.1) then
      call construct_full_from_lower(n, ptr, row, a_n, a_ne, a_ptr, a_row, &
         options, info%stat)
   else
      call construct_full_from_full(n, ptr, row, a_n, a_ne, a_ptr, a_row, &
         options, info%stat)
   endif
   if(info%stat.ne.0) go to 10

   ! Output summary of input matrix post-conversion
   if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) then
      ! Print out a_ptr and a_row
      write (options%unit_diagnostics,'(a8)') 'a_ptr = '
      write (options%unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n)
      write (options%unit_diagnostics,'(a8)') 'a_row = '
      write (options%unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne)
   else if (options%print_level.ge.1.and.options%unit_diagnostics.gt.0) then
      ! Print out first few entries of a_ptr and a_row
      write (options%unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n)) = '
      write (options%unit_diagnostics,'(5i15)') (a_ptr(i),i=1,min(5,a_n))
      write (options%unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne)) = '
      write (options%unit_diagnostics,'(5i15)') (a_row(i),i=1,min(5,a_ne))
   end if

   ! Allocate iperm and initialize to identity
   allocate (iperm(a_n),stat=info%stat)
   if (info%stat.ne.0) go to 10
   iperm(:) = (/ (i,i=1,a_n) /)

   ! Initialize all variables to not be in a seperator
   ! NB: To make code easier to handle, we always have a work_seps(:) array
   !     rather than optional arguments throughout (i.e. treat it as always
   !     present).
   allocate(work_seps(n), stat=info%stat)
   if (info%stat.ne.0) go to 10
   work_seps(:) = -1
   if(present(seps)) seps(:) = -1 ! Cover anything not sent to ND routine

   ! Remove any dense rows from matrix and modify iperm (if enabled)
   if (options%remove_dense_rows) &
      call remove_dense_rows(a_n, a_ne, a_ptr, a_row, iperm, options, info)

   ! Return if matrix is diagonal
   if (a_ne.eq.0) then ! (recall diagonal entries are not stored)
      call nd_print_diagnostic(1, options, ' ')
      call nd_print_diagnostic(1, options, 'Matrix is diagonal')

      info%nsuper = a_n
      info%nzsuper = 0

      goto 100 ! Skip straight to exit
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

   !!!!!!!!!!!!!!!!!!!!!!
   ! Main ordering
   !!!!!!!!!!!!!!!!!!!!!!

   info%nzsuper = a_ne_curr
   info%nsuper = a_n_curr

   if (options%amd_switch2.le.0 .or. &
         a_n_curr.le.max(2,max(options%amd_call,options%amd_switch1)) ) then
      ! Apply AMD to matrix
      ! Allocate work to have length 5*a_n_curr+a_ne_curr
      call nd_print_diagnostic(2, options, ' Form AMD ordering')
      lirn = a_ne_curr + a_n_curr
      allocate(amd_order_ptr(a_n_curr), amd_order_row(lirn), &
         amd_order_sep(a_n_curr), amd_order_perm(a_n_curr), &
         amd_order_work(7*a_n_curr), amd_order_iperm(a_n))
      if (info%stat.ne.0) go to 10

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

      call nd_print_diagnostic(2, options, ' Form ND ordering')

      if (nsvar+num_zero_row.ne.a_n) then
         ! Create shadow versions of iperm and seps that work on
         ! supervariables rather than variables
         allocate(work_iperm(a_n), stat=info%stat)
         if (info%stat.ne.0) go to 10
         work_iperm(1:a_n_curr) = (/ (i,i=1,a_n_curr) /)
      end if
      ! Allocate a workspace that can be reused at lower levels
      allocate (work(a_n+14*a_n_curr+a_ne_curr), stat=info%stat)
      if (info%stat.ne.0) go to 10

      use_multilevel = .true.
      sumweight = sum(a_weight(1:a_n_curr))
      lwork = 12*a_n_curr + sumweight + a_ne_curr
      if (nsvar+num_zero_row.eq.a_n) then
         call nd_nested_internal(a_n_curr, a_ne_curr, a_ptr(1:a_n_curr), &
            a_row(1:a_ne_curr), a_weight(1:a_n_curr), sumweight,         &
            iperm(1:a_n_curr), work(1:lwork),                            &
            work(lwork+1:lwork+a_n_curr),                                &
            work(lwork+a_n_curr+1:lwork+2*a_n_curr), 0, options, info,   &
            .false., use_multilevel, grid, work_seps)
      else
         ! Supervariables: use work_iperm(:) instead of perm(:)
         call nd_nested_internal(a_n_curr, a_ne_curr, a_ptr(1:a_n_curr), &
            a_row(1:a_ne_curr), a_weight(1:a_n_curr), sumweight,         &
            work_iperm, work(1:lwork), work(lwork+1:lwork+a_n_curr),     &
            work(lwork+a_n_curr+1:lwork+2*a_n_curr), 0, options, info,   &
            .false., use_multilevel, grid, work_seps)
      end if
      if (info%flag.lt.0) return

      if (grid%level.eq.1) call mg_grid_destroy(grid,info%flag)

      if (nsvar+num_zero_row.eq.a_n) then
         if (present(seps)) then
            ! Reorder seps
            do i = 1, a_n_curr
               j = iperm(i)
               work(j) = work_seps(i)
            end do
            seps(1:a_n_curr) = work(1:a_n_curr)
         end if
      else
         if (present(seps)) then
            ! Expand and reorder seps
            do i = 1, a_n_curr
               j = work_iperm(i)
               if (j.eq.1) then
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
            if (j.eq.1) then
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

   !!!!!!!!!!!!!!!!!!!!!!
   ! Postprocessing
   !!!!!!!!!!!!!!!!!!!!!!
   100 continue

   ! Create perm from iperm
   do i = 1, n
      j = iperm(i)
      perm(j) = i
   end do

   ! Output summary of results
   if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) then
      ! Print out perm
      write (options%unit_diagnostics,'(a8)') 'perm = '
      write (options%unit_diagnostics,'(5i15)') (perm(i),i=1,n)
   else if (options%print_level.ge.1.and.options%unit_diagnostics.gt.0) then
      ! Print out first few entries of perm
      write (options%unit_diagnostics,'(a21)') 'perm(1:min(5,n)) = '
      write (options%unit_diagnostics,'(5i15)') (perm(i),i=1,min(5,n))
   end if

   info%flag = 0
   call nd_print_diagnostic(1, options, ' nd_order: successful completion')
   return

   10 continue
   info%flag = ND_ERR_MEMORY_ALLOC
   call nd_print_error(info%flag, options, 'nd_order')
   return
end subroutine nd_order

!
! Main (recursive) routine for performing nested dissection
!
recursive subroutine nd_nested_internal(a_n, a_ne, a_ptr, a_row, &
      a_weight, sumweight, iperm, work, work_comp_n, work_comp_nz, level, &
      options, info, use_amd, use_multilevel, grid, seps)
   ! NB: Many array arguments are chopped up and reused as the matrix is
   ! partitioned.
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, dimension(a_n), intent(inout) :: a_ptr(a_n)
   integer, dimension(a_ne), intent(inout) :: a_row(a_ne)
   integer, intent(inout) :: a_weight(a_n) ! Weight of each column
   integer, intent(in) :: sumweight ! sum(a_weight(:))
   integer, intent(inout) :: iperm(a_n) ! On input, iperm(i) contains the
      ! row in the original matrix (when nd_nested was called) that
      ! row i in this sub problem maps to. On output, this is updated to
      ! reflect the computed permutation.
   integer, intent(out) :: work_comp_n(a_n)
   integer, intent(out) :: work_comp_nz(a_n)
   integer, intent(out) :: work(12*a_n+sumweight+a_ne) ! Used during the
      ! algorithm to reduce need for allocations. The output is discarded.
   integer, intent(in) :: level ! which level of nested dissection is this
   type (nd_options), intent(in) :: options
   type (nd_inform), intent(inout) :: info
   logical, intent(in) :: use_amd
   logical, intent(inout) :: use_multilevel
   type (nd_multigrid), intent(inout) :: grid
   integer, intent(inout) :: seps(a_n)
      ! seps(i) is -1 if vertex i of permuted submatrix is not in a
      ! separator; otherwise it is equal to l, where l is the nested
      ! dissection level at which it became part of the separator

   integer :: i, j, k, l, m, s
   integer :: unit_diagnostics ! unit on which to print diagnostics
   integer :: num_components ! Number of independent components found
   integer :: compwork, lwork
   integer :: sumweight_sub, maxdeg_max_component
   integer :: offset_ptr, offset_row
   integer :: a_n1, a_n2, a_ne1, a_ne2
   integer :: amd_order_irn, amd_order_ip, amd_order_sep, amd_order_perm, &
      amd_order_work, lirn
   logical :: printi, printd, use_amdi

   ! ---------------------------------------------
   ! Printing levels
   unit_diagnostics = options%unit_diagnostics
   printi = (options%print_level.eq.1 .and. unit_diagnostics.ge.0)
   printd = (options%print_level.ge.2 .and. unit_diagnostics.ge.0)
   use_amdi = .false.

   if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) then
      write (options%unit_diagnostics,'(a)') ' '
      write (options%unit_diagnostics,'(a,i6)') &
         'Nested dissection level ', level
   end if

   ! Check whether matrix is diagonal and act accordingly
   if (a_ne.eq.0) then ! Recall we don't store diagonal entries
      call nd_print_diagnostic(1, options, ' ')
      call nd_print_diagnostic(1, options, 'Submatrix is diagonal')
      return
   end if

   if (level.eq.0) then
      maxdeg_max_component = a_ne + 1 - a_ptr(a_n)
      do i = 1, a_n - 1
         if (maxdeg_max_component.lt.a_ptr(i+1)-a_ptr(i)) &
            maxdeg_max_component = a_ptr(i+1) - a_ptr(i)
      end do
   end if

   ! Check whether to stop recursion
   ! (i.e. if either max levels reached or matrix is smaller than nd_switch)
   if (level.ge.options%amd_switch2 .or. &
         a_n.le.max(2,options%amd_switch1) .or. &
         use_amd) &
      go to 10

   ! Find a partition into (B, W, S)
   lwork = 12*a_n + sumweight + a_ne
   call nd_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, level, &
      a_n1, a_n2, a_ne1, a_ne2, iperm, work(1:lwork), options, info,      &
      use_multilevel, grid)

   if (a_n1.eq.a_n) go to 10 ! Failed to find a partition


   if (a_n1.ne.0 .and. a_n2.ne.0 .and. a_n1+a_n2.eq.a_n) then ! i.e. |S| = 0
      ! matrix is reducible
      call nd_print_diagnostic(1, options, ' ')
      call nd_print_diagnostic(1, options, 'Matrix is reducible')
      compwork = 0
      ! Whilst B and W are independent, there may be more than two independent
      ! components, we want to find them all.
      ! NB: work array needs to be total length 5*a_n+a_ne
      call nd_find_indep_comps(a_n, a_ne, a_ptr, a_row, a_weight, iperm, &
         num_components, work_comp_n(1:a_n), work_comp_nz(1:a_n), &
         work(compwork+1:compwork+3*a_n+a_ne), options)
      if (num_components.eq.1) then
         ! This should never happen.
         info%flag = ND_ERR_INTERNAL
         call nd_print_diagnostic(0, options, 'Only found one component?!?!?')
         call nd_print_error(info%flag, options, 'nd_nested_internal')
         return
      end if
      k = level
      if (level.eq.0) info%num_components = num_components

      ! Apply the ND to each component - do not test for indep comps
      offset_ptr = a_n + 1
      offset_row = a_ne + 1
      i = num_components
      ! work from last component so that space at end of work_comp_n and
      ! work_comp_nz can be reused without worrying about overwriting
      ! important data
      j = a_n
      do while (i.ge.1)
         if (level.eq.0) use_multilevel = .true.
         l = work_comp_n(i)
         use_amdi = (level.eq.0 .and. l .le. options%amd_call)
         m = work_comp_nz(i)
         s = sum(a_weight(offset_ptr-l:offset_ptr-1))
         if (m.gt.0) then
            ! Matrix not diagonal
            call nd_nested_internal(l, m, a_ptr(offset_ptr-l:offset_ptr-1), &
               a_row(offset_row-m:offset_row-1),                            &
               a_weight(offset_ptr-l:offset_ptr-1), s,                      &
               iperm(offset_ptr-l:offset_ptr-1),                            &
               work(compwork+1:compwork+12*l+s+m), work_comp_n(j-l+1:j),    &
               work_comp_nz(j-l+1:j), k, options, info, use_amdi,           &
               use_multilevel, grid, seps(offset_ptr-l:offset_ptr-1) )
            if(info%flag.lt.0) return
         end if
         offset_ptr = offset_ptr - l
         offset_row = offset_row - m
         j = j - l
         i = i - 1
      end do
      return
   end if

   ! Otherwise, S is non-empty
   if (level.eq.0 .and. a_n.gt.info%n_max_component) then
      info%n_max_component = a_n
      info%nz_max_component = a_ne
      info%maxdeg_max_component = maxdeg_max_component
   end if

   seps(a_n1+a_n2+1:a_n) = level
   if (a_n1.gt.max(2,options%amd_switch1)) then
      sumweight_sub = sum(a_weight(1:a_n1))
      call nd_nested_internal(a_n1, a_ne1, a_ptr(1:a_n1), a_row(1:a_ne1), &
         a_weight(1:a_n1), sumweight_sub, iperm(1:a_n1),                  &
         work(1:12*a_n1+sumweight_sub+a_ne1), work_comp_n(1:a_n1),        &
         work_comp_nz(1:a_n1), level+1, options, info, use_amdi,          &
         use_multilevel, grid, seps(1:a_n1))
      if(info%flag.lt.0) return
   end if

   if (a_n2.gt.max(2,options%amd_switch1)) then
      if (a_n1.gt.max(2,options%amd_switch1)) then
         sumweight_sub = sum(a_weight(a_n1+1:a_n1+a_n2))
         call nd_nested_internal(a_n2, a_ne2, a_ptr(a_n1+1:a_n1+a_n2),  &
            a_row(a_ne1+1:a_ne1+a_ne2), a_weight(a_n1+1:a_n1+a_n2),     &
            sumweight_sub, iperm(a_n1+1:a_n1+a_n2),                     &
            work(1:12*a_n2+sumweight_sub+a_ne2), work_comp_n(1:a_n2),   &
            work_comp_nz(1:a_n2), level+1, options, info, use_amdi,     &
            use_multilevel, grid, seps=seps(a_n1+1:a_n1+a_n2) )
      else
         sumweight_sub = sum(a_weight(a_n1+1:a_n1+a_n2))
         call nd_nested_internal(a_n2, a_ne2, a_ptr(1:a_n2), a_row(1:a_ne2),&
            a_weight(a_n1+1:a_n1+a_n2), sumweight_sub,                      &
            iperm(a_n1+1:a_n1+a_n2), work(1:12*a_n2+sumweight_sub+a_ne2),   &
            work_comp_n(1:a_n2), work_comp_nz(1:a_n2), level+1, options,    &
            info, use_amdi, use_multilevel, grid,                           &
            seps=seps(a_n1+1:a_n1+a_n2) )
      end if
      if(info%flag.lt.0) return
   end if
   go to 20
   return

   ! No partition found or max number of levels have been reached
   10 continue
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
   call amd_order(a_n, a_ne, lirn, work(amd_order_ip+1:amd_order_ip+a_n),  &
      work(amd_order_irn+1:amd_order_irn+lirn),                            &
      work(amd_order_sep+1:amd_order_sep+a_n),                             &
      work(amd_order_perm+1:amd_order_perm+a_n),                           &
      work(amd_order_work+1:amd_order_work+7*a_n) )

   ! Extract perm from amd_order to apply to iperm
   do i = 1, a_n
      j = work(amd_order_perm+i)
      work(amd_order_work+i) = iperm(j)
   end do
   iperm(1:a_n) = work(amd_order_work+1:amd_order_work+a_n)

   20 continue
   info%flag = 0
   if (printi .or. printd) then
      call nd_print_message(info%flag,unit_diagnostics,'nd_nested_internal')
   end if
   return
end subroutine nd_nested_internal

!
! Finds and forms independent components in a matrix
!
subroutine nd_find_indep_comps(a_n, a_ne, a_ptr, a_row, a_weight, iperm, &
      comp_num, compsizes, compnzs, work, options)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(inout) :: a_ptr(a_n) ! On entry, input matrix. On exit,
      ! a_ptr(1:compsizes(1)) contains component 1;
      ! a_ptr(compsizes(1)+1:compsizes(1)+compsizes(2)) contains  compontent 2;
      ! etc.
   integer, intent(inout) :: a_row(a_ne) ! On entry, input matrix. On exit,
      ! a_row(1:compnzs(1)) contains component 1;
      ! a_ptr(compnzs(1)+1:compnzs(1)+compnzs(2)) contains compontent 2; etc.
   integer, intent(inout) :: a_weight(a_n) ! On entry, input matrix. On exit
      ! permuted as appropriate.
   integer, intent(inout) :: iperm(a_n) ! On input, iperm(i) contains the
      ! row in the original matrix (when nd_nested was called) that
      ! row i in this sub problem maps to. On output, this is updated to
      ! reflect the computed permutation.
   integer, intent(out) :: comp_num ! number independent components found
   integer, intent(out) :: compsizes(a_n) ! compsizes(i) will contain the
      ! size of component i
   integer, intent(out) :: compnzs(a_n) ! compnzs(i) will contain the number of
      ! nonzeros in component i
   integer, target, intent(out) :: work(3*a_n+a_ne)
   type (nd_options), intent(in) :: options

   integer :: root
   integer :: front_start, front_end ! start and end of front list
   integer :: num_assigned ! number of cols that have been assigned to a
      ! component
   integer :: i, j, u, v, front_pos
   integer, dimension(:), pointer :: ptr_work, row_work, mask, comp_rows

   call nd_print_diagnostic(1, options, ' ')
   call nd_print_diagnostic(1, options, 'Start independent component search')

   mask        => work(      1 : 1*a_n)
   comp_rows   => work(1*a_n+1 : 2*a_n)
   ptr_work    => work(2*a_n+1 : 3*a_n)
   row_work    => work(3*a_n+1 : 3*a_n + a_ne)

   ! Initialise mask
   mask(:) = 0

   !
   ! From each vertex (candidate root) explore fully its connected component.
   ! Then try next vertex (if not already explored) and so forth until all
   ! vertices have been tried as candidate roots.
   !
   num_assigned = 0
   comp_num = 0
   do root = 1, a_n
      if (mask(root).ne.0) cycle ! Vertex already visited
      ! Initialize new component
      comp_num = comp_num + 1
      front_start = num_assigned + 1
      front_end = front_start
      comp_rows(front_start) = root
      compsizes(comp_num) = 1
      mask(root) = compsizes(comp_num)
      num_assigned = num_assigned + 1

      ! Breadth-first search
      do while (front_end-front_start.ge.0) ! while list is non-empty
         do i = front_start, front_end
            v = comp_rows(i)
            do j = a_ptr(v), nd_get_ptr(v+1, a_n, a_ne, a_ptr) - 1
               ! pick a neighbour
               u = a_row(j)
               if (mask(u).ne.0) cycle ! Already added
               ! found unmasked vertex
               compsizes(comp_num) = compsizes(comp_num) + 1
               num_assigned = num_assigned + 1
               ! mask this vertex
               mask(u) = compsizes(comp_num)
               ! add vertex to component
               comp_rows(num_assigned) = u
            end do
         end do
         front_start = front_end + 1
         front_end = num_assigned
      end do
   end do

   if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) then
      write (options%unit_diagnostics,'(a)') ' '
      write (options%unit_diagnostics,'(i5,a)') comp_num, ' components found'
   end if

   if (comp_num.le.1) return ! No reordering to do

   ! Reorder matrix into block diagonal form
   call extract_matrices(a_n, a_ne, a_ptr, a_row, comp_num, compsizes, &
      comp_rows, compnzs, ptr_work, row_work, mask)
   a_ptr(:) = ptr_work(:)
   a_row(:) = row_work(:)

   ! Reorder iperm and a_weight
   do front_pos = 1, a_n
      j = comp_rows(front_pos)
      comp_rows(front_pos) = iperm(j)
      mask(front_pos) = a_weight(j)
   end do
   iperm(:) = comp_rows(:)
   a_weight(:) = mask(:)

   call nd_print_diagnostic(1, options, &
      ' nd_find_indep_comps: successful completion' &
      )
end subroutine nd_find_indep_comps

!
! Partition the matrix and if one (or more) of the generated submatrices is
! small enough, apply halo amd
!
subroutine nd_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, level, &
      a_n1, a_n2, a_ne1, a_ne2, iperm, work, options, info, use_multilevel, &
      grid)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(inout) :: a_ptr(a_n) ! On entry, input matrix. On exit,
      ! modified so parts are contigous
   integer, intent(inout) :: a_row(a_ne) ! On entry, input matrix. On exit,
      ! modified so parts are contigous
   integer, intent(inout) :: a_weight(a_n) ! As above
   integer, intent(in) :: sumweight ! Sum entries in a_weight.
   integer, intent(in) :: level ! Current nested dissection level
   integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
   integer, intent(out) :: a_ne1, a_ne2 ! no. nonzeros in two submatrices
   integer, intent(inout) :: iperm(a_n) ! On input, iperm(i) contains the
      ! row in the original matrix (when nd_nested was called) that
      ! row i in this sub problem maps to. On output, this is updated to
      ! reflect the computed permutation.
   logical, intent(inout) :: use_multilevel
   integer, target, intent(out) :: work(12*a_n+sumweight+a_ne)
   type (nd_options), intent(in) :: options
   type (nd_inform), intent(inout) :: info
   type (nd_multigrid), intent(inout) :: grid

   integer, dimension(:), pointer :: partition, partition2
   integer :: work_ptr ! pointer into work array
   integer :: a_ptr_sub_ptr ! pointer into work array
   integer :: a_row_sub_ptr ! pointer into work array
   integer :: i, j, k
   integer :: lwork
   integer :: a_weight_1, a_weight_2, a_weight_sep
   integer :: ref_method, ref_options
   integer :: a_n1_new, a_n2_new, a_weight_1_new, a_weight_2_new, &
      a_weight_sep_new
   real(wp) :: balance, balance_tol, tau_best, tau, band, depth
   logical :: imbal

   ! Internal options on what to do
   integer, parameter :: REFINE_MAXFLOW   = 0, &
                         REFINE_TRIM      = 1, &
                         REFINE_EDGE      = 2

   call nd_print_diagnostic(1, options, ' ')
   call nd_print_diagnostic(1, options, 'Start finding a partition')

   balance_tol = max(1.0_wp, options%balance)

   ! If matrix is full, then don't partition
   if (a_ne.eq.a_n*(a_n-1)) then
      a_n1 = a_n
      a_n2 = 0
      a_ne1 = a_ne
      a_ne2 = 0
      info%flag = 0
      call nd_print_diagnostic(1, options, 'No partition found')
      call nd_print_diagnostic(1, options, &
         ' nd_partition: successful completion' &
         )
      return
   end if

   ! Find the partition
   partition => work(1:a_n); work_ptr = a_n
   select case(options%coarse_partition_method)
   case(:ND_PARTITION_HALF_LEVEL_SET)
      call nd_half_level_set(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,    &
         level, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, partition,  &
         work(work_ptr+1:work_ptr+9*a_n+sumweight), options, band, depth,     &
         use_multilevel, info%flag)
   case(ND_PARTITION_LEVEL_SET:)
      call nd_level_set(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,         &
         level, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, partition,  &
         work(work_ptr+1:work_ptr+9*a_n+sumweight), options, band, depth,     &
         use_multilevel, info%flag)
   end select
   if(info%flag.ne.0) return ! it's all gone horribly wrong
   if (use_multilevel) then
      lwork = 9*a_n + sumweight
      call multilevel_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,  &
         partition, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, options, &
         info%flag, lwork, work(work_ptr+1:work_ptr+lwork), grid)
   end if

   ! If S is empty, return and caller will handle as special case
   if (a_n1+a_n2.eq.a_n) return

   if (level.eq.0) then
      if (   a_n.gt.info%n_max_component .or. &
            (a_n.eq.info%n_max_component .and. band.gt.info%band)) then
         info%band = band
         info%depth = depth
      end if
   end if

   ! If one or other partition is empty, give up
   if(a_n1.eq.0 .or. a_n2.eq.0) then
      info%flag = 0
      call nd_print_diagnostic(1, options, 'No partition found')
      call nd_print_diagnostic(1, options, &
         ' nd_partition: successful completion' &
         )
      return
   endif

   if ( .not. use_multilevel) then
      select case(options%refinement)
      case(:ND_REFINE_TRIM_FM_BOTH)
         ref_options = ND_REFINE_TRIM_FM_BOTH
      case(ND_REFINE_AUTO:)
         ref_options = ND_REFINE_TRIM_FM_AUTO
      case default
         ref_options = options%refinement
      end select

      select case (ref_options)
      case (ND_REFINE_TRIM_FM_BOTH)
         ref_method = REFINE_TRIM
      case (ND_REFINE_TRIM_FM_SMALLER)
         ref_method = REFINE_EDGE
      case (ND_REFINE_TRIM_FM_AUTO)
         balance = max(a_weight_1,a_weight_2) / &
            real(min(a_weight_1,a_weight_2)+a_weight_sep)
         if (balance.ge.balance_tol) then
            ref_method = REFINE_EDGE
         else
            ref_method = REFINE_TRIM
         end if
      case (ND_REFINE_MAXFLOW_BOTH)
         ref_method = REFINE_MAXFLOW
      case (ND_REFINE_MAXFLOW_SMALLER)
         ref_method = REFINE_EDGE
      case (ND_REFINE_MAXFLOW_AUTO)
         balance = max(a_weight_1,a_weight_2) / &
            real(min(a_weight_1,a_weight_2)+a_weight_sep)
         if (balance.ge.balance_tol) then
            ref_method = REFINE_EDGE
         else
            ref_method = REFINE_MAXFLOW
         end if
      end select
      if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) then
         write (options%unit_diagnostics,'(a)') 'Partition before refinement'
         write (options%unit_diagnostics,'(a,i10,a,i10,a,i10)') &
            'a_n1=', a_n1, ',  a_n2=', a_n2, ',  a_n_sep=', a_n - a_n1 - a_n2
      end if

      select case (ref_method)
      case (REFINE_MAXFLOW)
         call nd_refine_max_flow(a_n, a_ne, a_ptr, a_row, a_weight, a_n1,  &
            a_n2, a_weight_1, a_weight_2, a_weight_sep, partition,         &
            work(work_ptr+1:work_ptr+8), options)
      case (REFINE_TRIM)
         if (min(a_weight_1,a_weight_2)+a_weight_sep.lt. &
               max(a_weight_1,a_weight_2)) then
            call nd_refine_block_trim(a_n, a_ne, a_ptr, a_row, a_weight,    &
               sumweight, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, &
               partition, work(work_ptr+1:work_ptr+5*a_n),options)
         else
            call nd_refine_trim(a_n, a_ne, a_ptr, a_row, a_weight,          &
               sumweight, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, &
               partition, work(work_ptr+1:work_ptr+3*a_n), options)
         end if
      case (REFINE_EDGE)
         call nd_refine_edge(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, &
            a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep,              &
            partition, work(work_ptr+1:work_ptr+3*a_n), options)
      end select

      if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) then
         write (options%unit_diagnostics,'(a)') 'Partition after refinement'
         write (options%unit_diagnostics,'(a,i10,a,i10,a,i10)') &
            'a_n1=', a_n1, ',  a_n2=', a_n2, ',  a_n_sep=', a_n - a_n1 - a_n2
      end if

      if (options%max_improve_cycles.gt.0) then
         imbal = (balance_tol.le.real(sumweight-2))
         call cost_function(a_weight_1, a_weight_2, a_weight_sep, sumweight,&
            balance_tol, imbal, options%cost_function, tau_best)
         a_n1_new = a_n1
         a_n2_new = a_n2
         a_weight_1_new = a_weight_1
         a_weight_2_new = a_weight_2
         a_weight_sep_new = a_weight_sep
      end if

      partition2 => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
      partition2(:) = partition(:)

      k = options%max_improve_cycles
      do i = 1, k
         call expand_partition(a_n, a_ne, a_ptr, a_row, a_weight, a_n1_new, &
            a_n2_new, a_weight_1_new, a_weight_2_new, a_weight_sep_new,     &
            partition2, work(work_ptr+1:work_ptr+5*a_n))

         if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) then
            write (options%unit_diagnostics,'(a)') &
               'Partition sizes after expansion'
            write (options%unit_diagnostics,'(a,i10,a,i10,a,i10)') &
               'a_n1=', a_n1_new, ',  a_n2=', a_n2_new, &
               ',  a_n_sep=', a_n - a_n1_new - a_n2_new
         end if

         select case (ref_options)
         case (ND_REFINE_TRIM_FM_AUTO)
            balance = max(a_weight_1_new,a_weight_2_new) / &
               real(min(a_weight_1_new,a_weight_2_new)+a_weight_sep_new)
            if (balance.gt.balance_tol) then
               ref_method = REFINE_EDGE
            else
               ref_method = REFINE_TRIM
            end if
         case (ND_REFINE_MAXFLOW_AUTO)
            balance = max(a_weight_1_new,a_weight_2_new) / &
               real(min(a_weight_1_new,a_weight_2_new)+a_weight_sep_new)
            if (balance.gt.balance_tol) then
               ref_method = REFINE_EDGE
            else
               ref_method = REFINE_MAXFLOW
            end if
         end select

         select case (ref_method)
         case (REFINE_MAXFLOW)
            call nd_refine_max_flow(a_n, a_ne, a_ptr, a_row, a_weight,     &
               a_n1_new, a_n2_new, a_weight_1_new, a_weight_2_new,         &
               a_weight_sep_new, partition2, work(work_ptr+1:work_ptr+8),  &
               options)
         case (REFINE_TRIM)
            if (min(a_weight_1,a_weight_2)+a_weight_sep.lt. &
                  max(a_weight_1,a_weight_2)) then
               call nd_refine_block_trim(a_n, a_ne, a_ptr, a_row, a_weight, &
                  sumweight, a_n1_new, a_n2_new, a_weight_1_new,            &
                  a_weight_2_new, a_weight_sep_new, partition2,             &
                  work(work_ptr+1:work_ptr+5*a_n), options)
            else
               call nd_refine_trim(a_n, a_ne, a_ptr, a_row, a_weight,       &
                  sumweight, a_n1_new, a_n2_new, a_weight_1_new,            &
                  a_weight_2_new, a_weight_sep_new,                         &
                  partition2, work(work_ptr+1:work_ptr+3*a_n), options)
            end if
         case (REFINE_EDGE)
            call nd_refine_edge(a_n, a_ne, a_ptr, a_row, a_weight,   &
               sumweight, a_n1_new, a_n2_new, a_weight_1_new,        &
               a_weight_2_new, a_weight_sep_new, partition2,         &
               work(work_ptr+1:work_ptr+3*a_n), options)
         end select

         if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) then
            write (options%unit_diagnostics,'(a)') &
               'Partition sizes after refinement'
            write (options%unit_diagnostics,'(a,i10,a,i10,a,i10)') &
               'a_n1=', a_n1_new, ',  a_n2=', a_n2_new, &
               ',  a_n_sep=', a_n - a_n1_new - a_n2_new
         end if

         call cost_function(a_weight_1_new, a_weight_2_new,       &
            a_weight_sep_new, sumweight, balance_tol, imbal,      &
            options%cost_function, tau)
         if (tau.ge.tau_best) exit ! No improvement, stop
         tau_best = tau
         partition(:) = partition2(:)
         a_n1 = a_n1_new
         a_n2 = a_n2_new
         a_weight_1 = a_weight_1_new
         a_weight_2 = a_weight_2_new
         a_weight_sep = a_weight_sep_new
      end do

      call nd_refine_fm(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, a_n1, &
         a_n2, a_weight_1, a_weight_2, a_weight_sep, partition,             &
         work(work_ptr+1:work_ptr+8*a_n+sumweight), options)

      nullify(partition2); work_ptr = work_ptr - a_n
   end if

   if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) then
      write (options%unit_diagnostics,'(a)') 'Partition found:'
      write (options%unit_diagnostics,'(a,i10,a,i10,a,i10)') &
         'a_n1=', a_n1, ', a_n2=', a_n2, ', a_nsep=', a_n - a_n1 - a_n2
   end if

   if (  a_n1.le.max(2,options%amd_switch1) .and. &
         a_n2.le.max(2,options%amd_switch1) ) then
      ! apply halo amd to submatrics
      call amd_order_both(a_n, a_ne, a_ptr, a_row, a_n1, a_n2, &
         partition, iperm, a_weight, work(work_ptr+1:work_ptr+12*a_n+a_ne))
      a_ne1 = 0
      a_ne2 = 0
   else if (a_n1.le.max(2,options%amd_switch1)) then
      ! apply halo amd to [A1, B1'; B1, I] using two levels
      ! return A2 and apply ND to it
      call amd_order_one(a_n, a_ne, a_ptr, a_row, a_n1, a_n2, partition, &
         iperm, a_weight, a_ne2, work(work_ptr+1:work_ptr+12*a_n+a_ne))
      a_ne1 = 0
   else if (a_n2.le.max(2,options%amd_switch1)) then
      ! apply halo amd to [A2, B2'; B2, I] using two levels
      ! return A1 and apply ND to it
      call amd_order_one(a_n, a_ne, a_ptr, a_row, a_n1, a_n2, partition, &
         iperm, a_weight, a_ne1, work(work_ptr+1:work_ptr+12*a_n+a_ne))
      a_ne2 = 0
   else
      ! return A1 and A2 and apply ND to them
      a_ptr_sub_ptr = work_ptr + a_n ! length a_n
      a_row_sub_ptr = a_ptr_sub_ptr + a_n ! length a_ne

      call extract_both_matrices(a_n, a_ne, a_ptr, a_row, a_n1, a_n2,   &
         partition, a_ne1, a_ne2,                                       &
         work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2), a_ne,           &
         work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne),                      &
         work(work_ptr+1:work_ptr+a_n))

      ! Copy extracted matrices
      a_ptr(1:a_n1+a_n2) = work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2)
      a_row(1:a_ne1+a_ne2) = work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne1+a_ne2)

      ! Update iperm and a_weight
      do i = 1, a_n
         j = partition(i)
         work(a_ptr_sub_ptr+i) = iperm(j)
         work(work_ptr+i) = a_weight(j)
      end do
      do i = 1, a_n
         iperm(i) = work(a_ptr_sub_ptr+i)
         a_weight(i) = work(work_ptr+i)
      end do
   end if

   info%flag = 0
   call nd_print_diagnostic(1, options, ' nd_partition: successful completion')
end subroutine nd_partition








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
    if (i.eq.a_n) then
      l = a_ne
    else
      l = a_ptr(i+1) - 1
    end if
    do k = a_ptr(i), l
      m = a_row(k)
      if (work(mask+m).ne.0) then
        a_ptr_sub(j) = a_ptr_sub(j) + 1

      end if
    end do
  end do
  do j = a_n_part + 1, a_n_sub
    a_ptr_sub(j) = 0
    i = rows_sub(j)
    if (i.eq.a_n) then
      l = a_ne
    else
      l = a_ptr(i+1) - 1
    end if
    do k = a_ptr(i), l
      m = a_row(k)
      if (work(mask+m).gt.0) a_ptr_sub(j) = a_ptr_sub(j) + 1
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
    if (i.eq.a_n) then
      l = a_ne
    else
      l = a_ptr(i+1) - 1
    end if
    do k = a_ptr(i), l
      m = a_row(k)
      if (work(mask+m).ne.0) then
        a_ptr_sub(j) = a_ptr_sub(j) - 1
        p = a_ptr_sub(j)
        a_row_sub(p) = abs(work(mask+m))
      end if
    end do
  end do
  do j = a_n_part + 1, a_n_sub
    i = rows_sub(j)
    if (i.eq.a_n) then
      l = a_ne
    else
      l = a_ptr(i+1) - 1
    end if
    do k = a_ptr(i), l
      m = a_row(k)
      if (work(mask+m).gt.0) then
        a_ptr_sub(j) = a_ptr_sub(j) - 1
        p = a_ptr_sub(j)
        a_row_sub(p) = work(mask+m)
      end if
    end do
  end do
end subroutine extract_matrix

subroutine extract_matrices(n, ne_in, ptr_in, row_in, nparts, n_sub, &
      rows_sub, ne_out, ptr_out, row_out, work)
   ! Matrix to be partitioned
   integer, intent(in) :: n
   integer, intent(in) :: ne_in
   integer, dimension(n), intent(in) :: ptr_in
   integer, dimension(ne_in), intent(in) :: row_in
   ! Details of partitions
   integer, intent(in) :: nparts ! number of partitions
   integer, dimension(nparts), intent(in) :: n_sub ! dimensions of partitions
   integer, dimension(*), intent(in) :: rows_sub ! rows in each partition
   integer, dimension(nparts), intent(out) :: ne_out ! #entries in each part
   integer, intent(out) :: ptr_out(*) ! First n_sub(1) are for first matrix...
   integer, intent(out) :: row_out(*)
   integer, intent(out) :: work(n)

   integer :: i, j, cp, col, ins, part_start
   integer :: part, rsptr

   ! Loop over parts, extracting matrix
   rsptr = 0
   ins = 1
   do part = 1, nparts
      ! Set work(col) to index of col within new submatrix, or 0 if absent
      work(:) = 0
      do cp = rsptr+1, rsptr + n_sub(part)
         col = rows_sub(cp)
         work(col) = cp - rsptr
      end do
      ! Initialize refence start for ptr
      part_start = ins - 1
      ! Copy submatrix of marked variables
      do cp = rsptr+1, rsptr + n_sub(part)
         ptr_out(cp) = ins - part_start
         col = rows_sub(cp)
         do i = ptr_in(col), nd_get_ptr(col+1, n, ne_in, ptr_in)-1
            j = work( row_in(i) )
            if(j.eq.0) cycle ! Variable not in this partition
            row_out(ins) = j
            ins = ins + 1
         end do
      end do
      ne_out(part) = ins - part_start - 1
      ! Update for next partition
      rsptr = rsptr + n_sub(part)
   end do
end subroutine extract_matrices

subroutine extract_both_matrices(a_n, a_ne, a_ptr, a_row, a_n_part1, &
      a_n_part2, rows_sub, a_ne_sub1, a_ne_sub2, a_ptr_sub, len_a_row_sub, &
      a_row_sub, work)
   integer, intent(in) :: a_n ! order of matrix being partitioned
   integer, intent(in) :: a_ne ! no. entries in matrix being partitioned
   integer, intent(in) :: a_ptr(a_n) ! col ptrs for matrix being  partitioned
   integer, intent(in) :: a_row(a_ne) ! row indices for matrix being
      ! partitioned.
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

   integer, dimension(2) :: n_array, ne_array

   n_array(1) = a_n_part1
   n_array(2) = a_n_part2
   call extract_matrices(a_n, a_ne, a_ptr, a_row, 2, n_array, rows_sub, &
      ne_array, a_ptr_sub, a_row_sub, work)
   a_ne_sub1 = ne_array(1)
   a_ne_sub2 = ne_array(2)
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
  if (a_n1.lt.a_n2) then
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
  if (a_n1.lt.a_n2) then
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
  if (a_n1.lt.a_n2) then

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

subroutine multilevel_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, &
    partition, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, options, &
    info1, lwork, work, grid)

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

  type (nd_multigrid), intent(inout) :: grid ! the multilevel of
  ! graphs (matrices)

  integer :: i, j, k, inv1, inv2, ins
  integer :: mp
  integer :: mglevel_cur ! current level
  integer :: err, print_level ! printing
  logical :: lerr

  info1 = 0
  ! Set up printing
  if (options%print_level.lt.0) print_level = 0
  ! The default is options%print_level = 0
  if (options%print_level.eq.0) print_level = 1
  if (options%print_level.eq.1) print_level = 2
  if (options%print_level.gt.1) print_level = 3
  mp = options%unit_diagnostics
  if (mp.lt.0) print_level = 0
  ! Set error optionss
  lerr = options%unit_error .ge. 0 .and. print_level .gt. 0
  err = options%unit_error


  if (print_level.gt.1) then
    write (mp,'(a)') 'Start multilevel_partition:'
  end if

  ! construct the multigrid at this level

  if ( .not. allocated(grid%graph)) allocate (grid%graph)

  call nd_matrix_construct(grid%graph,a_n,a_n,a_ne,info1)
  if (info1.lt.0) then
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
  if (info1.lt.0) then
    if (lerr) call nd_print_message(info1,err, &
      ' multilevel_partition')
    return
  end if

  call nd_assoc(grid%row_wgt,a_n,info1)
  if (info1.lt.0) then
    if (lerr) call nd_print_message(info1,err, &
      ' multilevel_partition')
    return
  end if

  ! Initialise row weights
  grid%row_wgt(1:a_n) = a_weight(1:a_n)

  ! initialise mglevel_cur to the maximum number of levels
  ! allowed for this bisection
  mglevel_cur = options%stop_coarsening2
  call multilevel(grid,options,sumweight,mglevel_cur,mp,print_level, &
    lwork,work,info1)

  if (info1.ne.0) then
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
  if (info1.ne.0) then
    if (lerr) call nd_print_message(info1,err, &
      ' multilevel_partition')
    return
  end if

  ! deallocate (matrix%ptr,matrix%col,matrix%val,stat=st)
  ! if (st.ne.0) info1 = ND_ERR_MEMORY_DEALLOC
  if (info1.lt.0) then
    if (lerr) call nd_print_message(info1,err, &
      ' multilevel_partition')
    return
  end if

  if (print_level.gt.2) then
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
  real(wp) :: tau, balance_tol, tau_best
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
  if (grid%level.ge.mglevel_cur .or. grid%size.le.stop_coarsening1) then
    if (print_level.ge.2) call level_print(mp,'end of level ',grid%level)

    ! coarsest level in multilevel so compute separator
    a_ne = grid%graph%ptr(grid%graph%n+1) - 1
    call nd_coarse_partition(grid%graph%n,a_ne,grid%graph%ptr, &
      grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
      grid%part_div(2),grid%where,lwork,work,options,info)
    return
  end if

  ! Coarsest level not yet reached so carry on coarsening
  if (options%matching.eq.1) then
    lwk = grid%size
    call coarsen_hec(grid,lwk,work(1:lwk),info)
  else
    if (options%matching.gt.1) then
      lwk = 3*grid%size
      call coarsen_best(grid,lwk,work(1:lwk),info)
    else
      lwk = 2*grid%size
      call coarsen_cn(grid,lwk,work(1:lwk))
    end if
  end if
  if (info.lt.0) return

  cgrid => grid%coarse
  cnvtx = cgrid%size
  ! allocate coarse grid quantities
  call nd_assoc(cgrid%where,cnvtx,info)
  if (info.ne.0) then
    return
  end if

  call nd_assoc(cgrid%row_wgt,cnvtx,info)
  if (info.ne.0) then
    return
  end if

  ! see if the grid reduction is achieved, if not, set the allowed
  ! maximum level to current level and partition this level
  ! deallocate the coarse grid quantities that haves been allocated so
  ! far
  if (real(cgrid%size)/real(grid%size).gt.grid_rdc_fac_max .or. &
      real(cgrid%size)/real(grid%size).lt.grid_rdc_fac_min .or. &
      cgrid%size.lt.4) then

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
    if (info.lt.0) return

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
  if (info.lt.0) return

  ! check if matrix is full
  if (real(cgrid%graph%ptr(cgrid%graph%n+1)-1)/real(cgrid%graph%n).ge.real &
      (cgrid%graph%n-1)) then
    if (print_level.ge.2) then
      write (mp,'(a,i10,a)') 'at level ', grid%level, &
        ' further coarsening gives full matrix'
    end if

    ! set current grid level and recurse
    mglevel_cur = grid%level - 1
    call multilevel(grid,options,sumweight,mglevel_cur,mp,print_level, &
      lwork,work,info)
    if (info.lt.0) return

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
  if (cgrid%part_div(1).eq.0 .or. cgrid%part_div(2).eq.0) then
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
    if (info.lt.0) return

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
    if (fwhere(i).eq.ND_PART1_FLAG) then
      grid%part_div(1) = grid%part_div(1) + 1
    else
      if (fwhere(i).eq.ND_PART2_FLAG) then
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

  if (a_weight_sep.gt.0) then
    ! Do not refine if separable graph

    if (options%refinement.gt.6) then
      ref_options = 3
    else
      if (options%refinement.lt.1) then
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
          a_weight_2)+a_weight_sep).gt.max(real(1.0, &
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
          a_weight_2)+a_weight_sep).gt.max(real(1.0, &
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
      if (min(a_weight_1,a_weight_2)+a_weight_sep.lt. &
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

    if (options%max_improve_cycles.gt.0) then
      balance_tol = max(1.0_wp,options%balance)
      imbal = (balance_tol.le.real(sumweight-2))
      call cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight, &
        balance_tol,imbal,options%cost_function,tau_best)
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
            a_weight_1_new,a_weight_2_new)+a_weight_sep_new).gt.max(real( &
            1.0,wp),options%balance)) then
          ref_method = 2
        else
          ref_method = 1
        end if

      case (6)
        if (real(max(a_weight_1_new,a_weight_2_new))/real(min( &
            a_weight_1_new,a_weight_2_new)+a_weight_sep_new).gt.max(real( &
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
        if (min(a_weight_1,a_weight_2)+a_weight_sep.lt. &
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
        sumweight,balance_tol,imbal,options%cost_function,tau)
      if (tau.lt.tau_best) then
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

  if (info.lt.0) return

  if (print_level.eq.3) call level_print(mp,' after post smoothing ', &
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

  ! ---------------------------------------------
  ! Printing levels
  unit_diagnostics = options%unit_diagnostics
  printi = (options%print_level.eq.1 .and. unit_diagnostics.ge.0)
  printd = (options%print_level.ge.2 .and. unit_diagnostics.ge.0)

  if (printi .or. printd) then
    write (unit_diagnostics,'(a)') ' '
    write (unit_diagnostics,'(a)') 'Start finding a coarse partition'
    write (unit_diagnostics,'(a,i10,a,i10)') 'a_n=', a_n, ', a_ne=', &
      a_ne
  end if

  ! Find the partition
  if (options%coarse_partition_method.le.1) then
    partition_method = 1
  else
    partition_method = 2
  end if

  partition_ptr = 0 ! length a_n
  work_ptr = partition_ptr + a_n ! max length needed 9*a_n+a_ne

  allocate (work1(a_n),stat=st)
  if (st.ne.0) then
    info = ND_ERR_MEMORY_ALLOC
    return
  end if


  select case (partition_method)
  case (1)
    ! Half level set method
    use_multilevel = .false.
    call nd_half_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
      a_n2,a_weight_1,a_weight_2,a_weight_sep, &
      work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
      options,dummy,dummy1,use_multilevel,info)
    if(info.ne.0) return ! it's all gone horribly wrong

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
      options,dummy,dummy1,use_multilevel, info)
    if(info.ne.0) return ! it's all gone horribly wrong

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



  if (a_n1.ne.0 .and. a_n2.ne.0 .and. a_n.ge.3) then
    if (a_n1+a_n2.lt.a_n) then
      ! Refine the partition
      if (options%refinement.gt.6) then
        ref_options = 3
      else
        if (options%refinement.lt.1) then
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
            a_weight_2)+a_weight_sep).gt.max(real(1.0, &
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
            a_weight_2)+a_weight_sep).gt.max(real(1.0, &
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
        if (min(a_weight_1,a_weight_2)+a_weight_sep.lt. &
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
  if (st.ne.0) then
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

    if (grid%level.ne.1) then

      call multigrid_deallocate(grid,info)

    else

      call multigrid_deallocate_first(grid,info)

    end if

  else

    if (grid%level.ne.1) then

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
  if (info.ne.0) then
    return
  end if

  call nd_matrix_destruct(grid%p,info)
  if (info.ne.0) then
    return
  end if



  call nd_matrix_destruct(grid%r,info)
  if (info.ne.0) then
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
  if (ierr.ne.0) then
    info = ND_ERR_MEMORY_DEALLOC
    return
  end if

  call nd_matrix_destruct(grid%r,ierr)
  if (ierr.ne.0) then
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
    if (ierr.ne.0) then
      info = ND_ERR_MEMORY_DEALLOC
      return
    end if
  end if

  deallocate (grid%where,grid%row_wgt,stat=ierr)
  if (ierr.ne.0) info = ND_ERR_MEMORY_DEALLOC

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
  if (ierr.ne.0) then
    info = ND_ERR_MEMORY_DEALLOC
    return
  end if

  deallocate (matrix%val,stat=ierr)
  if (present(stat)) stat = ierr
  if (ierr.ne.0) then
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
  if (info.lt.0) then
    return
  end if

  call nd_alloc(p%col,ne,info)
  if (info.lt.0) then
    return
  end if

  call nd_alloc(p%val,ne,info)
  if (info.lt.0) then
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
    if (SIZE(v).lt.n) then
      deallocate (v,stat=st)
      if (st.lt.0) then
        info = ND_ERR_MEMORY_ALLOC
        return
      end if
    else
      return
    end if
  end if

  allocate (v(n),stat=st)
  if (st.lt.0) then
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
   if (st.ne.0) info = ND_ERR_MEMORY_ALLOC

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
    if (work(ptr_match+v).ne.unmatched) cycle
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
      if (work(ptr_match+u).eq.unmatched .and. maxwgt.lt.abs(graph%val(j))) &
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
    if (maxind.ne.v) then
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
    if (work(ptr_match+i).eq.i) then
      r%ptr(k) = j
      r%col(j) = i
      j = j + 1
      k = k + 1
    else
      if (work(ptr_match+i).gt.i) then
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
    if (k.eq.i) then
      p%col(p%ptr(i)) = j
      j = j + 1
    else
      if (k.gt.i) then
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
    if (work(ptr_match+v).ne.unmatched) cycle
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
      if (work(ptr_match+u).ne.unmatched) cycle
      num = 0
      do k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
        w = grid%graph%col(k)
        if (work(ptr_flag+w).eq.i) num = num + 1
      end do
      if (num.gt.max_neigh) then
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
    if (maxind.ne.v) then
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
    if (work(ptr_match+i).eq.i) then
      r%ptr(k) = j
      r%col(j) = i
      j = j + 1
      k = k + 1
    else
      if (work(ptr_match+i).gt.i) then
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
    if (k.eq.i) then
      p%col(p%ptr(i)) = j
      j = j + 1
    else
      if (k.gt.i) then
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
    if (work(ptr_match+v).ne.unmatched) cycle
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
      if (work(ptr_match+u).eq.unmatched .and. maxwgt.lt.abs(graph%val(j))) &
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
    if (work(ptr_match1+v).ne.unmatched) cycle
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
      if (work(ptr_match1+u).ne.unmatched) cycle
      num = 0
      do k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
        w = grid%graph%col(k)
        if (work(ptr_flag+w).eq.i) num = num + 1
      end do
      if (num.gt.max_neigh) then
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
    if (maxind.ne.v) then
      nz = nz + 1
    end if
  end do

  ! --------------------------------------------------------------------
  ! -
  if (cnvtx.le.cnvtx1) then
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
    if (matching(i).eq.i) then
      r%ptr(k) = j
      r%col(j) = i
      j = j + 1
      k = k + 1
    else
      if (matching(i).gt.i) then
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
    if (k.eq.i) then
      p%col(p%ptr(i)) = j
      j = j + 1
    else
      if (k.gt.i) then
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
  ! if (info65.lt.0) then
  ! info = info65
  ! return
  ! end if
  nvtx = matrix%n
  cnvtx = p%n

  ! get the size of the coarse matrix first
  call galerkin_graph_rap_size(nvtx,cnvtx,nz,p%ptr(nvtx+1)-1,p%col, &
    p%ptr,matrix%ptr(nvtx+1)-1,matrix%col,matrix%ptr,r%ptr(cnvtx+1)-1, &
    r%col,r%ptr,lwork,work(1:lwork))
  if (info.lt.0) return

  call nd_matrix_construct(cmatrix,cnvtx,cnvtx,nz,info)
  if (info.lt.0) then
    return
  end if

  call galerkin_graph_rap(nvtx,cnvtx,p%ptr(nvtx+1)-1,p%val,p%col,p%ptr, &
    matrix%ptr(nvtx+1)-1,matrix%val,matrix%col,matrix%ptr, &
    r%ptr(cnvtx+1)-1,r%val,r%col,r%ptr,nz,cmatrix%val,cmatrix%col, &
    cmatrix%ptr,lwork,work(1:lwork))
  if (info.lt.0) return

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
        if (work(ptr_mask+neighneigh).ne.i) then
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
        if (work(ptr_mask+neighneigh).ne.-i .and. neighneigh.ne.i) then
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
        if (nzz.eq.0) then
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
        if (neighneigh.eq.i) cycle
        nzz = work(ptr_mask+neighneigh)
        if (nzz.eq.0) then
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


end module spral_nd
