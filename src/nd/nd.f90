module spral_nd
   use spral_nd_hamd
   use spral_nd_maxflow
   use spral_nd_multilevel
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


end module spral_nd
