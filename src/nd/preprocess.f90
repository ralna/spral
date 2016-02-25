module spral_nd_preprocess
   use spral_nd_numaware
   use spral_nd_types
   use spral_nd_util
   implicit none

   private
   public :: compress_by_svar,            & ! detect and use supervariables
             construct_full_from_full,    & ! lwr+upr CSC -> lwr+upr internal
             construct_full_from_lower,   & ! lwr CSC -> lwr+upr internal
             remove_dense_rows              ! drop dense rows, upd perm

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for conversion from external to internal format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Constructs a full matrix (without diagonals) from one with only lower
! triangle stored (perhaps with diagonals)
!
subroutine construct_full_from_lower(n, ptr, row, flags, n_out, ne_out, &
      ptr_out, row_out, flags_out, flags_diag_out, options, st)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(:), allocatable, intent(in) :: flags
   integer, intent(out) :: n_out
   integer, intent(out) :: ne_out
   integer, dimension(:), allocatable, intent(out) :: ptr_out
   integer, dimension(:), allocatable, intent(out) :: row_out
   integer, dimension(:), allocatable, intent(out) :: flags_out
   integer, dimension(:), allocatable, intent(out) :: flags_diag_out
   type (nd_options), intent(in) :: options
   integer, intent(out) :: st

   integer :: i, j, k

   n_out = n
   if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
      write (options%unit_diagnostics,'(a,i10)') 'n = ', n_out

   ! Allocate space to store pointers for expanded matrix
   allocate (ptr_out(n),stat=st)
   if (st.ne.0) return

   ! Set ptr_out(j) to hold no. nonzeros in column j, without diagonal
   ptr_out(:) = 0
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         if (j.ne.i) then
            ptr_out(i) = ptr_out(i) + 1
            ptr_out(j) = ptr_out(j) + 1
         end if
      end do
   end do

   ! Set ptr_out(j) to point to where row indices will end in row_out
   do j = 2, n
      ptr_out(j) = ptr_out(j-1) + ptr_out(j)
   end do
   ne_out = ptr_out(n)

   if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
      write (options%unit_diagnostics,'(a,i10)') &
         'entries in expanded matrix with diags removed = ', ne_out

   ! Allocate space to store row indices of expanded matrix
   allocate (row_out(ne_out), stat=st)
   if (st.ne.0) return
   if (allocated(flags)) then
      allocate (flags_out(ne_out), flags_diag_out(n), stat=st)
      if (st.ne.0) return
      flags_diag_out(:) = FLAG_SMALL
   endif

   ! Fill row_out and ptr_out
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         if (j.ne.i) then
            row_out(ptr_out(i)) = j
            row_out(ptr_out(j)) = i
            if(allocated(flags)) then
               flags_out(ptr_out(j)) = flags(k)
               flags_out(ptr_out(i)) = flag_transpose(flags(k))
            endif
            ptr_out(i) = ptr_out(i) - 1
            ptr_out(j) = ptr_out(j) - 1
         else if(allocated(flags)) then
            flags_diag_out(i) = flags(k)
         end if
      end do
   end do

   ! Reset ptr_out to point to where column starts
   do j = 1, n
      ptr_out(j) = ptr_out(j) + 1
   end do
end subroutine construct_full_from_lower

!
! Constructs a new full matrix (without diagonals) in internal CSC format
! from user supplied matrix in standard CSC format (which may have diagonals)
!
subroutine construct_full_from_full(n, ptr, row, flags, n_out, ne_out, ptr_out,&
      row_out, flags_out, flags_diag_out, options, st)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(:), allocatable, intent(in) :: flags
   integer, intent(out) :: n_out
   integer, intent(out) :: ne_out
   integer, dimension(:), allocatable, intent(out) :: ptr_out
   integer, dimension(:), allocatable, intent(out) :: row_out
   integer, dimension(:), allocatable, intent(out) :: flags_out
   integer, dimension(:), allocatable, intent(out) :: flags_diag_out
   type (nd_options), intent(in) :: options
   integer, intent(out) :: st

   integer :: i, j, k, p
   integer :: ndiags

   ! Set the dimension of the expanded matrix
   n_out = n
   if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
      write (options%unit_diagnostics,'(a,i10)') 'n = ', n

   ! Work out how many diagonal entries need removing
   ndiags = 0
   do i = 1, n
      do j = ptr(i), ptr(i+1) - 1
         k = row(j)
         if (k.eq.i) ndiags = ndiags + 1
      end do
   end do
   ne_out = ptr(n+1) - 1 - ndiags

   ! Allocate space to store pointers and rows for expanded matrix
   allocate (ptr_out(n), row_out(ne_out), stat=st)
   if (st.ne.0) return
   if (allocated(flags)) then
      allocate (flags_out(ne_out), flags_diag_out(n), stat=st)
      if (st.ne.0) return
      flags_diag_out(:) = FLAG_SMALL
   endif

   if (ndiags.eq.0) then
      ! No diagonal entries so do direct copy
      ptr_out(1:n) = ptr(1:n)
      row_out(1:ne_out) = row(1:ne_out)
      flags_out(1:ne_out) = flags(1:ne_out)
   else
      ! Diagonal entries present
      k = 1
      do i = 1, n
         ptr_out(i) = k
         do p = ptr(i), ptr(i+1) - 1
            j = row(p)
            if (i.ne.j) then
               row_out(k) = j
               k = k + 1
               if(allocated(flags)) flags_out(k) = flags(p)
            else if(allocated(flags)) then
               flags_diag_out(i) = flags(p)
            end if
         end do
      end do
   end if
end subroutine construct_full_from_full

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for detecting and removing dense rows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Identifies and removes dense rows in place. Updates iperm to place dense
! rows at the back of the permutation.
!
subroutine remove_dense_rows(a_n, a_ne, a_ptr, a_row, iperm, options, info, &
      a_flags, a_flags_diag, a_match)
   integer, intent(inout) :: a_n
   integer, intent(inout) :: a_ne
   integer, dimension(a_n), intent(inout) :: a_ptr
   integer, dimension(a_ne), intent(inout) :: a_row
   integer, dimension(a_n), intent(inout) :: iperm
   type (nd_options), intent(in) :: options
   type (nd_inform), intent(inout) :: info
   integer, dimension(a_ne), optional, intent(inout) :: a_flags
   integer, dimension(a_n), optional, intent(inout) :: a_flags_diag
   integer, dimension(a_n), optional, intent(inout) :: a_match

   integer, dimension(:), allocatable :: deg, prev, next, dense
   integer :: ndense ! number of dense rows found
   integer :: max_deg ! maximum degree
   integer :: degree, i, j, k, l, l1, l2, m, m1
   integer :: a_n_out ! dimension of subproblem after dense rows removed
   integer :: a_ne_out ! no. nonzeros of subproblem after dense rows removed

   call nd_print_diagnostic(1, options, ' ')
   call nd_print_diagnostic(1, options, 'Find and remove dense rows')

   ! Set pointers into work array
   allocate(deg(a_n), prev(a_n), next(a_n), dense(a_n), &
      stat=info%stat)
   if(info%stat.ne.0) then
      info%flag = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info%flag, options, 'remove_dense_rows')
      return
   end if

   ! By the end of this loop dense(i) will be
   !   0 if row is not dense
   !  <0 otherwise.
   ! The larger the number, the earlier the row was determined to be dense.
   ndense = 0
   max_deg = 0
   deg(:) = 0

   ! Calculate degree of each row before anything removed
   do i = 1, a_n
      k = a_ptr(i)
      if (i.lt.a_n) then
         degree = a_ptr(i+1) - k
      else
         degree = a_ne - a_ptr(a_n) + 1
      end if
      dense(i) = degree
      if (degree.ne.0) then
         max_deg = max(max_deg,degree)
         call dense_add_to_list(a_n, next, prev, deg, i, degree)
      end if
   end do
   degree = max_deg
   a_n_out = a_n
   a_ne_out = a_ne

   do while (degree - real(a_ne_out)/a_n_out .ge. &
         40*(real(a_n_out-1)/a_n_out)*LOG(real(a_n_out)) .and. degree.gt.0)
      i = deg(degree)
      if(dense(i).ge.0) then
         ! Not already matched with a dense row we've removed
         ndense = ndense + 1
         dense(i) = -ndense
      endif
      if(present(a_flags)) then
         j = a_match(i)
         if(j.ne.-1 .and. j.ne.i) then
            ! Mark partner as dense too - but don't remove from lists
            ndense = ndense + 1
            dense(j) = -ndense
         endif
      endif
      call dense_remove_from_list(a_n, next, prev, deg, i, degree)
      ! update degrees of adjacent vertices
      do k = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
         j = a_row(k)
         if (dense(j).gt.0) then
            call dense_remove_from_list(a_n, next, prev, deg, j, dense(j))
            dense(j) = dense(j) - 1
            if (dense(j).gt.0) then
               call dense_add_to_list(a_n, next, prev, deg, j, dense(j))
            end if
         end if
      end do
      a_n_out = a_n_out - 1
      a_ne_out = a_ne_out - 2*degree
      if (deg(degree).eq.0) then
         ! Find next largest degree
         degree = degree - 1
         do
            if (degree.eq.0) exit
            if (deg(degree).gt.0) exit
            degree = degree - 1
         end do
      end if
   end do

   ! By the end of this loop dense(i) will be
   ! >=0 if row is not dense
   ! <0 otherwise. The larger the number, the earlier the row was
   ! was determined to be dense.
   ! !!!!

   if (ndense.gt.0) then
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) then
         write (options%unit_diagnostics,'(a)') ' '
         write (options%unit_diagnostics,'(i10,a)') &
            ndense, ' dense rows detected'
      end if
      info%dense = ndense

      a_n_out = 0
      j = a_n + 1
      do i = 1, a_n
         k = dense(i)
         if (k.ge.0) then
            a_n_out = a_n_out + 1
            dense(i) = a_n_out
            next(a_n_out) = i
         else
            next(j+k) = i
         end if
      end do

      k = 1
      j = 1

      do i = 1, a_n
         l1 = a_ptr(i)
         l2 = nd_get_ptr(i+1, a_n, a_ne, a_ptr) - 1
         if (dense(i).ge.0) then
            a_ptr(j) = k
            if(present(a_flags)) then
               a_flags_diag(j) = a_flags_diag(i)
               l = a_match(i)
               if(l.eq.-1) then
                  a_match(j) = -1
               else
                  if(dense(l).lt.1) then
                     a_match(j) = -1
                  else
                     a_match(j) = dense(l)
                  end if
               end if
            end if
            do l = l1, l2
               m = a_row(l)
               m1 = dense(m)
               if (m1.ge.0) then
                  a_row(k) = m1
                  if(present(a_flags)) &
                     a_flags(k) = a_flags(l)
                  k = k + 1
               end if
            end do
            j = j + 1
         end if
      end do
      a_ptr(j) = k

      if(present(a_flags)) &
         call make_a_match_symmetric(a_n_out, a_match(1:a_n_out))

      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) then
         write (options%unit_diagnostics,'(a11)') 'a_n_out = '
         write (options%unit_diagnostics,'(i15)') a_n_out
         write (options%unit_diagnostics,'(a11)') 'a_ne_out = '
         write (options%unit_diagnostics,'(i15)') a_ne_out
         if (options%print_level.ge.2) then
            ! Print out a_ptr and a_row in full
            write (options%unit_diagnostics,'(a8)') 'a_ptr = '
            write (options%unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n_out)
            write (options%unit_diagnostics,'(a8)') 'a_row = '
            write (options%unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne_out)
         else
            ! Print out first few entries of a_ptr and a_row
            write (options%unit_diagnostics,'(a21)') &
               'a_ptr(1:min(5,a_n_out)) = '
            write (options%unit_diagnostics,'(5i15)') a_ptr(1:min(5,a_n_out))
            write (options%unit_diagnostics,'(a21)') &
               'a_row(1:min(5,a_ne_out)) = '
            write (options%unit_diagnostics,'(5i15)') a_row(1:min(5,a_ne_out))
         endif
      end if
   else

      a_n_out = a_n
      a_ne_out = a_ne
      next(:) = (/ (i,i=1,a_n) /)
   end if

   do i = 1, a_n
      j = next(i)
      next(i) = iperm(j)
   end do

   do i = 1, a_n
      iperm(i) = next(i)
   end do

   info%flag = 0
   call nd_print_diagnostic(1, options, &
      'remove_dense_rows successful completion' &
      )
   if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) &
      write (options%unit_diagnostics,'(a,i10)') &
         ' No. dense rows removed: ', a_n - a_n_out

   a_n = a_n_out
   a_ne = a_ne_out

end subroutine remove_dense_rows

subroutine dense_remove_from_list(n, next, prev, deg, irm, ig)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: next
   integer, dimension(n), intent(inout) :: prev
   integer, dimension(n), intent(inout) :: deg
   integer, intent(in) :: irm, ig
   integer :: inext, ilast

   inext = next(irm)
   ilast = prev(irm)
   if (ilast.eq.0) then
      deg(ig) = inext
      if (inext.ne.0) prev(inext) = 0
   else
      next(ilast) = inext
      if (inext.ne.0) prev(inext) = ilast
   end if
end subroutine dense_remove_from_list

subroutine dense_add_to_list(n, next, prev, deg, irm, ig)
   integer, intent(in) :: n
   integer, intent(inout) :: next(n)
   integer, intent(inout) :: prev(n)
   integer, intent(inout) :: deg(n)
   integer, intent(in) :: irm, ig
   integer :: inext

   inext = deg(ig)
   deg(ig) = irm
   next(irm) = inext
   if (inext.ne.0) prev(inext) = irm
   prev(irm) = 0
end subroutine dense_add_to_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for detecting and exploiting supervariables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Detects supervariables and compresses graph in place
!
subroutine compress_by_svar(a_n, a_ne, a_ptr, a_row, a_weight, a_n_curr, &
      a_ne_curr, nsvar, svar, sinvp, num_zero_row, options, st)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, dimension(a_n), intent(inout) :: a_ptr
   integer, dimension(a_ne), intent(inout) :: a_row
   integer, dimension(a_n), intent(out) :: a_weight
   integer, intent(out) :: a_n_curr
   integer, intent(out) :: a_ne_curr
   integer, intent(out) :: nsvar
   integer, dimension(a_n), intent(out) :: svar
   integer, dimension(a_n), intent(out) :: sinvp
   integer, intent(out) :: num_zero_row
   type(nd_options) :: options
   integer, intent(out) :: st

   integer :: i, j, k
   integer :: nnz_rows ! number of non-zero rows
   integer, allocatable, dimension(:) :: ptr2, row2, perm

   allocate (ptr2(a_n+1), row2(a_ne), perm(a_n), stat=st)
   if (st.ne.0) return

   ! Construct simple identity permutation
   perm(:) = (/ (i,i=1,a_n) /)
   sinvp(:) = (/ (i,i=1,a_n) /)

   ! Identify supervariables
   nnz_rows = a_n
   call nd_supervars(nnz_rows, a_ne, a_ptr, a_row, perm, sinvp, nsvar, &
      svar, st)
   if (st.ne.0) return

   num_zero_row = a_n - nnz_rows
   if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) &
      write (options%unit_diagnostics,'(a,i10)') &
         'Number supervariables: ', nsvar + num_zero_row

   ! If there are no supervariables, don't bother compressing: return
   if (nsvar+num_zero_row.eq.a_n) then
      a_n_curr = a_n
      a_ne_curr = a_ne
      a_weight(:) = 1
      return
   end if

   ! Otherwise, compress the matrix
   call nd_compress_by_svar(a_n, a_ne, a_ptr, a_row, sinvp, nsvar, &
      svar, ptr2, row2, st)
   if (st.ne.0) return

   ! FIXME: what is happening below? can it be simplified?
   a_n_curr = nsvar
   if(options%reord.eq.1) then
      a_ne_curr = ptr2(nsvar+1)-1
      a_ptr(1:a_n_curr) = ptr2(1:a_n_curr)
      a_row(1:a_ne_curr) = row2(1:a_ne_curr)
   else
      ! Fill a_ptr removing any diagonal entries
      a_ptr(:) = 0

      ! Set a_ptr(j) to hold no. nonzeros in column j
      do j = 1, a_n_curr
         do k = ptr2(j), ptr2(j+1) - 1
            i = row2(k)
            if (j.lt.i) then
               a_ptr(i) = a_ptr(i) + 1
               a_ptr(j) = a_ptr(j) + 1
            end if
         end do
      end do

      ! Set a_ptr(j) to point to where row indices will end in a_row
      do j = 2, a_n_curr
         a_ptr(j) = a_ptr(j-1) + a_ptr(j)
      end do
      a_ne_curr = a_ptr(a_n_curr)
      ! Initialise all of a_row to 0
      a_row(1:a_ne_curr) = 0

      ! Fill a_row and a_ptr
      do j = 1, a_n_curr
         do k = ptr2(j), ptr2(j+1) - 1
            i = row2(k)
            if (j.lt.i) then
               a_row(a_ptr(i)) = j
               a_row(a_ptr(j)) = i
               a_ptr(i) = a_ptr(i) - 1
               a_ptr(j) = a_ptr(j) - 1
            end if
         end do
      end do

      ! Reset a_ptr to point to where column starts
      do j = 1, a_n_curr
         a_ptr(j) = a_ptr(j) + 1
      end do
   endif

   ! Initialise a_weight
   a_weight(1:a_n_curr) = svar(1:a_n_curr)
   a_weight(a_n_curr+1:a_n_curr+num_zero_row) = 1

   ! Add zero rows/cols to matrix
   a_ptr(a_n_curr+1:a_n_curr+num_zero_row) = a_ne_curr + 1
   a_n_curr = a_n_curr + num_zero_row

   ! set svar(:) such that svar(i) points to the end of the list of variables
   ! in sinvp for supervariable i
   do i = 2, nsvar
      svar(i) = svar(i) + svar(i-1)
   end do
   j = svar(nsvar)
   do i = 1, num_zero_row
      svar(nsvar+i) = j + 1
      j = j + 1
   end do
end subroutine compress_by_svar

!
! Find supervariables if any exist. Modifes n to reflect number of non-zero
! columns.
!
! This is a modified version of mc78_supervars().
!
subroutine nd_supervars(n,ne,ptr,row,perm,invp,nsvar,svar,st)
   integer, intent(inout) :: n ! Dimension of system. On exit num non-zero cols
   integer, intent(in) :: ne ! Number of entries
   integer, dimension(n), intent(in) :: ptr ! Column pointers
   integer, dimension(ne), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence.
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nsvar ! number of supervariables
   integer, dimension(n), intent(out) :: svar ! number of vars in each svar
   integer, intent(out) :: st

   logical :: full_rank ! flags if supervariable 1 has ever become empty.
      ! If it has not, then the varaibles in s.v. 1 are those that never occur
   integer :: i
   integer(long) :: ii
   integer :: j
   integer :: idx ! current index
   integer :: next_sv ! head of free sv linked list
   integer :: nsv ! new supervariable to move j to
   integer :: piv ! current pivot
   integer :: col ! current column of A
   integer :: sv ! current supervariable
   integer :: svc ! temporary holding supervariable count
   integer, dimension(:), allocatable :: sv_new ! Maps each supervariable to
      ! a new supervariable with which it is associated.
   integer, dimension(:), allocatable :: sv_seen ! Flags whether svariables have
      ! been seen in the current column. sv_seen(j) is set to col when svar j
      ! has been encountered.
   integer, dimension(:), allocatable :: sv_count ! number of variables in sv.

   allocate(sv_new(n+1), sv_seen(n+1), sv_count(n+1), stat=st)
   if (st.ne.0) return

   svar(:) = 1
   sv_count(1) = n
   sv_seen(1) = 0

   ! Setup linked list of free super variables
   next_sv = 2
   do i = 2, n
      sv_seen(i) = i + 1
   end do
   sv_seen(n+1) = -1

   ! Determine supervariables using modified Duff and Reid algorithm
   full_rank = .false.
   do col = 1, n
      if (nd_get_ptr(col+1,n,ne,ptr).ne.ptr(col)) then
         ! If column is not empty, add implicit diagonal entry
         j = col
         sv = svar(j)
         if (sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! MUST BE the first time that sv has been seen for this
            ! column, so just leave j in sv, and go to next variable.
            ! (Also there can be no other vars in this block pivot)
         else
            ! There is at least one other variable remaining in sv
            ! MUST BE first occurence of sv in the current row/column,
            ! so define a new supervariable and associate it with sv.
            sv_seen(sv) = col
            sv_new(sv) = next_sv
            nsv = next_sv
            next_sv = sv_seen(next_sv)
            sv_new(nsv) = nsv ! avoids problems with duplicates
            sv_seen(nsv) = col
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = 1
            ! This sv cannot be empty as initial sv_count was > 1
         end if
      end if
      do ii = ptr(col), nd_get_ptr(col+1, n, ne, ptr) - 1
         j = row(ii)
         sv = svar(j)
         if (sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! If so, and this is first time that sv has been seen for this
            ! column, then we can just leave j in sv, and go to next variable.
            if (sv_seen(sv).lt.col) cycle
            ! Otherwise, we have already defined a new supervariable associated
            ! with sv. Move j to this variable, then retire (now empty) sv.
            nsv = sv_new(sv)
            if (sv.eq.nsv) cycle
            svar(j) = nsv
            sv_count(nsv) = sv_count(nsv) + 1
            ! Old sv is now empty, add it to top of free stack
            sv_seen(sv) = next_sv
            next_sv = sv
         else
            ! There is at least one other variable remaining in sv
            if (sv_seen(sv).lt.col) then
               ! this is the first occurence of sv in the current row/column,
               ! so define a new supervariable and associate it with sv.
               sv_seen(sv) = col
               sv_new(sv) = next_sv
               sv_new(next_sv) = next_sv ! avoids problems with duplicates
               next_sv = sv_seen(next_sv)
               sv_count(sv_new(sv)) = 0
               sv_seen(sv_new(sv)) = col
            end if
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = sv_count(nsv) + 1
            ! This sv cannot be empty as sv_count was > 1
         end if
      end do
   end do


   ! Now modify pivot order such that all variables in each supervariable are
   ! consecutive. Do so by iterating over pivots in elimination order. If a
   ! pivot has not already been listed, then order that pivot followed by
   ! any other pivots in that supervariable.

   ! We will build a new inverse permutation in invp, and then find perm
   ! afterwards. First copy invp to perm:
   perm(:) = invp(:)
   ! Next we iterate over the pivots that have not been ordered already
   ! Note: as we begin, all entries of sv_seen are less than or equal to n+1
   ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable has been
   ! ordered.
   idx = 1
   nsvar = 0
   do piv = 1, n
      if (sv_seen(piv).gt.n+1) cycle ! already ordered
      ! Record information for supervariable
      sv = svar(perm(piv))
      if ( .not. full_rank .and. sv.eq.1) cycle ! Don't touch unused vars
      nsvar = nsvar + 1
      svc = sv_count(sv)
      sv_new(nsvar) = svc ! store # vars in s.v. to copy to svar later
      j = piv
      ! Find all variables that are members of sv and order them.
      do while (svc.gt.0)
         do j = j, n
            if (svar(perm(j)).eq.sv) exit
         end do
         sv_seen(j) = n + 2 ! flag as ordered
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
   end do
   ! Push unused variables to end - these are those vars still in s.v. 1
   if ( .not. full_rank) then
      svc = sv_count(1)
      ! Find all variables that are members of sv and order them.
      j = 1
      do while (svc.gt.0)
         do j = j, n
            if (svar(perm(j)).eq.1) exit
         end do
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
      n = n - sv_count(1)
   end if

   ! Recover perm as inverse of invp
   do piv = 1, n
      perm(invp(piv)) = piv
   end do
   ! sv_new has been used to store number of variables in each svar, copy
   ! into
   ! svar where it is returned.
   svar(1:nsvar) = sv_new(1:nsvar)
end subroutine nd_supervars


!
! This subroutine takes a set of supervariables and compresses the supplied
! matrix using them. It avoids adding diagonal entries as ND code doesn't need
! them.
!
subroutine nd_compress_by_svar(n, ne, ptr, row, invp, nsvar, svar, ptr2, row2, &
      st)
   integer, intent(in) :: n ! Dimension of system
   integer, intent(in) :: ne ! Number off-diagonal zeros in system
   integer, dimension(n), intent(in) :: ptr ! Column pointers
   integer, dimension(ne), intent(in) :: row ! Row indices
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar ! super variables of A
   integer, dimension(nsvar+1), intent(out) :: ptr2
   integer, dimension(ne), intent(out) :: row2
   integer, intent(out) :: st

   integer :: piv, svc, sv, col
   integer :: j, idx
   integer, dimension(:), allocatable :: flag, sv_map

   allocate(flag(nsvar),sv_map(n),stat=st)
   if (st.ne.0) return
   flag(:) = 0

   ! Setup sv_map
   piv = 1
   do svc = 1, nsvar
      do piv = piv, piv + svar(svc) - 1
         sv_map(invp(piv)) = svc
      end do
   end do

   piv = 1
   idx = 1
   do svc = 1, nsvar
      col = invp(piv)
      ptr2(svc) = idx
      do j = ptr(col), nd_get_ptr(col+1, n, ne, ptr) - 1
         sv = sv_map(row(j))
         if(sv.eq.svc) cycle ! Skip diagonals
         if (flag(sv).eq.piv) cycle ! Already dealt with this supervariable
         ! Add row entry for this sv
         row2(idx) = sv
         flag(sv) = piv
         idx = idx + 1
      end do
      piv = piv + svar(svc)
   end do
   ptr2(svc) = idx
end subroutine nd_compress_by_svar

end module spral_nd_preprocess
