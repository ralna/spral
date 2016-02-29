! COPYRIGHT (c) 2012-3 Science and Technology Facilities Council
! Authors: Jonathan Hogg and Jennifer Scott
! Note: This code is a heavily modified version of HSL_MC80

! Given a sparse symmetric  matrix A, this module provides routines to
! use a matching algorithm to compute an elimination
! order that is suitable for use with a sparse direct solver. 
! It optionally computes scaling factors.

module spral_match_order
   use spral_metis_wrapper, only : metis_order
   use spral_nd, only : nd_order, &
                        nd_options_type => nd_options, &
                        nd_inform_type => nd_inform
   use spral_scaling, only : hungarian_scale_sym, &
                             hungarian_options_type => hungarian_options, &
                             hungarian_inform_type => hungarian_inform
   implicit none

   private
   public :: match_order ! Find a matching-based ordering
   public :: mo_options, mo_inform

   integer, parameter :: wp = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)

   ! Error flags
   integer, parameter :: SUCCESS               = 0
   integer, parameter :: ERROR_ALLOCATION      = -1
   integer, parameter :: ERROR_A_N_OOR         = -2
   integer, parameter :: ERROR_SINGULAR        = -3
   integer, parameter :: ERROR_UNKNOWN         = -99

   ! warning flags
   integer, parameter :: WARNING_SINGULAR      = 1

   type mo_options
      integer :: match_method = 1 ! 1=mc64, 2=mc64-alt, 3=gupta, 4=gupta-zero
      integer :: order_method = 1 ! ordering method to use 1=metis, 2=nd
      type(nd_options_type) :: nd_options ! options for call to nd
      real(wp) :: u = 0.01 ! Threshold for match_method!=1
      logical :: abort_if_singular = .true.
      real(wp) :: small = 1e-18 ! Threshold for zero values match_method=4
   end type mo_options

   type mo_inform
      integer :: flag = 0
      integer :: st = 0
      integer :: nmatch = 0 ! Number of matched vars
      type(nd_inform_type) :: nd_inform
      type(hungarian_inform_type) :: hungarian_inform
   end type mo_inform

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! On input ptr, row , val hold the ** lower AND upper ** triangular 
! parts of the matrix.
! this reduces amount of copies of matrix required (so slightly
! more efficient on memory and does not need to expand supplied matrix)
!  
subroutine match_order(n, ptr, row, val, order, scale, options, inform, seps)
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   integer, dimension(:), intent(out) :: order ! order(i)  holds the position
      ! of variable i in the elimination order (pivot sequence). 
   real(wp), dimension(n), intent(out) :: scale ! returns the mc64 symmetric
      ! scaling
   type(mo_options), intent(in) :: options
   type(mo_inform), intent(out) :: inform
   integer, dimension(n), optional, intent(out) :: seps

   integer, dimension(:), allocatable :: match ! used to hold matching
   integer(long), dimension(:), allocatable :: ptr2 ! column pointers for
      ! expanded matrix. 
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix 
   real(wp), dimension(:), allocatable :: val2 ! entries of expanded matrix.

   integer :: ne
   type(hungarian_options_type) :: hoptions

   inform%flag = 0
   inform%st = 0

   ! check n has valid value
   if (n < 0) then
     inform%flag = ERROR_A_N_OOR
     return
   end if

   ! just return with no action if n = 0
   if (n.eq.0) return

   !
   ! Find Hungarian scaling. As log() is expensive, only use half matrix
   ! (even though it is expanded again internally to routine).
   !
   ne = ptr(n+1)-1
   allocate(ptr2(n+1), row2(ne), val2(ne), match(n), stat=inform%st)
   if (inform%st.ne.0) then
      inform%flag = ERROR_ALLOCATION
      return
   end if
   call full_to_half(n, ptr, row, val, ptr2, row2, val2)
   hoptions%scale_if_singular = (.not.options%abort_if_singular)
   call hungarian_scale_sym(n, ptr, row, val, scale, hoptions, &
      inform%hungarian_inform, match=match)
   select case(inform%hungarian_inform%flag)
   case(0)
      ! Success, do nothing
   case(1)
      ! Singular warning
      inform%flag = WARNING_SINGULAR
   case(-1)
      ! Allocation error
      inform%flag = ERROR_ALLOCATION
      inform%st = inform%hungarian_inform%stat
      return
   case(-2)
      ! Singular and options%abort_if_singular=.true.
      inform%flag = ERROR_SINGULAR
      return
   end select

   !
   ! Find a matching
   !
   select case(options%match_method)
   case(:1)
      ! Just symmetrise MC64 matching
      call symmetrise_matching(n, match, inform%st)
      if(inform%st.ne.0) then
         inform%flag = ERROR_SINGULAR
         return
      endif
   case(2)
      ! Symmetrise MC64 matching, scrapping unneeded matches
      call symmetrise_matching_relaxed(n, ptr, row, val, scale, match, &
         options%u, inform%st)
      if(inform%st.ne.0) then
         inform%flag = ERROR_SINGULAR
         return
      endif
   case(3:)
      ! Gupta-method, partner = arg max |a_ij| / nz_extra(i,j)
      ! If match_method=3, use 1x1 test with u. O'wise use zero test vs small
      call gupta_matching(n, ptr, row, val, scale, match, options, inform%st)
      if(inform%st.ne.0) then
         inform%flag = ERROR_SINGULAR
         return
      endif
   end select

   !
   ! Form compressed matrix, order it, and return uncompressed ordering
   !
   call compress_order(n,ptr,row,order,match,options,inform,seps=seps)
end subroutine match_order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Drops entries outwith lower part of matrix
subroutine full_to_half(n, ptr_in, row_in, val_in, ptr_out, row_out, val_out)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr_in
   integer, dimension(ptr_in(n+1)-1), intent(in) :: row_in
   real(wp), dimension(ptr_in(n+1)-1), intent(in) :: val_in
   integer(long), dimension(*), intent(out) :: ptr_out
   integer, dimension(*), intent(out) :: row_out
   real(wp), dimension(*), intent(out) :: val_out

   integer :: i, k
   integer(long) :: jj

   k = 1
   do i = 1, n
      ptr_out(i) = k
      do jj = ptr_in(i), ptr_in(i+1)-1
         if(row_in(jj).lt.i .or. row_in(jj).gt.n) cycle
         row_out(k) = row_in(jj)
         val_out(k) = val_in(jj)
         k = k + 1
      end do
   end do
end subroutine full_to_half

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Given an unsymmetric matching, constructs a symmetric one.
!
subroutine symmetrise_matching(n, match, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: match
   integer, intent(out) :: st

   integer :: i, j
   logical, dimension(:), allocatable :: work

   allocate(work(n), stat=st)
   if(st.ne.0) return

   work(:) = .false.
   do i = 1, n
      j = match(i)
      if(j.eq.-1 .or. j.eq.i) then
         ! Unmatched, or matched on diagonal
         match(i) = -1
         cycle
      endif
      if(work(i)) cycle ! i already matched to something else
      if(work(j)) then
         ! partner already matched; if to something else, mark i as unmatched
         if(match(j).ne.i) match(i) = -1
         cycle
      endif
      ! Otherwise, we are free to match i and j
      match(j) = i
      work(i) = .true.
      work(j) = .true.
   end do
end subroutine symmetrise_matching

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Given an unsymmetric matching, constructs a symmetric one.
!
subroutine symmetrise_matching_relaxed(n, ptr, row, val, scale, match, u, st)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   real(wp), dimension(n), intent(inout) :: scale
   integer, dimension(n), intent(inout) :: match
   real(wp), intent(in) :: u
   integer, intent(out) :: st

   integer :: i, j
   logical, dimension(:), allocatable :: work

   allocate(work(n), stat=st)
   if(st.ne.0) return

   work(:) = .false.
   do i = 1, n
      j = match(i)
      if(j.eq.-1 .or. j.eq.i) then
         ! Unmatched, or matched on diagonal
         match(i) = -1
         cycle
      endif
      if(work(i)) cycle ! i already matched to something else
      if(.not.needs_partner(i,ptr,row,val,scale,u)) then
         ! Doesn't actually need to be matched
         match(i) = -1
         cycle
      endif
      if(work(j)) then
         ! partner already matched; if to something else, mark i as unmatched
         if(match(j).ne.i) match(i) = -1
         cycle
      endif
      ! Otherwise, we are free to match i and j
      match(j) = i
      work(i) = .true.
      work(j) = .true.
   end do
end subroutine symmetrise_matching_relaxed

subroutine gupta_matching(n, ptr, row, val, scale, match, options, st)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   real(wp), dimension(n), intent(in) :: scale
   integer, dimension(n), intent(out) :: match
   type(mo_options), intent(in) :: options
   integer, intent(out) :: st

   integer :: i, k
   integer(long) :: jj, pp
   integer :: nz_extra, best_idx
   integer, dimension(:), allocatable :: pattern
   real(wp) :: score, best_score
   integer :: failed_to_pair

   allocate(pattern(n), stat=st)
   if(st.ne.0) return
   pattern(:) = 0

   failed_to_pair = 0
   match(:) = -1 ! Initially all unmatched
   do i = 1, n
      if(match(i).ne.-1) cycle ! Already matched
      if(options%match_method.eq.3) then
         if(.not.needs_partner(i,ptr,row,val,scale,options%u)) &
            cycle ! No match required
      else
         if(.not.is_diag_zero(i,ptr,row,val,scale,options%small)) &
            cycle ! No match required
      endif

      ! Create pattern(:) of column
      do jj = ptr(i), ptr(i+1)-1
         pattern(row(jj)) = i
      end do

      ! Build a list of candidate merge columns
      ! Score them based on number of extra entries merging will create
      best_idx = -1
      best_score = -huge(best_score)
      do jj = ptr(i), ptr(i+1)-1
         k = row(jj)
         if(match(k).gt.0) cycle ! k already matched
         if(k.eq.i) cycle ! no self matching
         ! Count number of extra entries from merging columns i and k
         nz_extra = int(ptr(i+1)-ptr(i)+1) + int(ptr(k+1)-ptr(k)+1)
         do pp = ptr(k), ptr(k+1)-1
            if(pattern(row(pp)) .ge. i) & ! Overlapping entry
               nz_extra = nz_extra - 2
         end do
         ! Calculate score and check if it is best yet
         score = abs(val(jj)*scale(k)) / nz_extra
         if(score.gt.best_score) then
            best_score = score
            best_idx = k
         endif
      end do
      if(best_idx.ne.-1) then
         match(i) = best_idx
         match(best_idx) = i
      else
         failed_to_pair = failed_to_pair + 1
      endif
   end do
   print *, "Failed to pair ", failed_to_pair
end subroutine gupta_matching

logical pure function is_diag_zero(idx, ptr, row, val, scale, small)
   integer, intent(in) :: idx
   integer, dimension(*), intent(in) :: ptr
   integer, dimension(*), intent(in) :: row
   real(wp), dimension(*), intent(in) :: val
   real(wp), dimension(*), intent(in) :: scale
   real(wp), intent(in) :: small

   integer :: j
   real(wp) :: v

   is_diag_zero = .true.
   do j = ptr(idx), ptr(idx+1)-1
      if(row(j).eq.idx) then
         v = scale(idx) * abs(val(j)) * scale(idx)
         is_diag_zero = (v.le.small)
         return
      endif
   end do
end function is_diag_zero

logical pure function needs_partner(idx, ptr, row, val, scale, u)
   integer, intent(in) :: idx
   integer, dimension(*), intent(in) :: ptr
   integer, dimension(*), intent(in) :: row
   real(wp), dimension(*), intent(in) :: val
   real(wp), dimension(*), intent(in) :: scale
   real(wp), intent(in) :: u

   integer :: j
   real(wp) :: dv, ov ! diagonal and off-diagonal values
   real(wp) :: v

   dv = 0.0
   ov = 0.0
   do j = ptr(idx), ptr(idx+1)-1
      v = abs(val(j)) * scale(row(j)) ! Technically *scale(i), but it cancels
      if(row(j).eq.idx) then
         dv = v
      else
         ov = max(ov, v)
      endif
   end do

   needs_partner = (dv .le. u*ov)
end function needs_partner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Split matching into 1- and 2-cycles only and then
! compress matrix and order using mc68.
!
! Input (ptr2, row2 , val2) holds the ** lower and upper triangles ** 
! of the matrix (with explicit zeros removed).
! Overwritten in the singular case
!
subroutine compress_order(n,ptr2,row2,order,match,options,inform,seps)
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr2 
   integer, dimension(:), intent(in) :: row2 
   integer, dimension(n), intent(out) :: order ! used to hold ordering
   integer, dimension(n), intent(in) :: match ! used to hold matching 
   type(mo_options), intent(in) :: options
   type(mo_inform), intent(inout) :: inform
   integer, dimension(n), optional, intent(out) :: seps

   integer, dimension(:), allocatable :: iwork ! work array
   integer, dimension(:), allocatable :: old_to_new, new_to_old
      ! holds mapping between original matrix indices and those in condensed
      ! matrix.
   integer, dimension(:), allocatable :: ptr3 ! column pointers for condensed 
      ! matrix.
   integer, dimension(:), allocatable :: row3 ! row indices for condensed 
      ! matrix.
   ! FIXME: Removes seps array as its just here for debug
   integer, dimension(:), allocatable :: seps2 ! Get list of seperators

   integer :: i, j, j1, j2, jj, k, krow, metis_flag
   integer(long) :: klong
   integer :: ncomp ! order of compressed matrix
   integer :: ncomp_matched ! order of compressed matrix (matched entries only)
   integer(long) :: ne ! number of non zeros
   integer, dimension(:), allocatable :: invp

   ne = ptr2(n+1) - 1
   allocate(ptr3(n+1), row3(ne), old_to_new(n), new_to_old(n), iwork(n), &
      stat=inform%st)
   if (inform%st.ne.0) return

   !
   ! Build maps for new numbering schemes
   !
   inform%nmatch = 0
   k = 1
   do i = 1,n
      j = match(i)
      if (j<i .and. j.gt.0) cycle
      old_to_new(i) = k
      new_to_old(k) = i ! note: new_to_old only maps to first of a pair
      if (j.gt.0) then
         old_to_new(j) = k   
         inform%nmatch = inform%nmatch + 1
      endif
      k = k + 1
   end do
   ncomp_matched = k-1

   !
   ! Produce a condensed version of the matrix for ordering.
   ! Hold pattern using ptr3 and row3.
   !
   ptr3(1) = 1
   iwork(:) = 0 ! Use to indicate if entry is in a paired column
   ncomp = 1
   jj = 1
   do i = 1, n
      j = match(i)
      if (j<i .and. j.gt.0) cycle ! already seen
      do klong = ptr2(i), ptr2(i+1)-1
         krow = old_to_new(row2(klong))
         if (iwork(krow).eq.i) cycle ! already added to column
         if (krow>ncomp_matched) cycle ! unmatched row not participating
         row3(jj) = krow
         jj = jj + 1
         iwork(krow) = i
      end do
      if (j.gt.0) then
         ! Also check column match(i)
         do klong = ptr2(j), ptr2(j+1)-1
            krow = old_to_new(row2(klong))
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

   ! Reorder compressed matrix
   select case(options%order_method)
   case(1) ! METIS
      allocate(invp(ncomp), stat=inform%st)
      if (inform%st.ne.0) return
      call metis_order(ncomp,ptr3,row3,order,invp,metis_flag,inform%st)
      select case(metis_flag)
      case(0)
         ! OK, do nothing
      case(-1)
         ! Allocation error
         inform%flag = ERROR_ALLOCATION
         return
      case default
         ! Unknown error, should never happen
         print *, "metis_order() returned unknown error ", metis_flag
         inform%flag = ERROR_UNKNOWN
         return
      end select
   case(2)
      ! Use spral_nd to order
      ! Note: No need for invp array
      allocate(seps2(n))
      call nd_order(0,ncomp,ptr3,row3,order,options%nd_options, &
         inform%nd_inform,seps=seps2)
      select case(inform%nd_inform%flag)
      case(0)
        ! OK, do nothing
      case default
         ! Unknown error, should never happen
         print *, "nd_order() returned unknown error ", inform%nd_inform%flag
         inform%flag = ERROR_UNKNOWN
         return
      end select
      print *, "ncomp = ", ncomp
      i = maxval(seps2(1:ncomp))
      print *, "max depth = ", i
      do i = -1, i
         print *, "sep lvl", i, " sz = ", count(seps2(1:ncomp).eq.i)
      end do
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
      if(present(seps)) seps(j) = seps2(i)
      k = k + 1
      if (match(j).gt.0) then
         j = match(j)
         order(j) = k
         if(present(seps)) seps(j) = seps2(i)
         k = k + 1
      end if
   end do

end subroutine compress_order

!**********************************************************************
end module spral_match_order
