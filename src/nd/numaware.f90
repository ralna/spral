module spral_nd_numaware
   use spral_nd_types
   use spral_scaling
   implicit none

   private
   public :: construct_aflags

   ! NB: We exploit the additive nature of the below when setting.
   integer, parameter :: FLAG_SMALL    = 0, &
                         FLAG_BIG_COL  = 1, &
                         FLAG_BIG_ROW  = 2, &
                         FLAG_BIG_BOTH = 3

   real(wp), parameter :: sing_tol = 100*epsilon(sing_tol) ! Singular tolerance

contains

subroutine construct_aflags(n, ptr, row, val, a_flags, options, inform)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   integer, dimension(:), allocatable, intent(inout) :: a_flags
   type(nd_options), intent(in) :: options
   type(nd_inform), intent(inout) :: inform

   integer :: ne
   integer :: i, j, k
   real(wp), dimension(:), allocatable :: scaling, a_val
   integer, dimension(:), allocatable :: a_match

   type(hungarian_options) :: hoptions
   type(hungarian_inform) :: hinform
   type(auction_options) :: aoptions
   type(auction_inform) :: ainform

   ne = ptr(n+1)-1
   allocate(scaling(n), a_match(n), stat=inform%stat)
   if (inform%stat.ne.0) return
   select case (options%scaling)
   case (0)
      scaling(:) = 1.0_wp
   case (1)
      call auction_scale_sym(n,ptr,row,val,scaling,aoptions,ainform,&
         match = a_match)
      if (ainform%flag.eq.-1) inform%stat = ainform%stat
   case (2)
      call hungarian_scale_sym(n,ptr,row,val,scaling,hoptions,hinform,&
          match = a_match)
      if (hinform%flag.eq.-1) inform%stat = hinform%stat
   end select
   if (inform%stat.ne.0) return
   call make_a_match_symmetric(n, a_match)

   ! Scale values
   allocate (a_val(ne), stat=inform%stat)
   if (inform%stat.ne.0) return
   ! Fill a_row and a_ptr (and a_val if val present)
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         a_val(k) = val(k)*scaling(j)*scaling(i)
      end do
   end do

   ! Find flags
   allocate (a_flags(ne), stat=inform%stat)
   if (inform%stat.ne.0) return
   select case(options%large_method)
   case(1)
      ! Relative method of 1x1 and 2x2 pivot tests
      call nd_set_a_flags(options%u, n, ptr, row, a_val, a_match, a_flags)
   case(2)
      ! Duff-Pralet absolute method
      call nd_set_a_flags_dp(options%dp_theta, options%dp_drop, n, &
         ptr, row, a_val, a_match, a_flags)
   end select
end subroutine construct_aflags

subroutine make_a_match_symmetric(n, match)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: match

   logical, dimension(:), allocatable :: seen

   integer :: i, j

   ! Mark all vertices as unseen
   allocate(seen(n))
   seen(:) = .false.

   ! Iterate over vertices and ensure we only have symmetric pairings
   do i = 1, n
      if(seen(i)) cycle ! Already part of a pair
      j = match(i)
      if(j.le.0) then
        ! Vertex is unmatched
        match(i) = -1
        cycle
      endif
      if(seen(j)) then
         match(i) = -1 ! Partner unavailable, set as unmatched
         cycle
      endif
      ! Break self-matches
      if(i.eq.j) then
         match(i) = -1
         cycle
      endif
      ! Mark i and j as seen and symmetrise match (NB often i = j)
      seen(i) = .true.
      seen(j) = .true.
      match(i) = j
      match(j) = i
   end do
end subroutine make_a_match_symmetric

! Should work with either half or full matrix
! We mark diagonal entries as big if they make a suitable pivot, or not if they
! don't.
subroutine nd_set_a_flags(u, n, ptr, row, val, match, a_flags)
   real(wp), intent(in) :: u ! Pivot tolerance, must be <= 0.5
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   integer, dimension(n), intent(in) :: match
   integer, dimension(ptr(n+1)-1), intent(out) :: a_flags

   real(wp), dimension(:), allocatable :: dval, cmax, cmax2

   integer :: i, k
   integer(long) :: jj
   logical :: bigrow, bigcol
   real(wp) :: v

   ! Find column maxima and diagonal entries
   allocate(dval(n), cmax(n), cmax2(n))
   dval(:) = 0.0 ! Diagonal values
   cmax(:) = 0.0 ! Column maxima (of off diagonal entries)
   cmax2(:) = 0.0 ! Runner up column maxima (of off diagonal entries)
   do i = 1, n
      do jj = ptr(i), ptr(i+1)-1
         k = row(jj)
         v = abs(val(jj))
         if(k.eq.i) then
            dval(i) = abs(val(jj))
         else
            if(v .ge. cmax(i)) then
               cmax2(i) = cmax(i)
               cmax(i) = v
            else if(v.gt.cmax2(i)) then
               cmax2(i) = v
            endif
            if(v .ge. cmax(k)) then
               cmax2(k) = cmax(k)
               cmax(k) = v
            else if(v.gt.cmax2(k)) then
               cmax2(k) = v
            endif
         endif
      end do
   end do
   !print *, "dval = ", dval(1:5)
   !print *, "cmax = ", cmax(1:5)

   ! Now fill out flag array
   do i = 1, n
      !print *, "===== COLUMN ", i, "======"
      do jj = ptr(i), ptr(i+1)-1
         k = row(jj)
         !print *, "ENTRY ", k, i, val(jj)
         if(k.eq.i) then
            ! Diagonal entry. Considered big if a good pivot, small otherwise
            bigcol = (dval(i).ge.u*cmax(i))
            !print *, "Diagonal test ", dval(i), " >= ", u, "*", cmax(i), "=", bigcol
            bigrow = bigcol
         else
            ! Off diagonal entry. Considered small if corresponding row/col is
            ! considered a good pivot, otherwise check how it works as a 2x2
            ! Check col
            if(dval(i).ge.u*cmax(i)) then
               bigcol = .false. ! Diagonal is sufficient
               !print *, "col diag sufficient"
            else
               bigcol = is_good_2x2(k, i, val(jj), dval, cmax, cmax2, u)
            endif
            ! Check row
            if(dval(k).ge.u*cmax(k)) then
               bigrow = .false. ! Diagonal is sufficient
               !print *, "row diag sufficient"
            else
               bigrow = is_good_2x2(i, k, val(jj), dval, cmax, cmax2, u)
            endif
         endif
         !print *, "RESULT ", bigrow, bigcol
         a_flags(jj) = FLAG_SMALL
         if(bigrow) a_flags(jj) = a_flags(jj) + FLAG_BIG_ROW
         if(bigcol) a_flags(jj) = a_flags(jj) + FLAG_BIG_COL
      end do
   end do

   ! Ensure any match recommended by MC64 is flagged as big
   do i = 1, n
      do jj = ptr(i), ptr(i+1)-1
         k = row(jj)
         if(match(i).eq.k .or. match(k).eq.i) then
            a_flags(jj) = FLAG_BIG_BOTH
         endif
      end do
   end do

   do i = 1, n
      do jj = ptr(i), ptr(i+1)-1
         k = row(jj)
         if(i.ne.301232 .and. k.ne.301232) cycle ! Don't care
         if(i.eq.k) print *, "Entry (", k, ",", i, ") is ", val(jj)
         if(a_flags(jj).eq.0) cycle ! Not big
         print *, "Entry (", k, ",", i, ") is big ", a_flags(jj), val(jj)
      end do
   end do

end subroutine nd_set_a_flags

! Should work with either half or full matrix
! We mark diagonal entries as big if they make a suitable pivot, or not if they
! don't.
subroutine nd_set_a_flags_dp(theta, drop, n, ptr, row, val, match, a_flags)
   real(wp), intent(in) :: theta ! Diagonal small tolerance
   real(wp), intent(in) :: drop ! Off-diagonal drop tolerance
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   integer, dimension(n), intent(in) :: match
   integer, dimension(ptr(n+1)-1), intent(out) :: a_flags

   integer :: i, k
   integer(long) :: jj

   ! Now fill out flag array
   do i = 1, n
      do jj = ptr(i), ptr(i+1)-1
         k = row(jj)
         a_flags(jj) = FLAG_SMALL
         if(k.eq.i) then
            ! Diagonal entry
            if(abs(val(jj)).ge.theta) &
               a_flags(jj) = FLAG_BIG_BOTH
         else
            ! Off diagonal entry
            if(abs(val(jj)).ge.drop) &
               a_flags(jj) = FLAG_BIG_BOTH
         endif
      end do
   end do

   ! Ensure any match recommended by MC64 is flagged as big
   do i = 1, n
      do jj = ptr(i), ptr(i+1)-1
         k = row(jj)
         if(match(i).eq.k .or. match(k).eq.i) then
            a_flags(jj) = FLAG_BIG_BOTH
         endif
      end do
   end do

end subroutine nd_set_a_flags_dp

logical function is_good_2x2(r, c, val, dval, cmax, cmax2, u)
   integer, intent(in) :: r ! row
   integer, intent(in) :: c ! column
   real(wp), intent(in) :: val ! value
   real(wp), dimension(*), intent(in) :: dval ! diagonal values
   real(wp), dimension(*), intent(in) :: cmax ! column maxima (off diagonal)
   real(wp), dimension(*), intent(in) :: cmax2 ! runner up cmax
   real(wp), intent(in) :: u ! pivot tolerance

   real(wp) :: detval, utest, v, cmr, cmc

   detval = abs(dval(c)*dval(r) - val**2)
   if(detval .le. sing_tol) then
      ! Pivot is singular, give up on it as pivot
      !print *, "2x2 test", r, c, "singular ", detval
      is_good_2x2 = .false.
      return
   endif
   ! Pivot is good as long as
   ! 1/det * | (  a_rr -a_cr ) | ( cmax(c) ) < u^-1 * ( 1 )
   !         | ( -a_cr  a_cc ) | ( cmax(r) )          ( 1 )
   v = abs(val)
   cmc = cmax(c)
   if(v.eq.cmc) cmc = cmax2(c)
   cmr = cmax(r)
   if(v.eq.cmr) cmr = cmax2(r)
   utest = detval / max( &
      dval(r)*cmc +  v     *cmr, &
       v     *cmc + dval(c)*cmr )
   is_good_2x2 = (utest.ge.u)
end function is_good_2x2



end module spral_nd_numaware
