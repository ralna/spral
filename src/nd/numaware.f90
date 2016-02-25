module spral_nd_numaware
   use spral_nd_types
   use spral_nd_util
   use spral_scaling
   implicit none

   private
   public :: construct_aflags, order_matched_pair, nd_trim_keep_valid, &
      nd_expand_to_valid, nd_match_order_sep, make_a_match_symmetric, &
      extract_compressed_matrix
   public :: check_matrix_sym

   real(wp), parameter :: sing_tol = 100*epsilon(sing_tol) ! Singular tolerance

   integer, parameter :: nd_both_flag = 99

contains

subroutine construct_aflags(n, ptr, row, val, a_match, a_flags, options, inform)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   integer, dimension(n), intent(out) :: a_match
   integer, dimension(:), allocatable, intent(inout) :: a_flags
   type(nd_options), intent(in) :: options
   type(nd_inform), intent(inout) :: inform

   integer :: ne
   integer :: i, j, k
   real(wp), dimension(:), allocatable :: scaling, a_val

   type(hungarian_options) :: hoptions
   type(hungarian_inform) :: hinform
   type(auction_options) :: aoptions
   type(auction_inform) :: ainform

   ne = ptr(n+1)-1
   allocate(scaling(n), stat=inform%stat)
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


subroutine order_matched_pair(j1, j2, a_n, a_ne, a_ptr, a_flags_diag, &
      partition)
   integer, intent(inout) :: j1
   integer, intent(inout) :: j2
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(in) :: a_ptr(a_n)
   integer, intent(in) :: a_flags_diag(a_n)
   integer, optional, intent(in) :: partition(*)

   integer :: j, p1, p2
   integer :: i1, i2

   if(j2.eq.-1) return ! No ordering to do

   ! In some cases indexing into a_flags_diag needs translating...
   if(present(partition)) then
      i1 = partition(j1)
      i2 = partition(j2)
   else
      i1 = j1
      i2 = j2
   endif

   ! If we have exactly one "large" diagonal, order it first
   if(a_flags_diag(i1).eq.3 .and. a_flags_diag(i2).ne.3) then
      ! order j1 first: i.e. do nothing 
      return
   else if(a_flags_diag(i1).ne.3 .and. a_flags_diag(i2).eq.3) then
      ! order j2 first
      j = j1
      j1 = j2
      j2 = j
      return
   endif

   ! Otherwise: either both have large diagonals, or neither does. Choose
   ! densest column to go first in hopes it encourages other column
   ! into the same supernode
   p1 = nd_get_ptr(j1+1, a_n, a_ne, a_ptr) - a_ptr(j1)
   p2 = nd_get_ptr(j2+1, a_n, a_ne, a_ptr) - a_ptr(j2)
   if(p2 > p1) then
      j = j1
      j1 = j2
      j2 = j
   endif
end subroutine order_matched_pair


SUBROUTINE nd_expand_to_valid(a_n,a_ne,a_ptr,a_row,a_weight,a_flags,&
    a_flags_diag,a_match,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,&
    partition,work)
  ! Expand the separator to ensure that no incoming red edges are cut
  INTEGER, INTENT (IN) :: a_n ! order of matrix
  INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
  INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
  ! position in a_row that entries for column i start.
  INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
  ! indices of the non-zero rows. Diagonal entries have been removed
  ! and the matrix expanded.
  INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
  ! the weight of column i
  INTEGER, INTENT (IN) :: a_flags(a_ne) ! big/small flags in
  ! matrix
  INTEGER, INTENT (IN) :: a_flags_diag(a_n) ! big/small flags
  INTEGER, INTENT (INOUT) :: a_match(a_n) ! if a_match present, 
  ! i is matched to j => a_match(i) = j, if i is not matched, 
  ! a_match(i)=i. This matching may be updated.
  INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
  INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
  INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
  ! hted
  ! size of partitions and separator
  INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end. This is updated to the new
  ! partition
  INTEGER, TARGET, INTENT (OUT) :: work(10*a_n) ! Work array

  ! Local variables
  integer, dimension(:), pointer :: part, backtrace, progress, visited
  integer :: i, ii, j, k, kk, p, q
  integer :: depth, vcnt

  ! Stats
  integer :: npullin, naug1, naug2 ! aug1 = augment fixes 1, aug2=fixes 2

  ! Init stats
  npullin = 0
  naug1 = 0
  naug2 = 0

  ! Partition workspace
  part      => work(1:a_n)
  backtrace => work(a_n+1:2*a_n)
  progress  => work(2*a_n+1:3*a_n)
  visited   => work(3*a_n+1:4*a_n)

  ! Set part(:) to hold flags to indicate what
  ! part of the partition the nodes are in
  ! Note: partition(:) is the output array, part(:) is a pseduo-inverse
  !       work array!
  part(partition(1:a_n1)) = nd_part1_flag
  part(partition(a_n1+1:a_n1+a_n2)) = nd_part2_flag
  part(partition(a_n1+a_n2+1:a_n)) = nd_sep_flag

  ! Initialise workspaces
  visited(:) = 0; vcnt = 0

  ! Loop over seperator. For each matching that is cut, seek an
  ! augmenting path to modify matching. If none exists, move partner
  ! into seperator too.
  do ii = a_n1+a_n2+1, a_n
    i = partition(ii)
    j = a_match(i)
    if(j.eq.-1) cycle ! Not matched
    if(part(j).eq.nd_sep_flag) cycle ! j in sep, so edge not cut
    !print *, "Considering ", i, j
    ! Check if diagonal entries are sufficiently large anyway
    if(a_flags_diag(j).eq.3) then
      ! Either:
      ! (a) Both diagonals are large, match not actually required; or
      ! (b) Diagonal in seperator is small, one in partition is large.
      !     We can live with this as seperator diagonal will be modified
      !     by Schur complement of a large entry if our conditions are
      !     enforced.
      ! In both cases, break match and move to next column.
      !print *, "-->Part idx is sufficiently large, break match"
      a_match(i) = -1
      a_match(j) = -1
      cycle
    endif
    ! Otherwise, we have a matched pair cut with a small diagonal in a
    ! partition. Seek an augmenting path through flagged entries that
    ! results in an uncut matching.
    ! Depth first search
    depth = 1
    backtrace(depth) = j
    progress(depth) = a_ptr(j)
    vcnt = vcnt + 1
    depth_first_search: &
    do
      k = backtrace(depth)
      ! Loop through arcs of k
      do kk = progress(depth), nd_get_ptr(k+1, a_n, a_ne, a_ptr)-1
        if(.not.is_big_in_col(a_flags(kk))) cycle ! Not big
        p = a_row(kk)
        if(part(p).eq.nd_sep_flag) cycle ! Arc leads into sep
        if(any(backtrace(1:depth).eq.p)) cycle ! Avoid cycles
        if(visited(p).ge.vcnt) cycle ! Already visited index
        visited(p) = vcnt ! mark as visited
        q = a_match(p)
        if(q.lt.0) exit depth_first_search ! Found an unmatched edge
        if(a_flags_diag(q).eq.3) then
          ! p is matched to q, but q doesn't actually require a match,
          ! so we can augment along it
          exit depth_first_search
        endif
        if(part(q).eq.nd_sep_flag) then
          ! q is in seperator, so (p,q) is a cut edge. We can make q
          ! unmatched and still be happy
          naug2 = naug2 + 1
          exit depth_first_search
        endif
        ! Store current status for k and descend to q
        progress(depth) = kk
        depth = depth + 1
        backtrace(depth) = q
        progress(depth) = a_ptr(q)
        cycle depth_first_search
      end do
      ! If we reach this point, no valid augmenting path was found:
      ! backtrack up one level
      depth = depth - 1
      if(depth.eq.0) exit depth_first_search ! No augmenting paths
      progress(depth) = progress(depth) + 1 ! Finished with that edge
    end do depth_first_search
    if(depth.ne.0) progress(depth) = kk ! Store progress if we suceeded
    if(depth.eq.0) then
      ! No augmenting path was found: move j into seperator
      npullin = npullin + 1
      !print *, "-->No augmentation found, moving ", j, " to seperator"
      if(part(j).eq.nd_part1_flag) then
        a_n1 = a_n1 - 1
        a_weight_1 = a_weight_1 - a_weight(j)
        a_weight_sep = a_weight_sep + a_weight(j)
      else
        a_n2 = a_n2 - 1
        a_weight_2 = a_weight_2 - a_weight(j)
        a_weight_sep = a_weight_sep + a_weight(j)
      endif
      part(j) = nd_sep_flag
    else
      ! Augment along path
      naug1 = naug1 + 1
      !print *, "-->Found an augmenting path"
      a_match(i) = -1 ! Mark seperator entry as unmatched
      do kk = 1, depth
        k = backtrace(kk)
        p = a_row(progress(kk))
        q = a_match(p)
        !print *, "----> ", p, q, " (set ", k, p, " to match)"
        a_match(k) = p
        a_match(p) = k
        if(q.gt.0 .and. p.ne.q) a_match(q) = -1
      end do
    endif
  end do
  
  ! FIXME: Remove below assertions
  if(a_n1 .ne. count(part(:).eq.nd_part1_flag)) then
    print *, "a_n1 and count disagree: ", a_n1, count(part(:).eq.nd_part1_flag)
    stop
 endif
  if(a_n2 .ne. count(part(:).eq.nd_part2_flag)) then
    print *, "a_n2 and count disagree: ", a_n2, count(part(:).eq.nd_part2_flag)
    stop
 endif

  ! Update partition to match modified part(:)
  i = 1
  j = i + a_n1
  k = j + a_n2
  DO p = 1, a_n
    SELECT CASE (part(p))
    CASE (nd_part1_flag)
      partition(i) = p
      i = i + 1
    CASE (nd_part2_flag)
      partition(j) = p
      j = j + 1
    CASE DEFAULT
      partition(k) = p
      k = k + 1
    END SELECT
  END DO

  naug1 = naug1 - naug2 ! Double counted some
  print *, "Post expand npullin, naug1, naug2 = ", npullin, naug1, naug2
  call check_match_valid_partition(a_n, a_ne, a_ptr, a_row, &
    a_flags, a_flags_diag, a_match, a_n1, a_n2, partition, work)
END SUBROUTINE nd_expand_to_valid

SUBROUTINE nd_trim_keep_valid(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
    a_flags,a_flags_diag,a_match,a_n1,a_n2,a_weight_1,a_weight_2,      &
    a_weight_sep,partition,work,options)
  ! Trim the separator to make minimal whilst making sure that no cut red
  ! edges are introduced
  INTEGER, INTENT (IN) :: a_n ! order of matrix
  INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
  INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
  ! position in a_row that entries for column i start.
  INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
  ! indices of the non-zero rows. Diagonal entries have been removed
  ! and the matrix expanded.
  INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
  ! the weight of column i
  INTEGER, INTENT (IN) :: sumweight ! Sum of entries in a_weight
  INTEGER, INTENT (IN) :: a_flags(a_ne) ! big/small flags in
  ! matrix
  INTEGER, INTENT (IN) :: a_flags_diag(a_n) ! big/small flags
  INTEGER, INTENT (INOUT) :: a_match(a_n) ! if a_match present, 
  ! i is matched to j => a_match(i) = j, if i is not matched, 
  ! a_match(i)=i. This matching may be updated.
  INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
  INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
  INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
  ! hted
  ! size of partitions and separator
  INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end. This is updated to the new
  ! partition
  INTEGER, TARGET, INTENT (OUT) :: work(10*a_n) ! Work array
  TYPE(nd_options), INTENT (IN) :: options

  ! Local variables
  INTEGER :: i, j, k, kk, m, ncand, j1, j2, j12, k1, k2, p
  INTEGER :: mvchoice
  LOGICAL :: next1,next2,imbal
  REAL(wp) :: tau1,tau2,ratio
  integer :: depth1, depth2, partner1, partner2, partner, vcnt
  integer, dimension(:), pointer :: part, mask, bnd, visited
  integer, dimension(:), pointer :: backtrace1, backtrace2
  integer, dimension(:), pointer :: progress1, progress2

  ! Stats
  integer :: nreject
  integer :: npair
  integer :: naug

  nreject = 0
  npair = 0
  naug = 0

  ! Partition workspace
  part       => work(      1:  a_n)
  mask       => work(  a_n+1:2*a_n)
  bnd        => work(2*a_n+1:3*a_n)
  backtrace1 => work(3*a_n+1:4*a_n)
  progress1  => work(4*a_n+1:5*a_n)
  backtrace2 => work(5*a_n+1:6*a_n)
  progress2  => work(6*a_n+1:7*a_n)
  visited    => work(7*a_n+1:8*a_n)

  ! Set part(:) to hold flags to indicate what
  ! part of the partition the nodes are in
  part(partition(1:a_n1)) = nd_part1_flag
  part(partition(a_n1+1:a_n1+a_n2)) = nd_part2_flag
  part(partition(a_n1+a_n2+1:a_n)) = nd_sep_flag

  ! Initialize visited array
  visited(:) = 0; vcnt = 0

  !
  ! Set mask(j) to non-zero if in separator according to
  !   nd_part1_flag if only on boundary of partition 1
  !   nd_part2_flag if only on boundary of partition 2
  !   nd_sep_flag if not on boundary
  !   nd_both_flag if on boundary of both partitions
  !
  bnd(:) = -1
  mask(:) = 0
  ncand = 0
  do i = a_n1+a_n2+1,a_n
     j = partition(i)
     next1 = .FALSE.
     next2 = .FALSE.
     do kk = a_ptr(j), nd_get_ptr(j+1, a_n, a_ne, a_ptr)-1
       m = a_row(kk)
       if(part(m).eq.nd_part1_flag) next1 = .true.
       if(part(m).eq.nd_part2_flag) next2 = .true.
     end do
     mask(j) = 1
     if(next1.and.next2) then
       ! Adjacent to both
       bnd(j) = nd_both_flag
       mask(j) = 0
     else if(next1) then
       ! Adjacent to 1 only
       !print *, "Set ", j, " as cand1"
       ncand = ncand + 1
       bnd(j) = nd_part1_flag
     else if(next2) then
       ! Adjacent to 2 only
       !print *, "Set ", j, " as cand2"
       ncand = ncand + 1
       bnd(j) = nd_part2_flag
     else
       ! Adjacent to neither (internal to seperator)
       !print *, "Set ", j, " as cand12"
       ncand = ncand + 1
       bnd(j) = nd_sep_flag
     endif
  end do
        
  ratio = max(1.0_wp, options%balance)
  imbal = (ratio.le.sumweight-2.0)

  ! Keep removing pivots from seperator until no more candidates
  do while(ncand.gt.0)
    !print *, "=== ncand ", ncand, "==="

    ! Search for entry j1 that can only be moved into partition 1, 
    ! j2 that can only be moved into partition 2; keep track of an 
    ! j12 that can be moved to either partition as backup.
    j1 = 0
    j2 = 0
    j12 = 0
    do j = 1,a_n
      if(mask(j).ne.1) cycle ! Not a candidate
      ! Can be moved
      select case(bnd(j))
      case(nd_part1_flag)
        if(j1.eq.0) j1 = j
      case(nd_part2_flag)
        if(j2.eq.0) j2 = j
      case(nd_sep_flag)
        if(j12.eq.0) j12 = j
      end select
      if(j1.gt.0 .and. j2.gt.0) exit ! Found some candidates
    end do

    ! If we failed to find a j1 or j2, then use j12 instead
    if(j1.eq.0) j1 = j12
    if(j2.eq.0) j2 = j12

    ! Ensure candidates are matching compliant
    if(.not.check_match_compliance(j1, nd_part1_flag, a_n, a_ne, a_ptr, &
        a_row, a_flags, a_flags_diag, a_match, part, bnd, partner1,     &
        vcnt, visited, depth1, backtrace1, progress1)) then
      ! j1 can't be moved into partition 1: stop it from being a
      ! candidate to do so then try again
      if(bnd(j1).eq.nd_part1_flag) then
         ! Can't be a candidate any more
         !print *, "Marking ", j1, " as bad candidate"
         nreject = nreject + 1
         ncand = ncand - 1
         mask(j1) = 0
      else ! bnd(j1).eq.nd_sep_flag
         ! Mark as only possible to go to partition 2
         !print *, "Marking ", j1, " as only valid for part2"
         bnd(j1) = nd_part2_flag
      endif
      cycle
    endif
    if(.not.check_match_compliance(j2, nd_part2_flag, a_n, a_ne, a_ptr, &
        a_row, a_flags, a_flags_diag, a_match, part, bnd, partner2,     &
        vcnt, visited, depth2, backtrace2, progress2)) then
      ! j2 can't be moved into partition 2: stop it from being a
      ! candidate to do so then try again
      if(bnd(j2).eq.nd_part2_flag) then
         ! Can't be a candidate any more
         !print *, "Marking ", j2, " as bad candidate"
         nreject = nreject + 1
         ncand = ncand - 1
         mask(j2) = 0
      else ! bnd(j2).eq.nd_sep_flag
         ! Mark as only possible to go to partition 1
         bnd(j2) = nd_part1_flag
         !print *, "Marking ", j2, " as only valid for part1"
      endif
      cycle
    endif

    if(j1.eq.0 .and. j2.eq.0) then
       print *, "No candidate? But ncand = ", ncand
       stop
    endif

    ! Make choice (NB at least one of j1 or j2 must be non-zero otherwise
    !              loop would have already ended)
    tau1 = huge(tau1)
    if(j1.ne.0) then
       k1 = a_weight(j1)
       if(partner1.ne.0) then
         if(part(partner1).eq.nd_sep_flag) &
            k1 = k1 + a_weight(partner1)
       endif
       call cost_function(a_weight_1+k1, a_weight_2, &
         a_weight_sep-k1, sumweight, ratio, imbal, options%cost_function, tau1)
    endif
    tau2 = huge(tau2)
    if(j2.ne.0) then
       k2 = a_weight(j2)
       if(partner2.ne.0) then
         if(part(partner2).eq.nd_sep_flag) &
            k2 = k2 + a_weight(partner2)
       endif
       call cost_function(a_weight_1, a_weight_2+k2, &
         a_weight_sep-k2, sumweight, ratio, imbal, options%cost_function, tau2)
    endif
    if (tau1.lt.tau2) then
      ! Move j1 (and partner if it exists) into partition 1
      !print *, "Moving ", j1, "(", partner1, ") into part1"
      mvchoice = 1
      i = j1
      partner = partner1
      part(i) = nd_part1_flag
      a_n1 = a_n1 + 1
      if(partner.ne.0) then
        part(partner) = nd_part1_flag
        a_n1 = a_n1 + 1
      endif
      a_weight_1 = a_weight_1 + k1
      a_weight_sep = a_weight_sep - k1
      call trim_update_matching(j1, depth1, backtrace1, progress1, &
        a_row, a_match, naug)
    else
      ! Move j2 (and partner if it exists) into position 2
      !print *, "Moving ", j2, "(", partner2, ") into part2"
      mvchoice = 2
      i = j2
      partner = partner2
      part(i) = nd_part2_flag
      a_n2 = a_n2 + 1
      if(partner.ne.0) then
        part(partner) = nd_part2_flag
        a_n2 = a_n2 + 1
      endif
      a_weight_2 = a_weight_2 + k2
      a_weight_sep = a_weight_sep - k2
      call trim_update_matching(j2, depth2, backtrace2, progress2, &
        a_row, a_match, naug)
    endif
    ncand = ncand - 1
    mask(i) = 0
    if(partner.ne.0) then
      if(mask(partner).eq.1) then
        mask(partner) = 0
        ncand = ncand - 1 ! Partner was a candidate as well
      endif
    endif

    if(partner.ne.0) npair = npair + 1
   
    ! Assess neighbours of i (and partner if it exists) to see if their
    ! classification needs updating
    call trim_update_neighbours(i, a_n, a_ne, a_ptr, a_row, part, &
      ncand, mask, bnd)
    if(partner.ne.0) &
      call trim_update_neighbours(partner, a_n, a_ne, a_ptr, a_row, &
        part, ncand, mask, bnd)
  end do

  ! FIXME: Remove below assertions
  if(a_n1 .ne. count(part(:).eq.nd_part1_flag)) then
    print *, "a_n1 and count disagree: ", a_n1, count(part(:).eq.nd_part1_flag)
    stop
  endif
  if(a_n2 .ne. count(part(:).eq.nd_part2_flag)) then
    print *, "a_n2 and count disagree: ", a_n2, count(part(:).eq.nd_part2_flag)
    stop
 endif

  i = 1
  j = i + a_n1
  k = j + a_n2
  DO p = 1, a_n
    SELECT CASE (part(p))
    CASE (nd_part1_flag)
      partition(i) = p
      i = i + 1
    CASE (nd_part2_flag)
      partition(j) = p
      j = j + 1
    CASE DEFAULT
      partition(k) = p
      k = k + 1
    END SELECT
  END DO

  print *, "Post trim: naug, npair, nreject = ", naug, npair, nreject
  call check_match_valid_partition(a_n, a_ne, a_ptr, a_row, &
    a_flags, a_flags_diag, a_match, a_n1, a_n2, partition, work)
END SUBROUTINE nd_trim_keep_valid

SUBROUTINE check_match_valid_partition(a_n, a_ne, a_ptr, a_row, &
    a_flags, a_flags_diag, a_match, a_n1, a_n2, partition, work)
  integer, intent(in) :: a_n
  integer, intent(in) :: a_ne
  integer, intent(in) :: a_ptr(a_n)
  integer, intent(in) :: a_row(a_ne)
  integer, intent(in) :: a_flags(a_ne)
  integer, intent(in) :: a_flags_diag(a_n)
  integer, intent(in) :: a_match(a_n)
  integer, intent(in) :: a_n1
  integer, intent(in) :: a_n2
  integer, intent(in) :: partition(a_n)
  integer, target, intent(out) :: work(10*a_n)

  integer :: i, j, k
  integer :: ii
  integer, dimension(:), pointer :: part, seen
  logical :: seen_match

  part       => work(0*a_n+1:1*a_n)
  seen       => work(1*a_n+1:2*a_n)

  seen(:) = 0
  do i = 1, a_n
    j = partition(i)
    if(seen(j).ne.0) then
      print *, "Repeated partition component ", j, " at ", seen(j), " and ", i
      stop
    endif
    seen(j) = i
  end do

  ! Set part(:) to hold flags to indicate what
  ! part of the partition the nodes are in
  if(a_n1+a_n2.gt.a_n) then
     print *, "a_n1+a_n2 (", a_n1+a_n2, ") > a_n (", a_n, ")"
     stop
  endif
  part(partition(1:a_n1)) = nd_part1_flag
  part(partition(a_n1+1:a_n1+a_n2)) = nd_part2_flag
  part(partition(a_n1+a_n2+1:a_n)) = nd_sep_flag

  ! Check seperation and matching conditions
  do i = 1, a_n
    if(part(i).eq.nd_sep_flag) cycle ! Seperator can connect to anything
    j = a_match(i)
    seen_match = .false.
    do ii = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
      k = a_row(ii)
      if(part(i).ne.part(k) .and. part(k).ne.nd_sep_flag) then
         print *, "Seperator does not seperate: Entry (", k, ", ", i, &
           ") with part(k) = ", part(k), " part(i) = ", part(i)
         stop
      endif
      if(k.eq.j) seen_match = .true.
      if(k.eq.j .and. .not.is_big_in_col(a_flags(ii)) .and. a_flags_diag(i).ne.3) then
        print *, "Found matching edge (", k, ",", i, ") but a_flags=", &
          a_flags(ii)
        stop
      endif
    end do
    if(a_flags_diag(i).ne.3 .and. a_match(i).eq.0) then
      print *, "Non-seperator index ", i, " is not matched but should be"
      stop
    endif
    if(a_match(i).ne.-1) then
      j = a_match(i)
      if(part(j).ne.part(i)) then
        print *, "Matched indices ", i, "(", part(i), ") and ", &
           j, " (", part(j), ") not in same part"
        stop
      endif
      if(a_match(j).ne.i) then
        print *, "a_match is unsymmetrical: a_match(", i, ") = ", &
           a_match(i), " a_match(", j, ") = ", a_match(j)
        stop
      endif
    endif
  end do
END SUBROUTINE check_match_valid_partition

LOGICAL ELEMENTAL FUNCTION is_big_in_col(flag)
  integer, intent(in) :: flag

  integer :: r

  r = iand(flag, FLAG_BIG_COL)

  is_big_in_col = (r.ne.0)
END FUNCTION is_big_in_col

! Function checks whether idx can move out of the seperator without
! causing trouble. If necessary an augmented path is returned in the
! arrays backtrace(:) and progress(:). If a match is found within the
! seperator, the partner (that must leave with idx) is the return value,
! otherwise 0 is returned.
LOGICAL FUNCTION check_match_compliance(idx, dest, a_n, a_ne, a_ptr,  &
    a_row, a_flags, a_flags_diag, a_match, part, bnd, partner, vcnt, &
    visited, depth, backtrace, progress)
  integer, intent(in) :: idx ! Index we want to move out of seperator
  integer, intent(in) :: dest ! nd_part1_flag or nd_part2_flag
  integer, intent(in) :: a_n
  integer, intent(in) :: a_ne
  integer, intent(in) :: a_ptr(a_n)
  integer, intent(in) :: a_row(a_ne)
  integer, intent(in) :: a_flags(a_ne)
  integer, intent(in) :: a_flags_diag(a_n)
  integer, intent(in) :: a_match(a_n)
  integer, intent(in) :: part(a_n)
  integer, intent(in) :: bnd(a_n)
  integer, intent(out) :: partner
  integer, intent(inout) :: vcnt
  integer, intent(inout) :: visited(a_n)
  integer, intent(out) :: depth
  integer, intent(out) :: backtrace(a_n)
  integer, intent(out) :: progress(a_n)

  integer :: j, k, kk, p, q

  !print *, "CHECKING idx ", idx

  ! Initialization and immediate return conditions
  check_match_compliance = .true. ! default to match found
  partner = 0 ! No partner
  depth = 0 ! No augmentation
  if(idx.eq.0) return ! Immediate return if no idx specified
  !if(a_flags_diag(idx).eq.3) print *, "->OK large diag"
  if(a_flags_diag(idx).eq.3) return ! No match required: large diagonal
  ! If we reach here, then idx needs a partner

  ! Look for a partner within seperator to leave with:
  ! must be adjacent only to dest, or internal to seperator
  do kk = a_ptr(idx), nd_get_ptr(idx+1, a_n, a_ne, a_ptr)-1
    if(.not.is_big_in_col(a_flags(kk))) cycle ! Not large entry
    j = a_row(kk)
    if(part(j).ne.nd_sep_flag) cycle ! j not in seperator
    if(bnd(j).eq.dest .or. bnd(j).eq.nd_sep_flag) then
       ! Good partner found in seperator, go with it
       depth = 1
       backtrace(1) = idx
       progress(1) = kk
       partner = j
       !print *, "->OK good partner in sep ", j
       return
    endif
  end do
  ! If we reach here, no good partner to leave seperator with

  ! Look instead for a partner in the destination partition, may need
  ! an augmenting path
  depth = 1
  backtrace(depth) = idx
  progress(depth) = a_ptr(idx)
  vcnt = vcnt + 1
  depth_first_search: &
  do
    k = backtrace(depth)
    ! Loop through arcs of k
    do kk = progress(depth), nd_get_ptr(k+1, a_n, a_ne, a_ptr)-1
      if(.not.is_big_in_col(a_flags(kk))) cycle ! Not big enough
      p = a_row(kk)
      if(visited(p).ge.vcnt) cycle ! Already visited
      if(any(backtrace(1:depth).eq.p)) cycle ! Avoid cycles
      visited(p) = vcnt ! Mark as visited
      if(part(p).eq.nd_sep_flag) then
        ! Arc leads into sep - check if suitable partner
        if(bnd(p).eq.dest) then
           ! p only connects to destination partition, so valid partner
           partner = p
           progress(depth) = kk
           !print *, "->OK found augmenting path back into sep at ", p
           return
        else
           ! p connects to both partitions, not a valid partner
           cycle
        endif
      endif
      q = a_match(p)
      if(q.lt.0) then
        ! p is an unmatched edge
        progress(depth) = kk
        !print *, "->OK found an augmenting path to unmatched node ", p
        return
      endif
      if(a_flags_diag(q).eq.3) then
        ! p is matched to q, but q doesn't actually require a match,
        ! so we can augment along it
        progress(depth) = kk
        !print *, "->OK found an augmenting path to unneedy node ", p, q
        return
      endif
      ! Store current status for k and descend to q
      progress(depth) = kk ! Store current progress
      depth = depth + 1
      backtrace(depth) = q
      progress(depth) = a_ptr(q)
      cycle depth_first_search
    end do
    ! If we reach this point, no valid augmenting path was found:
    ! backtrack up one level
    depth = depth - 1
    if(depth.eq.0) exit depth_first_search ! No augmenting paths
    progress(depth) = progress(depth) + 1 ! finished with that edge
  end do depth_first_search

  ! If we reach this point, no good partner exists
  !print *, "->Failed to make good"
  check_match_compliance = .false.
END FUNCTION check_match_compliance

subroutine check_matrix_sym(n, ne, ptr, row)
   integer, intent(in) :: n, ne, ptr(n), row(ne)

   integer :: i, j, k
   integer, dimension(:), allocatable :: ptr_trans, row_trans, work

   !print *, "IN"
   !do i = 1, n
   !   print *, i, ":", row(ptr(i):nd_get_ptr(i+1,n,ne,ptr)-1)
   !end do

   !
   ! Find A^T
   !
   allocate(ptr_trans(n+3), row_trans(ne))
   ! Count # entries at +2 offset
   ptr_trans(:) = 0
   do i = 1, n
      do j = ptr(i), nd_get_ptr(i+1, n, ne, ptr)-1
         k = row(j)
         ptr_trans(k+2) = ptr_trans(k+2) + 1
      end do
   end do
   ! Calculate row starts at offset +1
   ptr_trans(1:2) = 1
   do i = 1, n
      ptr_trans(i+2) = ptr_trans(i+1) + ptr_trans(i+2)
   end do
   ! Drop entries into place, counting up in ptr at offset +1
   do i = 1, n
      do j = ptr(i), nd_get_ptr(i+1, n, ne, ptr)-1
         k = row(j)
         row_trans(ptr_trans(k+1)) = i
         ptr_trans(k+1) = ptr_trans(k+1) + 1
      end do
   end do

   !print *, "TRANS"
   !do i = 1, n
   !   print *, i, ":", row_trans(ptr_trans(i):ptr_trans(i+1)-1)
   !end do

   !
   ! Compare A and A^T
   !
   allocate(work(n))
   work(:) = 0
   do i = 1, n
      work(row_trans(ptr_trans(i):ptr_trans(i+1)-1)) = i
      do j = ptr(i), nd_get_ptr(i+1, n, ne, ptr)-1
         k = row(j)
         if(work(k).lt.i) then
            print *, "Assymetry! entry (", k, ",", i, ")"
            stop
         endif
      end do
   end do

   print *, "MATRIX IS SYMMETRIC"

end subroutine check_matrix_sym

! This subroutine orders sep_list to try and keep 2x2 pivots together
SUBROUTINE nd_match_order_sep(a_n,a_ne,a_ptr,a_row,a_flags,a_flags_diag,&
    a_match,iperm,nsep,sep_list,work)
  INTEGER, INTENT (IN) :: a_n ! order of matrix
  INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
  INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
  ! position in a_row that entries for column i start.
  INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
  ! indices of the non-zero rows. Diagonal entries have been removed
  ! and the matrix expanded.
  INTEGER, INTENT (IN) :: a_flags(a_ne) ! big/small flags in
  ! matrix
  INTEGER, INTENT (IN) :: a_flags_diag(a_n) ! big/small flags
  INTEGER, INTENT (INOUT) :: a_match(a_n) ! if a_match present, 
  ! i is matched to j => a_match(i) = j, if i is not matched, 
  ! a_match(i)=i. This matching may be updated.
  INTEGER, INTENT (INOUT) :: iperm(a_n)
  INTEGER, INTENT (IN) :: nsep ! Size of seperator
  INTEGER, INTENT (INOUT) :: sep_list(nsep) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end. This is updated to the new
  ! partition
  INTEGER, TARGET, INTENT (OUT) :: work(10*a_n) ! Work array

  integer, dimension(:), pointer :: old_list, seen
  integer :: i, ii, j, jj, j1, j2
  integer :: ins ! insert location

  old_list => work(0*a_n+1:1*a_n)
  seen     => work(1*a_n+1:2*a_n)

  ! First, kill any unnecessary matchings;
  ! also mark in seen if vertex is in seperator
  seen(:) = 0
  do ii = 1, nsep
    i = sep_list(ii)
    seen(i) = 1
    j = a_match(i)
    if(j.eq.-1) cycle ! no matched
    if(a_flags_diag(i).ne.FLAG_BIG_BOTH) cycle ! i needs matching
    if(a_flags_diag(j).ne.FLAG_BIG_BOTH) cycle ! j needs matching
    ! FIXME: kill any self matchings!
    ! If we reach here neither i nor j actually needs this matching
    a_match(i) = -1
    a_match(j) = -1
  end do

  ! Second, match anything that requires it if we can: simple greedy alg
  do ii = 1, nsep
    i = sep_list(ii)
    if(a_match(i).ne.-1) cycle ! Already matched
    if(a_flags_diag(i).eq.FLAG_BIG_BOTH) cycle ! Doesn't need a match
    do jj = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
      j = a_row(i)
      if(seen(j).eq.0) cycle ! j not in seperator
      if(a_match(j).ne.-1) cycle ! j already matched
      if(.not.is_big_in_col(a_flags(jj))) cycle ! not a big entry
      ! Otherwise we can match on this
      a_match(i) = j
      a_match(j) = i
    end do
  end do

  ! Record an ordering: use seen to record if something already ordered
  ins = 1
  old_list(1:nsep) = sep_list(:)
  seen(:) = 0
  do ii = 1, nsep
    j1 = old_list(ii)
    if(seen(j1).ne.0) cycle ! Already ordered
    j2 = a_match(j1)
    if(j1.eq.j2) j2 = -1 ! FIXME: handle above instead
    call order_matched_pair(j1, j2, a_n, a_ne, a_ptr, a_flags_diag)
    ! Order j1
    sep_list(ins) = j1
    seen(j1) = 1
    ins = ins + 1
    ! Order j2 if it exists
    if(j2.eq.-1) cycle ! No partner
    sep_list(ins) = j2
    seen(j2) = 1
    ins = ins + 1
    ! Flag as a 2x2
    iperm(j1) = -iperm(j1)
    iperm(j2) = -iperm(j2)
  end do
END SUBROUTINE nd_match_order_sep

SUBROUTINE trim_update_matching(idx, depth, backtrace, progress, a_row, &
    a_match, naug)
  integer, intent(in) :: idx
  integer, intent(in) :: depth
  integer, intent(in) :: backtrace(depth)
  integer, intent(in) :: progress(depth)
  integer, intent(in) :: a_row(*)
  integer, intent(inout) :: a_match(*)
  integer, intent(inout) :: naug

  integer :: k, p, q
  integer :: kk

  if(depth.eq.0) then
    ! No augmentation needed, however we may be un an uneeded match
    q = a_match(idx)
    if(q.gt.0) then
      !print *, "Break unneeded match ", idx, q
      a_match(idx) = -1
      a_match(q) = -1
    endif
    return
  endif

  ! Break any existing matching of first index (that is in seperator and
  ! whose match (if any) must also be in seperator)
  p = backtrace(1)
  q = a_match(p)
  !if(q.gt.0) print *, "Breaking match ", p, q
  if(q.gt.0) a_match(q) = -1
  ! NB: We will overwrite a_match(p) wehn we augment along path

  ! Augment along path
  !print *, "Augmenting along path:"
  naug = naug + 1
  do kk = 1, depth
    k = backtrace(kk)
    p = a_row(progress(kk))
    q = a_match(p)
    !print *, p, "->", q, "(", k, "and", p, " now matched)"
    a_match(k) = p
    a_match(p) = k
    if(q.gt.0 .and. p.ne.q) a_match(q) = -1
  end do
END SUBROUTINE trim_update_matching

SUBROUTINE trim_update_neighbours(idx, a_n, a_ne, a_ptr, a_row, part, &
    ncand, mask, bnd)
  integer, intent(in) :: idx
  integer, intent(in) :: a_n
  integer, intent(in) :: a_ne
  integer, intent(in) :: a_ptr(a_n)
  integer, intent(in) :: a_row(a_ne)
  integer, intent(in) :: part(a_n)
  integer, intent(inout) :: ncand
  integer, intent(inout) :: mask(a_n)
  integer, intent(inout) :: bnd(a_n)

  integer :: j
  integer :: ii

  do ii = a_ptr(idx), nd_get_ptr(idx+1, a_n, a_ne, a_ptr)-1
    j = a_row(ii)
    if(part(j).ne.nd_sep_flag) cycle ! j not in seperator
    if(mask(j).ne.1) cycle ! j not a candidate: we don't care
    select case(bnd(j))
    case(nd_part1_flag)
      if(part(idx).eq.nd_part2_flag) then
        ! j now adjacent to both, so not a candidate
        bnd(j) = nd_both_flag
        mask(j) = 0
        ncand = ncand - 1
      endif
    case(nd_part2_flag)
      if(part(idx).eq.nd_part1_flag) then
        ! j now adjacent to both, so not a candidate
        bnd(j) = nd_both_flag
        mask(j) = 0
        ncand = ncand - 1
      endif
    case(nd_sep_flag)
       ! j now adjacent to wherever i is, but still a candidate
       if(part(idx).eq.nd_part1_flag) bnd(j) = nd_part1_flag
       if(part(idx).eq.nd_part2_flag) bnd(j) = nd_part2_flag
    end select
  end do
END SUBROUTINE trim_update_neighbours


SUBROUTINE extract_compressed_matrix(a_n,a_ne,a_ptr,a_row,a_match, &
    a_n_part,a_n_sep,rows_sub,a_n_sub,a_ne_sub,a_ptr_sub, &
    len_a_row_sub,a_row_sub,compress_invp,work)
  INTEGER, INTENT (IN) :: a_n ! order of matrix being partitioned
  INTEGER, INTENT (IN) :: a_ne ! no. entries in matrix being partitioned
  INTEGER, INTENT (IN) :: a_ptr(a_n) ! col ptrs for matrix being
  ! partitioned
  INTEGER, INTENT (IN) :: a_row(a_ne) ! row indices for matrix
  ! being partitioned.
  INTEGER, INTENT (IN) :: a_match(a_n) ! matching to compress with
  INTEGER, INTENT (IN) :: a_n_part ! no. rows in partition
  INTEGER, INTENT (IN) :: a_n_sep ! no. rows in partition
  INTEGER, INTENT (IN) :: rows_sub(a_n_part+a_n_sep) ! rows/cols of
  ! matrix
  ! to be extracted. Intersecting rows/cols of separator will be
  ! replaced
  ! by matrix of all zeros
  INTEGER, INTENT (OUT) :: a_n_sub ! no. columns in extracted matrix
  INTEGER, INTENT (OUT) :: a_ne_sub ! no. entries stored in extracted
  ! matrix
  INTEGER, INTENT (OUT) :: a_ptr_sub(a_n_part+a_n_sep) ! col ptrs for
  ! extracted matrix
  INTEGER, INTENT (IN) :: len_a_row_sub ! length of a_row_sub
  INTEGER, INTENT (OUT) :: a_row_sub(len_a_row_sub) ! row indices for
  ! extracted matrix
  INTEGER, INTENT (OUT) :: compress_invp(2*a_n)
  INTEGER, TARGET, INTENT (OUT) :: work(2*a_n)

  ! Local variables
  INTEGER :: i, ii, j, k, kk, p, ins, part_last, cnt
  integer, pointer, dimension(:) :: lperm, mark

  lperm => work(0*a_n+1:1*a_n)
  mark  => work(1*a_n+1:2*a_n)

  ! Set lperm
  a_n_sub = a_n_part + a_n_sep
  lperm(:) = 0
  ins = 1
  do ii = 1, a_n_part
    i = rows_sub(ii)
    j = a_match(i)
    if(j.eq.i) j = -1 ! Matches with self disregarded
    if(j.ne.-1) then
      if(lperm(j).ne.0) then
        ! Already handled partner
        a_n_sub = a_n_sub - 1 ! Not appearing in this matrix
        compress_invp(2*(lperm(j)-1)+2) = ii
        cycle
      endif
    endif
    ! This is the representative column of the pair
    lperm(i) = ins  ! +ive means representative
    if(j.ne.-1) lperm(j) = -ins ! -ive means not representative
    compress_invp(2*(ins-1)+1) = ii
    if(j.eq.-1) compress_invp(2*(ins-1)+2) = -1
    ins = ins + 1
  end do
  part_last = ins-1 ! Last local index in part
  do ii = a_n_part+1, a_n_part+a_n_sep
    i = rows_sub(ii)
    lperm(i) = ins
    compress_invp(2*(ins-1)+1) = ii
    compress_invp(2*(ins-1)+2) = -1
    ins = ins + 1
  end do

  ! Count number of entries in submatrix
  mark(:) = 0
  a_ptr_sub(1:a_n_sub) = 0
  DO ii = 1, a_n_part
    i = rows_sub(ii)
    if(lperm(i).lt.0) cycle ! Already handled this column
    mark(abs(lperm(i))) = ii ! Mark diagonal so we don't include it!
    cnt = 0
    DO kk = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
      k = abs(lperm(a_row(kk)))
      if(k.eq.0) cycle ! not in part or sep
      if(mark(k).ge.ii) cycle ! Already seen this local index
      mark(k) = ii
      cnt = cnt + 1
    END DO
    j = a_match(i)
    if(j.eq.i) j = -1 ! disregard self matches
    if(j.ne.-1) then
      DO kk = a_ptr(j), nd_get_ptr(j+1, a_n, a_ne, a_ptr)-1
        k = abs(lperm(a_row(kk)))
        if(k.eq.0) cycle ! not in part or sep
        if(mark(k).ge.ii) cycle ! Already seen this local index
        mark(k) = ii
        cnt = cnt + 1
      END DO
    endif
    a_ptr_sub(abs(lperm(i))) = cnt
  END DO
  DO ii = a_n_part + 1, a_n_part+a_n_sep
    i = rows_sub(ii)
    cnt = 0
    DO kk = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
      k = abs(lperm(a_row(kk)))
      if(k.eq.0) cycle ! not in part or sep
      if(k.gt.part_last) cycle ! not in part
      if(mark(k).ge.ii) cycle ! Already seen this local index
      mark(k) = ii
      cnt = cnt + 1
    END DO
    a_ptr_sub(abs(lperm(i))) = cnt
  END DO
  ! Set up column pointers
  a_ptr_sub(1) = a_ptr_sub(1) + 1
  DO j = 2, a_n_sub
    a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
  END DO
  a_ne_sub = a_ptr_sub(a_n_sub) - 1

  ! Form a_row_sub
  mark(:) = 0
  DO ii = 1, a_n_part
    i = rows_sub(ii)
    if(lperm(i).lt.0) cycle ! Already handled this column
    mark(abs(lperm(i))) = ii ! Mark diagonal so we don't include it!
    DO kk = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
      k = abs(lperm(a_row(kk)))
      if(k.eq.0) cycle ! not in part or sep
      if(mark(k).ge.ii) cycle ! Already seen this local index
      mark(k) = ii
      a_ptr_sub(abs(lperm(i))) = a_ptr_sub(abs(lperm(i))) - 1
      p = a_ptr_sub(abs(lperm(i)))
      a_row_sub(p) = k
    END DO
    j = a_match(i)
    if(j.eq.i) j = -1 ! disregard self matches
    if(j.eq.-1) cycle ! No partner
    DO kk = a_ptr(j), nd_get_ptr(j+1, a_n, a_ne, a_ptr)-1
      k = abs(lperm(a_row(kk)))
      if(k.eq.0) cycle ! not in part or sep
      if(mark(k).ge.ii) cycle ! Already seen this local index
      mark(k) = ii
      a_ptr_sub(abs(lperm(i))) = a_ptr_sub(abs(lperm(i))) - 1
      p = a_ptr_sub(abs(lperm(i)))
      a_row_sub(p) = k
    END DO
  END DO
  DO ii = a_n_part + 1, a_n_part+a_n_sep
    i = rows_sub(ii)
    DO kk = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
      k = abs(lperm(a_row(kk)))
      if(k.eq.0) cycle ! not in part or sep
      if(k.gt.part_last) cycle ! not in part
      if(mark(k).ge.ii) cycle ! Already seen this local index
      mark(k) = ii
      a_ptr_sub(abs(lperm(i))) = a_ptr_sub(abs(lperm(i))) - 1
      p = a_ptr_sub(abs(lperm(i)))
      a_row_sub(p) = k
    END DO
  END DO

  !print *, "Output:"
  !print *, "ptr = ", a_ptr_sub(1:a_n_sub)
  !print *, "ne = ", a_ne_sub
  !do i = 1, a_n_sub
  !  print *, i, ":", a_row_sub(a_ptr_sub(i):nd_get_ptr(i+1,a_n_sub,a_ne_sub,a_ptr_sub)-1)
  !end do
  !print *, "compress_invp = ", compress_invp(1:2*part_last)
  !stop
END SUBROUTINE extract_compressed_matrix


end module spral_nd_numaware
