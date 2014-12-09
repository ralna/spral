module spral_scaling
   use spral_matrix_util, only : half_to_full
   implicit none

   private
   ! Top level routines
   public :: auction_scale_sym,   & ! Symmetric scaling by Auction algorithm
             auction_scale_unsym, & ! Unsymmetric scaling by Auction algorithm
             equilib_scale_sym,   & ! Scaling by Equilibriation (MC77-like)
             hungarian_scale_sym    ! Scaling by Hungarian algorithm (MC64-like)
   ! Inner routines that allow calling internals
   public :: hungarian_match      ! Find a matching (no pre/post-processing)
   ! Data types
   public :: auction_options, auction_inform, &
             equilib_options, equilib_inform, &
             hungarian_options, hungarian_inform

   integer, parameter :: wp = kind(0d0)
   real(wp), parameter :: rinf = huge(rinf)

   type auction_options
      integer :: max_iterations = 30000
      integer :: max_unchanged(3) = (/ 10,   100, 100 /)
      real :: min_proportion(3) = (/ 0.90, 0.0, 0.0 /)
      real :: eps = 0.01
   end type auction_options

   type auction_inform
      integer :: flag = 0 ! success or failure
      integer :: stat = 0 ! Fortran stat value on memory allocation failure
      integer :: matched = 0 ! #matched rows/cols
      integer :: iterations = 0 ! #iterations
   end type auction_inform

   type equilib_options
      integer :: max_iterations = 10
      real :: tol = 1e-8
   end type equilib_options

   type equilib_inform
      integer :: flag
      integer :: stat
   end type equilib_inform

   type hungarian_options
      logical :: scale_if_singular = .false.
   end type hungarian_options

   type hungarian_inform
      integer :: flag
      integer :: stat
   end type hungarian_inform

   integer, parameter :: ERROR_ALLOCATION = -1
   integer, parameter :: ERROR_SINGULAR = -2

   integer, parameter :: WARNING_SINGULAR = 1

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Wrappers around scaling algorithms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**************************************************************
!
! Use matching-based scaling obtained using Hungarian algorithm
!
subroutine hungarian_scale_sym(n, ptr, row, val, scaling, options, inform, &
      match)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   type(hungarian_options), intent(in) :: options
   type(hungarian_inform), intent(out) :: inform
   integer, dimension(n), optional, intent(out) :: match

   integer :: flag
   integer, dimension(:), allocatable :: perm

   inform%flag = 0 ! Initialize to success
   
   if(present(match)) then
      call hungarian_wrapper(n,ptr,row,val,match,scaling,flag,inform%stat)
   else
      allocate(perm(n), stat=inform%stat)
      if (inform%stat .ne. 0) then
         inform%flag = ERROR_ALLOCATION
         return
      endif
      call hungarian_wrapper(n,ptr,row,val,perm,scaling,flag,inform%stat)
   endif
   select case(flag)
   case(0)
      ! success; do nothing
   case(1)
      ! structually singular matrix
      if(.not.options%scale_if_singular) then
         ! abort, set scaling to identity
         inform%flag = ERROR_SINGULAR
         scaling(:) = 1.0
         return
      else
         inform%flag = WARNING_SINGULAR
      end if
   case(ERROR_ALLOCATION)
      ! allocation error
      inform%flag = ERROR_ALLOCATION
      return
   end select

   scaling(1:n) = exp(scaling(1:n))

end subroutine hungarian_scale_sym

!**************************************************************
!
! Call auction algorithm to get a scaling, then symmetrize it
!
subroutine auction_scale_sym(n, ptr, row, val, scaling, options, inform, match)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   type(auction_options), intent(in) :: options
   type(auction_inform), intent(out) :: inform
   integer, dimension(n), optional, intent(out) :: match

   real(wp), dimension(:), allocatable :: rscaling, cscaling

   inform%flag = 0 ! Initialize to sucess

   ! Allocate memory
   allocate(rscaling(n), cscaling(n), stat=inform%stat)
   if(inform%stat.ne.0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif

   ! Call unsymmetric implementation with flag to expand half matrix
   call auction_match(.true., n, n, ptr, row, val, rscaling, cscaling, &
      options, inform, match)

   ! Average rscaling and cscaling to get symmetric scaling
   scaling(1:n) = exp( (rscaling(1:n)+cscaling(1:n))/2 )

end subroutine auction_scale_sym

!**************************************************************
!
! Call auction algorithm to get a scaling, then symmetrize it
!
subroutine auction_scale_unsym(m, n, ptr, row, val, rscaling, cscaling, &
      options, inform, match)
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(m), intent(out) :: rscaling
   real(wp), dimension(n), intent(out) :: cscaling
   type(auction_options), intent(in) :: options
   type(auction_inform), intent(out) :: inform
   integer, dimension(n), optional, intent(out) :: match

   inform%flag = 0 ! Initialize to sucess

   call auction_match(.false., m, n, ptr, row, val, rscaling, cscaling, &
      options, inform, match)

   rscaling(1:m) = exp(rscaling(1:m))
   cscaling(1:n) = exp(cscaling(1:n))

end subroutine auction_scale_unsym

!*******************************
!
! Call the infinity-norm equilibriation algorithm
!
subroutine equilib_scale_sym(n, ptr, row, val, scaling, options, inform)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   type(equilib_options), intent(in) :: options
   type(equilib_inform), intent(out) :: inform

   integer :: i

   real(wp), allocatable, dimension(:) :: scale_orig

   inform%flag = 0 ! Initialize to sucess

   allocate(scale_orig(n), stat=inform%stat)
   if(inform%stat.ne.0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif

   call inf_norm_equilib(n, ptr, row, val, scale_orig, options, inform)

   do i = 1, n
      scaling(i) = scale_orig(i)
   end do
   
end subroutine equilib_scale_sym

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inf-norm Equilibriation Algorithm Implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! We implement Algorithm 1 of:
! A Symmetry Preserving Algorithm for Matrix Scaling
! Philip Knight, Daniel Ruiz and Bora Ucar
! INRIA Research Report 7552 (Novemeber 2012)
! (This is similar to the algorithm used in MC77, but is a complete
!  reimplementation from the above paper to ensure it is 100% STFC
!  copyright and can be released as open source)
!
subroutine inf_norm_equilib(n, ptr, row, val, scaling, options, inform)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   real(wp), dimension(n), intent(out) :: scaling
   type(equilib_options), intent(in) :: options
   type(equilib_inform), intent(inout) :: inform

   integer :: itr, r, c, j
   real(wp) :: v
   real(wp), dimension(:), allocatable :: maxentry

   allocate(maxentry(n), stat=inform%stat)
   if(inform%stat.ne.0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif

   scaling(1:n) = 1.0
   do itr = 1, options%max_iterations
      ! Find maximum entry in each row and col
      ! Recall: matrix is symmetric, but we only have half
      maxentry(1:n) = 0.0
      do c = 1, n
         do j = ptr(c), ptr(c+1)-1
            r = row(j)
            v = abs(scaling(r) * val(j) * scaling(c))
            maxentry(r) = max(maxentry(r), v)
            maxentry(c) = max(maxentry(c), v)
         end do
      end do
      ! Update scaling (but beware empty cols)
      where(maxentry(1:n).gt.0) &
         scaling(1:n) = scaling(1:n) / sqrt(maxentry(1:n))
      ! Test convergence
      if(maxval(abs(1-maxentry(1:n))) .lt. options%tol) exit
   end do

end subroutine inf_norm_equilib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Hungarian Algorithm implementation (MC64)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**********************************************************************
!
! This subroutine wraps the core algorithm of hungarian_match(). It provides
! pre- and post-processing to transform a maximum product assignment to a
! minimum sum assignment problem (and back again). It also has post-processing
! to handle the case of a structurally singular matrix as per Duff and Pralet
! (though the efficacy of such an approach is disputed!)
!
! This code is adapted from HSL_MC64 v2.3.1
!
subroutine hungarian_wrapper(n,ptr,row,val,perm,scale,flag,stat)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(*), intent(in) :: row
   real(wp), dimension(*), intent(in) :: val
   integer, dimension(n), intent(out) :: perm
   real(wp), dimension(n), intent(out) :: scale
   integer, intent(out) :: flag
   integer, intent(out) :: stat

   integer, allocatable :: ptr2(:), row2(:), iw(:), new_to_old(:), &
      old_to_new(:), cperm(:)
   real(wp), allocatable :: val2(:), dw(:), cmax(:), cscale(:)
   real(wp) :: colmax
   integer :: i,j,k,ne,num,nn,j1,j2,jj
   real(wp), parameter :: zero = 0.0

   stat = 0
   ne = ptr(n+1)-1

   ! Reset ne for the expanded symmetric matrix
   ne = 2*ne

   ! Expand matrix, drop explicit zeroes and take log absolute values
   allocate (ptr2(n+1), row2(ne), val2(ne), &
             iw(5*n), dw(2*n), cmax(n), stat=stat)
   if (stat/=0) then
      flag = ERROR_ALLOCATION
      return
   endif

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if(val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
      ! Following log is seperated from above loop to expose expensive
      ! log operation to vectorization.
      val2(ptr2(i):k-1) = log(val2(ptr2(i):k-1))
   end do
   ptr2(n+1) = k
   call half_to_full(n, row2, ptr2, iw, a=val2)

   ! Compute column maximums
   do i = 1, n
      colmax = maxval(val2(ptr2(i):ptr2(i+1)-1))
      cmax(i) = colmax
      val2(ptr2(i):ptr2(i+1)-1) = colmax - val2(ptr2(i):ptr2(i+1)-1)
   end do

   call hungarian_match(n,ne,ptr2,row2,val2,perm,num,dw(1),dw(n+1),stat)
   if (stat.ne.0) then
      flag = ERROR_ALLOCATION
      return
   endif

   flag = 0

   if(num.eq.n) then ! Full rank
      ! Note that in this case m=n
      scale(1:n) = (dw(1:n)+dw(n+1:2*n)-cmax(1:n))/2
      return
   endif

   ! If we reach this point then structually rank deficient:
   ! Build a full rank submatrix and call matching on it.
   flag = 1 ! structually singular warning

   allocate(old_to_new(n),new_to_old(n),cperm(n),stat=stat)
   if (stat.ne.0) then
      flag = ERROR_ALLOCATION
      stat = stat
   end if

   j = num + 1
   k = 0
   do i = 1,n
      if (perm(i) < 0) then
         ! row i is not part of the matching
         old_to_new(i) = -j
         j = j + 1
      else
         k = k + 1
         ! old_to_new(i) holds the new index for variable i after
         ! removal of singular part and new_to_old(k) is the
         ! original index for k
         old_to_new(i) = k
         new_to_old(k) = i
      end if
   end do

   ! Overwrite ptr2, row2 and val2
   ne = 0
   k = 0
   ptr2(1) = 1
   j2 = 1
   do i = 1,n
      j1 = j2
      j2 = ptr2(i+1)
      ! skip over unmatched entries
      if (perm(i) < 0) cycle
      k = k + 1
      do j = j1,j2-1
         jj = row2(j)
         if (perm(jj) < 0) cycle
         ne = ne + 1
         row2(ne) = old_to_new(jj)
         val2(ne) = val2(j)
      end do
      ptr2(k+1) = ne + 1
   end do
   ! nn is order of non-singular part.
   nn = k
   call hungarian_match(nn,ne,ptr2,row2,val2,cperm,num,dw(1),dw(nn+1),stat)
   if (stat.ne.0) then
      flag = ERROR_ALLOCATION
      return
   endif

   do i = 1,n
      j = old_to_new(i)
      if (j < 0) then
         scale(i) = -huge(scale)
      else
         ! Note: we need to subtract col max using old matrix numbering
         scale(i) = (dw(j)+dw(nn+j)-cmax(i))/2
      end if
   end do

   perm(1:n) = -1
   do i = 1,nn
      j = cperm(i)
      perm(new_to_old(i)) = j
   end do

   do i = 1, n
      if(perm(i).eq.-1) then
         perm(i) = old_to_new(i)
      endif
   end do

   ! Apply Duff and Pralet correction to unmatched row scalings
   allocate(cscale(n), stat=stat)
   if(stat/=0) then
      flag = ERROR_ALLOCATION
      stat = stat
      return
   endif
   ! For columns i not in the matched set I, set
   !     s_i = 1 / (max_{k in I} | a_ik s_k |)
   ! with convention that 1/0 = 1
   cscale(1:n) = scale(1:n)
   do i = 1,n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         if(cscale(i).eq.-huge(scale).and.cscale(k).ne.-huge(scale)) then
            ! i not in I, k in I
            scale(i) = max(scale(i), log(abs(val(j)))+scale(k))
         endif
         if(cscale(k).eq.-huge(scale).and.cscale(i).ne.-huge(scale)) then
            ! k not in I, i in I
            scale(k) = max(scale(k), log(abs(val(j)))+scale(i))
         endif
      end do
   end do
   do i = 1,n
      if(cscale(i).ne.-huge(scale)) cycle ! matched part
      if(scale(i).eq.-huge(scale)) then
         scale(i) = zero
      else
         scale(i) = -scale(i)
      endif
   end do

end subroutine hungarian_wrapper

!**********************************************************************
!
! Provides the core Hungarian Algorithm implementation for solving the
! minimum sum assignment problem as per Duff and Koster.
!
! This code is adapted from MC64 v 1.6.0
!
subroutine hungarian_match(n,ne,ptr,row,val,iperm,num,u,d,st)
   integer, intent(in) :: n ! matrix dimension
   integer, intent(in) :: ne ! number of entries in matrix
   integer, intent(out) :: num ! cardinality of the matching
   integer, intent(in) :: ptr(n+1) ! column pointers
   integer, intent(in) :: row(ne) ! row pointers
   integer, intent(out) :: iperm(n) ! matching itself: row i is matched to
      ! column iperm(i).
   real(wp), intent(in) :: val(ne) ! value of the entry that corresponds to
      ! row(k). All values val(k) must be non-negative.
   real(wp), intent(out) :: u(n) ! u(i) is the reduced weight for row(i)
   real(wp), intent(out) :: d(n) ! d(i) is current shortest distance to row i
      ! from current column (d_i from Fig 4.1 of Duff and Koster paper)
   integer, intent(out) :: st

   integer, allocatable, dimension(:) :: jperm ! a(jperm(j)) is entry of A for
      ! matching in column j.
   integer, allocatable, dimension(:) :: out ! a(out(i)) is the new entry in a
      ! on which we match going along the short path back to original col.
   integer, allocatable, dimension(:) :: pr ! pr(i) is a pointer to the next
      ! column along the shortest path back to the original column
   integer, allocatable, dimension(:) :: Q ! q(1:qlen) forms a binary heap
      ! data structure sorted by d(q(i)) value. q(low:up) is a list of rows
      ! with equal d(i) which is lower or equal to smallest in the heap.
      ! q(up:n) is a list of already visited rows.
   integer, allocatable, dimension(:) :: L ! l(:) is an inverse of q(:)

   integer :: i,i0,ii,j,jj,jord,q0,qlen,jdum,isp,jsp
   integer :: k,k0,k1,k2,kk,kk1,kk2,up,low,lpos
   real(wp) :: csp,di,dmin,dnew,dq0,vj

   !
   ! Initialization
   !
   allocate(jperm(n), out(n), pr(n), q(n), l(n), stat=st)
   if(st.ne.0) return
   num = 0
   d(1:n) = 0.0
   iperm(1:n) = 0
   jperm(1:n) = 0
   pr(1:n) = ptr(1:n)
   l(1:n) = 0

   !
   ! Initialize u(i) using heurisitic to get an intial guess.
   ! Heuristic guaruntees that the generated partial matching is optimal
   ! on the restriction of the graph to the matched rows and columns.
   !
   u(1:n) = RINF
   do j = 1,n
      do k = ptr(j),ptr(j+1)-1
         i = row(k)
         if(val(k).gt.u(i)) cycle
         u(i) = val(k)
         iperm(i) = j
         l(i) = k
      end do
   end do
   do i = 1,n
      j = iperm(i)
      if(j.eq.0) cycle
      ! Row i is not empty
      iperm(i) = 0
      if(jperm(j).ne.0) cycle
      ! Don't choose cheap assignment from dense columns
      if(ptr(j+1)-ptr(j) .gt. n/10 .and. n.gt.50) cycle
      ! Assignment of column j to row i
      num = num + 1
      iperm(i) = j
      jperm(j) = l(i)
   end do
   if(num.eq.n) GO TO 1000
   ! Scan unassigned columns; improve assignment
   improve_assign: do j = 1,n
      ! jperm(j) ne 0 iff column j is already assigned
      if(jperm(j).ne.0) cycle
      k1 = ptr(j)
      k2 = ptr(j+1) - 1
      ! Continue only if column j is not empty
      if(k1.gt.k2) cycle
      i0 = row(k1)
      vj = val(k1) - u(i0)
      k0 = k1
      do k = k1+1,k2
         i = row(k)
         di = val(k) - u(i)
         if(di.gt.vj) cycle
         if(di.ge.vj .and. di.ne.RINF) then
            if(iperm(i).ne.0 .or. iperm(i0).eq.0) cycle
         endif
         vj = di
         i0 = i
         k0 = k
      end do
      d(j) = vj
      k = k0
      i = i0
      if(iperm(i).eq.0) then
         num = num + 1
         jperm(j) = k
         iperm(i) = j
         pr(j) = k + 1
         cycle
      endif
      do k = k0,k2
         i = row(k)
         if(val(k)-u(i).gt.vj) cycle
         jj = iperm(i)
         ! Scan remaining part of assigned column jj
         kk1 = pr(jj)
         kk2 = ptr(jj+1) - 1
         if(kk1.gt.kk2) cycle
         do kk = kk1,kk2
            ii = row(kk)
            if(iperm(ii).gt.0) cycle
            if(val(kk)-u(ii).le.d(jj)) then
               jperm(jj) = kk
               iperm(ii) = jj
               pr(jj) = kk + 1
               num = num + 1
               jperm(j) = k
               iperm(i) = j
               pr(j) = k + 1
               cycle improve_assign
            endif
         end do
         pr(jj) = kk2 + 1
      end do
      cycle
   end do improve_assign
   if(num.eq.n) go to 1000 ! If heurisitic got a complete matching, we're done

   !
   ! Repeatedly find augmenting paths until all columns are included in the
   ! matching. At every step the current matching is optimal on the restriction
   ! of the graph to currently matched rows and columns.
   !

   ! Main loop ... each pass round this loop is similar to Dijkstra's
   ! algorithm for solving the single source shortest path problem
   d(1:n) = RINF
   l(1:n) = 0
   isp=-1; jsp=-1 ! initalize to avoid "may be used unitialized" warning
   do jord = 1,n

      if(jperm(jord).ne.0) cycle
      ! jord is next unmatched column
      ! dmin is the length of shortest path in the tree
      dmin = RINF
      qlen = 0
      low = n + 1
      up = n + 1
      ! csp is the cost of the shortest augmenting path to unassigned row
      ! row(isp). The corresponding column index is jsp.
      csp = RINF
      ! Build shortest path tree starting from unassigned column (root) jord
      j = jord
      pr(j) = -1

      ! Scan column j
      do k = ptr(j),ptr(j+1)-1
         i = row(k)
         dnew = val(k) - u(i)
         if(dnew.ge.csp) cycle
         if(iperm(i).eq.0) then
            csp = dnew
            isp = k
            jsp = j
         else
            if(dnew.lt.dmin) dmin = dnew
            d(i) = dnew
            qlen = qlen + 1
            q(qlen) = k
         endif
      end do
      ! Initialize heap Q and Q2 with rows held in q(1:qlen)
      q0 = qlen
      qlen = 0
      do kk = 1,q0
         k = q(kk)
         i = row(k)
         if(csp.le.d(i)) then
            d(i) = RINF
            cycle
         endif
         if(d(i).le.dmin) then
            low = low - 1
            q(low) = i
            l(i) = low
         else
            qlen = qlen + 1
            l(i) = qlen
            call heap_update(i,n,Q,D,L)
         endif
         ! Update tree
         jj = iperm(i)
         out(jj) = k
         pr(jj) = j
      end do

      do jdum = 1,num

         ! If Q2 is empty, extract rows from Q
         if(low.eq.up) then
            if(qlen.eq.0) exit
            i = q(1) ! Peek at top of heap
            if(d(i).ge.csp) exit
            dmin = d(i)
            ! Extract all paths that have length dmin and store in q(low:up-1)
            do while(qlen.gt.0)
               i = q(1) ! Peek at top of heap
               if(d(i).gt.dmin) exit
               i = heap_pop(qlen,n,Q,D,L)
               low = low - 1
               q(low) = i
               l(i) = low
            end do
         endif
         ! q0 is row whose distance d(q0) to the root is smallest
         q0 = q(up-1)
         dq0 = d(q0)
         ! Exit loop if path to q0 is longer than the shortest augmenting path
         if(dq0.ge.csp) exit
         up = up - 1

         ! Scan column that matches with row q0
         j = iperm(q0)
         vj = dq0 - val(jperm(j)) + u(q0)
         do k = ptr(j),ptr(j+1)-1
            i = row(k)
            if(l(i).ge.up) cycle
            ! dnew is new cost
            dnew = vj + val(k)-u(i)
            ! Do not update d(i) if dnew ge cost of shortest path
            if(dnew.ge.csp) cycle
            if(iperm(i).eq.0) then
               ! Row i is unmatched; update shortest path info
               csp = dnew
               isp = k
               jsp = j
            else
               ! Row i is matched; do not update d(i) if dnew is larger
               di = d(i)
               if(di.le.dnew) cycle
               if(l(i).ge.low) cycle
               d(i) = dnew
               if(dnew.le.dmin) then
                  lpos = l(i)
                  if(lpos.ne.0) call heap_delete(lpos,qlen,n,Q,D,L)
                  low = low - 1
                  q(low) = i
                  l(i) = low
               else
                  if(l(i).eq.0) then
                     qlen = qlen + 1
                     l(i) = qlen
                  endif
                  call heap_update(i,n,Q,D,L) ! d(i) has changed
               endif
               ! Update tree
               jj = iperm(i)
               out(jj) = k
               pr(jj) = j
            endif
         end do
      end do

      ! If csp = RINF, no augmenting path is found
      if(csp.eq.RINF) GO TO 190
      ! Find augmenting path by tracing backward in pr; update iperm,jperm
      num = num + 1
      i = row(isp)
      iperm(i) = jsp
      jperm(jsp) = isp
      j = jsp
      do jdum = 1,num
         jj = pr(j)
         if(jj.eq.-1) exit
         k = out(j)
         i = row(k)
         iperm(i) = jj
         jperm(jj) = k
         j = jj
      end do

      ! Update U for rows in q(up:n)
      do kk = up,n
         i = q(kk)
         u(i) = u(i) + d(i) - csp
      end do
190   do kk = low,n
         i = q(kk)
         d(i) = RINF
         l(i) = 0
      end do
      do kk = 1,qlen
         i = q(kk)
         d(i) = RINF
         l(i) = 0
      end do

   end do ! End of main loop


   ! Set dual column variable in d(1:n)
1000 do j = 1,n
      k = jperm(j)
      if(k.ne.0) then
         d(j) = val(k) - u(row(k))
      else
         d(j) = 0.0
      endif
      if(iperm(j).eq.0) u(j) = 0.0
   end do

   ! Return if matrix has full structural rank
   if(num.eq.n) return

   ! Otherwise, matrix is structurally singular, complete iperm.
   ! jperm, out are work arrays
   jperm(1:n) = 0
   k = 0
   do i = 1,n
      if(iperm(i).eq.0) then
         k = k + 1
         out(k) = i
      else
         j = iperm(i)
         jperm(j) = i
      endif
   end do
   k = 0
   do j = 1,n
      if(jperm(j).ne.0) cycle
      k = k + 1
      jdum = out(k)
      iperm(jdum) = -j
   end do
end subroutine hungarian_match

!**********************************************************************
!
! Value associated with index i has decreased, update position in heap
! as approriate.
!
! This code is adapted from MC64 v 1.6.0
!
subroutine heap_update(idx,N,Q,val,L)
   integer, intent(in) :: idx
   integer, intent(in) :: N
   integer, intent(inout) :: Q(N)
   integer, intent(inout) :: L(N)
   real(wp), intent(in) :: val(N)

   integer :: pos,parent_pos,parent_idx
   real(wp) :: v

   ! Get current position of i in heap
   pos = L(idx)
   if(pos.le.1) then
      ! idx is already at root of heap, but set q as it may have only just
      ! been inserted.
      q(pos) = idx
      return
   endif

   ! Keep trying to move i towards root of heap until it can't go any further
   v = val(idx)
   do while(pos.gt.1) ! while not at root of heap
      parent_pos = pos / 2
      parent_idx = Q(parent_pos)
      ! If parent is better than idx, stop moving
      if(v.ge.val(parent_idx)) exit
      ! Otherwise, swap idx and parent
      Q(pos) = parent_idx
      L(parent_idx) = pos
      pos = parent_pos
   end do
   ! Finally set idx in the place it reached.
   Q(pos) = idx
   L(idx) = pos
end subroutine heap_update

!**********************************************************************
!
! The root node is deleted from the binary heap.
!
! This code is adapted from MC64 v 1.6.0
!
integer function heap_pop(QLEN,N,Q,val,L)
   integer, intent(inout) :: QLEN
   integer, intent(in) :: N
   integer, intent(inout) :: Q(N)
   integer, intent(inout) :: L(N)
   real(wp), intent(in) :: val(N)

   ! Return value is the old root of the heap
   heap_pop = q(1)

   ! Delete the root
   call heap_delete(1,QLEN,N,Q,val,L)

end function heap_pop

!**********************************************************************
!
! Delete element in poisition pos0 from the heap
!
! This code is adapted from MC64 v 1.6.0
!
subroutine heap_delete(pos0,QLEN,N,Q,D,L)
   integer :: pos0,QLEN,N
   integer :: Q(N),L(N)
   real(wp) :: D(N)

   integer :: idx,pos,parent,child,QK
   real(wp) :: DK,DR,v

   ! If we're trying to remove the last item, just delete it.
   if(QLEN.eq.pos0) then
      QLEN = QLEN - 1
      return
   endif

   ! Replace index in position pos0 with last item and fix heap property
   idx = Q(QLEN)
   v = D(idx)
   QLEN = QLEN - 1 ! shrink heap
   pos = pos0 ! pos is current position of node I in the tree

   ! Move up if appropriate
   if(pos.gt.1) then
      do
         parent = pos / 2
         QK = Q(parent)
         if(v.ge.D(QK)) exit
         Q(pos) = QK
         L(QK) = pos
         pos = parent
         if(pos.le.1) exit
      end do
   endif
   Q(pos) = idx
   L(idx) = pos
   if(pos.ne.pos0) return ! Item moved up, hence doesn't need to move down

   ! Otherwise, move item down
   do
      child = 2 * pos
      if(child.gt.QLEN) exit
      DK = D(Q(child))
      if(child.lt.QLEN) then
         DR = D(Q(child+1))
         if(DK.gt.DR) then
            child = child + 1
            DK = DR
         endif
      endif
      if(v.le.DK) exit
      QK = Q(child)
      Q(pos) = QK
      L(QK) = pos
      pos = child
   end do
   Q(pos) = idx
   L(idx) = pos

end subroutine heap_delete

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Auction Algorithm implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! An implementation of the auction algorithm to solve the assignment problem
! i.e. max_M sum_{(i,j)\in M} a_{ij}    where M is a matching.
! The dual variables u_i for row i and v_j for col j can be used to find
! a good scaling after postprocessing.
! We're aiming for:
! a_ij - u_i - v_j == 0    if (i,j) in M
! a_ij - u_i - v_j <= 0    otherwise
!
! Motivation of algorithm:
! Each unmatched column bids for its preferred row. Best bid wins.
! Prices (dual variables) are updated to reflect cost of using 2nd best instead
! for future auction rounds.
! i.e. Price of using entry (i,j) is a_ij - u_i
!
! In this implementation, only one column is bidding in each phase. This is
! largely equivalent to the motivation above but allows for faster progression
! as the same row can move through multiple partners during a single pass
! through the matrix.
!
subroutine auction_match_core(m, n, ptr, row, val, match, dualu, dualv, &
      options, inform)
   integer, intent(in) :: m
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   integer, dimension(n), intent(out) :: match
      ! match(j) = i => column j matched to row i
   real(wp), dimension(m), intent(out) :: dualu ! row dual variables
   real(wp), dimension(n), intent(inout) :: dualv ! col dual variables
   type(auction_options), intent(in) :: options
   type(auction_inform), intent(inout) :: inform

   integer, dimension(:), allocatable :: owner ! Inverse of match
   ! The list next(1:tail) is the search space of unmatched columns
   ! this is overwritten as we proceed such that next(1:insert) is the
   ! search space for the subsequent iteration.
   integer :: tail, insert
   integer, dimension(:), allocatable :: next
   integer :: unmatched ! Current number of unmatched cols

   integer :: itr, minmn
   integer :: i, j, k
   integer :: col, cptr, bestr
   real(wp) :: u, bestu, bestv

   real(wp) :: eps ! minimum improvement

   integer :: prev ! number of unmatched cols on previous iteration
   integer :: nunchanged ! number of iterations where #unmatched cols has been
      ! constant

   inform%flag = 0

   ! Allocate memory
   allocate(owner(m), next(n), stat=inform%stat)
   if(inform%stat.ne.0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif

   ! Set everything as unmatched
   minmn = min(m, n)
   unmatched = minmn
   match(1:n) = 0 ! 0 = unmatched, -1 = unmatched+ineligible
   owner(1:m) = 0
   dualu(1:m) = 0
   ! dualv is set for each column as it becomes matched, otherwise we use
   ! the value supplied on input (calculated as something sensible during
   ! preprocessing)

   ! Set up monitoring of progress
   prev = -1
   nunchanged = 0

   ! Initially all columns are unmatched
   tail = n
   do i = 1, n
      next(i) = i
   end do

   ! Iterate until we run out of unmatched buyers
   eps = options%eps
   do itr = 1, options%max_iterations
      if(unmatched.eq.0) exit ! nothing left to match
      ! Bookkeeping to determine number of iterations with no change
      if(unmatched.ne.prev) nunchanged = 0
      prev = unmatched
      nunchanged = nunchanged + 1
      ! Test if we satisfy termination conditions
      if(nunchanged                  .ge. options%max_unchanged(1) .and. &
         real(minmn-unmatched)/minmn .ge. options%min_proportion(1)) exit
      if(nunchanged                  .ge. options%max_unchanged(2) .and. &
         real(minmn-unmatched)/minmn .ge. options%min_proportion(2)) exit
      if(nunchanged                  .ge. options%max_unchanged(3) .and. &
         real(minmn-unmatched)/minmn .ge. options%min_proportion(3)) exit
      ! Update epsilon scaling
      eps = min(1.0_wp, eps+1.0_wp/(n+1))
      ! Now iterate over all unmatched entries listed in next(1:tail)
      ! As we progress, build list for next iteration in next(1:insert)
      insert = 0
      do cptr = 1, tail
         col = next(cptr)
         if(match(col).ne.0) cycle ! already matched or ineligible
         if(ptr(col).eq.ptr(col+1)) cycle ! empty col (only ever fails on first
            ! iteration - not put in next(1:insert) thereafter)
         ! Find best value of a_ij - u_i for current column
         ! This occurs for i=bestr with value bestu
         ! second best value stored as bestv
         j = ptr(col)
         bestr = row(j)
         bestu = val(j) - dualu(bestr)
         bestv = -huge(bestv)
         do j = ptr(col)+1, ptr(col+1)-1
            u = val(j) - dualu(row(j))
            if( u > bestu ) then
               bestv = bestu
               bestr = row(j)
               bestu = u
            elseif( u > bestv) then
               bestv = u
            endif
         end do
         if(bestv.eq.-huge(bestv)) bestv = 0.0 ! No second best
         ! Check if matching this column gives us a net benefit
         if(bestu > 0) then
            ! There is a net benefit, match column col to row bestr
            ! if bestr was previously matched to col k, unmatch it
            dualu(bestr) = dualu(bestr) + bestu - bestv + eps
            dualv(col) = bestv - eps ! satisfy a_ij - u_i - v_j = 0
            match(col) = bestr
            unmatched = unmatched - 1
            k = owner(bestr)
            owner(bestr) = col
            if(k.ne.0) then
               ! Mark column k as unmatched
               match(k) = 0 ! unmatched
               unmatched = unmatched + 1
               insert = insert + 1
               next(insert) = k
            endif
         else
            ! No net benefit, mark col as ineligible for future consideration
            match(col) = -1 ! ineligible
            unmatched = unmatched - 1
         endif
      end do
      tail = insert
   end do
   inform%iterations = itr-1

   ! We expect unmatched columns to have match(col) = 0
   where(match(:).eq.-1) match(:) = 0

end subroutine auction_match_core

! Find a scaling through a matching-based approach using the auction algorithm
! This subroutine actually performs pre/post-processing around the call to
! auction_match_core() to actually use the auction algorithm
!
! This consists of finding a2_ij = 2*maxentry - cmax + log(|a_ij|)
! where cmax is the log of the absolute maximum in the column
! and maxentry is the maximum value of cmax-log(|a_ij|) across entire matrix
! The cmax-log(|a_ij|) term converts from max product to max sum problem and
! normalises scaling across matrix. The 2*maxentry term biases the result
! towards a high cardinality solution.
!
! Lower triangle only as input (log(half)+half->full faster than log(full))
!
subroutine auction_match(expand,m,n,ptr,row,val,rscaling,cscaling,options, &
      inform,match)
   logical, intent(in) :: expand
   integer, intent(in) :: m
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(*), intent(in) :: row
   real(wp), dimension(*), intent(in) :: val
   type(auction_options), intent(in) :: options
   type(auction_inform), intent(inout) :: inform
   real(wp), dimension(m), intent(out) :: rscaling
   real(wp), dimension(n), intent(out) :: cscaling
   integer, dimension(n), optional, intent(out) :: match

   integer, allocatable :: ptr2(:), row2(:), iw(:), perm(:)
   real(wp), allocatable :: val2(:), cmax(:)
   real(wp) :: colmax
   integer :: i,j,k,ne
   real(wp), parameter :: zero = 0.0
   real(wp) :: maxentry

   inform%flag = 0
   inform%stat = 0

   ! Reset ne for the expanded symmetric matrix
   ne = ptr(n+1)-1
   ne = 2*ne

   ! Expand matrix, drop explicit zeroes and take log absolute values
   allocate (ptr2(n+1), row2(ne), val2(ne), cmax(n), stat=inform%stat)
   if(inform%stat.ne.0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if(val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
      ! Following log is seperated from above loop to expose expensive
      ! log operation to vectorization.
      val2(ptr2(i):k-1) = log(val2(ptr2(i):k-1))
   end do
   ptr2(n+1) = k
   if(expand) then
      if(m.ne.n) then
         ! Should never get this far with a non-square matrix
         inform%flag = -99
         return
      endif
      allocate (iw(5*n), cmax(n), stat=inform%stat)
      if(inform%stat.ne.0) then
         inform%flag = ERROR_ALLOCATION
         return
      endif
      call half_to_full(n, row2, ptr2, iw, a=val2)
   endif

   ! Compute column maximums
   do i = 1, n
      if(ptr2(i+1).le.ptr2(i)) then
         ! Empty col
         cmax(i) = 0.0
         cycle
      endif
      colmax = maxval(val2(ptr2(i):ptr2(i+1)-1))
      cmax(i) = colmax
      val2(ptr2(i):ptr2(i+1)-1) = colmax - val2(ptr2(i):ptr2(i+1)-1)
   end do

   maxentry = maxval(val2(1:ptr2(n+1)-1))
   ! Use 2*maxentry+1 to prefer high cardinality matchings (+1 avoids 0 cols)
   maxentry = 2*maxentry+1
   val2(1:ptr2(n+1)-1) = maxentry - val2(1:ptr2(n+1)-1)
   !cscaling(1:n) = maxentry - cmax(1:n) ! equivalent to scale=1.0 for unmatched
   !   ! cols that core algorithm doesn't visit
   cscaling(1:n) = - cmax(1:n) ! equivalent to scale=1.0 for unmatched
      ! cols that core algorithm doesn't visit
   if(present(match)) then
      call auction_match_core(m, n, ptr2, row2, val2, match, rscaling, &
         cscaling, options, inform)
      inform%matched = count(match.ne.0)
   else
      allocate (perm(n), stat=inform%stat)
      if(inform%stat.ne.0) then
         inform%flag = ERROR_ALLOCATION
         return
      endif
      call auction_match_core(m, n, ptr2, row2, val2, perm, rscaling, &
         cscaling, options, inform)
      inform%matched = count(perm.ne.0)
   endif
   if(inform%flag.ne.0) return

   ! Calculate an adjustment so row and col scaling similar orders of magnitude
   ! and undo pre processing
   if(present(match)) then
      call match_postproc(m, n, ptr, row, val, rscaling, cscaling, maxentry, &
         cmax, inform%matched, match, inform%flag, inform%stat)
   else
      call match_postproc(m, n, ptr, row, val, rscaling, cscaling, maxentry, &
         cmax, inform%matched, perm, inform%flag, inform%stat)
   endif
end subroutine auction_match

subroutine match_postproc(m, n, ptr, row, val, rscaling, cscaling, maxentry, &
      cmax, nmatch, match, flag, st)
   integer, intent(in) :: m
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   real(wp), dimension(m), intent(inout) :: rscaling
   real(wp), dimension(n), intent(inout) :: cscaling
   real(wp), intent(in) :: maxentry
   real(wp), dimension(n), intent(in) :: cmax
   integer, intent(in) :: nmatch
   integer, dimension(n), intent(in) :: match
   integer, intent(inout) :: flag
   integer, intent(inout) :: st

   integer :: i, j
   real(wp), dimension(:), allocatable :: rmax
   real(wp) :: ravg, cavg, adjust, colmax, v

   if(m.eq.n) then
      ! Square
      ! Just perform post-processing and magnitude adjustment
      ravg = maxentry - sum(rscaling(1:m)) / m
      cavg = (-sum(cmax(1:n)) - sum(cscaling(1:n))) /n
      adjust = (ravg - cavg) / 2
      rscaling(1:m) = -rscaling(1:m) + maxentry - adjust
      cscaling(1:n) = -cscaling(1:n) - cmax(1:n) + adjust
   elseif(m.lt.n) then
      ! More columns than rows
      ! First perform post-processing and magnitude adjustment based on match
      ravg = maxentry - sum(rscaling(1:m)) / m
      cavg = 0
      do i = 1, n
         if(match(i).eq.0) cycle
         cavg = cavg + cmax(i) + cscaling(i)
      end do
      cavg = -cavg / nmatch
      adjust = (ravg - cavg) / 2
      rscaling(1:m) = -rscaling(1:m) + maxentry - adjust
      cscaling(1:n) = -cscaling(1:n) - cmax(1:n) + adjust
      ! For each unmatched col, scale max entry to 1.0
      do i = 1, n
         if(match(i).ne.0) cycle
         colmax = 0.0
         do j = ptr(i), ptr(i+1)-1
            v = abs(val(j)) * exp( rscaling(row(j)) )
            colmax = max(colmax, v)
         end do
         if(colmax.eq.0.0) then
            cscaling(i) = 0.0
         else
            cscaling(i) = log(1/colmax)
         endif
      end do
   elseif(n.lt.m) then
      ! More rows than columns
      ! Allocate some workspace
      allocate (rmax(m), stat=st)
      if(st.ne.0) then
         flag = ERROR_ALLOCATION
         return
      endif
      ! First perform post-processing and magnitude adjustment based on match
      ! also record which rows have been matched
      ravg = 0
      do i = 1, n
         if(match(i).eq.0) cycle
         ravg = ravg + rscaling(match(i))
      end do
      ravg = maxentry - ravg / nmatch
      cavg = (-sum(cmax(1:n)) - sum(cscaling(1:n))) /n
      adjust = (ravg - cavg) / 2
      rscaling(1:m) = -rscaling(1:m) + maxentry - adjust
      cscaling(1:n) = -cscaling(1:n) - cmax(1:n) + adjust
      ! Find max column-scaled value in each row from unmatched cols
      rmax(:) = 0.0
      do i = 1, n
         if(match(i).ne.0) cycle
         do j = ptr(i), ptr(i+1)-1
            v = abs(val(j)) * exp( cscaling(i) )
            rmax(row(j)) = max(rmax(row(j)), v)
         end do
      end do
      ! Calculate scaling for each row, but overwrite with correct values for
      ! matched rows, then copy entire array over rscaling(:)
      where(rmax(1:m).eq.0.0)
         rmax(1:m) = 0.0
      elsewhere
         rmax(1:m) = log( 1 / rmax(1:m) )
      endwhere
      do i = 1, n
         if(match(i).eq.0) cycle
         rmax(match(i)) = rscaling(match(i))
      end do
      rscaling(1:m) = rmax(1:m)
   endif
end subroutine match_postproc

end module spral_scaling
