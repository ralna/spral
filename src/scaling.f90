module spral_scaling
   use spral_matrix_util, only : half_to_full
   implicit none

   private
   ! Top level routines
   public :: hungarian_scale_sym, & ! Scaling using Hungarian algorithm (MC64-like)
             auction_scale_sym,   & ! Scaling using Auction algorithm
             equilib_scale_sym      ! Scaling using Equilibriation (MC77-like)
   ! Inner routines that allow calling internals
   public :: hungarian_wrapper, & ! Find a matching (with pre/post-processing)
             hungarian_match      ! Find a matching (no pre/post-processing)
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
subroutine hungarian_scale_sym(n, ptr, row, val, scaling, options, inform)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   type(hungarian_options), intent(in) :: options
   type(hungarian_inform), intent(out) :: inform

   integer :: i, flag
   integer, dimension(:), allocatable :: perm
   real(wp), dimension(:), allocatable :: cscale

   inform%flag = 0 ! Initialize to success
   
   allocate(perm(2*n), cscale(n), stat=inform%stat)
   if (inform%stat .ne. 0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif
   call hungarian_wrapper(n,ptr,row,val,perm,cscale,flag,inform%stat)
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

   do i = 1, n
      scaling(i) = exp(cscale(i))
   end do

end subroutine hungarian_scale_sym

!**************************************************************
!
! Call auction algorithm to get a scaling, then symmetrize it
!
subroutine auction_scale_sym(n, ptr, row, val, scaling, options, inform)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   type(auction_options), intent(in) :: options
   type(auction_inform), intent(inout) :: inform

   integer :: i
   real(wp), dimension(:), allocatable :: cscale

   inform%flag = 0 ! Initialize to sucess

   allocate(cscale(n), stat=inform%stat)
   if (inform%stat .ne. 0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif

   call auction_match(n, ptr, row, val, cscale, options, inform)

   do i = 1, n
      scaling(i) = exp(cscale(i))
   end do

end subroutine auction_scale_sym

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
   integer, dimension(2*n), intent(out) :: perm
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
      ! Set column permutation and row permutation to be the same
      perm(n+1:n+n) = perm(1:n)
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

   perm(n+1:n+n) = perm(1:n)

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

   integer, allocatable, dimension(:) :: JPERM ! a(jperm(j)) is entry of A for
      ! matching in column j.
   integer, allocatable, dimension(:) :: OUT ! a(out(i)) is the new entry in a
      ! on which we match going along the short path back to original col.
   integer, allocatable, dimension(:) :: PR ! pr(i) is a pointer to the next
      ! column along the shortest path back to the original column
   integer, allocatable, dimension(:) :: Q ! q(1:qlen) forms a binary heap
      ! data structure sorted by d(q(i)) value. q(low:up) is a list of rows
      ! with equal d(i) which is lower or equal to smallest in the heap.
      ! q(up:n) is a list of already visited rows.
   integer, allocatable, dimension(:) :: L ! l(:) is an inverse of q(:)

   integer :: I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP
   integer :: K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
   real(wp) :: CSP,DI,DMIN,DNEW,DQ0,VJ

   !
   ! Initialization
   !
   allocate(jperm(n), out(n), pr(n), q(n), l(n), stat=st)
   if(st.ne.0) return
   NUM = 0
   D(1:N) = 0.0
   IPERM(1:N) = 0
   JPERM(1:N) = 0
   PR(1:N) = ptr(1:N)
   L(1:N) = 0

   !
   ! Initialize U(I) using heurisitic to get an intial guess.
   ! Heuristic guaruntees that the generated partial matching is optimal
   ! on the restriction of the graph to the matched rows and columns.
   !
   U(1:N) = RINF
   do J = 1,N
      do K = ptr(J),ptr(J+1)-1
         I = row(K)
         if(val(K).gt.U(I)) cycle
         U(I) = val(K)
         IPERM(I) = J
         L(I) = K
      end do
   end do
   do I = 1,N
      J = IPERM(I)
      if(J.eq.0) cycle
      ! Row I is not empty
      IPERM(I) = 0
      if(JPERM(J).ne.0) cycle
      ! Don't choose cheap assignment from dense columns
      if(ptr(J+1)-ptr(J) .gt. N/10 .and. N.gt.50) cycle
      ! Assignment of column J to row I
      NUM = NUM + 1
      IPERM(I) = J
      JPERM(J) = L(I)
   end do
   if(NUM.eq.N) GO TO 1000
   ! Scan unassigned columns; improve assignment
   improve_assign: do J = 1,N
      ! JPERM(J) ne 0 iff column J is already assigned
      if(JPERM(J).ne.0) cycle
      K1 = ptr(J)
      K2 = ptr(J+1) - 1
      ! Continue only if column J is not empty
      if(K1.gt.K2) cycle
      I0 = row(K1)
      VJ = val(K1) - U(I0)
      K0 = K1
      do K = K1+1,K2
         I = row(K)
         DI = val(K) - U(I)
         if(DI.gt.VJ) cycle
         if(DI.ge.VJ .and. DI.ne.RINF) then
            if(IPERM(I).ne.0 .or. IPERM(I0).eq.0) cycle
         endif
         VJ = DI
         I0 = I
         K0 = K
      end do
      D(J) = VJ
      K = K0
      I = I0
      if(IPERM(I).eq.0) then
         NUM = NUM + 1
         JPERM(J) = K
         IPERM(I) = J
         PR(J) = K + 1
         cycle
      endif
      do K = K0,K2
         I = row(K)
         if(val(K)-U(I).gt.VJ) cycle
         JJ = IPERM(I)
         ! Scan remaining part of assigned column JJ
         KK1 = PR(JJ)
         KK2 = ptr(JJ+1) - 1
         if(KK1.gt.KK2) cycle
         do KK = KK1,KK2
            II = row(KK)
            if(IPERM(II).gt.0) cycle
            if(val(KK)-U(II).le.D(JJ)) then
               JPERM(JJ) = KK
               IPERM(II) = JJ
               PR(JJ) = KK + 1
               NUM = NUM + 1
               JPERM(J) = K
               IPERM(I) = J
               PR(J) = K + 1
               cycle improve_assign
            endif
         end do
         PR(JJ) = KK2 + 1
      end do
      cycle
   end do improve_assign
   if(NUM.eq.N) GO TO 1000 ! If heurisitic got a complete matching, we're done

   !
   ! Repeatedly find augmenting paths until all columns are included in the
   ! matching. At every step the current matching is optimal on the restriction
   ! of the graph to currently matched rows and columns.
   !

   ! Main loop ... each pass round this loop is similar to Dijkstra's
   ! algorithm for solving the single source shortest path problem
   D(1:N) = RINF
   L(1:N) = 0
   ISP=-1; JSP=-1 ! initalize to avoid "may be used unitialized" warning
   do JORD = 1,N

      if(JPERM(JORD).ne.0) cycle
      ! JORD is next unmatched column
      ! DMIN is the length of shortest path in the tree
      DMIN = RINF
      QLEN = 0
      LOW = N + 1
      UP = N + 1
      ! CSP is the cost of the shortest augmenting path to unassigned row
      ! row(ISP). The corresponding column index is JSP.
      CSP = RINF
      ! Build shortest path tree starting from unassigned column (root) JORD
      J = JORD
      PR(J) = -1

      ! Scan column J
      do K = ptr(J),ptr(J+1)-1
         I = row(K)
         DNEW = val(K) - U(I)
         if(DNEW.ge.CSP) cycle
         if(IPERM(I).eq.0) then
            CSP = DNEW
            ISP = K
            JSP = J
         else
            if(DNEW.lt.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
         endif
      end do
      ! Initialize heap Q and Q2 with rows held in Q(1:QLEN)
      Q0 = QLEN
      QLEN = 0
      do KK = 1,Q0
         K = Q(KK)
         I = row(K)
         if(CSP.le.D(I)) then
            D(I) = RINF
            cycle
         endif
         if(D(I).le.DMIN) then
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
         else
            QLEN = QLEN + 1
            L(I) = QLEN
            call heap_update(I,N,Q,D,L)
         endif
         ! Update tree
         JJ = IPERM(I)
         OUT(JJ) = K
         PR(JJ) = J
      end do

      do JDUM = 1,NUM

         ! If Q2 is empty, extract rows from Q
         if(LOW.eq.UP) then
            if(QLEN.eq.0) exit
            I = Q(1) ! Peek at top of heap
            if(D(I).ge.CSP) exit
            DMIN = D(I)
            ! Extract all paths that have length DMIN and store in q(low:up-1)
            do while(QLEN.gt.0)
               I = Q(1) ! Peek at top of heap
               if(D(I).gt.DMIN) exit
               i = heap_pop(QLEN,N,Q,D,L)
               LOW = LOW - 1
               Q(LOW) = I
               L(I) = LOW
            end do
         endif
         ! Q0 is row whose distance D(Q0) to the root is smallest
         Q0 = Q(UP-1)
         DQ0 = D(Q0)
         ! Exit loop if path to Q0 is longer than the shortest augmenting path
         if(DQ0.ge.CSP) exit
         UP = UP - 1

         ! Scan column that matches with row Q0
         J = IPERM(Q0)
         VJ = DQ0 - val(JPERM(J)) + U(Q0)
         do K = ptr(J),ptr(J+1)-1
            I = row(K)
            if(L(I).ge.UP) cycle
            ! DNEW is new cost
            DNEW = VJ + val(K)-U(I)
            ! Do not update D(I) if DNEW ge cost of shortest path
            if(DNEW.ge.CSP) cycle
            if(IPERM(I).eq.0) then
               ! Row I is unmatched; update shortest path info
               CSP = DNEW
               ISP = K
               JSP = J
            else
               ! Row I is matched; do not update D(I) if DNEW is larger
               DI = D(I)
               if(DI.le.DNEW) cycle
               if(L(I).ge.LOW) cycle
               D(I) = DNEW
               if(DNEW.le.DMIN) then
                  LPOS = L(I)
                  if(LPOS.ne.0) call heap_delete(LPOS,QLEN,N,Q,D,L)
                  LOW = LOW - 1
                  Q(LOW) = I
                  L(I) = LOW
               else
                  if(L(I).eq.0) then
                     QLEN = QLEN + 1
                     L(I) = QLEN
                  endif
                  call heap_update(I,N,Q,D,L) ! D(I) has changed
               endif
               ! Update tree
               JJ = IPERM(I)
               OUT(JJ) = K
               PR(JJ) = J
            endif
         end do
      end do

      ! If CSP = RINF, no augmenting path is found
      if(CSP.eq.RINF) GO TO 190
      ! Find augmenting path by tracing backward in PR; update IPERM,JPERM
      NUM = NUM + 1
      I = row(ISP)
      IPERM(I) = JSP
      JPERM(JSP) = ISP
      J = JSP
      do JDUM = 1,NUM
         JJ = PR(J)
         if(JJ.eq.-1) exit
         K = OUT(J)
         I = row(K)
         IPERM(I) = JJ
         JPERM(JJ) = K
         J = JJ
      end do

      ! Update U for rows in Q(UP:N)
      do KK = UP,N
         I = Q(KK)
         U(I) = U(I) + D(I) - CSP
      end do
190   do KK = LOW,N
         I = Q(KK)
         D(I) = RINF
         L(I) = 0
      end do
      do KK = 1,QLEN
         I = Q(KK)
         D(I) = RINF
         L(I) = 0
      end do

   end do ! End of main loop


   ! Set dual column variable in D(1:N)
1000 do J = 1,N
      K = JPERM(J)
      if(K.ne.0) then
         D(J) = val(K) - U(row(K))
      else
         D(J) = 0.0
      endif
      if(IPERM(J).eq.0) U(J) = 0.0
   end do

   ! Return if matrix has full structural rank
   if(NUM.eq.N) return

   ! Otherwise, matrix is structurally singular, complete IPERM.
   ! JPERM, OUT are work arrays
   JPERM(1:N) = 0
   K = 0
   do I = 1,N
      if(IPERM(I).eq.0) then
         K = K + 1
         OUT(K) = I
      else
         J = IPERM(I)
         JPERM(J) = I
      endif
   end do
   K = 0
   do J = 1,N
      if(JPERM(J).ne.0) cycle
      K = K + 1
      JDUM = OUT(K)
      IPERM(JDUM) = -J
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
subroutine auction_match_core(n, ptr, row, val, match, dualu, dualv, options, &
      inform)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   integer, dimension(n), intent(out) :: match
      ! match(j) = i => column j matched to row i
   real(wp), dimension(n), intent(out) :: dualu ! row dual variables
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

   integer :: itr
   integer :: i, j, k
   integer :: col, cptr, bestr
   real(wp) :: u, bestu, bestv

   real(wp) :: eps ! minimum improvement

   integer :: prev ! number of unmatched cols on previous iteration
   integer :: nunchanged ! number of iterations where #unmatched cols has been
      ! constant

   inform%flag = 0

   ! Allocate memory
   allocate(owner(n), next(n), stat=inform%stat)
   if(inform%stat.ne.0) then
      inform%flag = ERROR_ALLOCATION
      return
   endif

   ! Set everything as unmatched
   unmatched = n
   match(1:n) = 0 ! 0 = unmatched, -1 = unmatched+ineligible
   owner(1:n) = 0
   dualu(1:n) = 0
   ! dualv is set for each column as it becomes matched, otherwise we use
   ! the value supplied on input (calculated as something sensible during
   ! preprocessing)

   ! Set up monitoring of progress
   prev = n
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
      if(nunchanged          .ge. options%max_unchanged(1) .and. &
         real(n-unmatched)/n .ge. options%min_proportion(1)) exit
      if(nunchanged          .ge. options%max_unchanged(2) .and. &
         real(n-unmatched)/n .ge. options%min_proportion(2)) exit
      if(nunchanged          .ge. options%max_unchanged(3) .and. &
         real(n-unmatched)/n .ge. options%min_proportion(3)) exit
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
subroutine auction_match(n,ptr,row,val,scale,options,inform)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(*), intent(in) :: row
   real(wp), dimension(*), intent(in) :: val
   type(auction_options), intent(in) :: options
   type(auction_inform), intent(inout) :: inform
   real(wp), dimension(n), intent(out) :: scale

   integer, allocatable :: ptr2(:), row2(:), iw(:), perm(:)
   real(wp), allocatable :: val2(:), dualu(:), dualv(:), cmax(:)
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
   allocate (ptr2(n+1), row2(ne), val2(ne), perm(n), &
             iw(5*n), dualu(n), dualv(n), cmax(n), stat=inform%stat)
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
   call half_to_full(n, row2, ptr2, iw, a=val2)

   ! Compute column maximums
   do i = 1, n
      colmax = maxval(val2(ptr2(i):ptr2(i+1)-1))
      cmax(i) = colmax
      val2(ptr2(i):ptr2(i+1)-1) = colmax - val2(ptr2(i):ptr2(i+1)-1)
   end do

   maxentry = maxval(val2(1:ptr2(n+1)-1))
   ! Use 2*maxentry to prefer high cardinality matchings
   maxentry = 2*maxentry
   val2(1:ptr2(n+1)-1) = maxentry - val2(1:ptr2(n+1)-1)
   dualv(1:n) = maxentry - cmax(1:n) ! equivalent to scale=1.0 for unmatched
      ! cols that core algorithm doesn't visit
   call auction_match_core(n, ptr2, row2, val2, perm, dualu, dualv, &
      options, inform)
   if(inform%flag.ne.0) return
   inform%matched = count(perm.ne.0)

   scale(1:n) = (maxentry-dualu(1:n)-dualv(1:n)-cmax(1:n))/2
end subroutine auction_match

end module spral_scaling
