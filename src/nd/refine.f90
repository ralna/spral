module spral_nd_refine
   use spral_nd_maxflow
   use spral_nd_types
   use spral_nd_util
   implicit none

   private
   public :: refine_partition, nd_refine_edge, nd_refine_fm, nd_refine_trim, expand_partition, &
      nd_refine_block_trim

contains

!
! Refines the current partition (seemingly in place to the caller)
!
subroutine refine_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,   &
      a_n1, a_n2, partition, a_weight_1, a_weight_2, a_weight_sep, options, &
      work)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, dimension(a_n), intent(in) :: a_ptr
   integer, dimension(a_ne), intent(in) :: a_row
   integer, dimension(a_n), intent(in) :: a_weight
   integer, intent(in) :: sumweight
   integer, intent(inout) :: a_n1
   integer, intent(inout) :: a_n2
   integer, dimension(a_n), intent(inout) :: partition
   integer, intent(inout) :: a_weight_1
   integer, intent(inout) :: a_weight_2
   integer, intent(inout) :: a_weight_sep
   type(nd_options), intent(in) :: options
   integer, dimension(9*a_n+sumweight), target, intent(out) :: work

   integer :: i, k
   integer :: a_n1_new, a_n2_new
   integer, dimension(:), pointer :: partition2
   integer :: a_weight_1_new, a_weight_2_new, a_weight_sep_new
   integer :: work_ptr
   integer :: ref_options, ref_method
   logical :: imbal
   real(wp) :: balance, balance_tol
   real(wp) :: tau, tau_best

   ! Internal options on what to do
   integer, parameter :: REFINE_MAXFLOW   = 0, &
                         REFINE_TRIM      = 1, &
                         REFINE_EDGE      = 2

   balance_tol = max(1.0_wp, options%balance)
   imbal = (balance_tol.le.real(sumweight-2))

   work_ptr = 0
   partition2 => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n

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
      call cost_function(a_weight_1, a_weight_2, a_weight_sep, sumweight,&
         balance_tol, imbal, options%cost_function, tau_best)
      a_n1_new = a_n1
      a_n2_new = a_n2
      a_weight_1_new = a_weight_1
      a_weight_2_new = a_weight_2
      a_weight_sep_new = a_weight_sep
   end if

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
end subroutine refine_partition

!
! Given a partition, refine the partition to improve the (weighted) value
! of the cost function. An edge separator is formed between the input
! separator and the larger partition, and this is then minimal using
! trimming or max flow
!
subroutine nd_refine_edge(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,  &
      a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, partition, work, &
      options)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(in) :: a_ptr(a_n)
   integer, intent(in) :: a_row(a_ne)
   integer, intent(in) :: a_weight(a_n)
   integer, intent(in) :: sumweight ! Sum of weights in a_weight
   integer, intent(inout) :: a_n1 ! Size of partition 1
   integer, intent(inout) :: a_n2 ! Size of partition 2
   integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weighted
      ! size of partitions and separator
   integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
      ! list of (local) indices in partition 1; next a_n2 entries
      ! contain list of (local) entries in partition 2; entries in
      ! separator are listed at the end. This is updated to the new partition
   integer, target, intent(out) :: work(3*a_n) ! Workspace
   type (nd_options), intent(in) :: options

   integer, dimension(:), pointer :: fm_flags

   ! Initialise fm_flags(:)
   fm_flags => work(1:a_n)
   call nd_convert_partition_flags(a_n, a_n1, a_n2, partition, ND_PART1_FLAG, &
      ND_PART2_FLAG, ND_SEP_FLAG, fm_flags)

   ! Create new separator by forming edge separator between input
   ! separator and largest of P1 and P2
   call nd_move_partition(a_n, a_ne, a_ptr, a_row, a_weight, a_n1, a_n2, &
      a_weight_1, a_weight_2, a_weight_sep, ND_PART1_FLAG, &
      ND_PART2_FLAG, ND_SEP_FLAG, fm_flags)

   ! Update partition
   call nd_convert_flags_partition(a_n, a_n1, a_n2, &
      fm_flags, ND_PART1_FLAG, ND_PART2_FLAG, partition(1:a_n))
   nullify(fm_flags) ! Workspace reused from now on

   select case (options%refinement)
   case(:ND_REFINE_TRIM_FM_AUTO)
      call nd_refine_trim(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, a_n1,  &
         a_n2, a_weight_1, a_weight_2, a_weight_sep, partition, work(1:3*a_n), &
         options)
   case(ND_REFINE_MAXFLOW_BOTH:)
      call nd_refine_max_flow(a_n, a_ne, a_ptr, a_row, a_weight, a_n1, a_n2, &
         a_weight_1, a_weight_2, a_weight_sep, partition, work(1:8), options)
   end select
end subroutine nd_refine_edge

!
! Given a partition, refine the partition using FM refinement.
! 
! This is just a wrapper for nd_refine_fm_inner that partitions memory and sets
! up various internal variables
subroutine nd_refine_fm(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, a_n1, &
      a_n2, a_weight_1, a_weight_2, a_weight_sep, partition, work, options)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, dimension(a_n), intent(in) :: a_ptr
   integer, dimension(a_ne), intent(in) :: a_row
   integer, dimension(a_n), intent(in) :: a_weight
   integer, intent(in) :: sumweight ! Sum of weights in a_weight
   integer, intent(inout) :: a_n1 ! Size of partition 1
   integer, intent(inout) :: a_n2 ! Size of partition 2
   integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weighted
      ! size of partitions and separator
   integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
      ! list of (local) indices in partition 1; next a_n2 entries
      ! contain list of (local) entries in partition 2; entries in
      ! separator are listed at the end. This is updated to the new
      ! partition
   integer, target, intent(out) :: work(8*a_n+sumweight) ! Workspace
   type (nd_options), intent(in) :: options

   integer :: work_ptr
   integer, dimension(:), pointer :: flags, ipart, next, last, gain1, gain2, &
      done, head, distance
   integer :: icut, mult ! Used within FM refinement
   integer :: band
   real(wp) :: balance_tol
   logical :: imbal ! Should we check for imbalance?

   if (options%refinement_band .lt. 1) return

   balance_tol = max(1.0_wp, options%balance)
   imbal = (balance_tol.le.real(sumweight-2))
   icut = min(sumweight-1, 3*(sumweight/a_n))
   icut = min(icut, 5*maxval(a_weight))
   mult = min(icut, sumweight - icut - 1)
   band = min(options%refinement_band, a_n)

   work_ptr = 0
   flags    => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
   ipart    => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
   next     => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
   last     => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
   gain1    => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
   gain2    => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
   done     => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n
   head     => work(work_ptr+1:work_ptr+icut+mult+1)
   work_ptr = work_ptr + icut+mult+1
   distance => work(work_ptr+1:work_ptr+a_n); work_ptr = work_ptr + a_n

   call nd_convert_partition_flags(a_n, a_n1, a_n2, partition, &
      ND_PART1_FLAG, ND_PART2_FLAG, ND_SEP_FLAG, flags)
   call nd_refine_fm_inner(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, icut, &
      mult, options%cost_function, a_n1, a_n2, a_weight_1, a_weight_2,         &
      a_weight_sep, band, balance_tol, flags, ipart, next, last, gain1, gain2, &
      done, head, distance)
   call nd_convert_flags_partition(a_n, a_n1, a_n2, flags, ND_PART1_FLAG, &
      ND_PART2_FLAG, partition)

end subroutine nd_refine_fm

!
! The subroutine nd_fm_refinement uses a version of the Fiduccia-Mattheyses
! refinement algorithm on a tripartite partitioing of the nodes of a
! graph where a node in the first partition is not connected to any node
! in the second partition and any path between nodes in partition 1 and
! partition 2 must go through a node in the cutset (partition 3).
! The intention of the algorithm is to reduce f(P1,P2,P3), where
! f(P1,P2,P3) = |P3|/(|P1||P2|) if min(|P1|,|P2|)/max(|P1|,|P2|) >= balance_tol
!             =  sumweight - 2 + max(|P1|,|P2|)/min(|P1|,|P2|), otherwise
! This is a banded version so only nodes a distance of at most band from the
! input separator can be moved into the new separator
!
subroutine nd_refine_fm_inner(n, a_ne, ptr, col, weight, sumweight, icut, &
      mult, costf, a_n1, a_n2, wnv1, wnv2, wns, band, balance_tol, flags, &
      ipart, next, last, gain1, gain2, done, head, dist)
   integer, intent(in) :: n
   integer, intent(in) :: a_ne
   integer, dimension(n), intent(in) :: ptr
   integer, dimension(a_ne), intent(in) :: col
   ! The array weight is used to hold a weight on the vertices indicating
   ! how many vertices from the finer graphs have been combined into the
   ! current coarse graph vertex.
   integer, dimension(n), intent(in) :: weight
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
   real(wp), intent(in) :: balance_tol ! balance_tol to determine whether
      ! partition is balanced

   ! flags holds a list of nodes stating which partition node i is in.
   ! The whole point of this routine is to return a revised partition
   ! with better properties.  Normally less nodes in the cutset while
   ! maintaining
   ! a balance between the number of nodes in the two components.
   ! flags(i) .eq. ND_PART1_FLAG : i is in partition 1
   ! flags(i) .eq. ND_PART2_FLAG : i is in partition 2
   ! flags(i) .eq. ND_SEP_FLAG   : i is in separator/cutset
   integer, intent(inout) :: flags(n)
   ! info holds parameters giving information about the performance of the
   ! subroutine
   integer, intent(out) :: ipart(n), next(n), last(n)
   integer, intent(out) :: gain1(n), gain2(n), done(n), head(-mult:icut)
   integer, intent(out) :: dist(n)

   ! Number nodes in each partition
   integer :: nv1, ns, nv2, inv1, inv2, ins
   ! Weighted nodes in each partition
   integer :: winv1, winv2, wins
   integer :: i, j, ii, jj, eye, k
   integer :: inn, outer
   integer :: move, ming, gain, old_gain, inext, idummy
   integer :: first, tail
   real(wp) :: eval, evalc, evalo, eval1, eval2
   logical :: imbal, found

   imbal = (balance_tol.le.real(sumweight-2))

   ! Set-up distance to hold (min) distance of node from separator. Only find
   ! distance for nodes within band distance from separator
   first = 0
   tail = 0
   do i = 1, n
      if (flags(i).eq.ND_SEP_FLAG) then
         dist(i) = 0
         if (first.eq.0) then
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

   do while (first.ne.0)
      j = dist(first)
      if (j.eq.band-1) then
         do i = ptr(first), nd_get_ptr(first+1, n, a_ne, ptr)-1
            if (dist(col(i)).eq.-2) then
               k = col(i)
               dist(k) = j + 1
            end if
         end do
      else
         do i = ptr(first), nd_get_ptr(first+1, n, a_ne, ptr)-1
            if (dist(col(i)).eq.-2) then
               k = col(i)
               dist(k) = j + 1
               next(tail) = k
               tail = k
            end if
         end do
      end if
      if (first.eq.tail) then
         first = 0
      else
         k = next(first)
         first = k
      end if
   end do
   next(1:n) = 0

   ! nv1,nv2,ns are the number of nodes in partitions 1, 2 and the cutset in
   ! the current partition
   ! inv1,inv2,ins,ipart are the equivalent quantities within the inner loop
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
   call cost_function(wnv1+1, wnv2+1, wns, sumweight, balance_tol, imbal, &
      costf, evalc)

   ! icut is set to limit search in inner loop .. may later be a parameter
   ! we allow gains of up to max(weight)*5

   head(-mult:icut) = 0

   ! Set up doubly linked list linking nodes with same gain and headers
   ! to starts (up to cut off value icut)

   ! Compute gains for nodes in cutset
   ming = sumweight
   do i = 1, n
      if (flags(i).eq.ND_SEP_FLAG) then
         ! Node i is in cutset
         ! gain1(i) is change to cutset size if node i is moved to partition 1.
         ! gain2(i) is change to cutset size if node i is moved to partition 2.
         ! Run through all neighbours of node i to see what loss/gain is if
         ! node i is moved

         call compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,flags)
         gain = max(-mult,min(gain1(i),gain2(i)))

         if (gain.lt.ming) ming = gain
         if (gain.gt.icut) cycle
         ! New node is put at head of list
         call add_to_list(n,mult,icut,next,last,head,i,gain)
      end if
   end do

   ! Initilialization finished.  Now perform F-M algorithm in two loops.
   ! In each inner loop we choose the best obtained and if the evaluation
   ! function for this is better than previous best we perform another inner
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
      inner: &
      do inn = 1, n
         ! Choose best eligible move
         found = .false.
         find_best_move: &
         do idummy = 1, n
            do gain = ming, icut
               if (head(gain).ne.0) exit
            end do
            if (gain.gt.icut) exit inner

            ! Now cycle through nodes of least gain
            ! Currently inefficient because of re-searching linked list
            inext = head(gain)
            do
               i = inext
               if (i.eq.0) cycle find_best_move
               ! Use node if it has not been considered already
               if (done(i).lt.outer) then
                  found = .true.
                  exit find_best_move
               end if
               inext = next(i)
            end do
         end do find_best_move
         if(.not.found) exit inner
         ! Node i has been selected as the best eligible node
         ! Set flag so only considered once in this pass
         done(i) = outer
         ! As i will not be chosen again in this pass, remove from list
         call remove_from_list(n, mult, icut, next, last, head, i, gain)
         ! Move the node to the appropriate partition and reset partition
         ! information
         ! We will try both weighted and unweighted

         if (wnv1.eq.0 .and. wnv2.gt.0) then
            ! Move node i to partition 1
            move = ND_PART1_FLAG
            inv1 = inv1 + 1
            winv1 = winv1 + weight(i)
            ins = ins - 1
            wins = wins - weight(i)
         else if (wnv2.eq.0 .and. wnv1.gt.0) then
            ! Move node i to partition 2
            move = ND_PART2_FLAG
            inv2 = inv2 + 1
            winv2 = winv2 + weight(i)
            ins = ins - 1
            wins = wins - weight(i)
         else
            call cost_function(winv1+weight(i), winv2+1-gain1(i)-weight(i), &
               wins+gain1(i)-1, sumweight, balance_tol, imbal, costf, eval1)
            call cost_function(winv1+1-gain2(i)-weight(i), winv2+weight(i), &
               wins+gain2(i)-1, sumweight, balance_tol, imbal, costf, eval2)
            if ((eval1.lt.eval2).or.((eval1.eq.eval2).and.(wnv1.lt.wnv2))) then
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
         do jj = ptr(i), nd_get_ptr(i+1, n, a_ne, ptr)-1
            j = col(jj)
            ! Check which partition node j is in and take appropriate action
            if (ipart(j).eq.move) cycle
            ! If node j is in cutset, update its gain value
            if (ipart(j).eq.ND_SEP_FLAG) then
               ! If it has already been chosen in this pass just skip it
               if (done(j).eq.outer .or. dist(j).eq.-2) cycle
               ! old_gain is present gain

               old_gain = max(-mult,min(gain1(j),gain2(j)))

               if (move.eq.ND_PART1_FLAG) gain2(j) = gain2(j) + weight(i)
               if (move.eq.ND_PART2_FLAG) gain1(j) = gain1(j) + weight(i)
               gain = max(-mult,min(gain1(j),gain2(j)))

               if (old_gain.eq.gain) cycle
               ! Remove from old list
               if (old_gain.le.icut) then
                  call remove_from_list(n, mult, icut, next, last, head, j, &
                     old_gain)
               end if
               ! gain has changed so move to new linked list if less than icut
               if (gain.le.icut) then
                  ! Reset ming if necessary
                  if (gain.lt.ming) ming = gain
                  call add_to_list(n,mult,icut,next,last,head,j,gain)
               end if
            end if
            if (ipart(j).eq.2-move) then
               ! We have a new node in the cutset.
               ipart(j) = ND_SEP_FLAG
               ! Compute gains for this new node in the cutset and place in
               ! linked list
               ! We intentionally did not do this earlier but we do now
               ! [maybe not since won't access this node again in this pass]
               ! We use done array to record this but probably not necessary
               ! as not put in head linked list so won't be accessed
               ! First check that it was not earlier moved from cutset
               if (done(j).ne.outer .and. dist(j).ne.-2) then
                  ! Compute gain
                  call compute_gain(n, a_ne, ptr, col, gain1, gain2, weight, &
                     j, ipart)
                  gain = max(-mult,min(gain1(j),gain2(j)))
                  ! !! Just added this
                  if (gain.lt.ming) ming = gain
                  ! Add to  list
                  if (gain.le.icut) then
                     call add_to_list(n,mult,icut,next,last,head,j,gain)
                  end if
               end if
               ! Update partition and gain of any nodes in cutset connected
               ! to node j
               ins = ins + 1
               wins = wins + weight(j)
               if (move.eq.ND_PART1_FLAG) then
                  inv2 = inv2 - 1
                  winv2 = winv2 - weight(j)
               end if
               if (move.eq.ND_PART2_FLAG) then
                  inv1 = inv1 - 1
                  winv1 = winv1 - weight(j)
               end if
               ! Check neighbours of j since any in cut set will have gain
               ! changed
               do ii = ptr(j), nd_get_ptr(j+1, n, a_ne, ptr)-1
                  eye = col(ii)
                  if (ipart(eye).ne.ND_SEP_FLAG) cycle
                  if (dist(eye).eq.-2) cycle
                  ! Neighbour is in cutset. Recompute gain and insert in 
                  ! linked list.
                  if (done(eye).eq.outer) cycle
                  ! old_gain is present gain
                  old_gain = max(-mult,min(gain1(eye),gain2(eye)))

                  if (move.eq.ND_PART1_FLAG) &
                     gain1(eye) = gain1(eye) - weight(j)
                  if (move.eq.ND_PART2_FLAG) &
                     gain2(eye) = gain2(eye) - weight(j)
                  ! gain is new gain
                  gain = max(-mult,min(gain1(eye),gain2(eye)))
                  if (old_gain.eq.gain) cycle
                  ! Remove from old list
                  if (old_gain.le.icut) &
                     call remove_from_list(n, mult, icut, next, last, head, &
                        eye, old_gain)
                  ! gain has changed so move to new linked list if < icut
                  if (gain.le.icut) then
                     ! Reset ming if necessary
                     if (gain.lt.ming) ming = gain
                     call add_to_list(n,mult,icut,next,last,head,eye,gain)
                  end if
               end do
            end if
            ! end of neighbours loop
         end do

         ! Evaluate new partition
         call cost_function(winv1+1, winv2+1, wins, sumweight, balance_tol, &
            imbal, costf, eval)
         ! Compare this with best so far in inner loop and store partition
         ! information if it is the best
         if (inv1*inv2.gt.0 .and. nv1*nv2.eq.0) then
            ! Might have to store gains and who is in the cutset
            evalc = eval
            nv1 = inv1
            nv2 = inv2
            ns = ins
            wnv1 = winv1
            wnv2 = winv2
            wns = wins
            flags = ipart
         else if (eval.lt.evalc .and. (inv1*inv2.gt.0)) then
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
      end do inner

      ! Leave loop if inner loop has not found better partition
      if (evalc.ge.(1.0-1.0/(log(real(sumweight))**2.3))*evalo) exit
      ! Otherwise we reset evalo and go back to inner loop
      evalo = evalc
      ! Recompute gains for this new partition
      ! Compute gains for nodes in cutset
      ! This is very inefficient but is in now to test functionality
      head(-mult:icut) = 0
      ming = icut + 1
      do i = 1, n
         if (flags(i).ne.ND_SEP_FLAG) cycle
         if (dist(i).eq.-2) cycle
         ! Node i is in cutset
         ! gain1(i) is change to cutset size if node i is moved to partition 1.
         ! gain2(i) is change to cutset size if node i is moved to partition 2.
         ! Run through all neighbours of node i to see what loss/gain is if node
         ! i is moved
         call compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,flags)
         ! Recalculate doubly linked list linking nodes with same gain and
         ! headers to starts (up to cut off value icut)
         ! Initialize array done that flags that a node has been considered in
         ! an inner loop pass
         gain = max(-mult,min(gain1(i),gain2(i)))
         ! gain = min(gain1(i),gain2(i))
         if (gain.gt.icut) cycle
         if (gain.lt.ming) ming = gain
         ! New node is put at head of list
         call add_to_list(n,mult,icut,next,last,head,i,gain)
      end do
      ! End of outer loop
   end do
   a_n1 = nv1
   a_n2 = nv2

end subroutine nd_refine_fm_inner

subroutine remove_from_list(n, mult, icut, next, last, head, irm, ig)
   integer, intent(in) :: n, mult, icut
   integer, intent(inout) :: next(n), last(n), head(-mult:icut)
   integer, intent(in) :: irm, ig
   integer :: inext, ilast

   inext = next(irm)
   ilast = last(irm)
   if (ilast.eq.0) then
      head(ig) = inext
      if (inext.ne.0) last(inext) = 0
   else
      next(ilast) = inext
      if (inext.ne.0) last(inext) = ilast
   end if
end subroutine remove_from_list

subroutine add_to_list(n, mult, icut, next, last, head, irm, ig)
   integer, intent(in) :: n, mult, icut ! order matrix
   integer, intent(inout) :: next(n), last(n), head(-mult:icut)
   integer, intent(in) :: irm, ig
   integer :: inext

   inext = head(ig)
   head(ig) = irm
   next(irm) = inext
   if (inext.ne.0) last(inext) = irm
   last(irm) = 0
end subroutine add_to_list

subroutine compute_gain(n, a_ne, ptr, col, gain1, gain2, weight, i, partit)
   integer, intent(in) :: n, a_ne, ptr(n), col(a_ne), weight(n)
   integer, intent(inout) :: gain1(n), gain2(n)
   integer, intent(in) :: i, partit(:)
   integer :: j, jj

   ! Initialize gain ... knowing node i will be removed from cutset
   ! The +1 is to give identical result to previous code when unit weights
   gain1(i) = -weight(i) + 1
   gain2(i) = -weight(i) + 1
   do jj = ptr(i), nd_get_ptr(i+1, n, a_ne, ptr)-1
      j = col(jj)
      ! Check which partition node j is in and adjust gain array
      ! appropriately
      if (partit(j).eq.ND_PART1_FLAG) then
         gain2(i) = gain2(i) + weight(j)
      end if
      if (partit(j).eq.ND_PART2_FLAG) then
         gain1(i) = gain1(i) + weight(j)
      end if
   end do
end subroutine compute_gain

!
! Given a partition, trim the partition to make it minimal
!
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
  real(wp) :: balance_tol

  balance_tol = max(1.0_wp,options%balance)
  imbal = (balance_tol.le.real(sumweight-2))

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
    if (j.lt.a_n) then
      k = a_ptr(j+1) - 1
    else
      k = a_ne
    end if
    do l = a_ptr(j), k
      m = a_row(l)
      if (work(work_part+m).eq.ND_PART1_FLAG) then
        next1 = .true.
      else if (work(work_part+m).eq.ND_PART2_FLAG) then
        next2 = .true.
      end if
    end do
    if ((next1 .and. .not. next2) .or. ( .not. next1 .and. .not. next2 &
        .and. a_n1.eq.0)) then
      ! Add to list 1
      if (head1.eq.0) then
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
        next2 .and. a_n2.eq.0)) then
      ! Add to list 2
      if (head2.eq.0) then
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

  do while (head1.gt.0 .or. head2.gt.0)
    if (head1.gt.0 .and. head2.gt.0) then
      w1 = a_weight(head1)
      w2 = a_weight(head2)
      call cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1, &
        sumweight,balance_tol,imbal,options%cost_function,t1)
      call cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2, &
        sumweight,balance_tol,imbal,options%cost_function,t2)

      if (t1.lt.t2) then
        go to 10
      else
        go to 20
      end if

    else if (head1.gt.0) then
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
    if (i.lt.a_n) then
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
        if (head1.eq.0) then
          head1 = j
        end if

      case (sep2)
        ! Remove from list 2
        p = work(work_prev+j)
        q = work(work_next+j)

        if (j.ne.head2 .and. j.ne.tail2) then
          work(work_prev+q) = p
          work(work_next+p) = q
          work(work_prev+j) = 0
          work(work_next+j) = 0
        else if (j.ne.head2 .and. j.eq.tail2) then
          work(work_next+p) = 0
          work(work_prev+j) = 0
          tail2 = p
        else if (j.ne.tail2 .and. j.eq.head2) then
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
    if (i.lt.a_n) then
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
        if (head2.eq.0) then
          head2 = j
        end if

      case (sep1)
        ! Remove from list 1
        p = work(work_prev+j)
        q = work(work_next+j)

        if (j.ne.head1 .and. j.ne.tail1) then
          work(work_prev+q) = p
          work(work_next+p) = q
          work(work_prev+j) = 0
          work(work_next+j) = 0
        else if (j.ne.head1 .and. j.eq.tail1) then
          work(work_next+p) = 0
          work(work_prev+j) = 0
          tail1 = p
        else if (j.ne.tail1 .and. j.eq.head1) then
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
    if (work(work_part+j).eq.ND_SEP_FLAG) then
      ! j is not on the boundary
      if (a_weight_1.lt.a_weight_2) then
        ! Move j into partition 1
        work(work_part+j) = ND_PART1_FLAG
        a_n1 = a_n1 + 1
        a_weight_1 = a_weight_1 + a_weight(j)
        a_weight_sep = a_weight_sep - a_weight(j)

        head1 = j
        tail1 = j
        do while (head1.gt.0)
          q = head1
          if (q.lt.a_n) then
            k = a_ptr(q+1) - 1
          else
            k = a_ne
          end if
          do l = a_ptr(q), k
            p = a_row(l)
            if (work(work_part+p).eq.ND_SEP_FLAG) then
              work(work_part+p) = ND_PART1_FLAG
              a_n1 = a_n1 + 1
              a_weight_1 = a_weight_1 + a_weight(p)
              a_weight_sep = a_weight_sep - a_weight(p)
              work(work_next+tail1) = p
              tail1 = p
            end if
          end do
          if (head1.eq.tail1) then
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

        do while (head2.gt.0)
          q = head2
          if (q.lt.a_n) then
            k = a_ptr(q+1) - 1
          else
            k = a_ne
          end if
          do l = a_ptr(q), k
            p = a_row(l)
            if (work(work_part+p).eq.ND_SEP_FLAG) then
              work(work_part+p) = ND_PART2_FLAG
              a_n2 = a_n2 + 1
              a_weight_2 = a_weight_2 + a_weight(p)
              a_weight_sep = a_weight_sep - a_weight(p)
              work(work_next+tail2) = p
              tail2 = p
            end if
          end do
          if (head2.eq.tail2) then
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
  real(wp) :: balance_tol

  balance_tol = max(1.0_wp,options%balance)
  imbal = (balance_tol.le.real(sumweight-2))

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
    if (j.lt.a_n) then
      k = a_ptr(j+1) - 1
    else
      k = a_ne
    end if
    do l = a_ptr(j), k
      m = a_row(l)
      if (work(work_part+m).eq.ND_PART1_FLAG) then
        next1 = .true.
      else if (work(work_part+m).eq.ND_PART2_FLAG) then
        next2 = .true.
      end if
    end do
    if (next1) then
      ! Add to list 1
      if (head1.eq.0) then
        head1 = j
      else
        work(work_next1+tail1) = j
      end if
      tail1 = j
      work(work_level1+j) = 1
    end if
    if (next2) then
      ! Add to list 2
      if (head2.eq.0) then
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
  do while (l1.gt.0)
    if (l1.lt.a_n) then
      k = a_ptr(l1+1) - 1
    else
      k = a_ne
    end if
    do l = a_ptr(l1), k
      m = a_row(l)
      if (work(work_part+m).eq.ND_SEP_FLAG .and. &
          work(work_level1+m).eq.0) then
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
  do while (l1.gt.0)
    if (l1.lt.a_n) then
      k = a_ptr(l1+1) - 1
    else
      k = a_ne
    end if
    do l = a_ptr(l1), k
      m = a_row(l)
      if (work(work_part+m).eq.ND_SEP_FLAG .and. &
          work(work_level2+m).eq.0) then
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
    if (work(work_level2+j).eq.0) then
      work(work_level2+j) = maxlevel2 + 1
    else
      if (work(work_level1+j).eq.0) then
        work(work_level1+j) = maxlevel1 + 1
      end if
    end if
  end do

  ! Trim the separator
  currlevel1 = 1
  currlevel2 = 1
  l1 = head1
  l2 = head2
  do while (currlevel1.le.maxlevel1 .or. currlevel2.le.maxlevel2)
    if (currlevel1.gt.maxlevel1) then
      t1 = huge(1.0_wp)
    else
      w1 = 0
      j = l1
      do while (work(work_level1+j).eq.currlevel1)
        if (work(work_level2+j).gt.currlevel2) then
          w1 = w1 + a_weight(j)
        end if
        j = work(work_next1+j)
        if (j.eq.0) exit
      end do
      if (w1.eq.0) then
        currlevel1 = currlevel1 + 1
        l1 = j
        cycle
      else
        call cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1, &
          sumweight,balance_tol,imbal,options%cost_function,t1)
      end if
    end if

    if (currlevel2.gt.maxlevel2) then
      t2 = huge(1.0_wp)
    else
      w2 = 0
      j = l2
      do while (work(work_level2+j).eq.currlevel2)
        if (work(work_level1+j).gt.currlevel1) then
          w2 = w2 + a_weight(j)
        end if
        j = work(work_next2+j)
        if (j.eq.0) exit
      end do
      if (w2.eq.0) then
        currlevel2 = currlevel2 + 1
        l2 = j
        cycle
      else
        call cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2, &
          sumweight,balance_tol,imbal,options%cost_function,t2)
      end if
    end if

    ! Add entries to relevant partition and update a_n1, a_n2 etc
    if (t1.lt.t2) then
      j = l1
      do while (work(work_level1+j).eq.currlevel1)
        if (work(work_level2+j).gt.currlevel2) then
          work(work_part+j) = ND_PART1_FLAG
          a_n1 = a_n1 + 1
        end if
        j = work(work_next1+j)
        if (j.eq.0) exit
      end do
      a_weight_1 = a_weight_1 + w1
      a_weight_sep = a_weight_sep - w1
      l1 = j
      if (j.eq.0) then
        currlevel1 = maxlevel1 + 1
      else
        currlevel1 = (work(work_level1+l1))
      end if

    else
      j = l2
      do while (work(work_level2+j).eq.currlevel2)
        if (work(work_level1+j).gt.currlevel1) then
          work(work_part+j) = ND_PART2_FLAG
          a_n2 = a_n2 + 1
        end if
        j = work(work_next2+j)
        if (j.eq.0) exit
      end do
      a_weight_2 = a_weight_2 + w2
      a_weight_sep = a_weight_sep - w2
      l2 = j
      if (j.eq.0) then
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
  if (a_weight_1.gt.a_weight_2) then
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
    if (flags(j).ne.flag_sep) cycle
    ! j is in initial separator
    if (j.eq.a_n) then
      k = a_ne
    else
      k = a_ptr(j+1) - 1
    end if
    stay_sep = .false.
    do l = a_ptr(j), k
      if (flags(a_row(l)).eq.part_no) then
        stay_sep = .true.
        if (part_no.eq.flag_1) then
          if (a_n1_temp.gt.1) then
            a_n1_temp = a_n1_temp - 1
            flags(a_row(l)) = -1
            a_weight_sep_temp = a_weight_sep_temp + a_weight(a_row(l))
            a_weight_1_temp = a_weight_1_temp - a_weight(a_row(l))
          end if
        else
          if (a_n2_temp.gt.1) then
            a_n2_temp = a_n2_temp - 1
            flags(a_row(l)) = -1
            a_weight_sep_temp = a_weight_sep_temp + a_weight(a_row(l))
            a_weight_2_temp = a_weight_2_temp - a_weight(a_row(l))
          end if
        end if
      else if (flags(a_row(l)).eq.-1) then
        stay_sep = .true.
      end if
    end do
    if ( .not. stay_sep) then
      if (part_no.eq.flag_1) then
        flags(j) = flag_2
      else
        flags(j) = flag_1
      end if
      a_weight_sep_temp = a_weight_sep_temp - a_weight(j)
      if (part_no.eq.flag_1) then
        a_n2_temp = a_n2_temp + 1
        a_weight_2_temp = a_weight_2_temp + a_weight(j)
      else
        a_n1_temp = a_n1_temp + 1
        a_weight_1_temp = a_weight_1_temp + a_weight(j)
      end if
    end if
  end do

  do j = 1, a_n
    if (flags(j).eq.-1) then
      flags(j) = flag_sep
    end if
  end do

  a_n1 = a_n1_temp
  a_n2 = a_n2_temp
  a_weight_1 = a_weight_1_temp
  a_weight_2 = a_weight_2_temp
  a_weight_sep = a_weight_sep_temp

end subroutine nd_move_partition


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
    if (j.eq.a_n) then
      k = a_ne
    else
      k = a_ptr(j+1) - 1
    end if
    do l = a_ptr(j), k
      m = a_row(l)
      if (work(work_part+m).eq.ND_PART1_FLAG .and. a_n1.gt.1) then
        ! if (side .eq. ND_PART1_FLAG .or. side .eq. ND_SEP_FLAG)
        ! then
        work(work_part+m) = ND_SEP_FLAG
        a_n1 = a_n1 - 1
        w = a_weight(m)
        a_weight_1 = a_weight_1 - w
        a_weight_sep = a_weight_sep + w
        ! end if
      else if (work(work_part+m).eq.ND_PART2_FLAG .and. a_n2.gt.1) then
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

end module spral_nd_refine
