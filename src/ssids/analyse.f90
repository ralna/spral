! (c) STFC 2010-2013
! Author: Jonathan Hogg
!
! Originally based on HSL_MA97 v2.2.0
module spral_ssids_analyse
   use, intrinsic :: iso_c_binding
   use spral_core_analyse, only : basic_analyse
   use spral_cuda, only : detect_gpu
   use spral_hw_topology, only : guess_topology, numa_region
   use spral_pgm, only : writePPM
   use spral_ssids_akeep, only : ssids_akeep
   use spral_ssids_cpu_subtree, only : construct_cpu_symbolic_subtree
   use spral_ssids_gpu_subtree, only : construct_gpu_symbolic_subtree
   use spral_ssids_datatypes
   use spral_ssids_inform, only : ssids_inform, ssids_print_flag
   implicit none

   private
   public :: analyse_phase,   & ! Calls core analyse and builds data strucutres
             check_order,     & ! Check order is a valid permutation
             expand_pattern,  & ! Specialised half->full matrix conversion
             expand_matrix      ! Specialised half->full matrix conversion

contains

!****************************************************************************

!
! Given lower triangular part of A held in row and ptr, expand to
! upper and lower triangular parts (pattern only). No checks.
!
! Note: we do not use half_to_full here to expand A since, if we did, we would
! need an extra copy of the lower triangle into the full structure before
! calling half_to_full
!
subroutine expand_pattern(n,nz,ptr,row,aptr,arow)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: nz
   integer, intent(in) :: ptr(n+1)
   integer, intent(in) :: row(nz)
   integer, intent(out) :: aptr(n+1)
   integer, intent(out) :: arow(2*nz)

   integer :: i,j,k

   ! Set aptr(j) to hold no. nonzeros in column j
   aptr(:) = 0
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         aptr(i) = aptr(i) + 1
         if (j.eq.i) cycle
         aptr(j) = aptr(j) + 1
      end do
   end do

   ! Set aptr(j) to point to where row indices will end in arow
   do j = 2, n
      aptr(j) = aptr(j-1) + aptr(j)
   end do
   aptr(n+1) = aptr(n) + 1

   ! Fill arow and aptr
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         arow(aptr(i)) = j
         aptr(i) = aptr(i) - 1
         if (j.eq.i) cycle
         arow(aptr(j)) = i
         aptr(j) = aptr(j) - 1
      end do
   end do
   do j = 1,n
      aptr(j) = aptr(j) + 1
   end do
end subroutine expand_pattern

!****************************************************************************
!
! Given lower triangular part of A held in row, val and ptr, expand to
! upper and lower triangular parts.

subroutine expand_matrix(n,nz,ptr,row,val,aptr,arow,aval)

   integer, intent(in)   :: n ! order of system
   integer, intent(in)   :: nz
   integer, intent(in)   :: ptr(n+1)
   integer, intent(in)   :: row(nz)
   real(wp), intent(in)  :: val(nz)
   integer, intent(out)  :: aptr(n+1)
   integer, intent(out)  :: arow(2*nz)
   real(wp), intent(out) :: aval(2*nz)

   integer :: i,j,k,ipos,jpos
   real(wp) :: atemp

   ! Set aptr(j) to hold no. nonzeros in column j
   aptr(:) = 0
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         aptr(i) = aptr(i) + 1
         if (j.eq.i) cycle
         aptr(j) = aptr(j) + 1
      end do
   end do

   ! Set aptr(j) to point to where row indices will end in arow
   do j = 2, n
      aptr(j) = aptr(j-1) + aptr(j)
   end do
   aptr(n+1) = aptr(n) + 1

   ! Fill arow, aval and aptr
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         atemp = val(k)
         ipos = aptr(i)
         arow(ipos) = j
         aval(ipos) = atemp
         aptr(i) = ipos - 1
         if (j.eq.i) cycle
         jpos = aptr(j)
         arow(jpos) = i
         aval(jpos) = atemp
         aptr(j) = jpos - 1
      end do
   end do
   do j = 1,n
      aptr(j) = aptr(j) + 1
   end do

end subroutine expand_matrix

!****************************************************************************
!
! This routine requires the LOWER triangular part of A
! to be held in CSC format.
! The user has supplied a pivot order and this routine checks it is OK
! and returns an error if not. Also sets perm, invp.
!
subroutine check_order(n, order, invp, akeep, options, inform)
    integer, intent(in) :: n ! order of system
    integer, intent(inout) :: order(:)
      ! If i is used to index a variable, |order(i)| must
      ! hold its position in the pivot sequence. If 1x1 pivot i required,
      ! the user must set order(i)>0. If a 2x2 pivot involving variables
      ! i and j is required, the user must set
      ! order(i)<0, order(j)<0 and |order(j)| = |order(i)|+1.
      ! If i is not used to index a variable, order(i) must be set to zero.
      ! !!!! In this version, signs are reset to positive value
   integer, intent(out) :: invp(n)
      ! Used to check order and then holds inverse of perm.
   type(ssids_akeep), intent(inout) :: akeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(inout) :: inform

   character(50)  :: context ! Procedure name (used when printing).

   integer :: i, j
   integer :: nout  ! stream for error messages

   context = 'ssids_analyse'
   nout = options%unit_error
   if (options%print_level < 0) nout = -1

   if (size(order) < n) then
      ! Order is too short
      inform%flag = SSIDS_ERROR_ORDER
      akeep%flag = inform%flag
      call ssids_print_flag(inform, nout, context)
      return
   end if

   ! initialise
   invp(:) = 0

   do i = 1,n
      order(i) = abs(order(i))
   end do
     
   ! Check user-supplied order and copy the absolute values to invp.
   ! Also add up number of variables that are not used (null rows)
   do i = 1, n
      j = order(i)
      if (j.le.0 .or. j.gt.n) exit ! Out of range entry
      if (invp(j) .ne. 0) exit ! Duplicate found
      invp(j) = i
   end do
   if (i-1 .ne. n) then
      inform%flag = SSIDS_ERROR_ORDER
      akeep%flag = inform%flag
      call ssids_print_flag(inform,nout,context)
      return
   end if
end subroutine check_order

!****************************************************************************

!
! Given an elimination tree, try and split it into at least min_npart parts
! of size at most max_flops.
!
! Parts are returned as contigous ranges of nodes. Part i consists of nodes
! part(i):part(i+1)-1
!
! FIXME: This really needs thought through in more detail for best performance
! it may be more sensible to come down from the top?
subroutine find_subtree_partition(nnodes, sptr, sparent, rptr, topology, &
      min_npart, max_flops, cpu_gpu_ratio, nparts, part, exec_loc, &
      contrib_ptr, contrib_idx, contrib_dest)
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   type(numa_region), dimension(:), intent(in) :: topology
   integer, intent(in) :: min_npart
   integer(long), intent(in) :: max_flops
   real, intent(in) :: cpu_gpu_ratio
   integer, intent(out) :: nparts
   integer, dimension(:), allocatable, intent(inout) :: part
   integer, dimension(:), allocatable, intent(inout) :: exec_loc ! 0=cpu, 1=gpu
   integer, dimension(:), allocatable, intent(inout) :: contrib_ptr
   integer, dimension(:), allocatable, intent(inout) :: contrib_idx
   integer, dimension(:), allocatable, intent(inout) :: contrib_dest

   integer :: i, j, k
   integer(long) :: jj, target_flops
   integer :: m, n, node
   integer(long), dimension(:), allocatable :: flops
   real :: cpu_flops, gpu_flops, total_flops
   integer :: st

   ! FIXME: stat parameters

   !print *, "min_npart = ", min_npart, " max_flops = ", max_flops
   ! Count flops below each node
   allocate(flops(nnodes+1))
   flops(:) = 0
   !print *, "There are ", nnodes, " nodes"
   do node = 1, nnodes
      m = int(rptr(node+1)-rptr(node))
      n = sptr(node+1)-sptr(node)
      do jj = m-n+1, m
         flops(node) = flops(node) + jj**2
      end do
      j = sparent(node)
      flops(j) = flops(j) + flops(node)
      !print *, "Node ", node, "parent", j, " flops ", flops(node)
   end do
   !print *, "Total flops ", flops(nnodes+1)

   ! Split into parts of at most target_flops. Each part is a subtree, not a
   ! subforest (unless it contains the artifical root node).
   target_flops = min(flops(nnodes+1) / min_npart, max_flops)
   !target_flops = flops(nnodes+1) ! FIXME: rm
   !print *, "Target flops ", target_flops
   total_flops = real(flops(nnodes+1))
   gpu_flops = 0.0
   cpu_flops = 0.0
   nparts = 0
   allocate(part(nnodes+1)) ! big enough to be safe FIXME: make smaller or alloc elsewhere?
   part(1) = 1
   node = 1
   do while(node.lt.nnodes+1)
      ! Head up from node until parent has too many flops
      j = node
      do while(j.lt.nnodes)
         if(flops(sparent(j)).gt.target_flops) exit ! Nodes node:j form a part
         j = sparent(j)
      end do
      ! Record node:j as a part
      nparts = nparts + 1
      part(nparts+1) = j+1
      node = j + 1
      if(j.eq.nnodes+1) exit ! done, root node is incorporated
      ! Remove subtree rooted at (node-1) from flops count of nodes above it
      do
         j = sparent(j)
         if(j.ge.nnodes+1) exit ! at top of tree
         flops(j) = flops(j) - flops(node-1)
      end do
   end do
   part(nparts+1) = nnodes+1 ! handle edge case so we don't access virtual root

   ! Figure out contribution blocks that are input to each part
   ! FIXME: consolidate all these deallocation by just calling free() at start of anal???
   deallocate(contrib_ptr, stat=st)
   deallocate(contrib_idx, stat=st)
   allocate(contrib_ptr(nparts+3), contrib_idx(nparts), contrib_dest(nparts))
   ! Count contributions at offset +2
   contrib_ptr(3:nparts+3) = 0
   do i = 1, nparts-1 ! by defn, last part has no parent
      j = sparent(part(i+1)-1) ! node index of parent
      if(j.gt.nnodes) cycle ! part is a root
      k = i+1 ! part index of j
      do while(j.ge.part(k+1))
         k = k + 1
      end do
      contrib_ptr(k+2) = contrib_ptr(k+2) + 1
   end do
   ! Figure out contrib_ptr starts at offset +1
   contrib_ptr(1:2) = 1
   do i = 1, nparts
      contrib_ptr(i+2) = contrib_ptr(i+1) + contrib_ptr(i+2)
   end do
   ! Drop sources into list
   do i = 1, nparts-1 ! by defn, last part has no parent
      j = sparent(part(i+1)-1) ! node index of parent
      if(j.gt.nnodes) then
         ! part is a root
         contrib_idx(i) = nparts+1
         cycle
      endif
      k = i+1 ! part index of j
      do while(j.ge.part(k+1))
         k = k + 1
      end do
      contrib_idx(i) = contrib_ptr(k+1)
      contrib_dest(contrib_idx(i)) = j
      contrib_ptr(k+1) = contrib_ptr(k+1) + 1
   end do
   contrib_idx(nparts) = nparts+1 ! last part must be a root

   ! Allocate subtrees to execution locations to try and get close to target
   ! cpu_gpu_ratio = flops(cpu) / flops(cpu+gpu)
   allocate(exec_loc(nparts)) !FIXME stat
   exec_loc(:) = -1 ! default to nowhere
   ! All subtrees with contribution must be on CPU
   do i = 1, nparts
      if(contrib_ptr(i).eq.contrib_ptr(i+1)) cycle ! no contrib to subtree
      exec_loc(i) = EXEC_LOC_CPU
      cpu_flops = cpu_flops + flops(part(i+1)-1)
   end do
   ! Handle any unallocated subtrees with greedy algorithm moving towards
   ! desired ratio
   do i = 1, nparts
      if(exec_loc(i).ne.-1) cycle ! already allocated
      if(cpu_flops / total_flops .lt. cpu_gpu_ratio) then
         exec_loc(i) = EXEC_LOC_CPU
         cpu_flops = cpu_flops + flops(part(i+1)-1)
      else
         exec_loc(i) = EXEC_LOC_GPU
         gpu_flops = gpu_flops + flops(part(i+1)-1)
      endif
   end do
end subroutine find_subtree_partition

!****************************************************************************

!
! This routine requires the LOWER and UPPER triangular parts of A
! to be held in CSC format using ptr2 and row2
! AND lower triangular part held using ptr and row.
!
! On exit from this routine, order is set to order
! input to factorization.
!
subroutine analyse_phase(n, ptr, row, ptr2, row2, order, invp, &
      akeep, options, inform, user_topology)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! col pointers (lower triangle) 
   integer, intent(in) :: row(ptr(n+1)-1) ! row indices (lower triangle)
   integer, intent(in) :: ptr2(n+1) ! col pointers (whole matrix)
   integer, intent(in) :: row2(ptr2(n+1)-1) ! row indices (whole matrix)
   integer, dimension(n), intent(inout) :: order
      !  On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(out) :: invp 
      ! Work array. Used to hold inverse of order but
      ! is NOT set to inverse for the final order that is returned.
   type(ssids_akeep), intent(inout) :: akeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(inout) :: inform
   type(numa_region), dimension(:), optional, intent(in) :: user_topology

   character(50)  :: context ! Procedure name (used when printing).
   integer, dimension(:), allocatable :: contrib_dest, exec_loc, level
   type(numa_region), dimension(:), allocatable :: topology

   real :: cpu_gpu_ratio
   integer :: nemin, flag
   integer :: blkm, blkn
   integer :: i, j
   integer :: nout, nout1 ! streams for errors and warnings
   integer :: nz ! ptr(n+1)-1
   integer :: st

   context = 'ssids_analyse'
   nout = options%unit_error
   if (options%print_level < 0) nout = -1
   nout1 = options%unit_warning
   if (options%print_level < 0) nout1 = -1
   st = 0

   ! Check nemin and set to default if out of range.
   nemin = options%nemin
   if(nemin < 1) nemin = nemin_default

   ! Perform basic analysis so we can figure out subtrees we want to construct
   call basic_analyse(n, ptr2, row2, order, akeep%nnodes, akeep%sptr, &
      akeep%sparent, akeep%rptr,akeep%rlist,                        &
      nemin, flag, inform%stat, inform%num_factor, inform%num_flops)
   select case(flag)
   case(0)
      ! Do nothing
   case(-1)
      ! Allocation error
      inform%flag = SSIDS_ERROR_ALLOCATION
      call ssids_print_flag(inform,nout,context)
      return
   case(1)
      ! Zero row/column.
      inform%flag = SSIDS_WARNING_ANAL_SINGULAR
   case default
      ! Should never reach here
      inform%flag = SSIDS_ERROR_UNKNOWN
   end select

   ! set invp to hold inverse of order
   do i = 1,n
      invp(order(i)) = i
   end do
   ! any unused variables are at the end and so can set order for them
   do j = akeep%sptr(akeep%nnodes+1), n
      i = invp(j)
      order(i) = 0
   end do

   ! Build map from A to L in nptr, nlist
   nz = ptr(n+1) - 1
   allocate(akeep%nptr(n+1), akeep%nlist(2,nz), stat=st)
   if (st .ne. 0) go to 100
   call build_map(n, ptr, row, order, invp, akeep%nnodes, akeep%sptr, &
      akeep%rptr, akeep%rlist, akeep%nptr, akeep%nlist, st)
   if (st .ne. 0) go to 100

   ! Sort out subtrees
   cpu_gpu_ratio = options%cpu_gpu_ratio
   if(.not.detect_gpu()) cpu_gpu_ratio = 1.0 ! Entirely on CPU
   call find_subtree_partition(akeep%nnodes, akeep%sptr, akeep%sparent, &
      akeep%rptr, akeep%topology, options%min_npart, options%max_flops_part, &
      cpu_gpu_ratio, akeep%nparts, akeep%part, exec_loc, &
      akeep%contrib_ptr, akeep%contrib_idx, contrib_dest)
   !print *, "invp = ", akeep%invp
   !print *, "sptr = ", akeep%sptr(1:akeep%nnodes+1)
   !print *, "sparent = ", akeep%sparent
   !print *, "Partition suggests ", akeep%nparts, " parts"
   !print *, "akeep%part = ", akeep%part(1:akeep%nparts+1)
   !print *, "exec_loc   = ", exec_loc(1:akeep%nparts)
   !print *, "parents = ", akeep%sparent(akeep%part(2:akeep%nparts+1)-1)
   !print *, "contrib_ptr = ", akeep%contrib_ptr(1:akeep%nparts+1)
   !print *, "contrib_idx = ", akeep%contrib_idx(1:akeep%nparts)
   !print *, "contrib_dest = ", &
   !   contrib_dest(1:akeep%contrib_ptr(akeep%nparts+1)-1)

   ! Construct symbolic subtrees
   if(allocated(akeep%subtree)) then
      do i = 1, size(akeep%subtree)
         if(associated(akeep%subtree(i)%ptr)) &
            deallocate(akeep%subtree(i)%ptr)
      end do
      deallocate(akeep%subtree)
   endif
   allocate(akeep%subtree(akeep%nparts))
   do i = 1, akeep%nparts
      select case(exec_loc(i))
      case(EXEC_LOC_CPU)
         akeep%subtree(i)%ptr => construct_cpu_symbolic_subtree(akeep%n, &
            akeep%part(i), akeep%part(i+1), akeep%sptr, akeep%sparent, &
            akeep%rptr, akeep%rlist, akeep%nptr, akeep%nlist, &
            contrib_dest(akeep%contrib_ptr(i):akeep%contrib_ptr(i+1)-1), &
            options)
      case(EXEC_LOC_GPU)
         akeep%subtree(i)%ptr => construct_gpu_symbolic_subtree(akeep%n, &
            akeep%part(i), akeep%part(i+1), akeep%sptr, akeep%sparent, &
            akeep%rptr, akeep%rlist, akeep%nptr, akeep%nlist, options)
      end select
   end do

   ! Info
   allocate(level(akeep%nnodes+1), stat=st)
   if (st .ne. 0) go to 100
   level(akeep%nnodes+1) = 0
   inform%maxfront = 0
   inform%maxdepth = 0
   do i = akeep%nnodes, 1, -1
      blkn = akeep%sptr(i+1) - akeep%sptr(i) 
      blkm = int(akeep%rptr(i+1) - akeep%rptr(i))
      level(i) = level(akeep%sparent(i)) + 1
      inform%maxfront = max(inform%maxfront, blkn)
      inform%maxdepth = max(inform%maxdepth, level(i))
   end do
   deallocate(level, stat=st)
   inform%matrix_rank = akeep%sptr(akeep%nnodes+1)-1
   inform%num_sup = akeep%nnodes

   ! Store copy of inform data in akeep
   akeep%flag = inform%flag
   akeep%matrix_dup = inform%matrix_dup
   akeep%matrix_missing_diag = inform%matrix_missing_diag
   akeep%matrix_outrange = inform%matrix_outrange
   akeep%maxdepth = inform%maxdepth
   akeep%num_sup = inform%num_sup
   akeep%num_flops = inform%num_flops

   return

   100 continue
   inform%stat = st
   if (inform%stat .ne. 0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      call ssids_print_flag(inform,nout,context)
   end if
   return
 
end subroutine analyse_phase

!****************************************************************************
!
! Build a map from A to nodes
! lcol( nlist(2,i) ) = val( nlist(1,i) )
! nptr defines start of each node in nlist
!
subroutine build_map(n, ptr, row, perm, invp, nnodes, sptr, rptr, rlist, &
      nptr, nlist, st)
   ! Original matrix A
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   ! Permutation and its inverse (some entries of perm may be negative to
   ! act as flags for 2x2 pivots, so need to use abs(perm))
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(in) :: invp
   ! Supernode partition of L
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   ! Row indices of L
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(in) :: rlist
   ! Output mapping
   integer, dimension(nnodes+1), intent(out) :: nptr
   integer, dimension(2, ptr(n+1)-1), intent(out) :: nlist
   ! Error check paramter
   integer, intent(out) :: st

   integer :: i, j, k, p
   integer(long) :: jj
   integer :: blkm
   integer :: col
   integer :: node
   integer, dimension(:), allocatable :: ptr2, row2, origin
   integer, dimension(:), allocatable :: map

   allocate(map(n), ptr2(n+3), row2(ptr(n+1)-1), origin(ptr(n+1)-1), stat=st)
   if(st.ne.0) return

   !
   ! Build transpose of A in ptr2, row2. Store original posn of entries in
   ! origin array.
   !
   ! Count number of entries in row i in ptr2(i+2). Don't include diagonals.
   ptr2(:) = 0
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         if (k.eq.i) cycle
         ptr2(k+2) = ptr2(k+2) + 1
      end do
   end do
   ! Work out row starts such that row i starts in posn ptr2(i+1)
   ptr2(1:2) = 1
   do i = 1, n
      ptr2(i+2) = ptr2(i+2) + ptr2(i+1)
   end do
   ! Drop entries into place
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         if (k.eq.i) cycle
         row2(ptr2(k+1)) = i
         origin(ptr2(k+1)) = j
         ptr2(k+1) = ptr2(k+1) + 1
      end do
   end do

   !
   ! Build nptr, nlist map
   !
   p = 1
   do node = 1, nnodes
      blkm = int(rptr(node+1) - rptr(node))
      nptr(node) = p

      ! Build map for node indices
      do jj = rptr(node), rptr(node+1)-1
         map(rlist(jj)) = int(jj-rptr(node)+1)
      end do

      ! Build nlist from A-lower transposed
      do j = sptr(node), sptr(node+1)-1
         col = invp(j)
         do i = ptr2(col), ptr2(col+1)-1
            k = abs(perm(row2(i))) ! row of L
            if (k<j) cycle
            nlist(2,p) = (j-sptr(node))*blkm + map(k)
            nlist(1,p) = origin(i)
            p = p + 1
         end do
      end do

      ! Build nlist from A-lower
      do j = sptr(node), sptr(node+1)-1
         col = invp(j)
         do i = ptr(col), ptr(col+1)-1
            k = abs(perm(row(i))) ! row of L
            if (k<j) cycle
            nlist(2,p) = (j-sptr(node))*blkm + map(k)
            nlist(1,p) = i
            p = p + 1
         end do
      end do
   end do
   nptr(nnodes+1) = p
   
end subroutine build_map

end module spral_ssids_analyse
