! (c) STFC 2010-2013
! Author: Jonathan Hogg
!
! Originally based on HSL_MA97 v2.2.0
module spral_ssids_analyse
  use, intrinsic :: iso_c_binding
  use spral_ssids_datatypes
  use spral_core_analyse, only : basic_analyse
  use spral_cuda, only : cudaMemcpy_h2d, cudaMalloc
  use spral_pgm, only : writePPM
  implicit none

  private
  public :: analyse_phase,   & ! Calls core analyse and builds data strucutres
            check_order,     & ! Check order is a valid permutation
            expand_pattern,  & ! Specialised half->full matrix conversion
            expand_pattern64,& ! Specialised half->full matrix conversion
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
  subroutine expand_pattern64(n,nz,ptr,row,aptr,arow)
    implicit none
    integer, intent(in) :: n ! order of system
    integer(long), intent(in) :: nz
    integer(long), intent(in) :: ptr(n+1)
    integer, intent(in) :: row(nz)
    integer(long), intent(out) :: aptr(n+1)
    integer, intent(out) :: arow(2*nz)

    integer :: i,j
    integer(long) :: k

    ! Set aptr(j) to hold no. nonzeros in column j
    aptr(:) = 0
    do j = 1, n
       do k = ptr(j), ptr(j+1) - 1
          i = row(k)
          aptr(i) = aptr(i) + 1
          if (j .eq. i) cycle
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
          if (j .eq. i) cycle
          arow(aptr(j)) = i
          aptr(j) = aptr(j) - 1
       end do
    end do
    do j = 1,n
       aptr(j) = aptr(j) + 1
    end do
  end subroutine expand_pattern64

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
    implicit none
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
          if (j .eq. i) cycle
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
          if (j .eq. i) cycle
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
    implicit none
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
          if (j .eq. i) cycle
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
          if (j .eq. i) cycle
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
    implicit none
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
    if (options%print_level .lt. 0) nout = -1

    if (size(order) .lt. n) then
       ! Order is too short
       inform%flag = SSIDS_ERROR_ORDER
       akeep%flag = inform%flag
       call ssids_print_flag(context,nout,inform%flag)
       return
    end if

    ! initialise
    invp(:) = 0

    do i = 1, n
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
       call ssids_print_flag(context,nout,inform%flag)
       return
    end if
  end subroutine check_order

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
       akeep, options, inform)
    implicit none
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

    character(50)  :: context ! Procedure name (used when printing).
    integer, dimension(:), allocatable :: child_next, child_head ! linked
      ! list for children, used to build akeep%child_ptr, akeep%child_list

    integer :: nemin, flag
    integer :: blkm, blkn
    integer :: i, j, k
    integer :: nout, nout1 ! streams for errors and warnings
    integer :: nz ! ptr(n+1)-1
    integer :: cuda_error, st

    context = 'ssids_analyse'
    nout = options%unit_error
    if (options%print_level .lt. 0) nout = -1
    nout1 = options%unit_warning
    if (options%print_level .lt. 0) nout1 = -1
    cuda_error = 0
    st = 0

    ! Check nemin and set to default if out of range.
    nemin = options%nemin
    if (nemin .lt. 1) nemin = nemin_default

    call basic_analyse(n, ptr2, row2, order, akeep%nnodes, akeep%sptr, &
         akeep%sparent, akeep%rptr,akeep%rlist,                        &
         nemin, flag, inform%stat, akeep%nfactor, inform%num_flops)
    select case (flag)
    case(0)
       ! Do nothing
    case(-1)
       ! Allocation error
       inform%flag = SSIDS_ERROR_ALLOCATION
       call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
       return
    case(1)
       ! Zero row/column.
       inform%flag = SSIDS_WARNING_ANAL_SINGULAR
    case default
       ! Should never reach here
       inform%flag = SSIDS_ERROR_UNKNOWN
    end select
    inform%num_factor = akeep%nfactor

    ! set invp to hold inverse of order
    do i = 1, n
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

    ! Build direct map for children
    allocate(akeep%rlist_direct(akeep%rptr(akeep%nnodes+1)-1), stat=st)
    if (st .ne. 0) go to 100
    call build_rlist_direct(n, akeep%nnodes, akeep%sparent, akeep%rptr, &
         akeep%rlist, akeep%rlist_direct, st)
    if (st .ne. 0) go to 100

    ! Find maxmn and setup levels
    allocate(akeep%level(akeep%nnodes+1), stat=st)
    if (st .ne. 0) go to 100

    akeep%maxmn = 0
    akeep%level(akeep%nnodes+1) = 0
    inform%maxfront = 0
    inform%maxdepth = 0
    do i = akeep%nnodes, 1, -1
       blkn = akeep%sptr(i+1) - akeep%sptr(i) 
       blkm = int(akeep%rptr(i+1) - akeep%rptr(i))
       akeep%maxmn = max(akeep%maxmn, blkm, blkn)
       akeep%level(i) = akeep%level(akeep%sparent(i)) + 1
       inform%maxfront = max(inform%maxfront, blkn)
       inform%maxdepth = max(inform%maxdepth, akeep%level(i))
    end do

    !call count_matrix_sizes(akeep%n, akeep%nnodes, akeep%level, &
    !   akeep%sptr, akeep%rptr)

    ! Setup child_ptr, child_next and calculate work per subtree
    allocate(child_next(akeep%nnodes+1), child_head(akeep%nnodes+1), &
         akeep%child_ptr(akeep%nnodes+2), akeep%child_list(akeep%nnodes), &
         akeep%subtree_work(akeep%nnodes+1), stat=st)
    if (st .ne. 0) go to 100
    child_head(:) = -1
    do i = akeep%nnodes, 1, -1 ! backwards so child list is in order
       blkn = akeep%sptr(i+1) - akeep%sptr(i) 
       blkm = int(akeep%rptr(i+1) - akeep%rptr(i))
       j = akeep%sparent(i)
       ! Add to parent's child linked list
       child_next(i) = child_head(j)
       child_head(j) = i
       ! Calculate extra work at this node
       akeep%subtree_work(i) = 0
       do k = blkm, blkm-blkn+1, -1
          akeep%subtree_work(i) = akeep%subtree_work(i) + k**2
       end do
    end do
    akeep%subtree_work(akeep%nnodes+1) = 0
    ! Add work up tree, build child_ptr and child_list
    akeep%child_ptr(1) = 1
    do i = 1, akeep%nnodes+1
       if (i .lt. akeep%nnodes+1) then
          j = akeep%sparent(i)
          akeep%subtree_work(j) = akeep%subtree_work(j) + akeep%subtree_work(i)
       end if
       j = child_head(i)
       akeep%child_ptr(i+1) = akeep%child_ptr(i)
       do while (j .ne. -1)
          akeep%child_list(akeep%child_ptr(i+1)) = j
          akeep%child_ptr(i+1) = akeep%child_ptr(i+1) + 1
          j = child_next(j)
       end do
    end do

    ! Sort out GPU-relevent data if needed
    call copy_analyse_data_to_device(2*(akeep%nptr(akeep%nnodes+1)-1), &
         akeep%nlist, akeep%rptr(akeep%nnodes+1)-1, akeep%rlist, &
         akeep%rlist_direct, akeep%gpu_nlist, akeep%gpu_rlist, &
         akeep%gpu_rlist_direct, cuda_error)
    if (cuda_error .ne. 0) goto 200

    ! Info
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
       call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
    end if

200 continue
    inform%cuda_error = cuda_error
    inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
    call ssids_print_flag(context,nout,inform%flag,cuda_error=inform%cuda_error)

  end subroutine analyse_phase

!****************************************************************************
!
! Build a direct mapping (assuming no delays) between a node's rlist and that
! of it's parent such that update can be done as:
! lcol_parent(rlist_direct(i)) = lcol(i)
!
  subroutine build_rlist_direct(n, nnodes, sparent, rptr, rlist, rlist_direct, st)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: nnodes
    integer, dimension(nnodes), intent(in) :: sparent
    integer(long), dimension(nnodes+1), intent(in) :: rptr
    integer, dimension(rptr(nnodes+1)-1), intent(in) :: rlist
    integer, dimension(rptr(nnodes+1)-1), intent(out) :: rlist_direct
    integer, intent(out) :: st

    integer :: node, parent
    integer(long) :: ii
    integer, dimension(:), allocatable :: map

    allocate(map(n), stat=st)
    if (st .ne. 0) return

    do node = 1, nnodes
       ! Build a map for parent
       parent = sparent(node)
       if (parent .gt. nnodes) cycle ! root of tree
       do ii = rptr(parent), rptr(parent+1)-1
          map(rlist(ii)) = int(ii-rptr(parent)+1)
       end do

       ! Build rlist_direct
       do ii = rptr(node), rptr(node+1)-1
          rlist_direct(ii) = map(rlist(ii))
       end do
    end do
  end subroutine build_rlist_direct

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
    if (st .ne. 0) return

    !
    ! Build transpose of A in ptr2, row2. Store original posn of entries in
    ! origin array.
    !
    ! Count number of entries in row i in ptr2(i+2). Don't include diagonals.
    ptr2(:) = 0
    do i = 1, n
       do j = ptr(i), ptr(i+1)-1
          k = row(j)
          if (k .eq. i) cycle
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
          if (k .eq. i) cycle
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
             if (k .lt. j) cycle
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
             if (k .lt. j) cycle
             nlist(2,p) = (j-sptr(node))*blkm + map(k)
             nlist(1,p) = i
             p = p + 1
          end do
       end do
    end do
    nptr(nnodes+1) = p

  end subroutine build_map

!****************************************************************************
!
! Copy data generated by analyse phase to GPU
! (Specifically, nlist, rlist and rlist_direct)
!
  subroutine copy_analyse_data_to_device(lnlist, nlist, lrlist, rlist, &
       rlist_direct, gpu_nlist, gpu_rlist, gpu_rlist_direct, cuda_error)
    implicit none
    integer, intent(in) :: lnlist
    integer(C_INT), dimension(lnlist), target, intent(in) :: nlist
    integer(long), intent(in) :: lrlist
    integer(C_INT), dimension(lrlist), target, intent(in) :: rlist
    integer(C_INT), dimension(lrlist), target, intent(in) :: rlist_direct
    type(C_PTR), intent(out) :: gpu_nlist
    type(C_PTR), intent(out) :: gpu_rlist
    type(C_PTR), intent(out) :: gpu_rlist_direct
    integer, intent(out) :: cuda_error ! Non-zero on error

    ! Copy nlist
    cuda_error = cudaMalloc(gpu_nlist, C_SIZEOF(nlist(1)) * lnlist)
    if (cuda_error .ne. 0) return
    cuda_error = cudaMemcpy_h2d(gpu_nlist, C_LOC(nlist), &
         C_SIZEOF(nlist(1)) * lnlist)
    if (cuda_error .ne. 0) return

    ! Copy rlist
    cuda_error = cudaMalloc(gpu_rlist, C_SIZEOF(rlist(1)) * lrlist)
    if (cuda_error .ne. 0) return
    cuda_error = cudaMemcpy_h2d(gpu_rlist, C_LOC(rlist), &
         C_SIZEOF(rlist(1)) * lrlist)
    if (cuda_error .ne. 0) return

    ! Copy rlist_direct
    cuda_error = cudaMalloc(gpu_rlist_direct, C_SIZEOF(rlist_direct(1)) * lrlist)
    if (cuda_error .ne. 0) return
    cuda_error = cudaMemcpy_h2d(gpu_rlist_direct, C_LOC(rlist_direct), &
         C_SIZEOF(rlist_direct(1)) * lrlist)
    if (cuda_error .ne. 0) return

  end subroutine copy_analyse_data_to_device

end module spral_ssids_analyse
