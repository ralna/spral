! COPYRIGHT (c) 2007-2013 Science & Technology Facilities Council
! Authors: Sue Thorne and Jonathan Hogg
! Origin: Heavily modified version of hsl_mc68
! 
module spral_metis_wrapper
  implicit none

  private
  public :: metis_order ! Calls metis on a symmetric matrix

  integer, parameter :: ERROR_ALLOC = -1
  integer, parameter :: ERROR_N_OOR = -2

contains

!
! Fortran wrapper around metis
!
  subroutine metis_order(n,ptr,row,perm,invp,flag,stat)
    implicit none
    integer, intent(in) :: n ! Must hold the number of rows in A
    integer, intent(in) :: ptr(n+1) ! ptr(j) holds position in row of start of
      ! row indices for column j. ptr(n)+1 must equal the number of entries
      ! stored + 1. Only the lower triangular entries are stored with no
      ! duplicates or out-of-range entries
    integer, intent(in) :: row(:) ! size at least ptr(n+1)-1
    integer, intent(out) :: perm(n) ! Holds elimination order on output
    integer, intent(out) :: invp(n) ! Holds inverse of elimination order on exit
    integer, intent(out) :: flag ! Return value
    integer, intent(out) :: stat ! Stat value on allocation failure

    ! ---------------------------------------------
    ! Local variables
    ! ---------------------------------------------
    integer, allocatable :: ptr2(:) ! copy of pointers which is later modified
    integer, allocatable :: row2(:) ! copy of row indices
    integer :: metis_opts(8) ! metis options array
    integer :: iwlen ! length of iw

    ! Initialise flag and stat
    flag = 0
    stat = 0

    !
    ! Check that restrictions are adhered to
    !
    if (n .lt. 1) then
       flag = ERROR_N_OOR
       return
    end if

    if (n .eq. 1) then
       ! Matrix is of order 1 so no need to call ordering subroutines
       perm(1) = 1
       return
    end if

    ! Set length of iw
    iwlen = 2*ptr(n+1) - 2

    ! Allocate arrays
    allocate (ptr2(n+1),row2(iwlen),stat=stat)
    if (stat .ne. 0) then
       flag = ERROR_ALLOC
       return
    end if

    ! Expand matrix, dropping diagonal entries
    call half_to_full_drop_diag(n, ptr, row, ptr2, row2)

    ! Carry out ordering
    metis_opts(1) = 0 ! MeTiS defaults
    call metis_nodend(n,ptr2,row2,1,metis_opts,invp,perm)
  end subroutine metis_order

  ! Convert a matrix in half storage to one in full storage.
  ! Drops any diagonal entries.
  subroutine half_to_full_drop_diag(n, ptr, row, ptr2, row2)
    implicit none
    integer, intent(in) :: n
    integer, dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(in) :: row
    integer, dimension(*), intent(out) :: ptr2
    integer, dimension(*), intent(out) :: row2

    integer :: i, j, k

    ! Set ptr2(j) to hold no. nonzeros in column j
    ptr2(1:n+1) = 0
    do j = 1, n
       do k = ptr(j), ptr(j+1) - 1
          i = row(k)
          if (j .ne. i) then
             ptr2(i) = ptr2(i) + 1
             ptr2(j) = ptr2(j) + 1
          end if
       end do
    end do

    ! Set ptr2(j) to point to where row indices will end in row2
    do j = 2, n
       ptr2(j) = ptr2(j-1) + ptr2(j)
    end do
    ptr2(n+1) = ptr2(n) + 1

    ! Fill ptr2 and row2
    do j = 1, n
       do k = ptr(j), ptr(j+1) - 1
          i = row(k)
          if (j .ne. i) then
             row2(ptr2(i)) = j
             row2(ptr2(j)) = i
             ptr2(i) = ptr2(i) - 1
             ptr2(j) = ptr2(j) - 1
          end if
       end do
    end do
    do j = 1, n
       ptr2(j) = ptr2(j) + 1
    end do

  end subroutine half_to_full_drop_diag

end module spral_metis_wrapper
