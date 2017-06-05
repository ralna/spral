! COPYRIGHT (c) 2007-2013 Science & Technology Facilities Council
! Authors: Sue Thorne and Jonathan Hogg
! Origin: Heavily modified version of hsl_mc68
! 
! This version is used to support METIS v5 (different API to v4)
!
module spral_metis_wrapper
  use, intrinsic :: iso_c_binding
  implicit none

  private
  public :: metis_order ! Calls metis on a symmetric matrix

  ! We use the C interface to METIS via the following interoperable interfaces
  interface
     subroutine METIS_SetDefaultOptions(options) &
          bind(C, name="METIS_SetDefaultOptions")
       use, intrinsic :: iso_c_binding
       implicit none
       integer(C_INT), dimension(*), intent(out) :: options
     end subroutine METIS_SetDefaultOptions
     integer(C_INT) function METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, &
          perm, iperm) bind(C, name="METIS_NodeND")
       use, intrinsic :: iso_c_binding
       implicit none
       integer(C_INT), intent(in) :: nvtxs
       integer(C_INT), dimension(*), intent(in) :: xadj, adjncy, options
       type(C_PTR), value :: vwgt
       integer(C_INT), dimension(*), intent(out) :: perm, iperm
     end function METIS_NodeND
  end interface
  ! Following array size based on #define in metis.h
  integer, parameter :: METIS_NOPTIONS = 40
  ! Following are based on enum in metis.h, adjusted to Fortran indexing.
  integer, parameter :: METIS_OPTION_CTYPE     =  3, &
                        METIS_OPTION_RTYPE     =  5, &
                        METIS_OPTION_DBGLVL    =  6, &
                        METIS_OPTION_NITER     =  7, &
                        METIS_OPTION_SEED      =  8, &
                        METIS_OPTION_NO2HOP    = 10, &
                        METIS_OPTION_COMPRESS  = 13, &
                        METIS_OPTION_CCORDER   = 14, &
                        METIS_OPTION_PFACTOR   = 15, &
                        METIS_OPTION_NSEPS     = 16, &
                        METIS_OPTION_UFACTOR   = 17, &
                        METIS_OPTION_NUMBERING = 18
  ! Following return codes based on enum in metis.h
  integer, parameter :: METIS_OK            =  1, &
                        METIS_ERROR_INPUT   = -2, &
                        METIS_ERROR_MEMORY  = -3, &
                        METIS_ERROR         = -4 

  ! Constants for this package
  integer, parameter :: ERROR_ALLOC = -1
  integer, parameter :: ERROR_N_OOR = -2
  integer, parameter :: ERROR_UNKNOWN = -999

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
    integer :: metis_opts(METIS_NOPTIONS) ! metis options array
    integer :: iwlen ! length of iw

    integer :: metis_flag

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
    call METIS_SetDefaultOptions(metis_opts)
    metis_opts(METIS_OPTION_NUMBERING) = 1 ! Fortran-style numbering
    metis_flag = METIS_NodeND(n, ptr2, row2, C_NULL_PTR, metis_opts, invp, perm)
    select case (metis_flag)
    case(METIS_OK)
       ! Everything OK, do nothing
    case(METIS_ERROR_MEMORY)
       ! Ran out of memory
       flag = ERROR_ALLOC
       stat = -99 ! As we don't know the real value, make up a distinctive one
       return
    case default
       ! Something else went wrong, report as unknown error
       ! (Should never encounter this code)
       print *, "Unknown metis error with code ", metis_flag
       flag = ERROR_UNKNOWN
    end select
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
