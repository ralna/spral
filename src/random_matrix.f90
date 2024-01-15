! Copyright (C) 2014 Science and Technology Facilities Council (STFC)
!
! This is a reimplementation of similar functionality to YM11 based on the
! description in the "Method" section of the YM11 documentation.
!
! FIXME: I don't think the positive definite case is implemented as per doc yet!
module spral_random_matrix
  use spral_random, only : random_state, random_integer, random_real
  use spral_matrix_util, only : SPRAL_MATRIX_UNSPECIFIED,          &
       SPRAL_MATRIX_REAL_RECT, SPRAL_MATRIX_REAL_UNSYM,            &
       SPRAL_MATRIX_REAL_SYM_PSDEF, SPRAL_MATRIX_REAL_SYM_INDEF,   &
       SPRAL_MATRIX_REAL_SKEW
  implicit none

  private
  public :: random_matrix_generate

  integer, parameter :: wp = kind(0d0)
  integer, parameter :: long = selected_int_kind(18)

  integer, parameter :: ERROR_ALLOCATION = -1, & ! Allocation failed
                        ERROR_MATRIX_TYPE= -2, & ! Bad matrix type
                        ERROR_ARG        = -3, & ! m, n or nnz < 1
                        ERROR_NONSQUARE  = -4, & ! m!=n contradicts matrix_type
                        ERROR_SINGULAR   = -5    ! request non-singular
                                                 ! but nnz<min(m,n)

  interface random_matrix_generate
     module procedure random_matrix_generate32, random_matrix_generate64
  end interface random_matrix_generate
contains

!
! Generate a random m x n matrix with nnz non-zeroes.
! User can additionally specify a symmetric matrix (requires m==n), forced
! non-singularity, and the sorting of entries within columns.
!
  subroutine random_matrix_generate32(state, matrix_type, m, n, nnz, ptr, row, &
       flag, stat, val, nonsingular, sort)
    implicit none
    type(random_state), intent(inout) :: state ! random generator to use
    integer, intent(in) :: matrix_type ! ignored except for symmetric/unsymmetric
      ! (in future will be used for complex sym vs hermitian at least)
    integer, intent(in) :: m ! number of rows
    integer, intent(in) :: n ! number of columns
    integer, intent(in) :: nnz ! number of entries
    integer, dimension(n+1), intent(out) :: ptr ! column pointers
    integer, dimension(nnz), intent(out) :: row ! row indices
    integer, intent(out) :: flag ! return code
    integer, optional, intent(out) :: stat ! allocate error code
    real(wp), dimension(nnz), optional, intent(out) :: val ! numerical values
    logical, optional, intent(in) :: nonsingular ! force matrix to be explicitly
      ! non-singular. If not present, treated as .false.
    logical, optional, intent(in) :: sort ! sort entries in columns by row index.
      ! If not present, treated as .false.

    integer(long), dimension(:), allocatable :: ptr64
    integer :: st

   ! Create temporary 64-bit version of ptr
    allocate(ptr64(n+1), stat=st)
    if (st .ne. 0) then
       flag = ERROR_ALLOCATION
       if (present(stat)) stat = st
       return
    end if

    ! Call 64-bit version
    call random_matrix_generate64(state, matrix_type, m, n, int(nnz,long), &
      ptr64, row, flag, stat=stat, val=val, nonsingular=nonsingular, sort=sort)

    ! ... and copy back to 32-bit ptr
    ptr(:) = int(ptr64(:))
  end subroutine random_matrix_generate32

!
! Generate a random m x n matrix with nnz non-zeroes.
! User can additionally specify a symmetric matrix (requires m==n), forced
! non-singularity, and the sorting of entries within columns.
!
! FIXME: This routine will be slow if we're asked for a (near) dense matrix
! In this case, we might be better served by finding holes or using the first
! part of random permutations
  subroutine random_matrix_generate64(state, matrix_type, m, n, nnz, ptr, row, &
       flag, stat, val, nonsingular, sort)
    implicit none
    type(random_state), intent(inout) :: state ! random generator to use
    integer, intent(in) :: matrix_type ! ignored except for symmetric/unsymmetric
      ! (in future will be used for complex sym vs hermitian at least)
    integer, intent(in) :: m ! number of rows
    integer, intent(in) :: n ! number of columns
    integer(long), intent(in) :: nnz ! number of entries
    integer(long), dimension(n+1), intent(out) :: ptr ! column pointers
    integer, dimension(nnz), intent(out) :: row ! row indices
    integer, intent(out) :: flag ! return code
    integer, optional, intent(out) :: stat ! allocate error code
    real(wp), dimension(nnz), optional, intent(out) :: val ! numerical values
    logical, optional, intent(in) :: nonsingular ! force matrix to be explicitly
      ! non-singular. If not present, treated as .false.
    logical, optional, intent(in) :: sort ! sort entries in columns by row index.
      ! If not present, treated as .false.

    integer :: i, j, k, minidx
    integer(long) :: ii, jj
    integer, dimension(:), allocatable :: cnt, rperm, cperm
    logical, dimension(:), allocatable :: rused
    logical :: lsymmetric, lnonsingular, lsort
    integer :: st

    ! Initialize return codes
    flag = 0
    if (present(stat)) stat = 0

    ! Generate local logical flags
    lnonsingular = .false.
    if (present(nonsingular)) lnonsingular = nonsingular
    lsort = .false.
    if (present(sort)) lsort = sort

    ! Handle matrix type
    select case (matrix_type)
    case(SPRAL_MATRIX_UNSPECIFIED, SPRAL_MATRIX_REAL_RECT)
       lsymmetric = .false.
    case(SPRAL_MATRIX_REAL_UNSYM)
       lsymmetric = .false.
       if (m .ne. n) then
          ! Matrix is not square - did user mean SPRAL_MATRIX_REAL_RECT?
          flag = ERROR_NONSQUARE
          return
       end if
    case(SPRAL_MATRIX_REAL_SYM_PSDEF, SPRAL_MATRIX_REAL_SYM_INDEF, &
         SPRAL_MATRIX_REAL_SKEW)
       lsymmetric = .true.
       if (m .ne. n) then
          ! Matrix is not square - did user mean SPRAL_MATRIX_REAL_RECT?
          flag = ERROR_NONSQUARE
          return
       end if
    case default
       ! COMPLEX or unknown matrix type
       flag = ERROR_MATRIX_TYPE
       return
    end select

    ! Check args
    if ((m .lt. 1) .or. (n .lt. 1) .or. (nnz .lt. 1)) then
       ! Args out of range
       flag = ERROR_ARG
       return
    end if
    if ((lsymmetric .and. (n*(n+1_long)/2 .lt. nnz)) .or. &
         ((.not. lsymmetric) .and. (m*(n+0_long) .lt. nnz))) then
       ! Too many non-zeroes for matrix
       flag = ERROR_ARG
       return
    end if
    if (lnonsingular .and. (nnz .lt. min(m,n))) then
       ! Requested a non-singular matrix, but not enough non-zeroes
       flag = ERROR_SINGULAR
       return
    end if

    ! Allocate non-zeroes to columns
    allocate(cnt(n), stat=st)
    if (st .ne. 0) goto 100
    cnt(:) = 0
    if (lsymmetric) then
       ! In symmetric case, structural non-singularity is guarunteed by adding
       ! the diagonal
       if (lnonsingular) then
          allocate(rperm(m), cperm(n), stat=st)
          if (st .ne. 0) goto 100
          ! To be consistent with unsymmetric case, we satisfy the following
          ! through identity permutations:
          ! If cperm(i)<=min(m,n) then that column has a structural non-zero
          ! in position rperm(cperm(i))
          do i = 1, n
             rperm(i) = i
             cperm(i) = i
          end do
          cnt(cperm(:)) = cnt(cperm(:)) + 1
       end if
       ! Generate column assignments of remaining entries
       ii = nnz; if(lnonsingular) ii = nnz - min(m,n) ! Allow for forced non-sing
       do ii = 1, ii
          j = random_sym_wt_integer(state, n)
          do while (cnt(j) .ge. (m-j+1))
             j = random_sym_wt_integer(state, n)
          end do
          cnt(j) = cnt(j) + 1
       end do
    else
       ! If we force (structural) non-singularity, generate locations and
       ! add to column counts
       if (lnonsingular) then
          allocate(rperm(m), cperm(n), stat=st)
          if (st .ne. 0) goto 100
          ! We generate random permutations of rows and columns
          ! We use the first min(m,n) of each permutation
          ! If cperm(i)<=min(m,n) then that column has a structural non-zero
          ! in position rperm(cperm(i))
          ! [i.e. rperm gives actual rows, cperm doesn't - its an inverse]
          call random_perm(state, m, rperm)
          call random_perm(state, n, cperm)
          where (cperm(:) .le. min(m,n)) cnt(:) = 1
       end if
       ! Generate column assignments of remaining entries
       ii = nnz; if(lnonsingular) ii = nnz - min(m,n) ! Allow for forced non-sing
       do ii = 1, ii
          j = random_integer(state, n)
          do while (cnt(j) .ge. m)
             j = random_integer(state, n)
          end do
          cnt(j) = cnt(j) + 1
       end do
       if (sum(cnt) .ne. nnz) stop
    end if

    ! Determine row values
    allocate(rused(m), stat=st)
    if (st .ne. 0) goto 100
    rused(:) = .false.
    ptr(1) = 1
    do i = 1, n
       ! Determine end of col
       ptr(i+1) = ptr(i) + cnt(i)
       jj = ptr(i)
       ! Add non-singular entry if required
       if (lnonsingular) then
          if (cperm(i) .le. min(m,n)) then
             k = rperm(cperm(i))
             row(jj) = k
             rused(k) = .true.
             jj = jj + 1
          end if
       end if
       ! Add normal entries
       minidx = 1
       if (lsymmetric) minidx = i
       do jj = jj, ptr(i+1)-1
          k = random_integer_in_range(state, minidx, m)
          do while (rused(k))
             k = random_integer_in_range(state, minidx, m)
          end do
          row(jj) = k
          rused(k) = .true.
       end do
       ! Reset rused(:)
       do jj = ptr(i), ptr(i+1)-1
          rused(row(jj)) = .false.
       end do
    end do

    ! Optionally, sort
    if (lsort) then
       call dbl_tr_sort(m, n, ptr, row, st)
       if (st .ne. 0) goto 100
    end if

    ! Determine values
    if (present(val)) then
       do jj = 1, ptr(n+1)-1
          val(jj) = random_real(state)
       end do
    end if

    return ! Normal return

100 continue
    ! Memory allocation failure
    flag = ERROR_ALLOCATION
    if (present(stat)) stat = st
    return
  end subroutine random_matrix_generate64

!
! Returns a random number in range [1,n] weighted by number of entries in
! lower half triangle
!
! Do this by only accepting randomly generated column with frequency
! proportional to number of entries in column
  integer function random_sym_wt_integer(state, n)
    implicit none
    type(random_state), intent(inout) :: state ! random state
    integer, intent(in) :: n ! dimension of matrix

    integer :: r1, r2

    r1 = random_integer(state, n)
    r2 = random_integer(state, n)
    do while (r2 .lt. r1)
       r1 = random_integer(state, n)
       r2 = random_integer(state, n)
    end do

    random_sym_wt_integer = r1
  end function random_sym_wt_integer

!
! Returns a random integer in range [minv,maxv] inclusive
!
  integer function random_integer_in_range(state, minv, maxv)
    implicit none
    type(random_state), intent(inout) :: state ! random state
    integer, intent(in) :: minv ! minimum value in range
    integer, intent(in) :: maxv ! maximum value in range

    random_integer_in_range = minv + random_integer(state, maxv-minv+1) - 1
  end function random_integer_in_range

!
! Returns a random permutation of length n in perm using Knuth shuffles
!
  subroutine random_perm(state, n, perm)
    implicit none
    type(random_state), intent(inout) :: state ! random generator to use
    integer, intent(in) :: n
    integer, dimension(n), intent(out) :: perm

    integer :: i, j, temp

    ! Initialize with identity
    do i = 1, n
       perm(i) = i
    end do

    ! Go through positions i=1:n-1
    do i = 1, n-1
       ! Swap perm(i) with perm(j), where j is random in [i:n]
       j = random_integer_in_range(state, i, n)
       temp = perm(i)
       perm(i) = perm(j)
       perm(j) = temp
    end do
  end subroutine random_perm

!
! Sort a matrix's columns to increase row order
!
! Uses a double transpose algorithm to do so
! Based on a modified version of subroutine from HSL_MC78
  subroutine dbl_tr_sort(m, n, ptr, row, st)
    implicit none
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer(long), dimension(n+1), intent(in) :: ptr
    integer, dimension(ptr(n+1)-1), intent(inout) :: row
    integer, intent(out) :: st

    integer :: i, j, node
    integer(long) :: ii, jj
    integer, dimension(:), allocatable :: col
    integer(long), dimension(:), allocatable :: ptr2, nptr

    allocate(ptr2(m+2), stat=st)
    if (st .ne. 0) return
    ptr2(:) = 0

    ! Count number of entries in each row. ptr2(i+2) = #entries in row i
    do node = 1, n
       do ii = ptr(node), ptr(node+1)-1
          j = row(ii)
          ptr2(j+2) = ptr2(j+2) + 1
       end do
    end do

    ! Determine row starts. ptr2(i+1) = start of row i
    ptr2(1:2) = 1
    do i = 1, m
       ptr2(i+2) = ptr2(i+1) + ptr2(i+2)
    end do

    jj = ptr2(m+2)-1 ! total number of entries
    allocate(col(jj), stat=st)
    if (st .ne. 0) return

    ! Now fill in col array
    do node = 1, n
       do ii = ptr(node), ptr(node+1)-1
          j = row(ii) ! row entry
          col( ptr2(j+1) ) = node
          ptr2(j+1) = ptr2(j+1) + 1
       end do
    end do

    ! Finally transpose back into nodes
    allocate(nptr(n), stat=st)
    if (st .ne. 0) return
    nptr(:) = ptr(1:n)
    do i = 1, m
       do jj = ptr2(i), ptr2(i+1)-1
          node = col(jj)
          row(nptr(node)) = i
          nptr(node) = nptr(node) + 1
       end do
    end do
  end subroutine dbl_tr_sort

end module spral_random_matrix
