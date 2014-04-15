program random_matrix
   use spral_matrix_util, only : SPRAL_MATRIX_UNSPECIFIED,       &
                                 SPRAL_MATRIX_REAL_RECT,       &
                                 SPRAL_MATRIX_REAL_SYM_PSDEF,  &
                                 SPRAL_MATRIX_REAL_SYM_INDEF,  &
                                 SPRAL_MATRIX_REAL_SKEW,       &
                                 SPRAL_MATRIX_CPLX_RECT
   use spral_random, only : random_state, random_integer, random_logical
   use spral_random_matrix, only : random_matrix_generate
   implicit none

   integer, parameter :: wp = kind(0d0)

   integer, parameter :: ERROR_ALLOCATION = -1, & ! Allocation failed
                         ERROR_MATRIX_TYPE= -2, & ! Bad matrix type
                         ERROR_ARG        = -3, & ! m, n or nnz < 1
                         ERROR_NONSQUARE  = -4, & ! m!=n contradicts matrix_type
                         ERROR_SINGULAR   = -5    ! request non-singular
                                                  ! but nnz<min(m,n)


   integer :: errors

   errors = 0

   call test_errors
   call test_random_symmetric
   call test_random_unsymmetric

   write(*,"(/a)") "================"
   if(errors.eq.0) then
      write(*,"(a)") "All tests passed"
   else
      write(*,"(a,i6,a)") "Encountered ", errors, " errors"
   endif
   write(*,"(a)") "================"

   if(errors.ne.0) stop 1 ! non-zero return code to OS on failure
contains

subroutine test_random_unsymmetric
   integer, parameter :: nprob = 100
   integer, parameter :: maxn = 10000
   integer, parameter :: maxnnz_factor = 10

   integer :: prblm
   integer :: matrix_type, m, n, nnz, flag
   integer, dimension(:), allocatable :: ptr, row
   real(wp), dimension(:), allocatable :: val
   type(random_state) :: state
   logical :: nonsingular, sort

   write(*,"(/a)") "==================================="
   write(*,"(a)")  "Testing random unsymmetric matrices"
   write(*,"(a)")  "==================================="

   allocate(ptr(maxn+1), row(maxnnz_factor*maxn), val(maxnnz_factor*maxn))

   matrix_type = SPRAL_MATRIX_UNSPECIFIED
   do prblm = 1, nprob
      if(prblm<10) then
         m = prblm
         n = prblm+1
      else
         m = random_integer(state, maxn)
         n = random_integer(state, maxn)
      endif
      nnz = random_integer(state, min((n+1)/2,maxnnz_factor)*n)
      nonsingular = random_logical(state)
      sort = random_logical(state)

      write(*, "(a,i5,a,i5,a,i5,a,i6,a,l1,l1,a)", advance="no") &
         " * no. ", prblm, " m = ", m, " n = ", n, " nnz = ", nnz, " flags = ",&
         nonsingular, sort, " ..."

      call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, &
         flag, val=val, nonsingular=nonsingular, sort=sort)
      if(nonsingular .and. nnz<min(m,n)) then
         call print_result(flag, ERROR_SINGULAR)
      else
         if(flag.ne.0) then
            write(*, "(a/a,i5)") "fail", "flag = ", flag
            errors = errors + 1
            cycle
         endif
         call chk_random_unsymmetric(m, n, nnz, ptr, row, val, nonsingular, &
            sort)
      endif

   end do

end subroutine test_random_unsymmetric

subroutine chk_random_unsymmetric(m, n, nnz, ptr, row, val, nonsingular, sort)
   integer, intent(in) :: m
   integer, intent(in) :: n
   integer, intent(in) :: nnz
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   logical, intent(in) :: nonsingular
   logical, intent(in) :: sort

   integer :: i, j
   integer, dimension(:), allocatable :: rcnt

   if(ptr(n+1)-1.ne.nnz) then
      write(*, "(a/a,2i6)") "fail", "ptr(n+1)-1 != nnz requested", &
         ptr(n+1)-1, nnz
      errors = errors + 1
      return
   endif

   do i = 1, n
      if(ptr(i+1).lt.ptr(i)) then
         write(*, "(a/a)") "fail", "ptr non-monotonic"
         errors = errors + 1
         return
      endif
      do j = ptr(i), ptr(i+1)-1
         if(row(j).lt.1 .or. row(j).gt.m) then
            write(*, "(a/a,i5,a,i5)") "fail", "col ", i, &
               " has out-of-range row index ", row(j)
            errors = errors + 1
            return
         endif
         if(sort .and. j.gt.ptr(i)) then
            if(row(j).le.row(j-1)) then
               write(*, "(a/a,i5,a,i5)") "fail", "col ", i, &
                  "has non-monotonic row numbers at j=", j
               errors = errors + 1
               return
            endif
         endif
         if(abs(val(j)).gt.1.0_wp) then
            write(*, "(a/a,i5,a,es12.4)") "fail", "val( ", j, &
               ") has out-of-range value ", val(j)
            errors = errors + 1
            return
         endif
      end do
   end do

   if(nonsingular) then
      ! Check at least min(m,n) rows and columns non-empty
      ! Note: to do this properly we'd need to find a maximal transveral,
      ! but probably not worth the effort in a test deck given we're 99%
      ! sure the code is correct
      allocate(rcnt(m))
      rcnt(:) = 0
      do i = 1, ptr(n+1)-1
         rcnt(row(i)) = rcnt(row(i)) + 1
      end do
      if(count(rcnt(:).ne.0) .lt. min(m,n)) then
         write(*, "(a/a,i5,a,i5)") "fail", "nonsingular requested but only ", &
            count(rcnt(:).ne.0), " rows have non-zeroes. min(m,n) = ", min(m,n)
         errors = errors + 1
         return
      endif
      j = 0
      do i = 1, n
         if(ptr(i+1)-ptr(i).gt.0) j = j + 1
      end do
      if(j .lt. min(m,n)) then
         write(*, "(a/a,i5,a,i5)") "fail", "nonsingular requested but only ", &
            j, " cols have non-zeroes. min(m,n) = ", min(m,n)
         errors = errors + 1
         return
      endif
   endif

   write(*, "(a)") "ok"
end subroutine chk_random_unsymmetric

subroutine test_random_symmetric
   integer, parameter :: nprob = 100
   integer, parameter :: maxn = 10000
   integer, parameter :: maxnnz_factor = 10

   integer :: prblm
   integer :: matrix_type, n, nnz, flag
   integer, dimension(:), allocatable :: ptr, row
   real(wp), dimension(:), allocatable :: val
   type(random_state) :: state
   logical :: nonsingular, sort

   write(*,"(/a)") "================================="
   write(*,"(a)")  "Testing random symmetric matrices"
   write(*,"(a)")  "================================="

   allocate(ptr(maxn+1), row(maxnnz_factor*maxn), val(maxnnz_factor*maxn))

   matrix_type = SPRAL_MATRIX_REAL_SYM_INDEF
   do prblm = 1, nprob
      if(prblm<10) then
         n = prblm
      else
         n = random_integer(state, maxn)
      endif
      nnz = random_integer(state, min((n+1)/2,maxnnz_factor)*n)
      nonsingular = random_logical(state)
      sort = random_logical(state)

      write(*, "(a,i5,a,i5,a,i6,a,l1,l1,a)", advance="no") &
         " * no. ", prblm, " n = ", n, " nnz = ", nnz, " flags = ", &
         nonsingular, sort, " ..."

      call random_matrix_generate(state, matrix_type, n, n, nnz, ptr, row, &
         flag, val=val, nonsingular=nonsingular, sort=sort)
      if(nonsingular .and. nnz<n) then
         call print_result(flag, ERROR_SINGULAR)
      else
         if(flag.ne.0) then
            write(*, "(a/a,i5)") "fail", "flag = ", flag
            errors = errors + 1
            cycle
         endif
         call chk_random_symmetric(n, nnz, ptr, row, val, nonsingular, sort)
      endif

   end do

end subroutine test_random_symmetric

subroutine chk_random_symmetric(n, nnz, ptr, row, val, nonsingular, sort)
   integer, intent(in) :: n
   integer, intent(in) :: nnz
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   real(wp), dimension(ptr(n+1)-1), intent(in) :: val
   logical, intent(in) :: nonsingular
   logical, intent(in) :: sort

   integer :: i, j
   logical :: dpresent

   if(ptr(n+1)-1.ne.nnz) then
      write(*, "(a/a)") "fail", "ptr(n+1)-1 != nnz requested"
      errors = errors + 1
      return
   endif

   do i = 1, n
      if(ptr(i+1).lt.ptr(i)) then
         write(*, "(a/a)") "fail", "ptr non-monotonic"
         errors = errors + 1
         return
      endif
      dpresent = .false.
      do j = ptr(i), ptr(i+1)-1
         if(row(j).eq.i) dpresent = .true.
         if(row(j).lt.i .or. row(j).gt.n) then
            write(*, "(a/a,i5,a,i5)") "fail", "col ", i, &
               "has out-of-range row index ", row(j)
            errors = errors + 1
            return
         endif
         if(sort .and. j.gt.ptr(i)) then
            if(row(j).le.row(j-1)) then
               write(*, "(a/a,i5,a,i5)") "fail", "col ", i, &
                  "has non-monotonic row numbers at j=", j
               errors = errors + 1
               return
            endif
         endif
         if(abs(val(j)).gt.1.0_wp) then
            write(*, "(a/a,i5,a,es12.4)") "fail", "val( ", j, &
               ") has out-of-range value ", val(j)
            errors = errors + 1
            return
         endif
      end do
      if(nonsingular .and. .not.dpresent) then
         ! In symmetric case nonsingular means diagonal is present
         write(*, "(a/a,i5)") "fail", "nonsingular requested but diagonal not &
            &present in column ", i
         errors = errors + 1
         return
      endif
   end do

   write(*, "(a)") "ok"
end subroutine chk_random_symmetric

subroutine test_errors
   integer :: matrix_type, m, n, nnz, flag
   integer, dimension(:), allocatable :: ptr, row
   type(random_state) :: state

   write(*,"(/a)") "======================"
   write(*,"(a)") "Testing errors:"
   write(*,"(a)") "======================"

   ! Set default values
   matrix_type = SPRAL_MATRIX_UNSPECIFIED
   m = 100
   n = 100
   nnz = 1000
   allocate(ptr(101), row(1000))

   !
   ! Test bad matrix types
   !
   write(*,"(a)",advance="no") " * Testing matrix_type = 100................."
   matrix_type = 100
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_MATRIX_TYPE)
   matrix_type = SPRAL_MATRIX_REAL_RECT ! restore

   write(*,"(a)",advance="no") " * Testing matrix_type = COMPLEX............."
   matrix_type = SPRAL_MATRIX_CPLX_RECT
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_MATRIX_TYPE)
   matrix_type = SPRAL_MATRIX_UNSPECIFIED ! restore

   !
   ! Test bad args (i.e. m,n,nnz < 1)
   !
   write(*,"(a)",advance="no") " * Testing m < 1............................."
   m = 0
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_ARG)
   m = 100 ! restore

   write(*,"(a)",advance="no") " * Testing n < 1............................."
   n = 0
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_ARG)
   n = 100 ! restore

   write(*,"(a)",advance="no") " * Testing nnz < 1..........................."
   nnz = 0
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_ARG)
   nnz = 1000 ! restore

   write(*,"(a)",advance="no") " * Testing nnz < m*n (unsym)................."
   nnz = m*n+1
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_ARG)
   nnz = 1000 ! restore

   write(*,"(a)",advance="no") " * Testing nnz < n*(n+1)/2 (sym)............."
   matrix_type = SPRAL_MATRIX_REAL_SYM_INDEF
   nnz = n*(n+1)/2 + 1
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_ARG)
   matrix_type = SPRAL_MATRIX_UNSPECIFIED; nnz = 1000 ! restore

   !
   ! Test non-square
   !
   matrix_type = SPRAL_MATRIX_REAL_SYM_PSDEF
   m = n+1
   write(*,"(a)",advance="no") " * Testing non-square + SYM_PSDEF............"
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_NONSQUARE)
   matrix_type = SPRAL_MATRIX_UNSPECIFIED; m = 100 ! restore

   matrix_type = SPRAL_MATRIX_REAL_SYM_INDEF
   m = n+1
   write(*,"(a)",advance="no") " * Testing non-square + SYM_INDEF............"
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_NONSQUARE)
   matrix_type = SPRAL_MATRIX_UNSPECIFIED; m = 100 ! restore

   matrix_type = SPRAL_MATRIX_REAL_SKEW
   m = n+1
   write(*,"(a)",advance="no") " * Testing non-square + SKEW................."
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag)
   call print_result(flag, ERROR_NONSQUARE)
   matrix_type = SPRAL_MATRIX_UNSPECIFIED; m = 100 ! restore

   !
   ! Test singular but insufficient nnz
   !
   nnz = n-1
   write(*,"(a)",advance="no") " * Testing non-singular but nnz too small...."
   call random_matrix_generate(state, matrix_type, m, n, nnz, ptr, row, flag, &
      nonsingular=.true.)
   call print_result(flag, ERROR_SINGULAR)
   nnz = 1000
end subroutine test_errors

subroutine print_result(actual, expected)
   integer :: actual
   integer :: expected

   if(actual.eq.expected) then
      write(*,"(a)") "ok"
   else
      write(*,"(a)") "fail"
      write(*,"(2(a,i4))") "returned ", actual, ", expected ", expected
      errors = errors + 1
   endif
end subroutine print_result

end program random_matrix
