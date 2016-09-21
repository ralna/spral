program rutherford_boeing_test
   use spral_matrix_util
   use spral_rutherford_boeing
   use spral_random
   use spral_random_matrix, only: random_matrix_generate
   implicit none

   integer, parameter :: long = selected_int_kind(18)
   integer, parameter :: wp = kind(0d0)

   character(len=*), parameter :: filename = "rb_test_matrix.rb"

   integer :: errors

   ! Return codes
   integer, parameter :: SUCCESS           =  0    ! No errors
   integer, parameter :: ERROR_BAD_FILE    = -1    ! Failed to open file
   integer, parameter :: ERROR_NOT_RB      = -2    ! Header not valid for RB
   integer, parameter :: ERROR_IO          = -3    ! Error return from io
   integer, parameter :: ERROR_TYPE        = -4    ! Tried to read bad type
   integer, parameter :: ERROR_ELT_ASM     = -5    ! Read elt as asm or v/v
   integer, parameter :: ERROR_MATRIX_TYPE = -6 ! Bad value of matrix_type
   integer, parameter :: ERROR_EXTRA_SPACE = -10   ! control%extra_space<1.0
   integer, parameter :: ERROR_LWR_UPR_FULL= -11   ! control%lwr_upr_full oor
   integer, parameter :: ERROR_VALUES      = -13   ! control%values oor
   integer, parameter :: ERROR_ALLOC       = -20   ! failed on allocate
   integer, parameter :: WARN_AUX_FILE     = 1     ! values in auxiliary file

   errors = 0

   call test_errors
   call test_warnings
   call test_special
   call test_random

   write(*, "(/a)") "=========================="
   write(*, "(a,i4)") "Total number of errors = ", errors

   if(errors.ne.0) stop 1 ! ERROR CODE for make check script

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Test error codes
subroutine test_errors()
   integer :: m, n, inform, iunit
   integer, dimension(:), allocatable :: ptr, row
   real(wp), dimension(:), allocatable :: val
   type(rb_read_options) :: read_options, default_read_options
   type(rb_write_options) :: write_options, default_write_options

   write(*, "(a)")
   write(*, "(a)") "====================="
   write(*, "(a)") "Testing error returns"
   write(*, "(a)") "====================="

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,"(/,a)") "rb_peek():"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Failure to open a non-existent file
   write(*,"(a)",advance="no") " * Failure to open file......................"
   call rb_peek("non_existent.rb", inform, m=m, n=n)
   call test_eq(inform, ERROR_BAD_FILE)

   ! Failure on i/o operation (e.g. file is too short)
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   close(iunit)
   write(*,"(a)",advance="no") " * Failure on i/o operation.................."
   call rb_peek(filename, inform, m=m, n=n)
   call test_eq(inform, ERROR_IO)

   ! Not an RB file 1
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a80/a80/a3,11x,4(1x,i13))") &
      "This is", "not an RB file", "xxx", 1, 2, 3, 4
   close(iunit)
   write(*,"(a)",advance="no") " * Not a valid RB file 1....................."
   call rb_peek(filename, inform, m=m, n=n)
   call test_eq(inform, ERROR_NOT_RB)

   ! Not an RB file 2
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a80/a80/a3,11x,4(1x,i13))") &
      "This is", "not an RB file", "rxx", 1, 2, 3, 4
   close(iunit)
   write(*,"(a)",advance="no") " * Not a valid RB file 2....................."
   call rb_peek(filename, inform, m=m, n=n)
   call test_eq(inform, ERROR_NOT_RB)

   ! Not an RB file 3
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a80/a80/a3,11x,4(1x,i13))") &
      "This is", "not an RB file", "rsx", 1, 2, 3, 4
   close(iunit)
   write(*,"(a)",advance="no") " * Not a valid RB file 3....................."
   call rb_peek(filename, inform, m=m, n=n)
   call test_eq(inform, ERROR_NOT_RB)

   ! Not an RB file 4
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 5, 1, 1, 3
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "rsa", 5, 5, 8, 1
   close(iunit)
   write(*,"(a)",advance="no") " * Not a valid RB file 4....................."
   call rb_peek(filename, inform, m=m, n=n)
   call test_eq(inform, ERROR_NOT_RB)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,"(/,a)") "rb_read():"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Failure to open a non-existent file
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Failure to open file......................"
   call rb_read("non_existent.rb", m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_BAD_FILE)

   ! Failure on i/o operation (e.g. file is too short)
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 5, 1, 1, 3
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "rsa", 5, 5, 8, 0
   write(iunit, "(a16, a16, a20)") "(40i2)", "(40i2)", "(3e24.16)"
   write(iunit, "(40i2)") 1, 3, 6, 8, 8, 9
   close(iunit)
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Failure on i/o operation.................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_IO)

   ! Not an RB file 1
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a80/a80/a3,11x,4(1x,i13))") &
      "This is", "not an RB file", "xxx", 1, 2, 3, 4
   close(iunit)
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Not a valid RB file 1....................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_NOT_RB)

   ! Not an RB file 2
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a80/a80/a3,11x,4(1x,i13))") &
      "This is", "not an RB file", "rxx", 1, 2, 3, 4
   close(iunit)
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Not a valid RB file 2....................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_NOT_RB)

   ! Not an RB file 3
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a80/a80/a3,11x,4(1x,i13))") &
      "This is", "not an RB file", "rsx", 1, 2, 3, 4
   close(iunit)
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Not a valid RB file 3....................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_NOT_RB)

   ! Not an RB file 4
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 5, 1, 1, 3
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "rsa", 5, 5, 8, 1
   write(iunit, "(a16, a16, a20)") "(40i2)", "(40i2)", "(3e24.16)"
   write(iunit, "(40i2)") 1, 3, 6, 8, 8, 9
   close(iunit)
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Not a valid RB file 4....................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_NOT_RB)

   ! Unsupported type: complex
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 5, 1, 1, 3
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "csa", 5, 5, 8, 0
   write(iunit, "(a16, a16, a20)") "(40i2)", "(40i2)", "(3e24.16)"
   write(iunit, "(40i2)") 1, 3, 6, 8, 8, 9
   close(iunit)
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Unsupported type: complex................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_TYPE)

   ! Unsupported type: element
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 5, 1, 1, 3
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "rse", 5, 5, 8, 0
   write(iunit, "(a16, a16, a20)") "(40i2)", "(40i2)", "(3e24.16)"
   write(iunit, "(40i2)") 1, 3, 6, 8, 8, 9
   close(iunit)
   read_options = default_read_options
   write(*,"(a)",advance="no") " * Unsupported type: element................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_ELT_ASM)

   ! options%extra_space < 1.0
   call write_simple_matrix()
   read_options = default_read_options
   read_options%extra_space = 0.9
   write(*,"(a)",advance="no") " * option%extra_space < 1.0.................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_EXTRA_SPACE)

   ! options%lwr_upr_full invalid value
   call write_simple_matrix()
   read_options = default_read_options
   read_options%lwr_upr_full = -1
   write(*,"(a)",advance="no") " * option%lwr_upr_full = -1.................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_LWR_UPR_FULL)

   ! options%lwr_upr_full invalid value
   call write_simple_matrix()
   read_options = default_read_options
   read_options%values = 100
   write(*,"(a)",advance="no") " * option%values = 100......................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, ERROR_VALUES)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,"(/,a)") "rb_write():"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Failure to open a non-existent file
   write(*,"(a)",advance="no") " * Failure to open file......................"
   write_options = default_write_options
   call get_simple_matrix(m,n,ptr,row,val)
   call rb_write("/does/not/exist/matrix.rb", SPRAL_MATRIX_REAL_SYM_INDEF, m, &
      n, ptr, row, write_options, inform)
   call test_eq(inform, ERROR_BAD_FILE)

   ! Invalid matrix_type: complex matrix_type to real call
   write(*,"(a)",advance="no") " * Complex matrix_type on call to real......."
   write_options = default_write_options
   call get_simple_matrix(m,n,ptr,row,val)
   call rb_write(filename, -1, m, n, ptr, row, write_options, inform, val=val)
   call test_eq(inform, ERROR_MATRIX_TYPE)

   ! Invalid matrix_type: matrix_type = 5 (undefined)
   write(*,"(a)",advance="no") " * matrix_type = 5 (undefined)..............."
   write_options = default_write_options
   call get_simple_matrix(m,n,ptr,row,val)
   call rb_write(filename, 5, m, n, ptr, row, write_options, inform)
   call test_eq(inform, ERROR_MATRIX_TYPE)

   ! Invalid matrix_type: matrix_type = 7 (out-of-range)
   write(*,"(a)",advance="no") " * matrix_type = 7 (out-of-range)............"
   write_options = default_write_options
   call get_simple_matrix(m,n,ptr,row,val)
   call rb_write(filename, 7, m, n, ptr, row, write_options, inform)
   call test_eq(inform, ERROR_MATRIX_TYPE)
end subroutine test_errors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Test warnings
subroutine test_warnings
   integer :: m, n, inform, iunit
   integer, dimension(:), allocatable :: ptr, row
   real(wp), dimension(:), allocatable :: val
   type(rb_read_options) :: read_options

   write(*, "(a)")
   write(*, "(a)") "======================="
   write(*, "(a)") "Testing warning returns"
   write(*, "(a)") "======================="

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,"(/,a)") "rb_read():"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Values in auxiliary file
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 2, 1, 1, 0
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "qsa", 5, 5, 8, 0
   write(iunit, "(a16, a16, a20)") "(40i2)", "(40i2)", "(3e24.16)"
   write(iunit, "(40i2)") 1, 3, 6, 8, 8, 9
   write(iunit, "(40i2)") 1, 2, 2, 3, 5, 3, 4, 5
   close(iunit)
   write(*,"(a)",advance="no") " * Values in auxiliary file.................."
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, WARN_AUX_FILE)
end subroutine test_warnings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Test special cases that aren't covered elsewhere
subroutine test_special
   integer :: m, n, inform, iunit
   integer, dimension(:), allocatable :: ptr, row
   real(wp), dimension(:), allocatable :: val
   type(rb_read_options) :: read_options, default_read_options

   write(*, "(a)")
   write(*, "(a)") "======================="
   write(*, "(a)") "Testing special"
   write(*, "(a)") "======================="

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,"(/,a)") "rb_read():"
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Integer data in file
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 2, 1, 1, 0
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "isa", 5, 5, 8, 0
   write(iunit, "(a16, a16, a20)") "(40i2)", "(40i2)", "(40i2)"
   write(iunit, "(40i2)") 1, 3, 6, 8, 8, 9
   write(iunit, "(40i2)") 1, 2, 2, 3, 5, 3, 4, 5
   write(iunit, "(40i2)") 1, 2, 3, 4, 5, 6, 7, 8
   close(iunit)
   write(*,"(a)",advance="no") " * Integer input............................."
   read_options = default_read_options
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, SUCCESS)

   ! Integer data in file, but only read pattern
   open(newunit=iunit,file=filename,status='replace')
   write(iunit,"(a72,a8)") "Matrix", "ID"
   write(iunit,"(i14, 1x, i13, 1x, i13, 1x, i13)") 2, 1, 1, 0
   write(iunit, "(a3, 11x, i14, 1x, i13, 1x, i13, 1x, i13)") &
      "isa", 5, 5, 8, 0
   write(iunit, "(a16, a16, a20)") "(40i2)", "(40i2)", "(40i2)"
   write(iunit, "(40i2)") 1, 3, 6, 8, 8, 9
   write(iunit, "(40i2)") 1, 2, 2, 3, 5, 3, 4, 5
   write(iunit, "(40i2)") 1, 2, 3, 4, 5, 6, 7, 8
   close(iunit)
   write(*,"(a)",advance="no") " * Integer input (read pattern only)........."
   read_options = default_read_options
   read_options%values = 1 ! pattern only
   call rb_read(filename, m, n, ptr, row, val, read_options, inform)
   call test_eq(inform, SUCCESS)
end subroutine test_special

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Test (actual.eq.expected).
!> If test passes, print "pass"
!> If test fails, print "fail" and values, and increment errors
subroutine test_eq(actual, expected)
   integer, intent(in) :: actual
   integer, intent(in) :: expected

   if(actual.eq.expected) then
      write(*, "(a)") "pass"
   else
      write(*, "(a)") "fail"
      write(*, "(2(a,i4))") "actual =", actual, " expected =", expected
      errors = errors + 1
   endif
end subroutine test_eq

subroutine get_simple_matrix(m, n, ptr, row, val)
   integer, intent(out) :: m, n
   integer, allocatable, intent(inout) :: ptr(:), row(:)
   real(wp), allocatable, intent(inout) :: val(:)

   integer :: st

   ! Ensure arrays are correct size
   deallocate(ptr, stat=st)
   deallocate(row, stat=st)
   deallocate(val, stat=st)
   allocate(ptr(6), row(8), val(8))

   ! Set values
   ! ( 2  1         )
   ! ( 1  4  1    8 )
   ! (    1  3  2   )
   ! (       2      )
   ! (    8       2 )
   m = 5
   n = 5
   ptr(1:n+1)        = (/ 1,        3,             6,      8,8,   9 /)
   row(1:ptr(n+1)-1) = (/ 1,   2,   2,   3,   5,   3,   4,   5   /)
   val(1:ptr(n+1)-1) = (/ 2.0, 1.0, 4.0, 1.0, 8.0, 3.0, 2.0, 2.0 /)
end subroutine get_simple_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Write simple matrix to file
subroutine write_simple_matrix()
   integer :: m, n, inform
   integer, allocatable :: ptr(:), row(:)
   real(wp), allocatable :: val(:)
   type(rb_write_options) :: options

   call get_simple_matrix(m, n, ptr, row, val)
   call rb_write(filename, SPRAL_MATRIX_REAL_SYM_INDEF, m, n, ptr, row, &
      options, inform, val=val)
   if(inform.ne.0) then
      print *, "write_simple_matrix: rb_write error. Aborting.", inform
      stop -2
   endif
end subroutine write_simple_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Tests writing then peeking at and reading random matrices
subroutine test_random()
   ! Parameters
   integer, parameter :: nprob = 100
   integer, parameter :: maxn = 1000
   integer, parameter :: maxnz = 1000000

   ! RNG
   type(random_state) :: state

   ! Matrix data (to write)
   character(len=72) :: title
   character(len=8) :: id
   integer :: matrix_type
   integer :: m, n
   integer(long), dimension(:), allocatable :: ptr64
   integer, dimension(:), allocatable :: ptr32, row
   real(wp), dimension(:), allocatable :: val
   logical :: pattern

   ! Matrix data (to read)
   character(len=72) :: title_in
   character(len=8) :: id_in
   character(len=3) :: type_code_in
   integer :: matrix_type_in
   integer :: m_in, n_in
   integer(long) :: nelt_in, nvar_in, nval_in
   integer(long), dimension(:), allocatable :: ptr64_in
   integer, dimension(:), allocatable :: ptr32_in, row_in
   real(wp), dimension(:), allocatable :: val_in

   ! Options
   type(rb_read_options) :: read_options, default_read_options
   type(rb_write_options) :: write_options

   ! Working variables
   integer :: problem, flag, st
   integer(long) :: nnz

   write(*, "(a)")
   write(*, "(a)") "======================="
   write(*, "(a)") "Testing random matrices"
   write(*, "(a)") "======================="

   allocate(ptr32(maxn+1), ptr64(maxn+1), row(maxnz), val(maxnz))

   do problem = 1, nprob
      ! Generate parameters
      matrix_type = random_integer(state, 5)
      if(matrix_type.eq.5) matrix_type = 6 ! matrix type +5 not defined
      n = random_integer(state, maxn)
      m = n
      if(matrix_type.eq.SPRAL_MATRIX_REAL_RECT) &
         m = random_integer(state, maxn)
      nnz = m*n/2 - max(m,n)
      nnz = max(m, n) + random_integer(state, nnz)
      pattern = random_logical(state)

      ! Write status line
      write(*,"(a,i5,a,i12,a,2i5,i8,a,i2,a)", advance="no") &
         "Test ", problem, " state ", random_get_seed(state), &
         " dimn =", m, n, nnz, " type =", matrix_type, "... "

      ! Generate random_matrix
      call random_matrix_generate(state, matrix_type, m, n, nnz, ptr64, row, &
         flag, val=val)
      if(flag.ne.0) then
         print *, "Bad return from random_matrix_generate(). Aborting", flag
         stop 2
      endif

      ! Write random matrix
      write(title, "(a,i5)") "Test matrix ", problem
      write(id, "(a2,i6)") "ID", problem
      if(random_logical(state)) then
         ! 32-bit ptr
         ptr32(1:n+1) = int( ptr64(1:n+1) )
         if(pattern) then
            call rb_write(filename, matrix_type, m, n, ptr32, row, &
               write_options, flag, title=title, identifier=id)
         else
            call rb_write(filename, matrix_type, m, n, ptr32, row, &
               write_options, flag, val=val, title=title, identifier=id)
         endif
      else
         ! 64-bit ptr
         if(pattern) then
            call rb_write(filename, matrix_type, m, n, ptr64, row, &
               write_options, flag, title=title, identifier=id)
         else
            call rb_write(filename, matrix_type, m, n, ptr64, row, &
               write_options, flag, val=val, title=title, identifier=id)
         endif
      endif
      if(flag.ne.0) then
         write(*, "(a,/,a,i3)") "fail", "rb_write() returned", flag
         errors = errors + 1
         cycle
      endif

      ! Peek random matrix
      call rb_peek(filename, flag, m=m_in, n=n_in, nelt=nelt_in, nvar=nvar_in, &
         nval=nval_in, matrix_type=matrix_type_in, type_code=type_code_in,     &
         title=title_in, identifier=id_in)
      if(flag.ne.0) then
         write(*, "(a,/,a,i3)") "fail", "rb_peek() returned", flag
         errors = errors + 1
         cycle
      endif
      ! Check data
      if(.not.check_data(m,m_in,n,n_in,matrix_type=matrix_type, &
         matrix_type2=matrix_type_in,title=title,title2=title_in, &
         id=id,id2=id_in,type_code=type_code_in,pattern=pattern, &
         advance="no")) cycle
      write(*, "(a)", advance="no") "... "

      ! Read random matrix
      read_options = default_read_options
      if(random_logical(state)) then
         ! 32-bit ptr
         call rb_read(filename, m_in, n_in, ptr32_in, row_in, val_in, &
            read_options, flag, matrix_type=matrix_type_in, &
            title=title_in, identifier=id_in)
         if(allocated(ptr32_in)) then
            deallocate(ptr64_in, stat=st)
            allocate(ptr64_in(n_in+1))
            ptr64_in(1:n_in+1) = ptr32_in(1:n_in+1)
         endif
      else
         ! 64-bit ptr
         call rb_read(filename, m_in, n_in, ptr64_in, row_in, val_in, &
            read_options, flag, matrix_type=matrix_type_in, &
            title=title_in, identifier=id_in)
      endif
      if(flag.ne.0) then
         write(*, "(a,/,a,i3)") "fail", "rb_read() returned", flag
         errors = errors + 1
         cycle
      endif

      ! Check data
      if(pattern) then
         if(.not.check_data(m,m_in,n,n_in,ptr=ptr64,ptr2=ptr64_in, &
            row=row,row2=row_in,matrix_type=matrix_type, &
            matrix_type2=matrix_type_in,title=title,title2=title_in, &
            id=id,id2=id_in)) cycle
         if(allocated(val_in)) then
            write(*, "(a,/,a)") &
               "fail", "Read pattern matrix, but val is allocated"
            errors = errors + 1
            cycle
         endif
      else
         if(.not.check_data(m,m_in,n,n_in,ptr=ptr64,ptr2=ptr64_in, &
            row=row,row2=row_in,val=val,val2=val_in,matrix_type=matrix_type, &
            matrix_type2=matrix_type_in,title=title,title2=title_in, &
            id=id,id2=id_in)) cycle
      endif

      ! Generate non-default read_options and test those too
      read_options%add_diagonal = random_logical(state)
      read_options%lwr_upr_full = random_integer(state, 3)
      read_options%values = random_integer(state, 7) - 3
      if(read_options%values.le.0) &
         read_options%values = read_options%values - 2
      if(random_logical(state)) then
         ! Pass state in
         call rb_read(filename, m_in, n_in, ptr64_in, row_in, val_in, &
            read_options, flag, title=title_in, identifier=id_in, state=state)
      else
         ! Don't pass state in
         call rb_read(filename, m_in, n_in, ptr64_in, row_in, val_in, &
            read_options, flag, title=title_in, identifier=id_in)
      endif
      if(flag.ne.0) then
         write(*, "(a,/,a,i3)") "fail", "rb_read() returned", flag
         errors = errors + 1
         cycle
      endif

   end do

end subroutine test_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Check if two versions of data match. Return .true. if they do, else .false.
!> Print "pass" or "fail"+details depending on whether it works
logical function check_data(m, m2, n, n2, ptr, ptr2, row, row2, val, val2, &
      matrix_type, matrix_type2, title, title2, id, id2, pattern, type_code, &
      advance)
   integer, intent(in) :: m, m2, n, n2, matrix_type, matrix_type2
   integer, optional, intent(in) :: row(:), row2(:)
   integer(long), optional, intent(in) :: ptr(:), ptr2(:)
   real(wp), optional, intent(in) :: val(:), val2(:)
   character(len=72), optional :: title, title2
   character(len=8), optional :: id, id2
   logical, optional :: pattern
   character(len=3), optional :: type_code
   character(len=*), optional :: advance

   integer :: mtexpect
   real(wp) :: vdiff

   ! Check we read the same data we wrote
   if(m .ne. m2) then
      write(*, "(a,/,a,2i8)") "fail", "m != m2", m, m2
      errors = errors + 1
      check_data = .false.
      return
   endif
   if(n .ne. n2) then
      write(*, "(a,/,a,2i8)") "fail", "n != n2", n, n2
      errors = errors + 1
      check_data = .false.
      return
   endif
   if(present(ptr)) then
      if(any(ptr(1:n+1) .ne. ptr2(1:n+1))) then
         write(*, "(a,/,a)") "fail", "ptr != ptr2"
         errors = errors + 1
         check_data = .false.
         return
      endif
   endif
   if(present(row)) then
      if(any(row(1:ptr(n+1)-1) .ne. row2(1:ptr(n+1)-1))) then
         write(*, "(a,/,a)") "fail", "row != row2"
         errors = errors + 1
         check_data = .false.
         return
      endif
   endif
   if(present(val)) then
      vdiff = maxval(abs( val(1:ptr(n+1)-1) - val2(1:ptr(n+1)-1) ))
      if(vdiff.gt.epsilon(1.0_wp)) then
         write(*, "(a,/,a,e10.2)") "fail", "val != val2", vdiff
         errors = errors + 1
         check_data = .false.
         return
      endif
   endif
   mtexpect = matrix_type
   if(matrix_type.eq.SPRAL_MATRIX_REAL_SYM_PSDEF) &
      mtexpect = SPRAL_MATRIX_REAL_SYM_INDEF
   if(mtexpect .ne. matrix_type2) then
      write(*, "(a,/,a,'''',a,'''',1x,'''',a,'''')") &
         "fail", "matrix_type != matrix_type2", matrix_type, matrix_type2
      errors = errors + 1
      check_data = .false.
      return
   endif
   if(title .ne. title2) then
      write(*, "(a,/,a,'''',a,'''',1x,'''',a,'''')") &
         "fail", "title != title2", title, title2
      errors = errors + 1
      check_data = .false.
      return
   endif
   if(id .ne. id2) then
      write(*, "(a,/,a,'''',a,'''',1x,'''',a,'''')") &
         "fail", "id != id2", id, id2
      errors = errors + 1
      check_data = .false.
      return
   endif
   if(present(pattern)) then
      if(pattern .and. type_code(1:1).ne."p") then
         write(*, "(a,/,a,a1)") &
            "fail", "pattern=true but type_code(1:1)=", type_code(1:1)
         errors = errors + 1
         check_data = .false.
         return
      endif
      if(.not.pattern .and. type_code(1:1).ne."r") then
         write(*, "(a,/,a,a1)") &
            "fail", "pattern=false but type_code(1:1)=", type_code(1:1)
         errors = errors + 1
         check_data = .false.
         return
      endif
   endif

   ! Otherwise ok
   check_data = .true.
   if(present(advance)) then
      write(*, "(a)", advance=advance) "pass"
   else
      write(*, "(a)") "pass"
   endif
end function check_data

end program rutherford_boeing_test
