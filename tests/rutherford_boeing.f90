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

   errors = 0

   call test_random

   write(*, "(/a)") "=========================="
   write(*, "(a,i4)") "Total number of errors = ", errors

   if(errors.ne.0) stop 1 ! ERROR CODE for make check script

contains

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

   ! Matrix data (to read)
   character(len=72) :: title_in
   character(len=8) :: id_in
   integer :: m_in, n_in
   integer(long), dimension(:), allocatable :: ptr64_in
   integer, dimension(:), allocatable :: ptr32_in, row_in
   real(wp), dimension(:), allocatable :: val_in

   ! Options
   type(rb_read_options) :: read_options
   type(rb_write_options) :: write_options

   ! Working variables
   integer :: problem, flag
   integer(long) :: nnz
   real(wp) :: vdiff

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

      ! Write status line
      write(*,"(a,i5,a,i12,a,2i5,i8,a,i2,a)", advance="no") &
         "Test ", problem, " state ", random_get_seed(state), &
         " dimn =", m, n, nnz, " type =", matrix_type, " ... "

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
      call rb_write(filename, matrix_type_to_sym(matrix_type), &
         m, n, ptr64, row, val, write_options, flag, title=title, identifier=id)
      if(flag.ne.0) then
         write(*, "(a,/,a,i3)") "fail", "rb_write() returned", flag
         errors = errors + 1
         cycle
      endif

      ! Read random matrix
      call rb_read(filename, m_in, n_in, ptr64_in, row_in, val_in, &
         read_options, flag, title=title_in, identifier=id_in)
      if(flag.ne.0) then
         write(*, "(a,/,a,i3)") "fail", "rb_read() returned", flag
         errors = errors + 1
         cycle
      endif

      ! Check we read the same data we wrote
      if(m .ne. m_in) then
         write(*, "(a,/,a,2i8)") "fail", "m != m_in", m, m_in
         errors = errors + 1
         cycle
      endif
      if(n .ne. n_in) then
         write(*, "(a,/,a,2i8)") "fail", "n != n_in", n, n_in
         errors = errors + 1
         cycle
      endif
      if(any(ptr64(1:n+1) .ne. ptr64_in(1:n+1))) then
         write(*, "(a,/,a)") "fail", "ptr != ptr_in"
         errors = errors + 1
         cycle
      endif
      if(any(row(1:ptr64(n+1)-1) .ne. row_in(1:ptr64(n+1)-1))) then
         write(*, "(a,/,a)") "fail", "row != row_in"
         errors = errors + 1
         cycle
      endif
      vdiff = maxval(abs( val(1:ptr64(n+1)-1) - val_in(1:ptr64(n+1)-1) ))
      if(vdiff.gt.epsilon(1.0_wp)) then
         write(*, "(a,/,a,e10.2)") "fail", "val != val_in", vdiff
         errors = errors + 1
         cycle
      endif
      if(title .ne. title_in) then
         write(*, "(a,/,a,'''',a,'''',1x,'''',a,'''')") &
            "fail", "title != title_in", title, title_in
         errors = errors + 1
         cycle
      endif
      if(id .ne. id_in) then
         write(*, "(a,/,a,'''',a,'''',1x,'''',a,'''')") &
            "fail", "id != id_in", id, id_in
         errors = errors + 1
         cycle
      endif

      ! Otherwise ok
      write(*, "(a)") "ok"
   end do

end subroutine test_random

character(len=1) function matrix_type_to_sym(matrix_type)
   integer, intent(in) :: matrix_type

   select case(matrix_type)
   case(SPRAL_MATRIX_UNSPECIFIED)
      matrix_type_to_sym = "r"
   case(SPRAL_MATRIX_REAL_RECT)
      matrix_type_to_sym = "r"
   case(SPRAL_MATRIX_REAL_UNSYM)
      matrix_type_to_sym = "u"
   case(SPRAL_MATRIX_REAL_SYM_PSDEF)
      matrix_type_to_sym = "s"
   case(SPRAL_MATRIX_REAL_SYM_INDEF)
      matrix_type_to_sym = "s"
   case(SPRAL_MATRIX_REAL_SKEW)
      matrix_type_to_sym = "z"
   end select
end function matrix_type_to_sym

end program rutherford_boeing_test
