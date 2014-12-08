!
! To convert from double:
! * Change wp
!
program main
   use spral_matrix_util, only : SPRAL_MATRIX_REAL_SYM_PSDEF, &
      SPRAL_MATRIX_REAL_SYM_INDEF
   use spral_random
   use spral_random_matrix, only : random_matrix_generate
   use spral_scaling
   implicit none

   integer, parameter :: wp = kind(0d0)
   real(wp), parameter :: err_tol = 5e-14

   integer :: errors

   type :: matrix_type
      integer :: n, ne
      integer, dimension(:), allocatable :: ptr, row
      real(wp), dimension(:), allocatable :: val
   end type matrix_type

   ! Error flags
   integer, parameter :: SCALING_SUCCESS                   = 0
   integer, parameter :: SCALING_ERROR_ALLOCATION          = -1

   errors = 0

   call test_auction_sym_random
   call test_equilib_sym_random
   call test_hungarian_sym_random

   write(*, "(/a)") "=========================="
   write(*, "(a,i4)") "Total number of errors = ", errors

   if(errors.ne.0) stop 1 ! ERROR CODE for make check script

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_auction_sym_random
   integer :: maxn = 1000
   integer :: maxnz =  1000000
   integer, parameter :: nprob = 100
   type(random_state) :: state

   type(matrix_type) :: a
   integer, allocatable, dimension(:) :: match, cnt
   real(wp), allocatable, dimension(:) :: scaling, rmax

   type(auction_options) :: options
   type(auction_inform) :: inform

   integer :: nza, prblm, i, j, nmatch
   real(wp) :: cmax, v

   write(*, "(a)")
   write(*, "(a)") "=================================================="
   write(*, "(a)") "Testing auction_scaling_sym() with random matrices"
   write(*, "(a)") "=================================================="

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz))
   allocate(scaling(maxn), match(maxn), rmax(maxn), cnt(maxn))

   prblm_loop: &
   do prblm = 1, nprob

      ! Generate parameters
      a%n = random_integer(state, maxn)
      if (prblm < 21) a%n = prblm ! check very small problems
      i = a%n**2/2 - a%n
      i = max(0,i)
      nza = a%n + random_integer(state, i)

      write(*, "(a, i3, a, i5, a, i7, a, i2, a)",advance="no") &
         " - no. ", prblm,  " n = ", a%n, " nza = ", nza, "..."

      if(nza.gt.maxnz .or. a%n.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5)") "n = ", a%n, " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      call gen_random_indef(a, nza, state)
      !print *, "n = ", a%n
      !do i = 1, a%n
      !   print *, "col ", i, ":", a%row(a%ptr(i):a%ptr(i+1)-1)
      !   print *, "                 :", a%val(a%ptr(i):a%ptr(i+1)-1)
      !end do

      !
      ! Call scaling
      !
      call auction_scale_sym(a%n, a%ptr, a%row, a%val, scaling, options, &
         inform, match=match)
      if(inform%flag .lt. 0) then
         write(*, "(a, i5)") "Returned inform%flag = ", inform%flag
         errors = errors + 1
         cycle prblm_loop
      endif
      !print *, "match = ", match(1:a%n)
      !print *, "scal = ", scaling(1:a%n)

      !
      ! Ensure most rows and columns are matched, and none more than once
      !
      nmatch = 0
      cnt(1:a%n) = 0
      do i = 1, a%n
         j = match(i)
         if(j.lt.0 .or. j.gt.a%n) then
            write(*, "(a, i5, a, i5)") "match(", i, ") = ", j
            errors = errors + 1
            cycle prblm_loop
         endif
         if(j.ne.0) then
            cnt(j) = cnt(j) + 1
            nmatch = nmatch + 1
         endif
      end do
      if(nmatch < 0.9*a%n) then
         write(*, "(a, i5, a, f4.1, a)") "Only matched ", nmatch, " pivots (", &
            (100.0*nmatch)/a%n, "%)"
         errors = errors + 1
         cycle prblm_loop
      endif
      if(any(cnt(1:a%n).gt.1)) then
         write(*, "(a)") "Mismatched row"
         errors = errors + 1
         cycle prblm_loop
      endif

      !
      ! Ensure all scaled entries are <= 1.0 and each row/col has an entry at 1
      !
      rmax(1:a%n) = 0.0
      do i = 1, a%n
         cmax = 0.0;
         do j = a%ptr(i), a%ptr(i+1)-1
            v = abs(scaling(i) * a%val(j) * scaling(a%row(j)))
            if(v >= 1.0+1.0) then
               write(*, "(a, es12.4)") "Scaled entry = ", v
               errors = errors + 1
               cycle prblm_loop
            endif
            cmax = max( cmax, v )
            rmax(a%row(j)) = max( rmax(a%row(j)), v )
         end do
         rmax(i) = max( rmax(i), cmax )
      end do

      do i = 1, a%n
         if(rmax(i) < 1.0-0.25) then
            write(*, "(a, i4, a, es12.4)") "rmax(", i, ") = ", rmax(i)
            errors = errors + 1
            cycle prblm_loop
         endif
      end do

      write(*, "(a)") "ok"

   end do prblm_loop

end subroutine test_auction_sym_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_equilib_sym_random
   integer :: maxn = 1000
   integer :: maxnz =  1000000
   integer, parameter :: nprob = 100
   type(random_state) :: state

   type(matrix_type) :: a
   real(wp), allocatable, dimension(:) :: scaling, rinf

   type(equilib_options) :: options
   type(equilib_inform) :: inform

   integer :: nza, prblm, i, j
   real(wp) :: cmax, v

   write(*, "(a)")
   write(*, "(a)") "=================================================="
   write(*, "(a)") "Testing equilib_scaling_sym() with random matrices"
   write(*, "(a)") "=================================================="

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz))
   allocate(scaling(maxn), rinf(maxn))

   prblm_loop: &
   do prblm = 1, nprob

      ! Generate parameters
      a%n = random_integer(state, maxn)
      if (prblm < 21) a%n = prblm ! check very small problems
      i = a%n**2/2 - a%n
      i = max(0,i)
      nza = a%n + random_integer(state, i)

      write(*, "(a, i3, a, i5, a, i7, a, i2, a)",advance="no") &
         " - no. ", prblm,  " n = ", a%n, " nza = ", nza, "..."

      if(nza.gt.maxnz .or. a%n.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5)") "n = ", a%n, " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      call gen_random_indef(a, nza, state)
      !print *, "n = ", a%n
      !do i = 1, a%n
      !   print *, "col ", i, ":", a%row(a%ptr(i):a%ptr(i+1)-1)
      !   print *, "                 :", a%val(a%ptr(i):a%ptr(i+1)-1)
      !end do

      !
      ! Call scaling
      !
      call equilib_scale_sym(a%n, a%ptr, a%row, a%val, scaling, options, inform)
      if(inform%flag .lt. 0) then
         write(*, "(a, i5)") "Returned inform%flag = ", inform%flag
         errors = errors + 1
         cycle prblm_loop
      endif
      !print *, "scal = ", scaling(1:a%n)

      !
      ! Ensure inf norm of all scaled rows  is close to 1.0
      !
      rinf(1:a%n) = 0.0
      do i = 1, a%n
         cmax = 0.0;
         do j = a%ptr(i), a%ptr(i+1)-1
            v = abs(scaling(i) * a%val(j) * scaling(a%row(j)))
            cmax = max( cmax, v )
            rinf(a%row(j)) = max( rinf(a%row(j)), v )
         end do
         rinf(i) = max( rinf(i), cmax )
      end do

      do i = 1, a%n
         if((1.0-rinf(i)) > 0.05) then
            write(*, "(a, i4, a, es12.4)") "rinf(", i, ") = ", rinf(i)
            errors = errors + 1
            cycle prblm_loop
         endif
      end do

      write(*, "(a)") "ok"

   end do prblm_loop

end subroutine test_equilib_sym_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_hungarian_sym_random
   integer :: maxn = 1000
   integer :: maxnz =  1000000
   integer, parameter :: nprob = 100
   type(random_state) :: state

   type(matrix_type) :: a
   integer, allocatable, dimension(:) :: match, cnt
   real(wp), allocatable, dimension(:) :: scaling, rmax

   type(hungarian_options) :: options
   type(hungarian_inform) :: inform

   integer :: nza, prblm, i, j
   real(wp) :: cmax, v

   write(*, "(a)")
   write(*, "(a)") "==================================================="
   write(*, "(a)") "Testing hungarian_scaling_sym() with random matrices"
   write(*, "(a)") "==================================================="

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz))
   allocate(scaling(maxn), match(maxn), rmax(maxn), cnt(maxn))

   prblm_loop: &
   do prblm = 1, nprob

      ! Generate parameters
      a%n = random_integer(state, maxn)
      if (prblm < 21) a%n = prblm ! check very small problems
      i = a%n**2/10 - a%n
      i = max(0,i)
      nza = a%n + random_integer(state, i)

      write(*, "(a, i3, a, i5, a, i7, a, i2, a)",advance="no") &
         " - no. ", prblm,  " n = ", a%n, " nza = ", nza, "..."

      if(nza.gt.maxnz .or. a%n.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5)") "n = ", a%n, " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      call gen_random_indef(a, nza, state)
      !print *, "n = ", a%n
      !do i = 1, a%n
      !   print *, "col ", i, ":", a%row(a%ptr(i):a%ptr(i+1)-1)
      !   print *, "                 :", a%val(a%ptr(i):a%ptr(i+1)-1)
      !end do

      !
      ! Call scaling
      !
      call hungarian_scale_sym(a%n, a%ptr, a%row, a%val, scaling, options, &
         inform, match=match)
      if(inform%flag .lt. 0) then
         write(*, "(a, i5)") "Returned inform%flag = ", inform%flag
         errors = errors + 1
         cycle prblm_loop
      endif
      !print *, "match = ", match(1:a%n)
      !print *, "scal = ", scaling(1:a%n)

      !
      ! Ensure each row and column are matched
      !
      cnt(1:a%n) = 0
      do i = 1, a%n
         j = match(i)
         if(j.lt.1 .or. j.gt.a%n) then
            write(*, "(a, i5, a, i5)") "match(", i, ") = ", j
            errors = errors + 1
            cycle prblm_loop
         endif
         cnt(j) = cnt(j) + 1
      end do
      if(any(cnt(1:a%n).ne.1)) then
         write(*, "(a)") "Mismatched row"
         errors = errors + 1
         cycle prblm_loop
      endif

      !
      ! Ensure all scaled entries are <= 1.0 and each row/col has an entry at 1
      !
      rmax(1:a%n) = 0.0
      do i = 1, a%n
         cmax = 0.0;
         do j = a%ptr(i), a%ptr(i+1)-1
            v = abs(scaling(i) * a%val(j) * scaling(a%row(j)))
            if(v >= 1.0+err_tol) then
               write(*, "(a, es12.4)") "Scaled entry = ", v
               errors = errors + 1
               cycle prblm_loop
            endif
            cmax = max( cmax, v )
            rmax(a%row(j)) = max( rmax(a%row(j)), v )
         end do
         rmax(i) = max( rmax(i), cmax )
      end do

      do i = 1, a%n
         if(rmax(i) < 1.0-err_tol) then
            write(*, "(a, i4, a, es12.4)") "rmax(", i, ") = ", rmax(i)
            errors = errors + 1
            cycle prblm_loop
         endif
      end do

      write(*, "(a)") "ok"

   end do prblm_loop

end subroutine test_hungarian_sym_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_random_indef(a, nza, state, zr)
   type(matrix_type), intent(inout) :: a
   integer, intent(in) :: nza
   type(random_state), intent(inout) :: state
   integer, optional, intent(in) :: zr ! if present, all entries in
     ! row zr are zero

   integer :: i, k, l, flag

   ! Generate a
   call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_INDEF, a%n, a%n, &
      nza, a%ptr, a%row, flag, val=a%val, nonsingular=.true., sort=.true.)
   if(flag.ne.0) print *, "Bad flag from random_matrix_generate()"

   if (present(zr)) then
      ! Scan along row
      do i = a%ptr(1), a%ptr(zr)-1
         if(a%row(i).eq.zr) a%val(i) = 0.0
      end do
      ! Scan along column
      do i = a%ptr(zr),a%ptr(zr+1)-1
         a%val(i) = 0.0
      end do
   elseif(a%n.gt.3) then
      ! Put some zeros on diagonal, observing first entry in column
      ! is always the diagonal after sorting
      ! but don't have all zeros in the col.
      l = random_integer(state,  a%n/2)
      do k = 1, a%n, max(1,l)
         if (a%ptr(k+1) > a%ptr(k) + 1) then
           i = a%ptr(k)
           a%val(i) = 0.0
         endif
      end do
      ! also make sure we have some large off diagonals
      do k = 1, a%n
         i = a%ptr(k+1) - 1
         a%val(i) = a%val(i)*1000
      end do
   endif

end subroutine gen_random_indef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_random_posdef(a, nza, state)
   type(matrix_type), intent(inout) :: a
   integer, intent(in) :: nza
   type(random_state), intent(inout) :: state

   integer :: i, j, k, flag
   real(wp) :: tempv

   ! Generate matrix
   call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_PSDEF, a%n, a%n, &
      nza, a%ptr, a%row, flag, val=a%val, nonsingular=.true., sort=.true.)
   if(flag.ne.0) print *, "Bad flag from random_matrix_generate()"

   ! Make a diagonally dominant, observing first entry in column
   ! is always the diagonal after sorting
   do k = 1, a%n
      tempv = 0.0
      do j = a%ptr(k)+1, a%ptr(k+1)-1
         tempv = tempv + abs(a%val(j))
         i = a%ptr(a%row(j))
         a%val(i) = a%val(i) + abs(a%val(j))
      end do
      i = a%ptr(k)
      a%val(i) = 1.0 + a%val(i) + tempv
   end do
end subroutine gen_random_posdef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_result(actual, expected, continued)
   integer :: actual
   integer :: expected
   logical, optional :: continued

   logical :: mycontinued

   mycontinued = .false.
   if(present(continued)) mycontinued = continued

   if(actual.eq.expected) then
      if(mycontinued) then
         write(*,"(a)", advance="no") "ok..."
      else
         write(*,"(a)") "ok"
      endif
      return
   endif

   write(*,"(a)") "fail"
   write(*,"(2(a,i4))") "returned ", actual, ", expected ", expected
   errors = errors + 1
end subroutine print_result


end program
