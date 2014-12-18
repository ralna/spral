!
! To convert from double:
! * Change wp
!
program main
   use spral_matrix_util, only : SPRAL_MATRIX_REAL_RECT, &
      SPRAL_MATRIX_REAL_SYM_INDEF
   use spral_random
   use spral_random_matrix, only : random_matrix_generate
   use spral_scaling
   implicit none

   integer, parameter :: wp = kind(0d0)
   real(wp), parameter :: err_tol = 5e-14

   integer :: errors

   type :: matrix_type
      integer :: m, n
      integer, dimension(:), allocatable :: ptr, row
      real(wp), dimension(:), allocatable :: val
   end type matrix_type

   ! Error flags
   integer, parameter :: SCALING_SUCCESS                   = 0
   integer, parameter :: SCALING_ERROR_ALLOCATION          = -1

   errors = 0

   call test_auction_sym_random
   call test_auction_unsym_random
   call test_equilib_sym_random
   call test_equilib_unsym_random
   call test_hungarian_sym_random
   call test_hungarian_unsym_random

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

   integer :: nza, prblm, i, j, k, nmatch
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

      call gen_random_sym(a, nza, state)
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
      match_check: &
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
            do k = a%ptr(j), a%ptr(j+1)-1
               if(a%row(k).eq.i) cycle match_check
            end do
            ! If we reach here, no entry in position (i,j)
            do k = a%ptr(i), a%ptr(i+1)-1
               if(a%row(k).eq.j) cycle match_check
            end do
            ! If we reach here, no entry in position (j,i) either
            write(*, "(a, i5, a, i5, a)") &
               "matched on (", i, ",", j, ") but no such entry"
            errors = errors + 1
            cycle prblm_loop
         endif
      end do match_check
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

subroutine test_auction_unsym_random
   integer :: maxn = 1000
   integer :: maxnz =  1000000
   integer, parameter :: nprob = 100

   real, parameter :: scale_tol = 2.0 ! How much above 1.0 can scaled value be?
   real, parameter :: max_tol = 0.01 ! How much less than 1.0 can col max be?
   real, parameter :: max_except = 0.05 ! proportion of bad entries allowed

   type(random_state) :: state

   type(matrix_type) :: a
   integer, allocatable, dimension(:) :: match, cnt
   real(wp), allocatable, dimension(:) :: rscaling, cscaling, rmax

   type(auction_options) :: options
   type(auction_inform) :: inform

   integer :: nza, prblm, i, j, k, rcnt, nmatch, nexcept
   real(wp) :: cmax, v


   write(*, "(a)")
   write(*, "(a)") "===================================================="
   write(*, "(a)") "Testing auction_scaling_unsym() with random matrices"
   write(*, "(a)") "===================================================="

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz))
   allocate(rscaling(maxn), cscaling(maxn), match(maxn), rmax(maxn), cnt(maxn))

   prblm_loop: &
   do prblm = 1, nprob

      ! Generate parameters
      a%n = random_integer(state, maxn)
      a%m = random_integer(state, maxn)
      if(random_integer(state, 2).eq.1) a%m = a%n ! 50% chance of unsym vs rect
      if (prblm < 21) then
         a%n = prblm ! check very small problems
         a%m = prblm ! check very small problems
      endif
      i = a%m*a%n/2 - max(a%m,a%n)
      i = max(0,i)
      nza = max(a%m,a%n) + random_integer(state, i)

      write(*, "(a, i3, a, i5, a, i5, a, i7, a, i2, a)",advance="no") &
         " - no. ", prblm,  " m = ", a%m, " n = ", a%n, " nza = ", nza, "..."

      !if(nza.gt.maxnz .or. a%m.gt.maxn .or. a%n.gt.maxn) then
      !   write(*, "(a)") "bad random matrix."
      !   write(*, "(a,i5,a,i5,a,i5)") "m = ", a%m, " n = ", a%n, " > maxn = ", maxn
      !   write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
      !   cycle
      !endif

      call gen_random_unsym(a, nza, state)
      !print *, "n = ", a%n
      !do i = 1, a%n
      !   print *, "col ", i, ":", a%row(a%ptr(i):a%ptr(i+1)-1)
      !   print *, "                 :", a%val(a%ptr(i):a%ptr(i+1)-1)
      !end do

      !
      ! Call scaling
      !
      call auction_scale_unsym(a%m, a%n, a%ptr, a%row, a%val, rscaling, &
         cscaling, options, inform, match=match)
      if(inform%flag .lt. 0) then
         write(*, "(a, i5)") "Returned inform%flag = ", inform%flag
         errors = errors + 1
         cycle prblm_loop
      endif
      !print *, "match = ", match(1:a%m)
      !print *, "rscal = ", rscaling(1:a%m)
      !print *, "cscal = ", cscaling(1:a%n)
      !print *
      !print "(a, 2es12.4)", "rscal = ", &
      !   minval(rscaling(1:a%m)), maxval(rscaling(1:a%m))
      !print "(a, 2es12.4)", "cscal = ", &
      !   minval(cscaling(1:a%n)), maxval(cscaling(1:a%n))

      !
      ! Ensure most rows and columns are matched, and none more than once
      !
      nmatch = 0
      cnt(1:a%n) = 0
      match_check: &
      do i = 1, a%m
         j = match(i)
         if(j.lt.0 .or. j.gt.a%n) then
            write(*, "(a, i5, a, i5)") "match(", i, ") = ", j
            errors = errors + 1
            cycle prblm_loop
         endif
         if(j.ne.0) then
            cnt(j) = cnt(j) + 1
            nmatch = nmatch + 1
            do k = a%ptr(j), a%ptr(j+1)-1
               if(a%row(k).eq.i) cycle match_check
            end do
            ! If we reach here, no entry in position (i,j)
            write(*, "(a, i5, a, i5, a)") &
               "matched on (", i, ",", j, ") but no such entry"
            errors = errors + 1
            cycle prblm_loop
         endif
      end do match_check
      if(nmatch < 0.9*min(a%m,a%n)) then
         write(*, "(a, i5, a, f4.1, a)") "Only matched ", nmatch, " pivots (", &
            (100.0*nmatch)/min(a%m,a%n), "%)"
         errors = errors + 1
         cycle prblm_loop
      endif
      if(any(cnt(1:a%m).gt.1)) then
         write(*, "(a)") "Mismatched row"
         errors = errors + 1
         cycle prblm_loop
      endif

      !
      ! Ensure all scaled entries are <= 1.0 and each row/col has an entry at 1
      !
      nexcept = 0
      rmax(1:a%m) = 0.0
      do i = 1, a%n
         cmax = 0.0;
         do j = a%ptr(i), a%ptr(i+1)-1
            v = abs(cscaling(i) * a%val(j) * rscaling(a%row(j)))
            if(v >= 1.0 + scale_tol) then
               nexcept = nexcept + 1
               !write(*, "(a, es12.4)") "Scaled entry = ", v
               !errors = errors + 1
               !cycle prblm_loop
            endif
            cmax = max( cmax, v )
            rmax(a%row(j)) = max( rmax(a%row(j)), v )
         end do
         if(a%ptr(i).ne.a%ptr(i+1) .and. cmax < 1.0 - max_tol) then
            nexcept = nexcept + 1
            !write(*, "(a, i4, a, es12.4)") "cmax(", i, ") = ", cmax
            !errors = errors + 1
            !cycle prblm_loop
         elseif(a%ptr(i).eq.a%ptr(i+1) .and. cscaling(i).ne. 1.0) then
            write(*, "(a, i4, a, es12.4, a, i5)") "cscaling(", i, &
               ") for empty col = ", cscaling(i), " match = ", match(i)
            errors = errors + 1
            cycle prblm_loop
         endif
      end do

      do i = 1, a%m
         if(rmax(i) < 1.0 - max_tol) then
            ! Check non-empty row before we complain
            rcnt = 0
            do j = 1, a%n
               do k = a%ptr(j), a%ptr(j+1)-1
                  if(a%row(k).eq.i) rcnt = rcnt + 1
               end do
            end do
            if(rcnt.ne.0) then
               nexcept = nexcept + 1
               !write(*, "(a, i4, a, es12.4)") "rmax(", i, ") = ", rmax(i)
               !errors = errors + 1
               !cycle prblm_loop
            endif
         endif
      end do
      if(nexcept > max_except*(a%ptr(a%n+1)-1)) then
         write(*, "(a, f12.2)") "Too many exceptional entries = ", &
            nexcept / (a%ptr(a%n+1)-1.0)
         errors = errors + 1
         cycle prblm_loop
      endif

      write(*, "(a)") "ok"

   end do prblm_loop

end subroutine test_auction_unsym_random

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

      call gen_random_sym(a, nza, state)
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

subroutine test_equilib_unsym_random
   integer :: maxn = 1000
   integer :: maxnz =  1000000
   integer, parameter :: nprob = 100
   type(random_state) :: state

   type(matrix_type) :: a
   real(wp), allocatable, dimension(:) :: rscaling, cscaling, rinf

   type(equilib_options) :: options
   type(equilib_inform) :: inform

   integer :: nza, prblm, i, j, k, rcnt
   real(wp) :: cmax, v

   write(*, "(a)")
   write(*, "(a)") "===================================================="
   write(*, "(a)") "Testing equilib_scaling_unsym() with random matrices"
   write(*, "(a)") "===================================================="

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz))
   allocate(rscaling(maxn), cscaling(maxn), rinf(maxn))

   prblm_loop: &
   do prblm = 1, nprob

      ! Generate parameters
      a%n = random_integer(state, maxn)
      a%m = random_integer(state, maxn)
      if(random_integer(state, 2).eq.1) a%m = a%n ! 50% chance of unsym vs rect
      if (prblm < 21) then
         a%n = prblm ! check very small problems
         a%m = prblm ! check very small problems
      endif
      i = a%m*a%n/2 - max(a%m,a%n)
      i = max(0,i)
      nza = max(a%m,a%n) + random_integer(state, i)

      write(*, "(a, i3, a, i5, a, i5, a, i7, a, i2, a)",advance="no") &
         " - no. ", prblm,  " m = ", a%m, " n = ", a%n, " nza = ", nza, "..."

      if(nza.gt.maxnz .or. a%n.gt.maxn .or. a%m.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5,a,i5)") "m = ", a%m, "n = ", a%n, &
            " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      call gen_random_unsym(a, nza, state)
      !print *, "n = ", a%n
      !do i = 1, a%n
      !   print *, "col ", i, ":", a%row(a%ptr(i):a%ptr(i+1)-1)
      !   print *, "                 :", a%val(a%ptr(i):a%ptr(i+1)-1)
      !end do

      !
      ! Call scaling
      !
      call equilib_scale_unsym(a%m, a%n, a%ptr, a%row, a%val, rscaling, &
         cscaling, options, inform)
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
         if(a%ptr(i).eq.a%ptr(i+1)) cycle ! Empty column
         cmax = 0.0;
         do j = a%ptr(i), a%ptr(i+1)-1
            v = abs(cscaling(i) * a%val(j) * rscaling(a%row(j)))
            cmax = max( cmax, v )
            rinf(a%row(j)) = max( rinf(a%row(j)), v )
         end do
         if(1.0-cmax > 0.05) then
            write(*, "(a, i4, a, es12.4)") "cinf(", i, ") = ", cmax
            errors = errors + 1
            cycle prblm_loop
         endif
      end do

      do i = 1, a%n
         if((1.0-rinf(i)) > 0.05) then
            ! Check non-empty row before we complain
            rcnt = 0
            do j = 1, a%n
               do k = a%ptr(j), a%ptr(j+1)-1
                  if(a%row(k).eq.i) rcnt = rcnt + 1
               end do
            end do
            if(rcnt > 0) then
               write(*, "(a, i4, a, es12.4)") "rinf(", i, ") = ", rinf(i)
               errors = errors + 1
               cycle prblm_loop
            endif
         endif
      end do

      write(*, "(a)") "ok"

   end do prblm_loop

end subroutine test_equilib_unsym_random

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

   integer :: nza, prblm, i, j, k
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

      call gen_random_sym(a, nza, state)
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
      match_check: &
      do i = 1, a%n
         j = match(i)
         if(j.lt.1 .or. j.gt.a%n) then
            write(*, "(a, i5, a, i5)") "match(", i, ") = ", j
            errors = errors + 1
            cycle prblm_loop
         endif
         cnt(j) = cnt(j) + 1
         do k = a%ptr(j), a%ptr(j+1)-1
            if(a%row(k).eq.i) cycle match_check
         end do
         ! If we reach here, no entry in position (i,j)
         do k = a%ptr(i), a%ptr(i+1)-1
            if(a%row(k).eq.j) cycle match_check
         end do
         ! If we reach here, no entry in position (j,i) either
         write(*, "(a, i5, a, i5, a)") &
            "matched on (", i, ",", j, ") but no such entry"
         errors = errors + 1
         cycle prblm_loop
      end do match_check
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_hungarian_unsym_random
   integer :: maxn = 1000
   integer :: maxnz =  1000000
   integer, parameter :: nprob = 100
   type(random_state) :: state

   type(matrix_type) :: a
   integer, allocatable, dimension(:) :: match, cnt
   real(wp), allocatable, dimension(:) :: rscaling, cscaling, rmax

   type(hungarian_options) :: options
   type(hungarian_inform) :: inform

   integer :: nza, prblm, i, j, k, nmatch, rcnt
   real(wp) :: cmax, v

   write(*, "(a)")
   write(*, "(a)") "======================================================"
   write(*, "(a)") "Testing hungarian_scaling_unsym() with random matrices"
   write(*, "(a)") "======================================================"

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz))
   allocate(rscaling(maxn), cscaling(maxn), match(maxn), rmax(maxn), cnt(maxn))

   prblm_loop: &
   do prblm = 1, nprob

      ! Generate parameters
      a%n = random_integer(state, maxn)
      a%m = random_integer(state, maxn)
      if(random_integer(state, 2).eq.1) a%m = a%n ! 50% chance of unsym vs rect
      if (prblm < 21) then
         a%n = prblm ! check very small problems
         a%m = prblm ! check very small problems
      endif
      i = a%m*a%n/2 - max(a%m,a%n)
      i = max(0,i)
      nza = max(a%m,a%n) + random_integer(state, i)

      write(*, "(a, i3, a, i5, a, i5, a, i7, a, i2, a)",advance="no") &
         " - no. ", prblm,  " m = ", a%m, " n = ", a%n, " nza = ", nza, "..."

      if(nza.gt.maxnz .or. a%n.gt.maxn .or. a%m.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5,a,i5)") "m = ", a%m, "n = ", a%n, &
            " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      call gen_random_unsym(a, nza, state)
      !print *, "n = ", a%n
      !do i = 1, a%n
      !   print *, "col ", i, ":", a%row(a%ptr(i):a%ptr(i+1)-1)
      !   print *, "                 :", a%val(a%ptr(i):a%ptr(i+1)-1)
      !end do

      !
      ! Call scaling
      !
      call hungarian_scale_unsym(a%m, a%n, a%ptr, a%row, a%val, rscaling, &
         cscaling, options, inform, match=match)
      if(inform%flag .lt. 0) then
         write(*, "(a, i5)") "Returned inform%flag = ", inform%flag
         errors = errors + 1
         cycle prblm_loop
      endif
      !print *, "match = ", match(1:a%n)
      !print *, "rscal = ", rscaling(1:a%m)
      !print *, "cscal = ", cscaling(1:a%n)

      !
      ! Ensure each row and column are matched [and on an entry that exists]
      !
      nmatch = 0
      cnt(1:a%n) = 0
      match_check: &
      do i = 1, a%m
         j = match(i)
         if(j.lt.0 .or. j.gt.a%n) then
            write(*, "(a, i5, a, i5)") "match(", i, ") = ", j
            errors = errors + 1
            cycle prblm_loop
         endif
         if(j.ne.0) then
            cnt(j) = cnt(j) + 1
            nmatch = nmatch + 1
            do k = a%ptr(j), a%ptr(j+1)-1
               if(a%row(k).eq.i) cycle match_check
            end do
            ! If we reach here, no entry in position (i,j)
            write(*, "(a, i5, a, i5, a)") &
               "matched on (", i, ",", j, ") but no such entry"
            errors = errors + 1
            cycle prblm_loop
         endif
      end do match_check
      if(nmatch .ne. min(a%m,a%n)) then
         write(*, "(a,i5,a,i5,a,i5,a)") &
            "Only matched ", nmatch, " in ", a%m, "x", a%n, " matrix"
         errors = errors + 1
         cycle prblm_loop
      endif
      if(any(cnt(1:a%m).gt.1)) then
         write(*, "(a)") "Mismatched row"
         errors = errors + 1
         cycle prblm_loop
      endif

      !
      ! Ensure all scaled entries are <= 1.0 and each row/col has an entry at 1
      !
      rmax(1:a%m) = 0.0
      do i = 1, a%n
         cmax = 0.0;
         do j = a%ptr(i), a%ptr(i+1)-1
            v = abs(cscaling(i) * a%val(j) * rscaling(a%row(j)))
            if(v >= 1.0+err_tol) then
               write(*, "(a, es12.4)") "Scaled entry = ", v
               errors = errors + 1
               cycle prblm_loop
            endif
            cmax = max( cmax, v )
            rmax(a%row(j)) = max( rmax(a%row(j)), v )
         end do
         if(cmax < 1.0 - err_tol .and. a%ptr(i).ne.a%ptr(i+1)) then
            write(*, "(a, i4, a, es12.4)") "cmax(", i, ") = ", cmax
            errors = errors + 1
            cycle prblm_loop
         endif
      end do

      do i = 1, a%m
         if(rmax(i) < 1.0-err_tol) then
            ! Check non-empty row before we complain
            rcnt = 0
            do j = 1, a%n
               do k = a%ptr(j), a%ptr(j+1)-1
                  if(a%row(k).eq.i) rcnt = rcnt + 1
               end do
            end do
            if(rcnt > 0) then
               write(*, "(a, i4, a, es12.4)") "rmax(", i, ") = ", rmax(i)
               errors = errors + 1
               cycle prblm_loop
            endif
         endif
      end do

      write(*, "(a)") "ok"

   end do prblm_loop

end subroutine test_hungarian_unsym_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_random_sym(a, nza, state, zr)
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

end subroutine gen_random_sym

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_random_unsym(a, nza, state)
   type(matrix_type), intent(inout) :: a
   integer, intent(in) :: nza
   type(random_state), intent(inout) :: state

   integer :: k, flag

   ! Generate a
   call random_matrix_generate(state, SPRAL_MATRIX_REAL_RECT, a%m, a%n, &
      nza, a%ptr, a%row, flag, val=a%val, nonsingular=.true., sort=.true.)
   if(flag.ne.0) print *, "Bad flag from random_matrix_generate()"

   ! make sure we have some large entries
   k = 1
   do while(k.le.a%ptr(a%n+1)-1)
      a%val(k) = 1000 * a%val(k)
      k = k + random_integer(state, 5)
   end do

end subroutine gen_random_unsym

end program
