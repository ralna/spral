! FIXME: Actually test stuff predicted is actually happening. :(

program main
   use spral_nd
   implicit none

   integer, parameter :: we_unit = 10
   integer, parameter :: dl_unit = 11

   integer :: nerror

   open(unit=we_unit, file='we.out', status="replace")
   open(unit=dl_unit, file='dl.out', status="replace")

   nerror = 0

   call test_errors
   call test_input
   call test_dense
   call test_halflevelset
   call test_levelset
   call test_shift
   call test_multilevel
   call test_amd
   call test_misc
   call test_refine
   call test_fm

   write (*,'()')
   if (nerror.eq.0) then
     write (*,'(a)') 'All tests successfully completed'
   else
     write (*,'(a, i4, a)') 'There were ', nerror, ' errors'
   end if

   close (we_unit)
   close (dl_unit)

contains

subroutine setup_options(options)
   type (nd_options), intent(out) :: options

   type (nd_options) :: default_options

   options = default_options

   options%unit_error = we_unit
   options%unit_diagnostics = dl_unit
end subroutine setup_options

subroutine check_flag(actual, expect, advance)
   integer, intent(in) :: actual
   integer, intent(in) :: expect
   character(len=*), optional, intent(in) :: advance

   character(len=3) :: myadvance

   myadvance = "yes"
   if(present(advance)) myadvance = advance

   if(actual.eq.expect) then
      write (*, '(a)', advance=myadvance) "ok."
   else
      nerror = nerror + 1
      write (*, '(a)', advance=myadvance) "fail."
      write (*, '(a,i3,a,i3,a)') &
         "info%flag = ", actual, " (expected ", expect, ")"
   endif
end subroutine check_flag

subroutine check_success(n, perm, info, seps, expect_nsuper, expect_nzsuper)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   type(nd_inform), intent(in) :: info
   integer, dimension(n), optional, intent(in) :: seps
   integer, optional, intent(in) :: expect_nsuper
   integer, optional, intent(in) :: expect_nzsuper

   integer :: i
   integer, dimension(:), allocatable :: check

   ! Check flag first
   if(info%flag .ne. 0) then
      nerror = nerror + 1
      write (*, '(a)') "fail."
      write (*, '(a,i3,a,i3,a)') &
         "info%flag = ", info%flag, " (expected ", 0, ")"
      return
   endif

   ! Check we have a permutation
   allocate(check(n))
   check(:) = 0
   do i = 1, n
      if(perm(i).lt.1 .or. perm(i).gt.n) then
         nerror = nerror + 1
         write (*, '(a)') "fail."
         write (*, '(a,i3,a,i3,a,i3,a)') &
            "perm(", i, ") = ", perm(i), "(out of range 1:", n, ")"
         return
      endif
      check(perm(i)) = check(perm(i)) + 1
   end do
   if(any(check(1:n).ne.1)) then
      nerror = nerror + 1
      write (*, '(a)') "fail."
      write (*, '(a)') "failed to return a valid permutation"
      return
   endif

   ! Check seps is sensible (if present)
   if(present(seps)) then
      if(any(seps.lt.-1)) then
         nerror = nerror + 1
         write (*, '(a)') "fail."
         write (*, '(a)') "seps(:) contains values < -1"
         return
      endif
   endif

   ! Check nsuper (if present)
   if(present(expect_nsuper)) then
      if(info%nsuper.ne.expect_nsuper) then
         nerror = nerror + 1
         write (*, '(a)') "fail."
         write (*, '(a,i3,a,i3,a)') &
            "info%nsuper = ", info%nsuper, " (expected ", expect_nsuper, ")"
         return
      endif
   endif

   ! Check nzsuper (if present)
   if(present(expect_nzsuper)) then
      if(info%nzsuper.ne.expect_nzsuper) then
         nerror = nerror + 1
         write (*, '(a)') "fail."
         write (*, '(a,i3,a,i3,a)') &
            "info%nzsuper = ", info%nzsuper, " (expected ", expect_nzsuper, ")"
         return
      endif
   endif

   ! Otherwise, we're good!
   write (*, '(a)') "ok."
end subroutine check_success

subroutine test_errors
   integer :: n, ne
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Error Returns              *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test error flag n<1
   ! --------------------------------------
   n = 0
   ne = 0
   allocate (ptr(n+1), row(ne), perm(n), seps(n))
   ptr(1) = 1
   call setup_options(options)
   options%amd_switch1 = 2

   write (*, '(a)', advance="no") ' * Testing n=0, lower entry...'
   call nd_order(0, n, ptr, row, perm, options, info, seps)
   call check_flag(info%flag, -3)

   write (*, '(a)', advance="no") ' * Testing n=0, lower+upper entry...'
   call nd_order(1, n, ptr, row, perm, options, info, seps)
   call check_flag(info%flag, -3)
end subroutine test_errors

! *****************************************************************

subroutine test_input
   integer :: n, ne, i, j
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   logical :: corr
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Simple Inputs              *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test diagonal matrix
   ! --------------------------------------
   n = 10
   ne = 0
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = 1
   call setup_options(options)
   options%amd_switch1 = 2

   write (*, '(a)', advance="no") &
      ' * Testing diagonal matrix, lower entry.........'
   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_flag(info%flag, 0, advance="no")
   corr = .true.
   do j = 1, n
      corr = (corr .and. (perm(j)==j))
   end do
   if (corr) then
      write (*, '(a)') "ok."
   else
      nerror = nerror + 1
      write (*, '(a)') "fail."
      write (*, '(a, 10i5)') 'perm(:) = ', perm(:)
   endif

   write (*, '(a)', advance="no") &
      ' * Testing diagonal matrix, lower+upper entry...'
   call nd_order(1,n,ptr,row,perm,options,info,seps)
   call check_flag(info%flag, 0, advance="no")
   corr = .true.
   do j = 1, n
      corr = (corr .and. (perm(j)==j))
   end do
   if (corr) then
      write (*, '(a)') "ok."
   else
      nerror = nerror + 1
      write (*, '(a)') "fail."
      write (*, '(a, 10i5)') 'perm(:) = ', perm(:)
   endif

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expansion of matrix from lower triangle input
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing matrix expand, lower entry...........'

   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%print_level = 2

   n = 10
   ne = 9
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n-1) = (/ (i,i=1,n-1) /)
   ptr(n:n+1) = n
   row(1:ne) = (/ (i,i=2,n) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=18)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expansion of matrix from lower triangle format and that diags are
   ! removed
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing matrix expand w diag, lower entry....'

   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2

   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n) = (/ (2*i-1,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
     10 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=18)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expansion of matrix from lower and upper triangle input with no
   ! diags
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing matrix expand, lower+upper entry.....'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%print_level = 2

   n = 10
   ne = 18
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n) = (/ (2*(i-1),i=2,n) /)
   ptr(n+1) = ne + 1
   row(1) = 2
   do i = 2, n - 1
     row(ptr(i)) = i - 1
     row(ptr(i)+1) = i + 1
   end do
   row(ptr(n)) = n - 1

   call nd_order(1,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=18)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expansion of matrix from lower and upper triangle input with diag
   ! entry
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing matrix expand w diag, lower entry....'

   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%print_level = 2

   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n) = (/ (2*(i-1),i=2,n) /)
   ptr(n+1) = ne + 1
   row(1) = 2
   do i = 2, n - 1
     row(ptr(i)) = i - 1
     row(ptr(i)+1) = i + 1
   end do
   row(ptr(n)) = n - 1
   row(ne) = n

   call nd_order(1,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=18)

   deallocate (ptr,row,perm,seps)

end subroutine test_input


! *****************************************************************


subroutine test_dense
   integer :: n, ne, i
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Dense Inputs               *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test dense row removal - row 1 dense
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing single dense row, no removal.........'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 20
   options%print_level = 2
   options%remove_dense_rows = .false.

   n = 800
   ne = n
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n+1) = ne + 1
   row(1:ne) = (/ (i,i=1,n) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
     if (info%flag>=0 .and. info%nsuper==800 .and. info%nzsuper==1598) &
   call check_success(n, perm, info, seps=seps, expect_nsuper=800, &
      expect_nzsuper=1598)

   deallocate (ptr,row,perm,seps)

   write (*, '(a)', advance="no") &
      ' * Testing single dense row, with removal.......'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 20
   options%print_level = 2

   n = 800
   ne = n
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n+1) = ne + 1
   row(1:ne) = (/ (i,i=1,n) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nsuper=799, &
      expect_nzsuper=0)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test dense row removal - first row dense
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing subdiag + first row dense............'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 20

   n = 800
   ne = 2*n - 3
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n) = (/ (n+i-2,i=2,n) /)
   ptr(n+1) = ne + 1
   row(1:ptr(2)-1) = (/ (i,i=2,n) /)
   row(ptr(2):ptr(n)-1) = (/ (i+1,i=2,n-1) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nsuper=799, &
      expect_nzsuper=1596)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test dense row removal - last row dense
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing last row dense.......................'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 20

   n = 800
   ne = n - 1
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n-1) = (/ (i,i=1,n-1) /)
   ptr(n:n+1) = ne + 1
   row(1:ne) = n

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nsuper=799, &
      expect_nzsuper=0)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test dense row removal - first two rows dense but row 2 has max degree
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Testing two rows dense, row 2 max degree.....'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 20

   n = 800
   ne = n - 3 + n - 2
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2) = n - 2
   ptr(3:n+1) = ne + 1
   row(ptr(1):ptr(2)-1) = (/ (i+3,i=1,n-3) /)
   row(ptr(2):ne) = (/ (i+2,i=1,n-2) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nsuper=798, &
      expect_nzsuper=0)

   deallocate (ptr,row,perm,seps)

end subroutine test_dense


! *****************************************************************


subroutine test_halflevelset
   integer :: n, ne
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Half Level-set Method      *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test Half level-set (Ashcraft) method
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 1.......................................'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%print_level = 2

   n = 10
   ne = 14
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
   row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test Ashcraft method
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2.......................................'
   call setup_options(options)
   options%amd_switch1 = 2
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2

   n = 10
   ne = 15
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 4, 6, 8, 9, 10, 11, 13, 15, 16, 16 /)
   row(1:ne) = (/ 2, 3, 4, 3, 5, 4, 8, 6, 7, 9, 8, 10, 9, 10, 10 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test Ashcraft method
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 3.......................................'
   call setup_options(options)
   options%amd_switch1 = 2
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2

   n = 13
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 4, 6, 8, 10, 12, 13, 13, 14, 14, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 4, 11, 5, 6, 7, 8, 9, 10, 8, 10, 12, 13, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test Ashcraft method
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 4.......................................'
   call setup_options(options)
   options%amd_switch1 = 2
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%balance = 20.0

   n = 13
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 4, 6, 8, 10, 12, 13, 13, 14, 14, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 4, 11, 5, 6, 7, 8, 9, 10, 8, 10, 12, 13, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

 end subroutine test_halflevelset


! *****************************************************************


subroutine test_levelset
   integer :: n, ne
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Full Level-set Method      *"
   write (*, '(a)') "****************************************"


   ! --------------------------------------
   ! Test one-sided levelset method
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 1.......................................'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%balance = 20.0
   options%print_level = 2

   n = 10
   ne = 14
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
   row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test one-sided levelset method
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2.......................................'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%balance = 1.5

   n = 10
   ne = 14
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
   row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

end subroutine test_levelset

! *****************************************************************


subroutine test_shift
   integer :: n, ne
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Shift Refinement           *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test DM-style refinement with 2-sided partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM-style, 2-sided part, refinement=2....'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%print_level = 2
   options%refinement = 2

   n = 8
   ne = 11
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 4, 5, 6, 8, 9, 11, 12, 12 /)
   row(1:ne) = (/ 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 8 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)


   write (*, '(a)', advance="no") &
      ' * Test DM-style, 2-sided part, refinement=5....'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%print_level = 2
   options%refinement = 5

   n = 8
   ne = 11
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 4, 5, 6, 8, 9, 11, 12, 12 /)
   row(1:ne) = (/ 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 8 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 2-sided, test 1..............'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%refinement = 2
   options%amd_call = 2

   n = 21
   ne = 37
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 13, 15, 16, 20, 22, 24, 26, 28, 29, &
     29, 32, 33, 35, 37, 38, 38 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 16, &
     17, 11, 13, 12, 13, 13, 15, 14, 15, 15, 17, 18, 19, 18, 19, 20, 20, &
     21, 21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 1-sided partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 1-sided, test 1..............'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%refinement = 2
   options%amd_call = 2
   options%balance = 33.0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))

   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
      27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
      13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
      24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 2-sided, test 2..............'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%refinement = 2
   options%amd_call = 2
   options%balance = 12.0

   n = 7
   ne = 6
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 6, 7, 7 /)
   row(1:ne) = (/ 7, 7, 7, 7, 7, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 1-sided partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 1-sided, test 2..............'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%refinement = 5
   options%amd_call = 2
   options%balance = 12.0

   n = 7
   ne = 6
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 6, 7, 7 /)
   row(1:ne) = (/ 7, 7, 7, 7, 7, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 1-sided partition with balanced partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 1-sided, balanced part.......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%refinement = 5
   options%balance = 2.0
   options%amd_call = 2
   options%print_level = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
     27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
     13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
     24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 2-sided, balanced part, T1...'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%refinement = 5
   options%balance = 2.0
   options%amd_call = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
     27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
     13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
     24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 2-sided, balanced part, T2...'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%refinement = 2
   options%balance = 1.05
   options%refinement_band = n
   options%amd_call = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)


   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 2-sided, balanced part, T3...'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%refinement = 2
   options%balance = 1.30
   options%refinement_band = n
   options%amd_call = 2

   n = 13
   ne = 15
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 6, 7, 8, 9, 11, 11, 13, 14, 15, 16, 16 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test DM refine, 2-sided, balanced part, T4...'
   call setup_options(options)
   options%amd_switch1 = 5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%refinement = 2
   options%balance = 1.30
   options%refinement_band = n
   options%amd_call = 2

   n = 13
   ne = 15
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 16 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 7, 7, 8, 9, 10, 11, 11, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

end subroutine test_shift

! *****************************************************************


subroutine test_multilevel
   integer :: n, ne, i, j, k
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Multilevel                 *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test multilevel, 2-sided partitioning
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, 2-sided......................'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%print_level = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! -------------------------------------------------------
   ! Test multilevel, 1-sided partitioning, common neighbours
   ! -------------------------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, 1-sided, matching=0..........'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%print_level = 2
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! -------------------------------------------------------
   ! Test multilevel, 1-sided partitioning, heavy edge
   ! -------------------------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, 1-sided, matching=1..........'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%print_level = 2
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! -------------------------------------------------------
   ! Test multilevel, 1-sided partitioning, heavy edge (matching=2)
   ! -------------------------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, 1-sided, matching=2..........'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%print_level = 2
   options%matching = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel coarsening
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel Coarsening, Test 1............'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%print_level = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel coarsening, common neighbours
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel Coarsening, matching=0........'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%matching = 0
   options%print_level = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel coarsening, heavy-edge
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel Coarsening, matching=1........'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%matching = 1
   options%print_level = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel coarsening
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel Coarsening, Test 2............'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%min_reduction = 0.4
   options%matching = 1
   options%print_level = 2
   options%matching = 2

   n = 12
   ne = 18
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 4, 6, 7, 9, 10, 12, 13, 16, 18, 19, 19 /)
   row(1:ne) = (/ 2, 9, 9, 4, 10, 10, 6, 11, 11, 8, 12, 12, 10, 11, 12, 11, &
     12, 12 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=0, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 0
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=0, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 0
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=1, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 1
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=1, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 1
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=2, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 2
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=2, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 2
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=3, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 3
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=3, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 3
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=4, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 4
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=4, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 4
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=5, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 5
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=5, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 5
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=6, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 6
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=6, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 6
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=7, matching=0......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 0
   options%matching = 0

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel refinement
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel refinement=7, matching=1......'
   call setup_options(options)
   options%amd_switch1 = 5
   options%balance = 1.05
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 5
   options%refinement = 0
   options%matching = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic multilevel-no multilevel check
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, No Multilevel Check..........'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%balance = 8.0
   options%partition_method = 2
   options%coarse_partition_method = 2
   options%ml_call = 10
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7

   n = 200
   ne = 199
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n) = (/ (i,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ (i,i=2,n) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic multilevel-no multilevel check
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, No Multilevel Check, CPM=1...'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%balance = 8.0
   options%partition_method = 2
   options%ml_call = 10
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7
   options%coarse_partition_method = 1

   n = 200
   ne = 7*(n-7)+21
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   j = 1
   do i = 1,n
      ptr(i) = j
      do k = i+1, min(i+7,n)
         row(j) = k
         j = j + 1
      end do
   end do
   ptr(n+1) = j

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic multilevel-no multilevel check
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, No Multilevel Check, CPM=2...'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%balance = 8.0
   options%partition_method = 2
   options%ml_call = 10
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7
   options%coarse_partition_method = 2

   n = 200
   ne = 7*(n-7)+21
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   j = 1
   do i = 1,n
      ptr(i) = j
      do k = i+1, min(i+7,n)
         row(j) = k
         j = j + 1
      end do
   end do
   ptr(n+1) = j

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel reduction ratio
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, Reduction Ratio, CPM=1.......'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%balance = 2.0
   options%partition_method = 1
   options%ml_call = 10
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7
   options%find_supervariables = .false.
   options%print_level = 2
   options%coarse_partition_method = 1

   n = 20
   ne = n-1
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n+1) = ne + 1
   row(1:ne) = (/ (i,i=2,n) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test multilevel reduction ratio
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test Multilevel, Reduction Ratio, CPM=2.......'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%balance = 2.0
   options%partition_method = 1
   options%ml_call = 10
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7
   options%find_supervariables = .false.
   options%print_level = 2
   options%coarse_partition_method = 2

   n = 20
   ne = n-1
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n+1) = ne + 1
   row(1:ne) = (/ (i,i=2,n) /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for full galerkin_graph_rap coverage
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test galerkin_graph_rap().....................'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 1
   options%print_level = 0
   options%refinement = 1

   n = 16
   ne = 42
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   j = 1
   do i = 1, 16
     ptr(i) = j
     if (mod(i,4)/=0) then
       row(j) = i + 1
       j = j + 1
     end if
     if (i<13) then
       if (mod(i,4)==1) then
         row(j) = i + 4
         j = j + 1
         row(j) = i + 5
         j = j + 1
       else
         if (mod(i,4)==0) then
           row(j) = i + 3
           j = j + 1
           row(j) = i + 4
           j = j + 1
         else
           row(j) = i + 3
           j = j + 1
           row(j) = i + 4
           j = j + 1
           row(j) = i + 5
           j = j + 1
         end if
       end if
     end if
   end do
   ptr(n+1) = j

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

end subroutine test_multilevel

! *****************************************************************


subroutine test_amd
   integer :: n, ne
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing AMD                        *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Call AMD with no nested
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test AMD, non nested, Test 1..................'
   call setup_options(options)
   options%amd_call = 40
   options%print_level = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
     27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
     13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
     24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Call AMD with no nested
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test AMD, non nested, Test 2..................'
   call setup_options(options)
   options%amd_call = 40
   options%print_level = 2

   n = 31
   ne = 1
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   ptr(2:n+1) = 2
   row(1:ne) = (/ 2 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

end subroutine test_amd

! *****************************************************************

subroutine test_refine
   integer :: n, ne
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Refine                     *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 1........................................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 0

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2, refinement=1..........................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 1

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2, refinement=2..........................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 2

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2, refinement=3..........................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 3

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2, refinement=4..........................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 4

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2, refinement=5..........................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 5

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2, refinement=6..........................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 6

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for refinement optionss
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test 2, refinement=7..........................'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%refinement = 7

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 1, refine=1....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement_band = -1
   options%refinement = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 1, refine=2....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement_band = -1
   options%refinement = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 1, refine=3....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement_band = -1
   options%refinement = 3

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 1, refine=4....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement_band = -1
   options%refinement = 4

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 1, refine=5....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement_band = -1
   options%refinement = 5

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 1, refine=6....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement_band = -1
   options%refinement = 6

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 1, refine=7....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 2
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement_band = -1
   options%refinement = 7

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 2, refine=1....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 1

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 2, refine=2....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 2

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 2, refine=3....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 3

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 2, refine=4....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 4

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 2, refine=5....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 5

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 2, refine=6....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 6

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 2, refine=7....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 3, refine=1....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 1

   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 3, refine=2....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 2

   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 3, refine=3....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 3

   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 3, refine=4....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 4

   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 3, refine=5....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 5

   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 3, refine=6....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 6

   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 3, refine=7....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 1
   options%ml_call = 4
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 7

   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for full refine_trim coverage
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test refine_trim(), Test 1....................'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 1
   options%print_level = 0
   options%refinement = 1

   n = 16
   ne = 22
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 7, 8, 9, 12, 14, 16, 17, 18, 20, 21, 22, 23, &
     23 /)
   row(1:ne) = (/ 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13, 13, 14, &
     14, 15, 15, 16, 16 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for full refine_trim coverage
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test refine_trim(), Test 2....................'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 1
   options%print_level = 2
   options%refinement = 1

   n = 16
   ne = 22
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 14, 15, 16, 18, 19, 20, 22, &
     23, 23 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 10, 11, 12, 13, &
     14, 15, 15, 16, 16 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for full refine_trim coverage
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test refine_trim(), Test 3....................'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%min_reduction = 0.5
   options%balance = 1.5
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 1
   options%print_level = 2
   options%refinement = 1

   n = 25
   ne = 44
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 5, 8, 10, 12, 14, 16, 18, 20, 21, 22, 23, 29, 30, &
     30, 32, 34, 36, 38, 39, 41, 43, 44, 45, 45 /)
   row(1:ne) = (/ 2, 3, 4, 5, 4, 6, 7, 5, 8, 9, 10, 7, 13, 8, 13, 9, 13, &
     10, 13, 13, 12, 13, 14, 16, 17, 18, 19, 20, 15, 17, 21, 18, 21, 19, &
     22, 20, 23, 23, 22, 24, 23, 24, 24, 25 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for not checking balance
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test for not checking balance.................'
   call setup_options(options)
   options%amd_switch1 = 4
   options%stop_coarsening1 = 3
   options%min_reduction = 0.5
   options%balance = real(n+1)
   options%partition_method = 1
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 1
   options%print_level = 0
   options%refinement = 7

   n = 25
   ne = 44
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 5, 8, 10, 12, 14, 16, 18, 20, 21, 22, 23, 29, 30, &
     30, 32, 34, 36, 38, 39, 41, 43, 44, 45, 45 /)
   row(1:ne) = (/ 2, 3, 4, 5, 4, 6, 7, 5, 8, 9, 10, 7, 13, 8, 13, 9, 13, &
     10, 13, 13, 12, 13, 14, 16, 17, 18, 19, 20, 15, 17, 21, 18, 21, 19, &
     22, 20, 23, 23, 22, 24, 23, 24, 24, 25 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 4, refine=1....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 1

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 4, refine=2....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 2

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 4, refine=3....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 3

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 4, refine=4....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 4

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 4, refine=5....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 5

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 4, refine=6....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 6

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 4, refine=7....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 7

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 5, refine=1....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 10.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 1

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 5, refine=2....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 10.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 2

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 5, refine=3....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 10.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 3

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 5, refine=4....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 10.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 4

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 5, refine=5....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 10.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 5

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 5, refine=6....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 10.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 6

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 5, refine=7....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 10.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%stop_coarsening1 = 2
   options%max_improve_cycles = 2
   options%refinement = 7

   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 6, refine=1....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 1

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 6, refine=2....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 2

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 6, refine=3....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 3

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 6, refine=4....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 4

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 6, refine=5....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 5

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 6, refine=6....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 6

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 6, refine=7....'
   call setup_options(options)
   options%amd_switch1 = 3
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7

   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 7, refine=1....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 1

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 7, refine=2....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 2

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 7, refine=3....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 3

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 7, refine=4....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 4

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 7, refine=5....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 5

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 7, refine=6....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 6

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test exand-refind cycles, test 7, refine=7....'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 8.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%refinement = 7

   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T1, refine=1..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 1

   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T1, refine=2..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 2

   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T1, refine=3..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 3

   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T1, refine=4..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 4

   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T1, refine=5..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 5

   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T1, refine=6..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 6

   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T1, refine=7..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 7

   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T2, refine=1..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 1

   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T2, refine=2..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 2

   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T2, refine=3..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 3

   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T2, refine=4..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 4

   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T2, refine=5..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 5

   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T2, refine=6..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 6

   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test automatic choice of shift, T2, refine=7..'
   call setup_options(options)
   options%amd_switch1 = 4
   options%balance = 1.01
   options%partition_method = 1
   options%coarse_partition_method = 2
   options%refinement_band = 0
   options%stop_coarsening1 = 3
   options%stop_coarsening2 = 3
   options%amd_call = 2
   options%max_improve_cycles = 2
   options%print_level = 2
   options%refinement = 7

   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

end subroutine test_refine

! *****************************************************************

subroutine test_misc
   integer :: n, ne, i, j, k
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: options
   type (nd_inform) :: info

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing Miscellaneous              *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test diagonal submatrix - not diag!!!!
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test diagonal *submatrix*.....................'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%print_level = 1

   n = 10
   ne = 10
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n-2) = (/ (i,i=1,n-2) /)
   ptr(n-1) = ptr(n-2) + 2
   ptr(n:n+1) = ne + 1
   row(1:ne-2) = n
   row(ne-1) = n - 1
   row(ne) = n

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=16)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test independent component detection
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test independent component detection, CPM=1...'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%print_level = 2
   options%coarse_partition_method = 1

   n = 5
   ne = 6
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 6, 7 /)
   row(1:ne) = (/ 2, 1, 4, 3, 5, 4 /)

   call nd_order(1,n,ptr,row,perm,options,info)
   call check_success(n, perm, info)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test independent component detection
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test independent component detection, CPM=2...'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%print_level = 2
   options%coarse_partition_method = 2

   n = 5
   ne = 6
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 6, 7 /)
   row(1:ne) = (/ 2, 1, 4, 3, 5, 4 /)

   call nd_order(1,n,ptr,row,perm,options,info)
   call check_success(n, perm, info)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for supervariable detection turned off
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test supervariable detection disabled.........'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%find_supervariables = .false.

   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n) = (/ (2*i-1,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
     10 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=18)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for full matrix
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test full matrix..............................'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%find_supervariables = .false.
   options%print_level = 2

   n = 10
   ne = 45
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1) = 1
   do i = 2, n
     ptr(i) = ptr(i-1) + n - i + 1
   end do
   ptr(n+1) = ne + 1

   j = 1
   do i = 1, n - 1
     do k = i + 1, n
       row(j) = k
       j = j + 1
     end do
   end do

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=90)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test cost function cost2
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test cost function 2..........................'
   call setup_options(options)
   options%amd_switch1 = 2
   options%amd_call = 2
   options%find_supervariables = .false.
   options%cost_function = 2

   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n) = (/ (2*i-1,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
     10 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps, expect_nzsuper=18)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test for no separator return, partition 2 larger than amd_switch1 but
   ! partition 1 equal to amd_switch1
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test no sep, p2>amd_switch1, p1==amd_switch1..'
   call setup_options(options)
   options%amd_switch1 = 6
   options%amd_call = 2
   options%find_supervariables = .false.
   options%partition_method = 0
   options%print_level = 2

   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   call nd_order(0,n,ptr,row,perm,options,info)
   call check_success(n, perm, info, expect_nzsuper=32)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test matrix extraction involving end row
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test matrix extraction w end row, Test 1......'
   call setup_options(options)
   options%amd_switch1 = 7
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 0
   options%print_level = 2
   options%refinement = 7

   n = 14
   ne = 26
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 6, 8, 10, 13, 17, 19, 23, 24, 26, 27, 27, 27, &
     27 /)
   row(1:ne) = (/ 2, 3, 3, 4, 5, 5, 14, 5, 13, 6, 13, 14, 8, 9, 13, 14, 10, &
     13, 9, 10, 11, 13, 11, 11, 12, 12 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test matrix extraction involving end row
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test matrix extraction w end row, Test 2......'
   call setup_options(options)
   options%amd_switch1 = 7
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 0
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 0
   options%print_level = 2
   options%refinement = 7

   n = 16
   ne = 31
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1, 3, 6, 8, 10, 13, 14, 16, 21, 25, 27, 29, 31, 31, 31, &
     32, 32 /)
   row(1:ne) = (/ 2, 3, 3, 4, 5, 5, 6, 5, 15, 6, 15, 16, 16, 11, 15, 9, 11, &
     12, 15, 16, 10, 12, 13, 16, 16, 13, 12, 14, 13, 14, 16 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)

   ! --------------------------------------
   ! Test no balanced partition found so switch to multilevel
   ! --------------------------------------
   write (*, '(a)', advance="no") &
      ' * Test no balanced partition so switch to ml....'
   call setup_options(options)
   options%amd_switch1 = 7
   options%stop_coarsening1 = 2
   options%min_reduction = 0.5
   options%balance = 1.0
   options%partition_method = 2
   options%coarse_partition_method = 1
   options%amd_call = 2
   options%max_improve_cycles = 0
   options%print_level = 2
   options%refinement = 7

   n = 18
   ne = 41
   allocate (ptr(n+1),row(ne),perm(n),seps(n))
   ptr(1:n+1) = (/ 1,2,3,8,10,12,15,19,21,24,28,30,33,37,39,40,42,42,42 /)
   row(1:ne) = (/ 2,3,4,5,6,7,8,6,7,7,8,7,9,10,8,9,10,11,10,11,10,12,13,11,&
        12,13,14,13,14,13,15,16,14,15,16,17,16,17,16,17,18 /)

   call nd_order(0,n,ptr,row,perm,options,info,seps)
   call check_success(n, perm, info, seps=seps)

   deallocate (ptr,row,perm,seps)
end subroutine test_misc

! *****************************************************************

subroutine test_fm
   integer :: n, ne, an1,an2,swgt,aw1,aw2,aws
   integer, allocatable, dimension(:) :: ptr,row,wgt,work,part
   type (nd_options) :: options

   write (*, '()')
   write (*, '(a)') "****************************************"
   write (*, '(a)') "*   Testing FM                         *"
   write (*, '(a)') "****************************************"

   ! --------------------------------------
   ! Test fm with empty partition 2
   ! --------------------------------------
   write (*, '(a)') &
      ' * Test FM with empty part 2.....................'
   call setup_options(options)
   options%print_level = 1

   n = 11
   ne = 14
   swgt = n
   allocate (ptr(n),row(ne),wgt(n),work(8*n+swgt),part(n))
   ptr(1:n) = (/ 1,2,4,6,8,9,11,12,13,14,15  /)
   row(1:ne) = (/ 2,3,4,5,6,6,7,8,8,9,9,10,10,11  /)
   wgt(:)=1
   an1 = 3
   an2 = 0
   aw1 = an1
   aw2 = an2
   aws = swgt - aw1-aw2
   part(1:n) = (/ 1,2,4,8,10,11,3,5,6,7,9 /)

   call nd_refine_fm(n,ne,ptr,row,wgt,swgt,an1,an2,aw1,aw2,aws,part,work,&
      options)
   ! FIXME: No check routine

   deallocate (ptr,row,wgt,work,part)

   ! --------------------------------------
   ! Test fm with empty partition 1
   ! --------------------------------------
   write (*, '(a)') &
      ' * Test FM with empty part 1.....................'
   call setup_options(options)
   options%print_level = 1

   n = 11
   ne = 14
   swgt = n
   allocate (ptr(n),row(ne),wgt(n),work(8*n+swgt),part(n))
   ptr(1:n) = (/ 1,2,4,6,8,9,11,12,13,14,15  /)
   row(1:ne) = (/ 2,3,4,5,6,6,7,8,8,9,9,10,10,11  /)
   wgt(:)=1
   an1 = 0
   an2 = 3
   aw1 = an1
   aw2 = an2
   aws = swgt - aw1-aw2
   part(1:n) = (/ 1,2,4,8,10,11,3,5,6,7,9 /)

   call nd_refine_fm(n,ne,ptr,row,wgt,swgt,an1,an2,aw1,aw2,aws,part,work,&
      options)
   ! FIXME: No check routine

   deallocate (ptr,row,wgt,work,part)

end subroutine test_fm

end program main
