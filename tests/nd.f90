program main
   use spral_nd
   implicit none

   integer, parameter :: we_unit = 10
   integer, parameter :: dl_unit = 11

   integer :: nerror
   logical :: ok

   open(unit=we_unit, file='we.out', status="replace")
   open(unit=dl_unit, file='dl.out', status="replace")

   nerror = 0

   call test_errors
   call test_input(ok)
   if(.not. ok) nerror = nerror + 1
   call test_dense(ok)
   if(.not. ok) nerror = nerror + 1
   call test_ashcraft(ok)
   if(.not. ok) nerror = nerror + 1
   call test_levelset(ok)
   if(.not. ok) nerror = nerror + 1
   call test_shift(ok)
   if(.not. ok) nerror = nerror + 1
   call test_multigrid(ok)
   if(.not. ok) nerror = nerror + 1
   call test_amd(ok)
   if(.not. ok) nerror = nerror + 1
   call test_misc(ok)
   if(.not. ok) nerror = nerror + 1
   call test_refine(ok)
   if(.not. ok) nerror = nerror + 1
   call test_fm(ok)
   if(.not. ok) nerror = nerror + 1

   if (nerror.eq.0) then
     write (*,'(a)') 'All tests successfully completed'
   else
     write (*,'(a, i4, a)') 'Failed ', nerror, ' tests'
   end if

   close (we_unit)
   close (dl_unit)

contains

   subroutine reset_control(control_orig,control_reset)


     type (nd_options), INTENT (IN) :: control_orig

     type (nd_options), INTENT (OUT) :: control_reset

     control_reset%print_level = control_orig%print_level
     control_reset%unit_diagnostics = control_orig%unit_diagnostics
     control_reset%unit_error = control_orig%unit_error
     control_reset%amd_call = control_orig%amd_call
     control_reset%amd_switch1 = control_orig%amd_switch1
     control_reset%amd_switch2 = control_orig%amd_switch2
     control_reset%cost_function = control_orig%cost_function
     control_reset%partition_method = control_orig%partition_method
     control_reset%matching = control_orig%matching
     control_reset%coarse_partition_method = control_orig% &
       coarse_partition_method
     control_reset%refinement = control_orig%refinement
     control_reset%refinement_band = control_orig%refinement_band
     control_reset%remove_dense_rows = control_orig%remove_dense_rows
     control_reset%stop_coarsening2 = control_orig%stop_coarsening2
     control_reset%stop_coarsening1 = control_orig%stop_coarsening1
     control_reset%ml_call = control_orig%ml_call
     control_reset%min_reduction = control_orig%min_reduction
     control_reset%max_reduction = control_orig%max_reduction
     control_reset%balance = control_orig%balance
     control_reset%max_improve_cycles = control_orig%max_improve_cycles
     control_reset%find_supervariables = control_orig%find_supervariables


   end subroutine reset_control

   subroutine testpr(n,perm)
     ! .. Scalar Arguments ..
     integer n
     ! ..
     ! .. Array Arguments ..
     integer perm(n)
     ! ..
     ! .. Local Scalars ..
     integer i, ip
     ! .. Local Arrays
     integer, allocatable, dimension (:) :: iw
     ! ..
     ! Initialize array used to test for valid permutation.
     allocate (iw(n))
     iw(1:n) = 1
     do i = 1, n
       ip = abs(perm(i))
       if (iw(ip)==0) then
         write (68,*) i, ip
         go to 10
       else
         iw(ip) = 0
       end if
     end do
     deallocate (iw)
     return
10      write (6,'(A)') '**** Error in permutation'
     deallocate (iw)
   end subroutine testpr

   subroutine ndef(n,nz,a,ip,ind,w,sym)
     ! This subroutine augments the diagonal to make the matrix pos def
     ! Then large entries put in off-diagonal to make it a little not so.
     ! .. Parameters ..
     doUBLE PRECISION zero, one
     PARAMETER (zero=0.0D0,one=1.0D0)
     ! ..
     ! .. Scalar Arguments ..
     integer n, nz
     logical sym
     ! ..
     ! .. Array Arguments ..
     doUBLE PRECISION a(nz), w(n)
     integer ind(nz), ip(*)
     ! ..
     ! .. Local Scalars ..
     doUBLE PRECISION rmax
     integer i, i1, i2, idiag, ii, ioff, j, maxoff, mprint, numoff
     ! ..
     ! .. Intrinsic Functions ..
     INTRINSIC abs, max
     ! ..
     numoff = 0
     ! !! MAXOFF was 1 .. now 10 to try to have more 2 x 2 pivots
     maxoff = 10
     mprint = 6
     do i = 1, n
       w(i) = zero
     end do
     do j = 1, n
       rmax = zero
       if (sym) rmax = w(j)
       idiag = 0
       ioff = 0
       i1 = ip(j)
       i2 = ip(j+1) - 1
       if (i2<i1) then
         go to 30
       else
         do ii = i1, i2
           i = ind(ii)
           rmax = rmax + abs(a(ii))
           w(i) = w(i) + abs(a(ii))
           if (i==j) idiag = ii
           if (i/=j) ioff = ii
         end do
         if (idiag==0) then
           go to 30
         else
           a(idiag) = rmax + one
           if (ioff/=0 .and. numoff<maxoff) then
             a(ioff) = 1.1*a(idiag)
             ! !! added
             a(idiag) = zero
             numoff = numoff + 1
           end if
         end if
       end if
     end do
     if ( .not. sym) then
       do j = 1, n
         i1 = ip(j)
         i2 = ip(j+1) - 1
         do ii = i1, i2
           i = ind(ii)
           if (i==j) go to 10
         end do
         go to 20
10          a(ii) = max(a(ii),w(i)+one)
20          continue
       end do
     end if
     go to 40
30      write (mprint,fmt=90000) j
     nz = -1
40      return
90000   FORMAT (' NO DIAGONAL ENTRY IN COLUMN ',I8,/, &
       ' SO CANNOT MAKE MATRIX DIAGONALLY doMINANT')
   end subroutine ndef

   subroutine pdef(n,nz,a,ip,ind,w,sym)
     ! THIS subroutine AUGMENTS THE DIAGONAL TO MAKE THE MATRIX POS DEF.
     ! .. Parameters ..
     doUBLE PRECISION zero, one
     PARAMETER (zero=0.0D0,one=1.0D0)
     ! ..
     ! .. Scalar Arguments ..
     integer n, nz
     logical sym
     ! ..
     ! .. Array Arguments ..
     doUBLE PRECISION a(nz), w(n)
     integer ind(nz), ip(*)
     ! ..
     ! .. Local Scalars ..
     doUBLE PRECISION rmax
     integer i, i1, i2, idiag, ii, j, mprint
     ! ..
     ! .. Intrinsic Functions ..
     INTRINSIC abs, max
     ! ..
     mprint = 6
     do i = 1, n
       w(i) = zero
     end do
     do j = 1, n
       rmax = zero
       if (sym) rmax = w(j)
       idiag = 0
       i1 = ip(j)
       i2 = ip(j+1) - 1
       if (i2<i1) then
         go to 30
       else
         do ii = i1, i2
           i = ind(ii)
           if (i/=j) then
             rmax = rmax + abs(a(ii))
             w(i) = w(i) + abs(a(ii))
           end if
           if (i==j) idiag = ii
         end do
         if (idiag==0) then
           go to 30
         else
           a(idiag) = rmax + one
         end if
       end if
     end do
     if ( .not. sym) then
       do j = 1, n
         i1 = ip(j)
         i2 = ip(j+1) - 1
         do ii = i1, i2
           i = ind(ii)
           if (i==j) go to 10
         end do
         go to 20
10          a(ii) = max(a(ii),w(i)+one)
20          continue
       end do
     end if
     go to 40
30      write (mprint,fmt=90000) j
     nz = -1
40      return
90000   FORMAT (' NO DIAGONAL ENTRY IN COLUMN ',I8,/, &
       ' SO CANNOT MAKE MATRIX DIAGONALLY doMINANT')
   end subroutine pdef

subroutine check_flag(actual, expect)
   integer, intent(in) :: actual
   integer, intent(in) :: expect

   if(actual.eq.expect) then
      write (*, '(a)') "ok"
   else
      nerror = nerror + 1
      write (*, '(a)') "fail"
      write (*, '(a,i3,a,i3,a)') &
         "info%flag = ", actual, " (expected ", expect, ")"
   endif
end subroutine check_flag

subroutine test_errors
   integer :: n, ne
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control
   type (nd_inform) :: info

   control%unit_error = we_unit
   control%unit_diagnostics = dl_unit

   ! --------------------------------------
   ! Test error flag n<1
   ! --------------------------------------
   n = 0
   ne = 0
   allocate (ptr(n+1), row(ne), perm(n), seps(n))
   ptr(1) = 1
   control%amd_switch1 = 2

   write (*, '(a)', advance="no") 'Testing n=0, lower entry...'
   call nd_order(0, n, ptr, row, perm, control, info, seps)
   call check_flag(info%flag, -3)

   write (*, '(a)', advance="no") 'Testing n=0, lower+upper entry...'
   call nd_order(1, n, ptr, row, perm, control, info, seps)
   call check_flag(info%flag, -3)
end subroutine test_errors

! *****************************************************************

  subroutine test_input(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i, j
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   logical :: corr
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info


   ok = .true.

   test_count = 0
   ! --------------------------------------
   ! Test diagonal matrix
   ! --------------------------------------
   testi_count = 0
   write (6,'(a1)') ' '
   write (6,'(a)') '**** Test inputs: diagonal matrix ****'
   n = 10
   ne = 0
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = 1

   control%amd_switch1 = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     corr = .true.
     do j = 1, n
       corr = (corr .and. (perm(j)==j))
     end do
     if (info%flag>=0 .and. corr) testi_count = testi_count + 1
     call nd_order(1,n,ptr,row,perm,control,info,seps)
     corr = .true.
     do j = 1, n
       corr = (corr .and. (perm(j)==j))
     end do
     if (info%flag>=0 .and. corr) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20

   if (testi_count/=2) then
     write (6,'(a29)') 'Code failure in test section '
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test expansion of matrix from lower triangle input
   ! --------------------------------------
   test = 1
   write (6,'(a1)') ' '
   write (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
   n = 10
   ne = 9
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n-1) = (/ (i,i=1,n-1) /)
   ptr(n:n+1) = n
   row(1:ne) = (/ (i,i=2,n) /)

   control%amd_switch1 = 2
   control%amd_call = 2
   do i = 0, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==18) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test expansion of matrix from lower triangle format and that diags are
   ! removed
   ! --------------------------------------

   test = 2
   write (6,'(a1)') ' '
   write (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n) = (/ (2*i-1,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
     10 /)

   control%amd_switch1 = 2
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==18) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test expansion of matrix from lower and upper triangle input with no
   ! diags
   ! --------------------------------------
   test = 3
   write (6,'(a1)') ' '
   write (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
   n = 10
   ne = 18
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1) = 1
   ptr(2:n) = (/ (2*(i-1),i=2,n) /)
   ptr(n+1) = ne + 1
   row(1) = 2
   do i = 2, n - 1
     row(ptr(i)) = i - 1
     row(ptr(i)+1) = i + 1
   end do
   row(ptr(n)) = n - 1

   control%amd_switch1 = 2
   control%amd_call = 2
   do i = 2, 2
     control%print_level = i
     call nd_order(1,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==18) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test expansion of matrix from lower and upper triangle input with diag
   ! entry
   ! --------------------------------------
   test = 4
   write (6,'(a1)') ' '
   write (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

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

   control%amd_switch1 = 2
   control%amd_call = 2
   do i = 0, 2
     control%print_level = i
     call nd_order(1,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==18) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if





   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

  end subroutine test_input


! *****************************************************************


  subroutine test_dense(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0
   ! --------------------------------------
   ! Test dense row removal - row 1 dense
   ! --------------------------------------
   test = 1
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
   n = 800
   ne = n
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1) = 1
   ptr(2:n+1) = ne + 1
   row(1:ne) = (/ (i,i=1,n) /)

   control%amd_switch1 = 2
   control%amd_call = 20
   do i = 0, 2
     control%print_level = i
     control%remove_dense_rows = .true.
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nsuper==799 .and. info%nzsuper==0) &
       testi_count = testi_count + 1
     control%remove_dense_rows = .false.
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nsuper==800 .and. info%nzsuper==1598) &
       testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=6) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test dense row removal - first row dense
   ! --------------------------------------
   test = 2
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
   n = 800
   ne = 2*n - 3
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1) = 1
   ptr(2:n) = (/ (n+i-2,i=2,n) /)
   ptr(n+1) = ne + 1
   row(1:ptr(2)-1) = (/ (i,i=2,n) /)
   row(ptr(2):ptr(n)-1) = (/ (i+1,i=2,n-1) /)

   control%amd_switch1 = 2
   control%amd_call = 20
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nsuper==799 .and. info%nzsuper==1596) &
       testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test dense row removal - last row dense
   ! --------------------------------------
   test = 3
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
   n = 800
   ne = n - 1
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n-1) = (/ (i,i=1,n-1) /)
   ptr(n:n+1) = ne + 1
   row(1:ne) = n

   control%amd_switch1 = 2
   control%amd_call = 20
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nsuper==799 .and. info%nzsuper==0) &
       testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test dense row removal - first two rows dense but row 2 has max degree
   ! --------------------------------------
   test = 4
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
   n = 800
   ne = n - 3 + n - 2
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1) = 1
   ptr(2) = n - 2
   ptr(3:n+1) = ne + 1
   row(ptr(1):ptr(2)-1) = (/ (i+3,i=1,n-3) /)
   row(ptr(2):ne) = (/ (i+2,i=1,n-2) /)

   control%amd_switch1 = 2
   control%amd_call = 20
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nsuper==798 .and. info%nzsuper==0) &
       testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_dense


! *****************************************************************


  subroutine test_ashcraft(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0


   ! --------------------------------------
   ! Test Ashcraft method
   ! --------------------------------------

   test = 1
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
   n = 10
   ne = 14
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
   row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   do i = 0, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test Ashcraft method
   ! --------------------------------------

   test = 2
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
   n = 10
   ne = 15
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 4, 6, 8, 9, 10, 11, 13, 15, 16, 16 /)
   row(1:ne) = (/ 2, 3, 4, 3, 5, 4, 8, 6, 7, 9, 8, 10, 9, 10, 10 /)

   control%amd_switch1 = 2
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test Ashcraft method
   ! --------------------------------------

   test = 3
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
   n = 13
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 4, 6, 8, 10, 12, 13, 13, 14, 14, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 4, 11, 5, 6, 7, 8, 9, 10, 8, 10, 12, 13, 13 /)

   control%amd_switch1 = 2
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test Ashcraft method
   ! --------------------------------------

   test = 4
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
   n = 13
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 4, 6, 8, 10, 12, 13, 13, 14, 14, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 4, 11, 5, 6, 7, 8, 9, 10, 8, 10, 12, 13, 13 /)

   control%amd_switch1 = 2
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%balance = 20.0
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_ashcraft


! *****************************************************************


  subroutine test_levelset(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0


   ! --------------------------------------
   ! Test one-sided levelset method
   ! --------------------------------------

   test = 1
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Level set test ', test, ' ***'
   n = 10
   ne = 14
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
   row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 2
   control%amd_call = 2
   control%balance = 20.0
   do i = 0, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test one-sided levelset method
   ! --------------------------------------

   test = 2
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Level set test ', test, ' ***'
   n = 10
   ne = 14
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
   row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 2
   control%balance = 1.5
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_levelset

! *****************************************************************


  subroutine test_shift(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0

  ! --------------------------------------
   ! Test DM-style refinement with 2-sided partition
   ! --------------------------------------

   test = 1
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 8
   ne = 11
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 4, 5, 6, 8, 9, 11, 12, 12 /)
   row(1:ne) = (/ 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 8 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   do i = 2, 2
     control%print_level = i
     control%refinement = 2
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
     control%refinement = 5
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
     control%print_level = i
     control%refinement = 2
     call nd_order(0,n,ptr,row,perm,control,info)
     if (info%flag>=0) testi_count = testi_count + 1
     control%refinement = 5
     call nd_order(0,n,ptr,row,perm,control,info)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=4) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition
   ! --------------------------------------

   test = 2
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 21
   ne = 37
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 13, 15, 16, 20, 22, 24, 26, 28, 29, &
     29, 32, 33, 35, 37, 38, 38 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 16, &
     17, 11, 13, 12, 13, 13, 15, 14, 15, 15, 17, 18, 19, 18, 19, 20, 20, &
     21, 21 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%refinement = 2
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test DM refinement with 1-sided partition
   ! --------------------------------------

   test = 3
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
     27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
     13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
     24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 2
   control%refinement = 2
   control%amd_call = 2
   control%balance = 33.0
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition
   ! --------------------------------------

   test = 4
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 7
   ne = 6
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 6, 7, 7 /)
   row(1:ne) = (/ 7, 7, 7, 7, 7, 7 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 2
   control%refinement = 2
   control%amd_call = 2
   control%balance = 12.0
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test DM refinement with 1-sided partition
   ! --------------------------------------

   test = 5
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 7
   ne = 6
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 6, 7, 7 /)
   row(1:ne) = (/ 7, 7, 7, 7, 7, 7 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 2
   control%refinement = 5
   control%amd_call = 2
   control%balance = 12.0
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test DM refinement with 1-sided partition with balanced partition
   ! --------------------------------------

   test = 6
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
     27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
     13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
     24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 2
   control%refinement = 5
   control%balance = 2.0
   control%amd_call = 2
   do i = 0, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------

   test = 7
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
     27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
     13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
     24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%refinement = 5
   control%balance = 2.0
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------

   test = 8
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%refinement = 2
   control%balance = 1.05
   control%refinement_band = n
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------

   test = 9
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 13
   ne = 15
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 3, 6, 7, 8, 9, 11, 11, 13, 14, 15, 16, 16 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%refinement = 2
   control%balance = 1.30
   control%refinement_band = n
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test DM refinement with 2-sided partition with balanced partition
   ! --------------------------------------

   test = 10
   write (6,'(a1)') ' '
   write (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
   n = 13
   ne = 15
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 16 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 7, 7, 8, 9, 10, 11, 11, 11, 12, 13 /)

   control%amd_switch1 = 5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%refinement = 2
   control%balance = 1.30
   control%refinement_band = n
   control%amd_call = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if





   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_shift

! *****************************************************************


  subroutine test_multigrid(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i, j, k
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0

   ! --------------------------------------
   ! Test multigrid, 2-sided partitioning
   ! --------------------------------------

   test = 1
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 5
   control%balance = 1.05
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 2
   control%amd_call = 2
   control%stop_coarsening1 = 5
   do i = 0, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test multigrid, 1-sided partitioning
   ! --------------------------------------

   test = 2
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 5
   control%balance = 1.05
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 5
   do i = 0, 2
     do j = 0, 2
       control%print_level = i
       control%matching = j
       call nd_order(0,n,ptr,row,perm,control,info,seps)
       if (info%flag>=0) testi_count = testi_count + 1
     end do
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=9) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test multigrid coarsening
   ! --------------------------------------

   test = 3
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 5
   control%balance = 1.05
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 2
   control%amd_call = 2
   control%stop_coarsening1 = 5
   do i = 0, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test multigrid coarsening
   ! --------------------------------------

   test = 4
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 5
   control%balance = 1.05
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 5
   do i = 0, 2
     do j = 0, 2
       control%print_level = i
       control%matching = j
       call nd_order(0,n,ptr,row,perm,control,info,seps)
       if (info%flag>=0) testi_count = testi_count + 1
     end do
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=9) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test multigrid coarsening
   ! --------------------------------------

   test = 5
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 12
   ne = 18
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 4, 6, 7, 9, 10, 12, 13, 16, 18, 19, 19 /)
   row(1:ne) = (/ 2, 9, 9, 4, 10, 10, 6, 11, 11, 8, 12, 12, 10, 11, 12, 11, &
     12, 12 /)

   control%amd_switch1 = 5
   control%balance = 1.05
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 5
   control%min_reduction = 0.4
   control%matching = 1
   do i = 0, 2
     control%print_level = i
     control%matching = j
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=3) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test multigrid refinement
   ! --------------------------------------

   test = 6
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 5
   control%balance = 1.05
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 5
   do i = 0, 7
     do j = 0, 2
       control%refinement = i
       control%matching = j
       call nd_order(0,n,ptr,row,perm,control,info,seps)
       if (info%flag>=0) testi_count = testi_count + 1
     end do
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=24) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test no. multigrid - automatic choice refinement
   ! --------------------------------------

   test = 7
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   control%amd_switch1 = 4
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 4
   control%print_level = 2
   do i = 0, 7
     do j = 0, 2
       control%refinement = 3
       control%matching = j
       call nd_order(0,n,ptr,row,perm,control,info,seps)
       if (info%flag>=0) testi_count = testi_count + 1
     end do
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=24) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test multigrid - automatic choice refinement
   ! --------------------------------------

   test = 8
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   control%amd_switch1 = 4
   control%balance = 1.0
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 2
   do i = 0, 7
     do j = 0, 2
       control%refinement = i
       control%matching = j
       call nd_order(0,n,ptr,row,perm,control,info,seps)
       if (info%flag>=0) testi_count = testi_count + 1
     end do
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=24) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if




   ! --------------------------------------
   ! Test automatic multilevel-no multilevel check
   ! --------------------------------------

   test = 9
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 200
   ne = 199
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10


   ptr(1:n) = (/ (i,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ (i,i=2,n) /)

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%balance = 8.0
   control%partition_method = 2
   control%coarse_partition_method = 2
   control%ml_call = 10
   control%amd_call = 2
   control%max_improve_cycles = 2
   control%refinement = 7
   do i = 0, 0
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test automatic multilevel-no multilevel check
   ! --------------------------------------

   test = 10
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 200
   ne = 7*(n-7)+21
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   j = 1
   do i = 1,n
    ptr(i) = j
    do k = i+1, min(i+7,n)
      row(j) = k
      j = j + 1
    end do
   end do
   ptr(n+1) = j

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%balance = 8.0
   control%partition_method = 2
   control%ml_call = 10
   control%amd_call = 2
   control%max_improve_cycles = 2
   control%refinement = 7
   do i = 1, 2
     control%coarse_partition_method = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=2) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test multigrid reduction ratio
   ! --------------------------------------

   test = 11
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 20
   ne = n-1
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10


   ptr(1) = 1

   ptr(2:n+1) = ne + 1

   row(1:ne) = (/ (i,i=2,n) /)

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%balance = 2.0
   control%partition_method = 1
   control%ml_call = 10
   control%amd_call = 2
   control%max_improve_cycles = 2
   control%refinement = 7
   control%find_supervariables = .false.
   control%print_level = 2
   do i = 1, 2
     control%coarse_partition_method = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=2) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test for full galerkin_graph_rap coverage
   ! --------------------------------------

   test = 12
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
   n = 16
   ne = 42
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

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

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%min_reduction = 0.5
   control%balance = 1.0
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 1
   control%print_level = 0
   do i = 0, 0
     control%refinement = 1
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_multigrid

! *****************************************************************


  subroutine test_amd(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0

   ! --------------------------------------
   ! Call AMD with no nested
   ! --------------------------------------

   test = 1
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****AMD test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
     27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
     13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
     24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)


   control%amd_call = 40
   do i = 2, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Call AMD with no nested
   ! --------------------------------------

   test = 2
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****AMD test ', test, '****'
   n = 31
   ne = 1
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1) = 1
   ptr(2:n+1) = 2
   row(1:ne) = (/ 2 /)

   control%amd_call = 40
   do i = 2, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_amd

! *****************************************************************


  subroutine test_refine(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i, j
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0

   ! --------------------------------------
   ! Test for refinement controls
   ! --------------------------------------

   test = 1
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   control%amd_switch1 = 6
   control%amd_call = 2
   control%find_supervariables = .false.
   control%partition_method = 0
   control%refinement = 0
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info)
     if (info%flag>=0 .and. info%nzsuper==32) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test for refinement controls
   ! --------------------------------------

   test = 2
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   control%amd_switch1 = 6
   control%amd_call = 2
   control%find_supervariables = .false.
   control%partition_method = 0
   control%refinement = 0
   do i = 0, 0
     control%print_level = i
     do j = 1, 7
       control%refinement = j
       call nd_order(0,n,ptr,row,perm,control,info)
       if (info%flag>=0 .and. info%nzsuper==32) testi_count = testi_count + &
         1
     end do
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=7) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if





   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------

   test = 3
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 4
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 2
   control%amd_call = 2
   control%max_improve_cycles = 2
   control%refinement_band = -1
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------

   test = 4
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   control%amd_switch1 = 3
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------

   test = 5
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 13
   ne = 27
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
     7, 7, 8, 9, 10, 11, 12, 13 /)

   control%amd_switch1 = 3
   control%stop_coarsening1 = 2
   control%min_reduction = 0.5
   control%balance = 1.0
   control%partition_method = 1
   control%ml_call = 4
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 2
   control%print_level = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test for full refine_trim coverage
   ! --------------------------------------

   test = 6
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 16
   ne = 22
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 6, 7, 8, 9, 12, 14, 16, 17, 18, 20, 21, 22, 23, &
     23 /)
   row(1:ne) = (/ 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13, 13, 14, &
     14, 15, 15, 16, 16 /)

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%min_reduction = 0.5
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 1
   control%print_level = 0
   do i = 0, 0
     control%refinement = 1
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test for full refine_trim coverage
   ! --------------------------------------

   test = 7
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 16
   ne = 22
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 14, 15, 16, 18, 19, 20, 22, &
     23, 23 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 10, 11, 12, 13, &
     14, 15, 15, 16, 16 /)

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%min_reduction = 0.5
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 1
   control%print_level = 2
   do i = 0, 0
     control%refinement = 1
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test for full refine_trim coverage
   ! --------------------------------------

   test = 8
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 25
   ne = 44
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 5, 8, 10, 12, 14, 16, 18, 20, 21, 22, 23, 29, 30, &
     30, 32, 34, 36, 38, 39, 41, 43, 44, 45, 45 /)
   row(1:ne) = (/ 2, 3, 4, 5, 4, 6, 7, 5, 8, 9, 10, 7, 13, 8, 13, 9, 13, &
     10, 13, 13, 12, 13, 14, 16, 17, 18, 19, 20, 15, 17, 21, 18, 21, 19, &
     22, 20, 23, 23, 22, 24, 23, 24, 24, 25 /)

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%min_reduction = 0.5
   control%balance = 1.5
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 1
   control%print_level = 2
   do i = 0, 0
     control%refinement = 1
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test for not checking balance
   ! --------------------------------------

   test = 9
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 25
   ne = 44
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 5, 8, 10, 12, 14, 16, 18, 20, 21, 22, 23, 29, 30, &
     30, 32, 34, 36, 38, 39, 41, 43, 44, 45, 45 /)
   row(1:ne) = (/ 2, 3, 4, 5, 4, 6, 7, 5, 8, 9, 10, 7, 13, 8, 13, 9, 13, &
     10, 13, 13, 12, 13, 14, 16, 17, 18, 19, 20, 15, 17, 21, 18, 21, 19, &
     22, 20, 23, 23, 22, 24, 23, 24, 24, 25 /)

   control%amd_switch1 = 4
   control%stop_coarsening1 = 3
   control%min_reduction = 0.5
   control%balance = real(n+1)
   control%partition_method = 1
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 1
   control%print_level = 0
   do i = 0, 0
     control%refinement = 7
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------

   test = 10
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   control%amd_switch1 = 4
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 2
   control%max_improve_cycles = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------

   test = 11
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 9
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
   row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
     9 /)

   control%amd_switch1 = 4
   control%balance = 10.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%stop_coarsening1 = 2
   control%max_improve_cycles = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------

   test = 12
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 7
   ne = 12
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
   row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

   control%amd_switch1 = 3
   control%balance = 8.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test expand-refine cycles
   ! --------------------------------------

   test = 13
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 31
   ne = 49
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
     27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
   row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
     14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
     24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

   control%amd_switch1 = 4
   control%balance = 8.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------

   test = 14
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 29
   ne = 68
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
         51, 55,59,62,64,66,67,67,68,69,69,69 /)
   row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
            14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
            19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29,&
            29,29,29,29 /)

   control%amd_switch1 = 4
   control%balance = 1.01
   control%partition_method = 1
   control%coarse_partition_method = 2
   control%refinement_band = 0
   control%stop_coarsening1 = 3
   control%stop_coarsening2 = 3
   control%amd_call = 2
   control%max_improve_cycles = 2
   control%print_level = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test automatic choice of shift within expand-refine cycle
   ! --------------------------------------

   test = 15
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
   n = 21
   ne = 34
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                    35,35 /)
   row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
        16,17,16,18,17,18,18,19,20,20,21,21 /)

   control%amd_switch1 = 4
   control%balance = 1.01
   control%partition_method = 1
   control%coarse_partition_method = 2
   control%refinement_band = 0
   control%stop_coarsening1 = 3
   control%stop_coarsening2 = 3
   control%amd_call = 2
   control%max_improve_cycles = 2
   control%print_level = 2
   do i = 0, 7
     control%refinement = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=8) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_refine




! *****************************************************************


  subroutine test_misc(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i, j, k
   integer, allocatable, dimension(:) :: ptr, row, perm, seps
   type (nd_options) :: control, control_orig
   type (nd_inform) :: info

   ok = .true.
   testi_count = 0
   test_count = 0
   ! --------------------------------------
   ! Test diagonal submatrix - not diag!!!!
   ! --------------------------------------
   test = 1
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 10
   ne = 10
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n-2) = (/ (i,i=1,n-2) /)
   ptr(n-1) = ptr(n-2) + 2
   ptr(n:n+1) = ne + 1
   row(1:ne-2) = n
   row(ne-1) = n - 1
   row(ne) = n

   control%amd_switch1 = 2
   control%amd_call = 2
   do i = 1, 1
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==16) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if




   ! --------------------------------------
   ! Test independent component detection
   ! --------------------------------------
   test = 2
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 5
   ne = 6
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 3, 4, 6, 7 /)
   row(1:ne) = (/ 2, 1, 4, 3, 5, 4 /)

   control%amd_switch1 = 2
   control%amd_call = 2
   do i = 0, 2
     control%print_level = i
     control%coarse_partition_method = 2
     call nd_order(1,n,ptr,row,perm,control,info)
     if (info%flag>=0) testi_count = testi_count + 1

     control%coarse_partition_method = 1
     call nd_order(1,n,ptr,row,perm,control,info)
     if (info%flag>=0) testi_count = testi_count + 1

   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=6) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if






   ! --------------------------------------
   ! Test for supervariable detection turned off
   ! --------------------------------------

   test = 3
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n) = (/ (2*i-1,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
     10 /)

   control%amd_switch1 = 2
   control%amd_call = 2
   control%find_supervariables = .false.
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==18) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test for full matrix
   ! --------------------------------------

   test = 4
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 10
   ne = 45
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

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

   control%amd_switch1 = 2
   control%amd_call = 2
   control%find_supervariables = .false.
   do i = 2, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==90) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test cost function cost2
   ! --------------------------------------

   test = 5
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 10
   ne = 19
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n) = (/ (2*i-1,i=1,n) /)
   ptr(n+1) = ne + 1
   row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
     10 /)

   control%amd_switch1 = 2
   control%amd_call = 2
   control%find_supervariables = .false.
   control%cost_function = 2
   do i = 0, 0
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0 .and. info%nzsuper==18) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test for no separator return, partition 2 larger than amd_switch1 but
   ! partition
   ! 1 equal to amd_switch1
   ! --------------------------------------

   test = 6
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 14
   ne = 16
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
   row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

   control%amd_switch1 = 6
   control%amd_call = 2
   control%find_supervariables = .false.
   control%partition_method = 0
   do i = 2, 2
     control%print_level = i
     call nd_order(0,n,ptr,row,perm,control,info)
     if (info%flag>=0 .and. info%nzsuper==32) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if






   ! --------------------------------------
   ! Test matrix extraction involving end row
   ! --------------------------------------

   test = 7
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 14
   ne = 26
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 6, 8, 10, 13, 17, 19, 23, 24, 26, 27, 27, 27, &
     27 /)
   row(1:ne) = (/ 2, 3, 3, 4, 5, 5, 14, 5, 13, 6, 13, 14, 8, 9, 13, 14, 10, &
     13, 9, 10, 11, 13, 11, 11, 12, 12 /)

   control%amd_switch1 = 7
   control%stop_coarsening1 = 2
   control%min_reduction = 0.5
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 0
   control%print_level = 2
   do i = 0, 0
     control%refinement = 7
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if



   ! --------------------------------------
   ! Test matrix extraction involving end row
   ! --------------------------------------

   test = 8
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 16
   ne = 31
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1, 3, 6, 8, 10, 13, 14, 16, 21, 25, 27, 29, 31, 31, 31, &
     32, 32 /)
   row(1:ne) = (/ 2, 3, 3, 4, 5, 5, 6, 5, 15, 6, 15, 16, 16, 11, 15, 9, 11, &
     12, 15, 16, 10, 12, 13, 16, 16, 13, 12, 14, 13, 14, 16 /)

   control%amd_switch1 = 7
   control%stop_coarsening1 = 2
   control%min_reduction = 0.5
   control%balance = 1.0
   control%partition_method = 0
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 0
   control%print_level = 2
   do i = 0, 0
     control%refinement = 7
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test no balanced partition found so switch to multilevel
   ! --------------------------------------

   test = 9
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****Misc test ', test, '****'
   n = 18
   ne = 41
   allocate (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n+1) = (/ 1,2,3,8,10,12,15,19,21,24,28,30,33,37,39,40,42,42,42 /)
   row(1:ne) = (/ 2,3,4,5,6,7,8,6,7,7,8,7,9,10,8,9,10,11,10,11,10,12,13,11,&
        12,13,14,13,14,13,15,16,14,15,16,17,16,17,16,17,18 /)

   control%amd_switch1 = 7
   control%stop_coarsening1 = 2
   control%min_reduction = 0.5
   control%balance = 1.0
   control%partition_method = 2
   control%coarse_partition_method = 1
   control%amd_call = 2
   control%max_improve_cycles = 0
   control%print_level = 2
   do i = 0, 0
     control%refinement = 7
     call nd_order(0,n,ptr,row,perm,control,info,seps)
     if (info%flag>=0) testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,perm,seps,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_misc




! *****************************************************************


  subroutine test_fm(ok)
   logical, intent(out) :: ok

   integer :: test_count, testi_count, test
   integer :: n, ne, st, i, an1,an2,swgt,aw1,aw2,aws
   integer, allocatable, dimension(:) :: ptr,row,wgt,work,part
   type (nd_options) :: control, control_orig

   ok = .true.
   testi_count = 0
   test_count = 0

   ! --------------------------------------
   ! Test fm with empty partition 2
   ! --------------------------------------
   test = 1
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****FM test ', test, '****'
   n = 11
   ne = 14
   swgt = n
   allocate (ptr(n),row(ne),wgt(n),work(8*n+swgt),part(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n) = (/ 1,2,4,6,8,9,11,12,13,14,15  /)
   row(1:ne) = (/ 2,3,4,5,6,6,7,8,8,9,9,10,10,11  /)
   wgt(:)=1
   an1 = 3
   an2 = 0
   aw1 = an1
   aw2 = an2
   aws = swgt - aw1-aw2
   part(1:n) = (/ 1,2,4,8,10,11,3,5,6,7,9 /)

   do i = 1, 1
     control%print_level = i
     call nd_refine_fm(n,ne,ptr,row,wgt,swgt,an1,an2,aw1,aw2,aws,part,work,&
       control)
     testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,wgt,work,part,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if

   ! --------------------------------------
   ! Test fm with empty partition 1
   ! --------------------------------------
   test = 2
   write (6,'(a1)') ' '
   write (6,'(a20,i5,a4)') '****FM test ', test, '****'
   n = 11
   ne = 14
   swgt = n
   allocate (ptr(n),row(ne),wgt(n),work(8*n+swgt),part(n),STAT=st)
   if (st/=0) go to 10

   ptr(1:n) = (/ 1,2,4,6,8,9,11,12,13,14,15  /)
   row(1:ne) = (/ 2,3,4,5,6,6,7,8,8,9,9,10,10,11  /)
   wgt(:)=1
   an1 = 0
   an2 = 3
   aw1 = an1
   aw2 = an2
   aws = swgt - aw1-aw2
   part(1:n) = (/ 1,2,4,8,10,11,3,5,6,7,9 /)

   do i = 1, 1
     control%print_level = i
     call nd_refine_fm(n,ne,ptr,row,wgt,swgt,an1,an2,aw1,aw2,aws,part,work,&
       control)
     testi_count = testi_count + 1
   end do
   call reset_control(control_orig,control)

   deallocate (ptr,row,wgt,work,part,STAT=st)
   if (st/=0) go to 20
   if (testi_count/=1) then
     write (6,'(a29,i5)') 'Code failure in test section ', test
     ok = .false.
     return
   else
     test_count = test_count + 1
     testi_count = 0
   end if


   return

10    write (*,'(a,i4)') 'Allocation error during test section ', test
   return

20    write (*,'(a,i4)') 'Deallocation error during test section ', test

 end subroutine test_fm
end program main
