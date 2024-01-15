! Copyright (c) 2010, 2012 David Fong, Michael Saunders, Stanford University
! Copyright (c) 2014, 2015, 2016 Science and Technology Facilites Council (STFC)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! Test code. Originally developed by Michael Saunders <saunders@stanford.edu>.
! Extended and redesigned to provide more comprehensive testing
! (as required by SPRAL) and for reverse communication interface.
!
program spral_lsmr_test
  use spral_lsmr
  use spral_random
  use spral_random_matrix, only : random_matrix_generate
  implicit none

  integer(4),  parameter :: ip = kind( 1 )
  integer(4),  parameter :: wp = kind( 1.0d+0 )
  real(wp),    parameter :: zero = 0.0_wp, one = 1.0_wp, two = 2.0_wp
  real(wp),    parameter :: power = 0.5_wp
  real(wp),    parameter :: tol2 = 1.0d-7
  real(wp),    parameter :: tol1 = 1.0d-9

  ! termination flags
   integer(ip), parameter :: lsmr_stop_x0              = 0
   integer(ip), parameter :: lsmr_stop_compatible      = 1
   integer(ip), parameter :: lsmr_stop_LS_atol         = 2
   integer(ip), parameter :: lsmr_stop_ill             = 3
   integer(ip), parameter :: lsmr_stop_Ax              = 4
   integer(ip), parameter :: lsmr_stop_LS              = 5
   integer(ip), parameter :: lmsr_stop_condA           = 6
   integer(ip), parameter :: lsmr_stop_itnlim          = 7
   !integer(ip), parameter :: lsmr_stop_allocation      = 8
   !integer(ip), parameter :: lsmr_stop_deallocation    = 9
   integer(ip), parameter :: lsmr_stop_m_oor           = 10

  !---------------------------------------------------------------------
  ! Michael's tests:
  ! This program calls lsmrtest(...) to generate a series of test problems
  ! Ax = b or Ax ~= b and solve them with LSMR.
  ! The matrix A is m x n.  It is defined by routines in lsmrReverse_TestModule.
  !
  ! 23 Sep 2007: First version of lsmrTestProgram.f90.
  ! 24 Oct 2007: Use real(wp) instead of compiler option -r8.
  ! 23 Nov 2015: Reverse communication variant
  !---------------------------------------------------------------------

  real(wp), allocatable :: d(:), hy(:), hz(:) ! These define A = Y D Z.
  real(wp), allocatable :: wm(:), wn(:)       ! Work vectors.


  ! Local variables
  integer(ip)  :: flag
  integer(ip)  :: m,n,nbar,ndamp,nduplc,nfail,npower
  integer(ip)  :: localSize  ! No. of vectors involved in local
                             ! reorthogonalization (>= 0)
  integer(ip) :: prob
  real(wp)    :: damp

  integer(ip) :: nout, mout, pout



   ! Switch off printing by setting nout, mout and pout to <0
   ! If these streams are positive, they must all be different.
   nout = 6
   mout = -10
   pout = -11

  nfail = 0
  prob = 0

  ! this file is for tests with printing on
  if (pout.gt.0) open(pout,file='LSMR_print_output',status="replace")

  ! the output in LSMR_Mtests.output can be directly compared with the output
  ! from Michael's code.
  if (mout.gt.0) open(mout,file='LSMR_Mtests.output',status='replace')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (nout.gt.0) then
      write(nout,"(a)") "======================"
      write(nout,"(a)") "Running Michaels tests:"
      write(nout,"(a)") "======================"
    end if

  nbar   = 100  ! 1000
  nduplc =   4  ! 40


  m = 2*nbar        ! Over-determined systems
  n = nbar
  localSize = 0

  do ndamp = 2,6
     npower = ndamp
     damp   = 0.0
     if (ndamp > 2) damp   = 10.0**(-ndamp)
     call lsmrtest(m,n,nduplc,npower,damp,localSize,flag,mout)
     if (flag.lt.0) then
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' On exit from lsmrtest flag = ', flag
        nfail = nfail + 1
     else if (flag.eq.1) then
        prob = prob + 1
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' LSMR  appears to be successful for problem ',prob
     end if
  end do

  localsize = 10    ! Repeat last test with local reorthogonalization
  call lsmrtest(m,n,nduplc,npower,damp,localSize,flag,mout)
     if (flag.lt.0) then
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' On exit from lsmrtest flag = ', flag
        nfail = nfail + 1
     else if (flag.eq.1) then
        prob = prob + 1
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' LSMR  appears to be successful for problem ',prob
     end if


  m = nbar          ! Square systems
  n = nbar
  localSize = 0

  do ndamp = 2,6
     npower = ndamp
     damp   = 0.0
     if (ndamp > 2) damp   = 10.0**(-ndamp-6)
     call lsmrtest(m,n,nduplc,npower,damp,localSize,flag,mout)
     if (flag.lt.0) then
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' On exit from lsmrtest flag = ', flag
        nfail = nfail + 1
     else if (flag.eq.1) then
        prob = prob + 1
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' LSMR  appears to be successful for problem ',prob
     end if
  end do

  localsize = 10    ! Repeat last test with local reorthogonalization
  call lsmrtest(m,n,nduplc,npower,damp,localSize,flag,mout)
     if (flag.lt.0) then
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' On exit from lsmrtest flag = ', flag
        nfail = nfail + 1
     else if (flag.eq.1) then
        prob = prob + 1
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' LSMR  appears to be successful for problem ',prob
     end if


  m = nbar          ! Under-determined systems
  n = 2*nbar
  localSize = 0

  do ndamp = 2,6
     npower = ndamp
     damp   = 0.0
     if (ndamp > 2) damp   = 10.0**(-ndamp-6)
     call lsmrtest(m,n,nduplc,npower,damp,localSize,flag,mout)
     if (flag.lt.0) then
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' On exit from lsmrtest flag = ', flag
        nfail = nfail + 1
     else if (flag.eq.1) then
        prob = prob + 1
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' LSMR  appears to be successful for problem ',prob
     end if
  end do

  localsize = 10    ! Repeat last test with local reorthogonalization
  call lsmrtest(m,n,nduplc,npower,damp,localSize,flag,mout)
     if (flag.lt.0) then
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' On exit from lsmrtest flag = ', flag
        nfail = nfail + 1
     else if (flag.eq.1) then
        prob = prob + 1
        if (nout.gt.0) write (nout,'(a,i3)') &
        ' LSMR  appears to be successful for problem ',prob
     end if

     if (mout.gt.0) write (mout,'(//a,i5)') &
         " At end of Michael's tests number of failures = ",nfail

     if (nout.gt.0) write (nout,'(//a,i5)') &
         " At end of Michael's tests number of failures = ",nfail


  if (mout.gt.0) close(mout)

  !!! End Michael's tests !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  call test_special

  call test_random

  if (pout.gt.0) close(pout)
  if (nfail.gt.0) stop -1

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine test_special

   integer(ip), parameter :: maxm =  40 ! largest problem size (rows)
   integer(ip), parameter :: maxn =  30 ! largest problem size (cols)
   integer(ip), parameter :: nprob = 3  ! number of error tests
   integer(ip), parameter :: mprob = 10 ! number of other tests
   type (random_state) :: state

   integer(ip), allocatable :: ptr(:),row(:)
   real(wp),    allocatable :: val(:), work(:)
   real(wp)    :: a(maxm,maxn), b(maxm), x(maxn), y(maxn), rhs(maxm,1)
   real(wp)    :: u(maxm), v(maxn) ! reverse communication arrays
   real(wp)    :: w(maxn) , z(maxm)

    type ( lsmr_keep )    :: keep
    type ( lsmr_options ) :: options
    type ( lsmr_inform )  :: inform

   integer(ip)  :: action
   integer(ip) :: flag
   integer(ip) :: i,j
   integer(ip) :: lwork
   integer(ip) :: m1,n,n1,ne
   integer(ip) :: prblm
   integer(ip) :: stat

   logical :: ldamp

!  real(wp) :: normAPr
   real(wp) :: normAr, normr, normx, norm_rhs, normb, normAb, ratio
   real(wp) :: dnrm2
   real(wp) :: lambda

    if (nout.gt.0) then
      write(nout,"(a)") "======================"
      write(nout,"(a)") "Testing bad arguments:"
      write(nout,"(a)") "======================"
    end if


 main: do prblm = 1,nprob+mprob

    if (nout.gt.0) then
      if (prblm.eq.1) then
         write(nout,"(a)") "======================"
         write(nout,"(a)") "Testing bad arguments:"
         write(nout,"(a)") "======================"

      else if (prblm.eq.nprob+1) then
         write(nout,"(a)") "======================================"
         write(nout,"(a)") "Testing some options with printing on:"
         write(nout,"(a)") "======================================"
     end if
   end if

         call matrix_gen(state, maxm, m, maxn, n, ne, ptr, row, val, b, a, &
            flag)

         if (flag /= 0) then
           if (nout.gt.0) write(nout,'(a,i3)') &
         ' Unexpected error when generating matrix. flag=',flag
           nfail = nfail + 1;   cycle main
         end if

         ! run with printing on (thus some tests are not for errors but
         ! just to get the print statements executed)
         call set_controls(options)
         options%unit_error = pout
         options%unit_diagnostics = pout

         action = 0
         flag = 0
         u(1:m) = b(1:m)
         lambda = 1.0d-08
         ldamp = .false.

         if (prblm.eq.1) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing m=0....................."
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing m=0................."
            m1 = 0
            options%ctest = 1 ! use this to get some printing tested
            call lsmr_solve(action, m1, n, u, v, y, keep, options, inform, &
              damp = lambda)
            call print_result(nout, nfail, inform%flag, lsmr_stop_m_oor)
            cycle
         end if

         if (prblm.eq.2) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing n=0....................."
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing n=0................."
            n1 = 0
            options%ctest = 1 ! use this to get some printing tested
            call lsmr_solve(action, m, n1, u, v, y, keep, options, inform)
            call print_result(nout, nfail, inform%flag, lsmr_stop_m_oor)
            cycle
         end if

         if (prblm.eq.3) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing itnlim too small........"
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing itnlim too small...."
            options%itnlim = 1
         end if

         if (prblm.eq.nprob+1) then
            if (nout.gt.0) &
            write(nout,"(a)") " * Testing options%ctest=1........"
            if (pout.gt.0) &
            write(pout,"(/a)") " * Testing options%ctest=1....."
            options%ctest = 1
            options%print_freq_itn = 1
            options%print_freq_head = 2
         end if

         if (prblm.eq.nprob+2) then
            if (nout.gt.0) &
            write(nout,"(a)") " * Testing options%ctest=2........"
            if (pout.gt.0) &
            write(pout,"(/a)") " * Testing options%ctest=2....."
            options%ctest = 2
            options%print_freq_itn = 1
            options%print_freq_head = 1
         end if

         if (prblm.eq.nprob+3) then
            if (nout.gt.0) &
            write(nout,"(a)") " * Testing options%ctest=2........"
            if (pout.gt.0) &
            write(pout,"(/a)") " * Testing options%ctest=2....."
            options%ctest = 2
            options%print_freq_itn = 1
            options%print_freq_head = 1
            ldamp = .true.
         end if

         if (prblm.eq.nprob+4) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing options%ctest=3........"
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing options%ctest=3....."
            options%ctest = 3
            options%print_freq_itn = 1
            options%print_freq_head = 1
            ldamp = .true.
         end if

         if (prblm.eq.nprob+5) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing options%ctest=3........"
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing options%ctest=3....."
            options%ctest = 3
            options%print_freq_itn = 1
            options%print_freq_head = 1
         end if

         if (prblm.eq.nprob+6) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing b=0...................."
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing b=0................."
            b(1:m) = zero
            u(1:m) = b(1:m)
            options%ctest = 3
         end if

         if (prblm.eq.nprob+9) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing b=0...................."
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing b=0................."
            b(1:m) = zero
            u(1:m) = b(1:m)
         end if

         if (prblm.eq.nprob+10) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing compatible............."
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing compatible.........."
            x(1:n) = one
            do i = 1,m
              b(i) = zero
              do j = 1,n
                 b(i) = b(i) + a(i,j)
              end do
            end do
            u(1:m) = b(1:m)
         end if

         ! compute w = A'*b
         call matrix_mult_trans(m,n,ptr,row,val,b,w,one,zero)
         normAb = dnrm2(n,w,1)
         normb  = dnrm2(m,b,1)

         do

            if (ldamp) then
               call lsmr_solve(action, m, n, u, v, y, keep, options, inform, &
                  damp=lambda)
            else
               call lsmr_solve(action, m, n, u, v, y, keep, options, inform)
            end if

            if (action.eq.0) then
              ! we are done. check whether we have an unexpected error.

              if (prblm.eq.3) then
                 call print_result(nout, nfail, inform%flag, lsmr_stop_itnlim)
                 flag = -1
                 exit
              else if (inform%flag.ge.7) then
                nfail = nfail + 1
                flag = -1
                exit
              end if

              if (prblm.eq.nprob+4 .or. prblm.eq.nprob+5) then
                 ! expect to have LS solution
                 call print_result(nout, nfail, inform%flag, lsmr_stop_LS_atol,&
                   expected1=lsmr_stop_compatible)
              end if

              if (prblm.eq.nprob+7 .or. prblm.eq.nprob+8) then
                 ! test with A'b=0 should have x=0
                 call print_result(nout, nfail, inform%flag, lsmr_stop_x0)
                 flag = -1
                 exit
              end if

              if (prblm.eq.nprob+6 .or. prblm.eq.nprob+9) then
                 ! test with b=0 should have x=0
                 call print_result(nout, nfail, inform%flag, lsmr_stop_x0)
                 flag = -1
                 exit
              end if

              if (prblm.eq.nprob+10) then
                 ! test with compatible
                 call print_result(nout, nfail, inform%flag,&
                   lsmr_stop_compatible)
              end if

              ! recover solution
              x(1:n) = y(1:n)

              exit

           !!!!!!!!!!!!!!!!!!!!!!!!!
           else if (action.eq.1) then

                 ! Compute v = v + A'*u without altering u
                  call matrix_mult_trans(m,n,ptr,row,val,u,v,one,one)

         ! artifically set v = 0 to test out some code.
         if (prblm.eq.nprob+7) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing A'b=0.................."
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing A'b=0..............."
            v(1:n) = zero
         end if

         if (prblm.eq.nprob+8) then
            if (nout.gt.0) &
            write(nout,"(a)",advance="no") " * Testing A'b=0.................."
            if (pout.gt.0) &
            write(pout,"(/a)",advance="no") " * Testing A'b=0..............."
            v(1:n) = zero
            options%ctest = 1
         end if

           !!!!!!!!!!!!!!!!!!!!!!!!!
           else if (action.eq.2) then

                 ! Compute u = u + A*v  without altering v
                 call matrix_mult(m,n,ptr,row,val,v,u,one,one)

           !!!!!!!!!!!!!!!!!!!!!!!!!

           else if (action.eq.3) then

              ! test for convergence (only occurs if options%ctest = 1 or 2)

              if (options%ctest.eq.2) then
                 ! we can use the data in inform for testing
                 normr = inform%normr
                 normAr = inform%normAPr
                 ratio = (normAr/normr)/(normAb/normb)
                 if (normr/normb < tol1 .or. ratio < tol2) then
                    x(1:n) = y(1:n)
                    exit
                 end if
              else
                 ! ctest = 1. We have to do the testing.
                 ! recover solution x = Py
                 ! test x for convergence
                 x(1:n) = y(1:n)

                 ! compute residuals z = r = b - A*x
                 z(1:m) = b(1:m)
                 call matrix_mult(m,n,ptr,row,val,x,z,-one,one)
                 normr = dnrm2(m,z,1)

                ! compute w = A'*r = A'*z
                call matrix_mult_trans(m,n,ptr,row,val,z,w,one,zero)
                normAr = dnrm2(n,w,1)

                ratio = (normAr/normr)/(normAb/normb)
                ! write (*,*) normr,normr/normb,ratio

                if (normr/normb < tol1) then
                   exit
                else if (ratio < tol2) then
                   exit
                end if

              end if

           !!!!!!!!!!!!!!!!!!!!!!!!!
           else

              ! we should not get here!
              if (nout.gt.0) &
              write (nout,'(a,i3)') ' We have a bug. action = ',action
              flag = -1
              nfail = nfail + 1
              exit
           end if

        end do

        if (flag.eq.-1) cycle

        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! compute residuals u = r = b - A*x
        u(1:m) = b(1:m)
        call matrix_mult(m,n,ptr,row,val,x,u,-one,one)
        normr = dnrm2(m,u,1)

        v(1:n) = zero
        ! compute v = A'*r = A'*u
        call matrix_mult_trans(m,n,ptr,row,val,u,v,one,zero)

     !   normAPr = dnrm2(n,v,1)

        ! crude check on normr and normAPr against values returned by LSMR
     !   if (options%ctest.eq.2 .or. options%ctest.eq.3) then
     !     if (abs(normr-inform%normr).gt.1000*epsilon(one) .or.  &
     !        abs(normAPr-inform%normAPr).gt.1000*epsilon(one) ) then
     !        nfail = nfail + 1
     !        write (nout,'(a)') &
     !          ' Something wrong when checking normr and normAPr'
      ! write (*,*)normr,inform%normr,normAPr,inform%normAPr
      ! write (*,*) abs(normr-inform%normr),abs(normAPr-inform%normAPr)
     !     end if
     !   end if
        !!!!!!!!!!!!!!!!!!!!!!!!!

        ! use Lapack to check the solution
        ! First find size of workspace
        lwork = -1
        allocate (work(1))
        call dgels( 'n', m, n, 1, a, maxm, rhs, maxm, work, lwork, flag)

        if (flag /= 0) then
           write (nout,'(a,i3)') &
         ' Unexpected error from dgels (initial call). info = ',flag
           nfail = nfail + 1;   cycle main
        end if

        ! use output to allocate required space and run dgels
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))
        rhs(1:m,1) = b(1:m)

        call dgels( 'n', m, n, 1, a, maxm, rhs, maxm, work, lwork, flag)

        if (flag /= 0) then
           if (nout.gt.0) &
           write (nout,'(a,i3)') ' Unexpected error from dgels. info = ',flag
           nfail = nfail + 1;   cycle main
        end if

        u(1:m) = b(1:m)
        call matrix_mult(m,n,ptr,row,val,rhs,u,-one,one)

        ! the computed solution from dgels is in rhs.
        ! the LSMR computed solution is in x.
        normx = dnrm2(n,x,1)
        norm_rhs = dnrm2(n,rhs,1)
        if (abs(normx-norm_rhs)/normx.gt.10*sqrt(epsilon(one)) ) then
           nfail = nfail + 1
           write (nout,'(a)') ' Something wrong when checking norm of solution'
           ! note: if options%atol and options%btol too large (ie stopping
           ! criteria for LSMR not stringent enough) may find we get this
           ! error. If so, try reducing the tolerances in stopping
           ! criteria to get more accurate answer.
         !   write (*,*) abs(normx-norm_rhs)/normx, normx, norm_rhs,m,n
        end if


       ! write (*,*) rhs(1:min(10,m),1)
       ! write (*,*) x(1:min(10,m))

        deallocate (work)

      end do main

      ! call for deallocating components of keep
      call LSMR_free (keep, stat)

      if (stat.ne.0) then
         nfail = nfail + 1
         if (nout.gt.0) &
         write (nout,'(a)') ' Non zero stat parameter returned from LSMR_free'
      end if

      if (nout.gt.0) write (nout,'(//a,i5)') &
         ' At end of test_special number of failures = ',nfail

   end subroutine test_special

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine test_random

   integer(ip), parameter :: maxm =  50!200 ! largest problem size (rows)
   integer(ip), parameter :: maxn =  40!100 ! largest problem size (cols)
   integer(ip), parameter :: nprob = 200 ! number of error-free runs

   type (random_state) :: state

   integer(ip), allocatable :: ptr(:),row(:)
   real(wp),    allocatable :: val(:), work(:)
   real(wp)    :: a(maxm,maxn), b(maxm), x(maxn), y(maxn), rhs(maxm,1)
   real(wp)    :: u(maxm), v(maxn) ! reverse communication arrays
   real(wp)    :: d(maxn), w(maxn) , z(maxm)

    type ( lsmr_keep )    :: keep
    type ( lsmr_options ) :: options
    type ( lsmr_inform )  :: inform

   integer(ip)  :: action
   integer(ip) :: flag
   integer(ip) :: i
   integer(ip) :: lwork
   integer(ip) :: n,ne
   integer(ip) :: prblm
   integer(ip) :: stat

   logical :: precond

   real(wp) :: normAr, normr, normx, norm_rhs, normb, normAb, ratio
   real(wp) :: dnrm2
   real(wp) :: lambda

    if (nout.gt.0) then
      write(nout, "(a)")
      write(nout, "(a)") "======================="
      write(nout, "(a)") "Testing random matrices"
      write(nout, "(a)") "======================="
    end if

 main: do prblm = 1,nprob

         call matrix_gen(state, maxm, m, maxn, n, ne, ptr, row, val, b, a, &
            flag)

         if (flag /= 0) then
           if (nout.gt.0) write(nout,'(a,i3)') &
         ' Unexpected error when generating matrix. flag=', flag
           nfail = nfail + 1;   cycle main
         end if

         if (nout.gt.0) write(nout, "(/a, i3, 3(1x, a, i8),a)") " * number ", &
         prblm, "m = ", m, "n = ", n, "ne = ", ne, " ..."

         call set_controls(options)

         ! switch off printing
         options%unit_diagnostics = -1
         options%unit_error       = -1

         i = random_integer(state, 10)
         options%ctest = i
         if (i.gt.3) options%ctest = 3 ! Michael's stopping criteria
         ! choose these small values to get results that are close to
         ! lapack results
      !   options%ctest = 1
         options%atol = 10d-12
         options%btol = 10d-12

         ! Choose localsize
         options%localSize = random_integer(state, 100)
         if (options%localSize.gt.min(m,n)) options%localSize = 0

         ! decide whether to use preconditioning
         precond = .true.
         i = random_integer(state, 10)
         if (i.lt.5) precond = .false.
         if (precond) then
            if (m.ge.n) then
              ! compute norms of cols of A (for diagonal preconditioning)
               call norm_colA (n, ptr, val, d)
            else
              ! for underdetermined problems, use P = 2*I (so that
              ! we get a solution we can compare with lapack dgels)
              d(1:n)= two
            end if
         end if

         lambda = zero
         if (prblm/2*2.eq.prblm) lambda = 1.0d-08

         action = 0
         flag = 0
         u(1:m) = b(1:m)

         ! compute w = A'*b
         call matrix_mult_trans(m,n,ptr,row,val,b,w,one,zero)
         normAb = dnrm2(n,w,1)
         normb  = dnrm2(m,b,1)

         do
            if (lambda.ne.zero) then
               call lsmr_solve(action, m, n, u, v, y, keep, options, inform, &
                 damp=lambda)
            else
               call lsmr_solve(action, m, n, u, v, y, keep, options, inform)
            end if

           if (action.eq.0) then
              ! we are done. check whether we have an unexpected error.
              if (inform%flag.ge.7) then
                nfail = nfail + 1
                flag = -1
                exit
              end if

              ! recover solution x = Py
              if (precond) then
                 do i = 1,n
                    x(i) = y(i)*d(i)
                 end do
              else
                 x(1:n) = y(1:n)
              end if
              exit

           !!!!!!!!!!!!!!!!!!!!!!!!!
           else if (action.eq.1) then

              if (precond) then
                 ! Compute v = v + P'A'*u without altering u
                 ! first w = A'*u then v = v + P'*w
                  call matrix_mult_trans(m,n,ptr,row,val,u,w,one,zero)

                 ! Diagonal P
                 do i = 1,n
                    v(i) = v(i) + w(i)*d(i)
                 end do
              else
                 ! Compute v = v + A'*u without altering u
                  call matrix_mult_trans(m,n,ptr,row,val,u,v,one,one)
              end if

           !!!!!!!!!!!!!!!!!!!!!!!!!
           else if (action.eq.2) then

              if (precond) then
                 ! Compute u = u + AP*v  without altering v.
                 ! first w = P*v then u = u + A*w
                ! Diagonal P
                 do i = 1,n
                    w(i) = v(i)*d(i)
                 end do
                 call matrix_mult(m,n,ptr,row,val,w,u,one,one)
              else
                 ! Compute u = u + A*v  without altering v
                 call matrix_mult(m,n,ptr,row,val,v,u,one,one)
              end if

           !!!!!!!!!!!!!!!!!!!!!!!!!

           else if (action.eq.3) then

              ! test for convergence (only occurs if options%ctest = 1 or 2)

              if (options%ctest.eq.2 .and. .not.precond) then
                 ! we can use the data in inform for testing
                 normr = inform%normr
                 normAr = inform%normAPr
                 ratio = (normAr/normr)/(normAb/normb)
                 if (normr/normb < tol1 .or. ratio < tol2) then
                    if (precond) then
                       do i = 1,n
                          x(i) = y(i)*d(i)
                       end do
                    else
                       x(1:n) = y(1:n)
                    end if
                    exit
                 end if
              else
                 ! ctest = 1. We have to do the testing.
                 ! recover solution x = Py
                 ! test x for convergence
                 if (precond) then
                    do i = 1,n
                       x(i) = y(i)*d(i)
                    end do
                 else
                    x(1:n) = y(1:n)
                 end if

                 ! compute residuals z = r = b - A*x
                 z(1:m) = b(1:m)
                 call matrix_mult(m,n,ptr,row,val,x,z,-one,one)
                 normr = dnrm2(m,z,1)

                ! compute w = A'*r = A'*z
                call matrix_mult_trans(m,n,ptr,row,val,z,w,one,zero)
                normAr = dnrm2(n,w,1)

                ratio = (normAr/normr)/(normAb/normb)
                ! write (*,*) normr,normr/normb,ratio

                if (normr/normb < tol1) then
                   exit
                else if (ratio < tol2) then
                   exit
                end if

              end if

           !!!!!!!!!!!!!!!!!!!!!!!!!
           else

              ! we should not get here!
              if (nout.gt.0) write (nout,'(a,i3)') &
              ' We have a bug. action = ',action
              flag = -1
              nfail = nfail + 1
              exit
           end if

        end do

        if (flag.eq.-1) cycle

        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! compute residuals u = r = b - A*x
        u(1:m) = b(1:m)
        call matrix_mult(m,n,ptr,row,val,x,u,-one,one)
        normr = dnrm2(m,u,1)

        v(1:n) = zero
        if (precond) then
           ! compute v = P'A'*r = P'A'*u
           call matrix_mult_trans(m,n,ptr,row,val,u,v,one,zero)
           ! Diagonal P
            do i = 1,n
               v(i) = v(i)*d(i)
            end do
        else
           ! compute v = A'*r = A'*u
           call matrix_mult_trans(m,n,ptr,row,val,u,v,one,zero)
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!

        ! use Lapack to check the solution
        ! First find size of workspace
        lwork = -1
        allocate (work(1))
        call dgels( 'n', m, n, 1, a, maxm, rhs, maxm, work, lwork, flag)

        if (flag /= 0) then
           if (nout.gt.0) write (nout,'(a,i3)') &
         ' Unexpected error from dgels (initial call). info = ',flag
           nfail = nfail + 1;   cycle main
        end if

        ! use output to allocate required space and run dgels
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))
        rhs(1:m,1) = b(1:m)

        call dgels( 'n', m, n, 1, a, maxm, rhs, maxm, work, lwork, flag)

        if (flag /= 0) then
           if (nout.gt.0) &
           write (nout,'(a,i3)') ' Unexpected error from dgels. info = ',flag
           nfail = nfail + 1;   cycle main
        end if

        u(1:m) = b(1:m)
        call matrix_mult(m,n,ptr,row,val,rhs,u,-one,one)

        ! the computed solution from dgels is in rhs.
        ! the LSMR computed solution is in x.
        normx = dnrm2(n,x,1)
        norm_rhs = dnrm2(n,rhs,1)
        if (abs(normx-norm_rhs)/normx.gt.10*sqrt(epsilon(one)) ) then
           nfail = nfail + 1
           if (nout.gt.0) &
           write (nout,'(a)') ' Something wrong when checking norm of solution'

           ! note: if options%atol and options%btol too large (ie stopping
           ! criteria for LSMR not stringent enough) may find we get this
           ! error. If so, try reducing the tolerances in stopping
           ! criteria to get more accurate answer.
         !   write (*,*) abs(normx-norm_rhs)/normx, normx, norm_rhs,m,n
        end if


       ! write (*,*) rhs(1:min(10,m),1)
       ! write (*,*) x(1:min(10,m))

        deallocate (work)

      end do main

      ! call for deallocating components of keep
      call LSMR_free (keep, stat)

      if (stat.ne.0) then
         nfail = nfail + 1
         if (nout.gt.0) &
         write (nout,'(a)') ' Non zero stat parameter returned from LSMR_free'
      end if

      if (nout.gt.0) write (nout,'(//a,i5)') &
         ' At end of test_random number of failures = ',nfail

    end subroutine test_random


!********************************************
      ! set default control values

      subroutine set_controls(options)

      type(lsmr_options), intent(inout) :: options

        options%atol             = sqrt(epsilon(one))
        options%btol             = sqrt(epsilon(one))
        options%conlim           = 1/(10*sqrt(epsilon(one)))
        options%ctest            = 3
        options%itnlim           = -1
        options%localSize        = 0
        options%print_freq_head  = 20
        options%print_freq_itn   = 10
        options%unit_diagnostics = 6
        options%unit_error       = 6

      end subroutine set_controls

!********************************************
      ! Generate random rectangular matrix
      ! Stored in CSC format using ptr, row, val.
      ! Also held as dense matrix in a(1:m,1:n)

   subroutine matrix_gen(state, maxm, m, maxn, n, ne, ptr, row, val, b, a, &
      flag)

      type (random_state), intent(inout) :: state
      integer(ip), intent(in)  :: maxm ! max row dim. of matrix
      integer(ip), intent(out) :: m    ! row dim. of matrix
      integer(ip), intent(in)  :: maxn ! max col. dim. of matrix
      integer(ip), intent(out) :: n    ! col. dim. of matrix
      integer(ip), intent(out) :: ne   ! entries in matrix
      integer(ip), intent(out) :: flag
      integer(ip), intent(out), dimension(:), allocatable :: ptr
      integer(ip), intent(out), dimension(:), allocatable :: row
      real(wp),    intent(out)                            :: a(maxm,maxn)
      real(wp),    intent(out)                            :: b(maxm)
      real(wp),    intent(out), dimension(:), allocatable :: val
      integer(ip) :: i,j,nzin,matrix_type,st

      flag = 0
      m = random_integer(state, maxm)
      m = max(1,m)

      n = random_integer(state, maxn)
      n = max(1,n)

      nzin = random_integer(state, m*n/3)
      nzin = max(3*m,nzin)
      nzin = min(n*m,nzin)
      if (m == 1 .and. n == 1) nzin = 1



      deallocate(row,stat=st)
      deallocate(ptr,stat=st)
      deallocate(val,stat=st)

      allocate(row(nzin),val(nzin),ptr(n+1),stat=st)
      if (st /= 0) then
        flag = -1
        return
      end if

      matrix_type = 1
      call random_matrix_generate(state, matrix_type, m, n, &
      nzin, ptr, row, flag, val=val, nonsingular=.true., sort=.true.)

      if (flag.lt.0) return

      ne = ptr(n+1) - 1

      ! copy into dense matrix a
      a = zero
      do j = 1,n
         do i = ptr(j),ptr(j+1)-1
            a(row(i),j) = val(i)
         end do
      end do

      ! generate a random right-hand side vector b

      do i = 1,m
         b(i) = random_real(state)
      end do

    end subroutine matrix_gen

!**************************************************************

   subroutine norm_colA (n, ptr, val, d)

      ! computes 2-norm of columns of A (held in CSC format)
      ! and set d(j) = 1/sqrt(2-norm(j)). Gives diagonal scaling.

      integer(ip), intent(in) :: n
      integer(ip), intent(in) :: ptr(n+1)
      real(wp), intent(in)    :: val(:)
      real(wp), intent(out)   :: d(n)

      integer(ip) :: i,j
      real(wp) :: sum

      do j = 1,n
         sum = zero
         do i = ptr(j),ptr(j+1)-1
            sum = sum + val(i)**2
         end do
      ! may have null cols so make sure we don't have too small a scaling factor
         if (sum.lt.sqrt(epsilon(one))) sum = one
         d(j) = one/sqrt(sum)
      end do

   end subroutine norm_colA

!**************************************************************

      ! Takes b and computes v = beta*v + alpha* A^T * u (A in CSC format)

      subroutine matrix_mult_trans(m,n,ptr,row,val,u,v,alpha,beta)

      integer,  intent(in) :: m,n
      integer,  intent(in) :: ptr(n+1),row(:)
      real(wp), intent(in) :: val(:),u(m)
      real(wp), intent(inout) :: v(n) !v must be defined on input if
       ! beta/=0.0; otherwise v does not need to be defined
      real(wp), intent(in)  :: alpha,beta

      integer:: i,j,k
      real(wp) :: sum

      if (beta.ne.zero) then
         do j = 1,n
            sum = zero
            do k = ptr(j),ptr(j+1)-1
               i = row(k)
               sum = sum + val(k)*u(i)
            end do
            v(j) = beta*v(j) + alpha*sum
         end do
      else
         do j = 1,n
            sum = zero
            do k = ptr(j),ptr(j+1)-1
               i = row(k)
               sum = sum + val(k)*u(i)
            end do
            v(j) = alpha*sum
         end do
      end if

      end subroutine matrix_mult_trans

!**************************************************************

      ! Takes b and computes u = beta*u + alpha*A * v (A in CSC format)

      subroutine matrix_mult(m,n,ptr,row,val,v,u,alpha,beta)

      integer,  intent(in) :: m,n
      integer,  intent(in) :: ptr(n+1),row(:)
      real(wp), intent(in) :: val(:),v(n)
      real(wp), intent(inout) :: u(m)
      real(wp), intent(in)  :: alpha,beta

      integer:: i,j,k
      real(wp) :: temp

         call dscal(m,beta,u,1)

         do j = 1,n
            temp = alpha*v(j)
            do k = ptr(j),ptr(j+1)-1
                i = row(k)
                u(i) = u(i) + val(k)*temp
            end do
         end do


      end subroutine matrix_mult

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lsmrtest(m,n,nduplc,npower,damp,localSize,flag,mout)

    integer(ip),  intent(in) :: m, n, nduplc, npower, &
                                localSize,        & ! Local reorthogonalization
                                mout
    integer(ip), intent(out) :: flag ! error flag. set to <0 if error found
    real(wp), intent(in)     :: damp

    !-------------------------------------------------------------------
    ! This is an example driver routine for running LSMR.
    ! It generates a test problem, solves it, and examines the results.
    !
    ! 27 Sep 2007: f90 version of lsqrtest.
    ! 17 Jul 2010: LSMR version derived from LSQR equivalent.
    !              localSize specifies that each bidiagonalization vector v
    !              should be reorthogonalized wrto the last "localSize" v's.
    !              localSize >= 0.
    ! 18 March 2016 : added flag for checking (removed stop statement)
    !------------------------------------------------------------------------

    intrinsic       :: epsilon

    ! Local arrays and variables
    real(wp), allocatable :: b(:), w1(:), w2(:), x(:), xtrue(:)
    real(wp), allocatable :: u(:), v(:)  ! reverse communication vectors.
    integer(ip)     :: j, minmn, nprint, stat
    real(wp)        :: condA, eps, normr,  norme, etol, normw, normx

    real(wp) :: dnrm2

    integer(ip)     :: action
    type ( lsmr_keep )    :: keep
    type ( lsmr_options ) :: options
    type ( lsmr_inform )  :: inform

    ! Local constants
    character(len=*),  &
              parameter :: line = '----------------------------------'

    flag = 0
    eps    = epsilon(eps)
    if (mout.gt.0) then
     if (eps > 1e-9) then
       write(nout,*) ' '
       write(nout,*) 'WARNING: '
       write(nout,*) 'MACHINE PRECISION EPS =', eps, '  SEEMS TO BE INADEQUATE'
       write(nout,*) ' '
     end if
    end if

    allocate ( b(m), x(n), xtrue(n) )

    !------------------------------------------------------------------------
    ! Generate the specified test problem
    ! and check that Aprod1, Aprod2 form y + Ax and x + A'y consistently.
    !------------------------------------------------------------------------
    do j = 1,n                         ! Set the desired solution xtrue.
       xtrue(j) = 0.1*j                ! For least-squares problems, this is it.
    end do                             ! If m<n, lstp will alter it.

    minmn  = min(m,n)                       ! Allocate arrays.
    allocate( d(minmn), hy(m), hz(n) )      ! Vectors defining A = Y D Z.
    allocate( wm(m), wn(n) )                ! Work vectors for Aprod1, Aprod2.
    allocate( u(m), v(n) )                  ! Vectors for reverse communication.

    ! Generate test problem. If m<n, xtrue is altered.
    call lstp (m,n,nduplc,npower,damp,xtrue,b,condA,normr,hy,hz,d,wm,wn)

    if (mout.gt.0) write(mout,1000) line,line,m,n,nduplc,npower,damp, &
                     condA,normr,line,line

    ! Check Aprod1, Aprod2. Use u,v as work arrays.
    allocate( w1(m), w2(n) )
    call Acheck(m,n,Aprod1,Aprod2,mout,flag, u, v, w1, w2)
    if (flag > 0) then
       if (mout.gt.0) write(mout,'(a)') 'Check tol in subroutine Acheck'
       flag = -1
       return
    end if
    deallocate (w1)
    deallocate (w2)

    !------------------------------------------------------------------------
    ! Set input parameters for LSMR
    ! and solve the problem defined by Aprod1, Aprod2, b, damp.
    ! Next line would ask for standard errors only if they are well-defined.
    !------------------------------------------------------------------------
    options%atol   = eps**0.99                      ! Asks for high accuracy.
    options%btol   = options%atol
    options%conlim = 1000.0 * condA
    options%itnlim = 4*(m + n + 50)

    options%unit_diagnostics = mout
    options%unit_error       = mout
    options%localsize = localsize


    action = 0
    u(1:m) = b(1:m)
    options%ctest = 3 ! Michael's stopping criteria
    do
      call lsmr_solve(action, m, n,  u, v, x, keep, &
         options, inform, damp=damp)

      if (action.eq.0) then
        ! we are done
        exit
      else if (action.eq.1) then
        ! Compute v = v + A'*u without altering u
          call Aprod2(m,n,v,u)

      else if (action.eq.2) then
         ! Compute u = u + A*v  without altering v
          call Aprod1(m,n,v,u)
         ! test for convergence if options%ctest = 1 or 2
      else if (action.eq.3) then

      else
         ! we should not get here!
         if (mout.gt.0) write(mout,'(a,i3)') ' We have a bug. action = ',action
         flag = -2
         go to 50
      end if

    end do

    ! call for deallocating components of keep
    call LSMR_free (keep, stat)

    if (stat.ne.0) write(mout,'(a)') &
     ' Non zero stat parameter returned from LSMR_free'

    call xcheck( m, n, Aprod1, Aprod2, b, damp, x, &    ! Check x
                 inform%normAP, mout, flag, u, v )

    nprint = min(m,n,8)
    if (mout.gt.0) &
       write(mout,2500) (j, x(j), j=1,nprint)    ! Print some of the solution

    wn     = x - xtrue                           ! Print a clue about whether
    normw  = dnrm2 (n,wn,1)                      ! the solution looks OK.
    normx  = dnrm2 (n,xtrue,1)
    norme  = normw/(one + normx)
    etol   = 1.0e-3
    if (norme <= etol) then
       if (mout.gt.0) write(mout,3000) norme
       flag = 1
    else
       if (mout.gt.0) write(mout,3100) norme
       flag = -3
    end if

 50 continue

    deallocate( wn )
    deallocate( wm )
    deallocate( hz )
    deallocate( hy )
    deallocate( d )
    deallocate (u)
    deallocate (v)
    deallocate ( b )
    deallocate ( x )
    deallocate ( xtrue )

    return

 1000 format(1p &
      // 1x, 2a &
      /  ' Least-Squares Test Problem      P(', 4i5, e12.2, ' )'      &
      /  ' Condition no. =', e12.4, '     Residual function =', e17.9 &
      /  1x, 2a)
 2500 format(/' Solution  x:'          / 4(i6, g14.6))
 3000 format(1p &
      /  ' LSMR  appears to be successful.'                           &
      /  ' Relative error in  x  =', e10.2)
 3100 format(1p &
      /  ' LSMR  appears to have failed.  '                           &
      /  ' Relative error in  x  =', e10.2)

  end subroutine lsmrtest

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Hprod (n,z,x)

    integer(ip),  intent(in)    :: n
    real(wp), intent(in)    :: z(n)
    real(wp), intent(inout) :: x(n)

    !-------------------------------------------------------------------
    ! Hprod  applies a Householder transformation stored in z
    ! to return x = (I - 2*z*z')*x.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    integer(ip)             :: i
    real(wp)                :: s

    s = zero
    do i = 1,n
       s = z(i)*x(i) + s
    end do

    s = s + s
    x = x - s*z

  end subroutine Hprod

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod1(m,n,x,y)

    integer(ip),  intent(in)    :: m,n
    real(wp), intent(in)        :: x(n)
    real(wp), intent(inout)     :: y(m)

    !-------------------------------------------------------------------
    ! Aprod1 computes y = y + A*x without altering x,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic           :: min
    integer(ip)         :: minmn

    minmn = min(m,n)
    wn    = x
    call Hprod (n,hz,wn)
    wm(1:minmn) = d(1:minmn)*wn(1:minmn)
    wm(n+1:m)   = zero
    call Hprod (m,hy,wm)
    y = y + wm

  end subroutine Aprod1

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod2(m,n,x,y)

    integer(ip),  intent(in)    :: m,n
    real(wp), intent(inout)     :: x(n)
    real(wp), intent(in)        :: y(m)

    !-------------------------------------------------------------------
    ! Aprod2 computes x = x + A'*y without altering y,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic               :: min
    integer(ip)             :: minmn

    minmn = min(m,n)
    wm    = y
    call Hprod (m,hy,wm)
    wn(1:minmn) = d(1:minmn)*wm(1:minmn)
    wn(m+1:n)   = zero
    call Hprod (n,hz,wn)
    x = x + wn

  end subroutine Aprod2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lstp (m,n,nduplc,npower,damp,x,b,condA,normr,hy,hz,d,wm,wn)

    integer(ip),  intent(in)    :: m, n, nduplc, npower
    real(wp), intent(in)        :: damp
    real(wp), intent(inout)     :: x(n)
    real(wp), intent(out)       :: b(m)
    real(wp), intent(out)       :: condA, normr
    real(wp), intent(out)       :: hy(m), hz(n), d(min(m,n)),wm(m),wn(n)

    !-------------------------------------------------------------------
    ! lstp  generates a sparse least-squares test problem of the form
    !           (   A    )*x = ( b )
    !           ( damp*I )     ( 0 )
    ! for solution by LSMR, or a sparse underdetermined system
    !            Ax + damp*s = b
    ! for solution by LSMR or CRAIG.  The matrix A is m by n and is
    ! constructed in the form  A = Y*D*Z,  where D is an m by n
    ! diagonal matrix, and Y and Z are Householder transformations.
    !
    ! m and n may contain any positive values.
    ! If m >= n  or  damp = 0, the true solution is x as given.
    ! Otherwise, x is modified to contain the true solution.
    !
    ! 1982---1991: Various versions implemented.
    ! 06 Feb 1992: lstp generalized to allow any m and n.
    ! 07 Sep 2007: Line by line translation for Fortran 90 compilers
    !              by Eric Badel <badel@nancy.inra.fr>.
    ! 23 Sep 2007: Fortran 90 version with modules.
    !-------------------------------------------------------------------

    intrinsic           :: min, cos, sin, sqrt
    integer(ip)         :: i, j, minmn
    real(wp)            :: alfa, beta, dampsq, fourpi, t
    real(wp) :: dnrm2

    !-------------------------------------------------------------------
    ! Make two vectors of norm 1.0 for the Householder transformations.
    ! fourpi  need not be exact.
    !-------------------------------------------------------------------
    minmn  = min(m,n)
    dampsq = damp**2
    fourpi = 4.0 * 3.141592
    alfa   = fourpi / m
    beta   = fourpi / n

    do i = 1,m
       hy(i) = sin( alfa*i )
    end do

    do i = 1,n
       hz(i) = cos( beta*i )
    end do

    alfa   = dnrm2 ( m, hy, 1 )
    beta   = dnrm2 ( n, hz, 1 )
    call dscal ( m, (- one/alfa), hy, 1 )
    call dscal ( n, (- one/beta), hz, 1 )

    !-------------------------------------------------------------------
    ! Set the diagonal matrix D.  These are the singular values of A.
    !-------------------------------------------------------------------
    do i = 1,minmn
       j    = (i-1+nduplc) / nduplc
       t    =  j*nduplc
       t    =  t / minmn
       d(i) =  t**npower
    end do

    condA  = (d(minmn)**2 + dampsq) / (d(1)**2 + dampsq)
    condA  = sqrt( condA )

    !-------------------------------------------------------------------
    ! If m >=n, the input x will be the true solution.
    ! If m < n, reset x = different true solution.
    ! It must be of the form x = Z(w1) for some w1.
    !                             (0 )
    ! We get w1 from the top of w = Zx.
    !-------------------------------------------------------------------
    wn = x
    call Hprod (n,hz,wn)
    if (m < n) then
       wn(m+1:n) = zero
       call Hprod (n,hz,wn)
       x = wn
    end if

    ! Let r = Y rbar and xbar = Z x.
    ! Solve D r1bar = damp^2 x1bar, where r1bar is in wm(1:minmn).

    wm(1:minmn)   = dampsq * wn(1:minmn)/d(1:minmn)

    wm(minmn+1:m) = one      ! Set r2bar to be anything (empty if m <= n).
    call Hprod (m,hy,wm)     ! Form r = Y rbar  in wm.

    normr  = dnrm2 (m,wm,1)  ! Compute norm(r)
    b      = wm              ! and  b = r + Ax.
    call Aprod1(m,n,x,b)

  end subroutine lstp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Acheck( m, n, Aprod1, Aprod2, nout, inform, w, v, y, x )

    integer(ip), intent(in)    :: m, n   ! No. of rows and cols of A
    integer(ip), intent(in)    :: nout   ! Output file number
    integer(ip), intent(out)   :: inform ! = 0 if Aprod1, Aprod2 seem ok
                                         ! = 1 otherwise
    real(wp), intent(out)      :: v(n), w(m), x(n), y(m) ! work arrays
    interface
       subroutine Aprod1(m,n,x,y)                   ! y := y + A*x
         integer(4),  parameter :: ip = kind( 1 )
         integer(4),  parameter :: wp = kind( 1.0D+0 )
         integer(ip), intent(in)    :: m,n
         real(wp),    intent(in)    :: x(n)
         real(wp),    intent(inout) :: y(m)
       end subroutine Aprod1

       subroutine Aprod2(m,n,x,y)                   ! x := x + A'*y
         integer(4),  parameter :: ip = kind( 1 )
         integer(4),  parameter :: wp = kind( 1.0D+0 )
         integer(ip), intent(in)    :: m,n
         real(wp),    intent(inout) :: x(n)
         real(wp),    intent(in)    :: y(m)
       end subroutine Aprod2
    end interface

    !-------------------------------------------------------------------
    ! One-liner: Acheck checks Aprod1 and Aprod2 for LSMR.
    !
    ! Purpose:   Acheck tests the user subroutines Aprod1 and Aprod2
    !   called by LSMR.  For some m x n matrix A,
    !   Aprod1 computes y := y + A*x  from given x,y without altering x,
    !   Aprod2 computes x := x + A'*y from given x,y without altering y.
    !   Acheck tries to verify that A and A' refer to the same matrix.
    !
    ! Method:    We cook up some unlikely vectors x and y of unit length
    !   and test if  y'(y + Ax)  =  x'(x + A'y).
    !
    ! Parameter Constants:
    !   Param   Type   Description
    !   power   real   eps**power is the tolerance for judging if
    !                  y'(y + Ax) = x'(x + A'y) to sufficient accuracy.
    !                  power should be in the range (0.25, 0.9) say.
    !                  For example, power = 0.75 means we are happy
    !                  if 3/4 of the available digits agree.
    !                  power = 0.5 seems a reasonable requirement
    !                  (asking for half the digits to agree).
    !
    ! History:
    ! 04 Sep 1991  Initial design and code.
    !              Michael Saunders, Dept of Operations Research,
    !              Stanford University.
    ! 10 Feb 1992  Aprod added as parameter.
    !              tol defined via power.
    ! 10 Feb 1992: Acheck revised and xcheck implemented.
    ! 27 May 1993: Acheck and xcheck kept separate from test problems.
    ! 23 Sep 2007: Acheck implemented as part of this f90 module.
    !-------------------------------------------------------------------

    intrinsic           :: abs, epsilon, sqrt

    ! Local variables
    integer(ip)         :: i, j
    real(wp)            :: alfa, beta, eps, t, test1, test2, test3, tol
    real(wp)            :: ddot, dnrm2

    eps    = epsilon(eps)
    tol    = eps**power
    if (nout > 0) write(nout,1000)

    !===================================================================
    ! Cook up some unlikely vectors x and y of unit length.
    !===================================================================
    t = one
    do j=1,n
       t    = t + one
       x(j) = sqrt(t)
    end do

    t = one
    do i=1,m
       t    = t + one
       y(i) = one/sqrt(t)
    end do

    alfa = dnrm2 (n,x,1)
    beta = dnrm2 (m,y,1)
    call dscal (n, (one/alfa), x, 1)
    call dscal (m, (one/beta), y, 1)

    !===================================================================
    ! Test if y'(y + Ax) = x'(x + A'y).
    !===================================================================
    w(1:m) = y(1:m)               ! First set  w = y + Ax,  v = x + A'y.
    v(1:n) = x(1:n)
    call Aprod1(m,n,x,w)
    call Aprod2(m,n,v,y)

    alfa   = ddot  (m,y,1,w,1)    ! Now set    alfa = y'w,  beta = x'v.
    beta   = ddot  (n,x,1,v,1)
    test1  = abs(alfa - beta)
    test2  = one + abs(alfa) + abs(beta)
    test3  = test1 / test2

    if (test3 <= tol) then        ! See if alfa and beta are essentially
       inform = 0                 ! the same.
       if (nout > 0) write(nout,1010) test3
    else
       inform = 1
       if (nout > 0) write(nout,1020) test3
    end if

    return

 1000 format(//' Enter Acheck.')
 1010 format(1p, &
        ' Aprod1, Aprod2 seem OK.  Relative error =', e10.1)
 1020 format(1p, &
        ' Aprod1, Aprod2 seem incorrect.  Relative error =', e10.1)

  end subroutine Acheck

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine xcheck( m, n, Aprod1, Aprod2, b, damp, x, &
                     Anorm, nout,                      &
                     inform, r, v )

    integer(ip),  intent(in)    :: m, n     ! No. of rows and cols of A
    integer(ip),  intent(in)    :: nout     ! Output file number
    integer(ip),  intent(out)   :: inform   ! = 0 if b = 0 and x = 0.
                                        ! = 1 2 or 3 if x seems to
                                        ! solve systems 1 2 or 3 below.
    real(wp), intent(in)    :: Anorm    ! An estimate of norm(A) or
                                        ! norm( A, delta*I ) if delta > 0.
                                        ! Provided by LSMR.
    real(wp), intent(in)    :: damp     ! Defines problem 3 below.
    real(wp), intent(in)    :: b(m)     ! The right-hand side of Ax ~= b.
    real(wp), intent(in)    :: x(n)     ! The given solution estimate.
    real(wp), intent(out)   :: r(m),v(n) ! work arrays

    interface
       subroutine Aprod1(m,n,x,y)                   ! y := y + A*x
  integer(4),  parameter :: ip = kind( 1 )
  integer(4),  parameter :: wp = kind( 1.0D+0 )
         integer(ip),  intent(in)    :: m,n
         real(wp), intent(in)        :: x(n)
         real(wp), intent(inout)     :: y(m)
       end subroutine Aprod1

       subroutine Aprod2(m,n,x,y)                   ! x := x + A'*y
  integer(4),  parameter :: ip = kind( 1 )
  integer(4),  parameter :: wp = kind( 1.0D+0 )
         integer(ip),  intent(in)    :: m,n
         real(wp), intent(inout)     :: x(n)
         real(wp), intent(in)        :: y(m)
       end subroutine Aprod2
    end interface

    !-------------------------------------------------------------------
    ! One-liner: xcheck tests if x solves a certain least-squares problem.
    !
    ! Purpose:   xcheck computes residuals and norms associated with
    ! the vector x and the least-squares problem solved by LSMR.
    ! It determines whether x seems to be a solution to any of three
    ! possible systems:  1. Ax = b
    !                    2. min norm(Ax - b)
    !                    3. min norm(Ax - b)^2 + damp^2 * norm(x)^2.
    !
    ! History:
    ! 07 Feb 1992  Initial design and code.
    !              Michael Saunders, Dept of Operations Research,
    !              Stanford University.
    ! 23 Sep 2007: xcheck implemented as part of this f90 module.
    !-------------------------------------------------------------------

    intrinsic           :: epsilon

    ! Local variables
    real(wp)            :: bnorm,dampsq,eps,rho1,rho2,sigma1,sigma2, &
                           test1,test2,test3,tol,snorm,xnorm,xsnorm

    real(wp)            :: dnrm2

    eps    = epsilon(eps)
    dampsq = damp**2
    tol    = eps**power

    r(1:m) = -b(1:m)        ! Compute the residual r = b - Ax
    call Aprod1(m,n,x,r)    ! via  r = -b + Ax,
    r(1:m) = -r(1:m)        !      r = -r.

    v(1:n) = zero           ! Compute v = A'r
    call Aprod2(m,n,v,r)    ! via  v = 0,  v = v + A'r.

    bnorm  = dnrm2 (m,b,1)  ! Compute the norms of b, x, r, v.
    xnorm  = dnrm2 (n,x,1)
    rho1   = dnrm2 (m,r,1)
    sigma1 = dnrm2 (n,v,1)

    if (nout > 0) write(nout,2200) damp, xnorm, rho1, sigma1

    if (damp == zero) then
       rho2   = rho1
       sigma2 = sigma1
    else
       v(1:n) = v(1:n) - dampsq*x(1:n)  ! v = A'r - damp**2 x.
       rho2   = sqrt(rho1**2 + dampsq*xnorm**2)
       sigma2 = dnrm2 (n,v,1)
       snorm  = rho1/damp
       xsnorm = rho2/damp
       if (nout > 0) write(nout,2300) snorm, xsnorm, rho2, sigma2
    end if

    !-------------------------------------------------------------------
    ! See if x seems to solve Ax = b  or  min norm(Ax - b)
    ! or the damped least-squares system.
    !-------------------------------------------------------------------
    if (bnorm == zero  .and.  xnorm == zero) then
       inform = 0
       test1  = zero
       test2  = zero
       test3  = zero
    else
       inform = 4
       test1  = rho1 / (bnorm + Anorm*xnorm)
       test2  = zero
       if (rho1  > zero) test2  = sigma1 / (Anorm*rho1)
       test3  = test2
       if (rho2  > zero) test3  = sigma2 / (Anorm*rho2)

       if (test3 <= tol) inform = 3
       if (test2 <= tol) inform = 2
       if (test1 <= tol) inform = 1
    end if

    if (nout > 0) write(nout,3000) inform,tol,test1,test2,test3
    return

 2200 format(1p                                             &
      // ' Enter xcheck.     Does x solve Ax = b, etc?'     &
      /  '    damp            =', e10.3                     &
      /  '    norm(x)         =', e10.3                     &
      /  '    norm(r)         =', e15.8, ' = rho1'          &
      /  '    norm(A''r)       =',e10.3, '      = sigma1')
 2300 format(1p/  '    norm(s)         =', e10.3            &
      /  '    norm(x,s)       =', e10.3                     &
      /  '    norm(rbar)      =', e15.8, ' = rho2'          &
      /  '    norm(Abar''rbar) =',e10.3, '      = sigma2')
 3000 format(1p/  '    inform          =', i2               &
      /  '    tol             =', e10.3                     &
      /  '    test1           =', e10.3, ' (Ax = b)'        &
      /  '    test2           =', e10.3, ' (least-squares)' &
      /  '    test3           =', e10.3, ' (damped least-squares)')

  end subroutine xcheck

!**************************************************************

subroutine print_result(nout, nfail, actual, expected, expected1)
   integer(ip), intent(in)    :: nout
   integer(ip), intent(inout) :: nfail
   integer(ip) :: actual
   integer(ip) :: expected
   integer(ip), optional :: expected1

   if(actual.eq.expected) then
      if (nout.gt.0) write(nout,"(a)") "ok"
      return
   endif

   if (present(expected1)) then
     if(actual.eq.expected1) then
       if (nout.gt.0) write(nout,"(a)") "ok"
       return
     endif
   endif

   if (nout.gt.0) write(nout,"(a)") "fail"
   if (nout.gt.0) &
   write(nout,"(2(a,i4))") "returned ", actual, ", expected ", expected
   nfail = nfail + 1

end subroutine print_result

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end program spral_lsmr_test
