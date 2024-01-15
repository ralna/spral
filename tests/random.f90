!
! Note: This isn't intended to be a comprehensive PRNG test suite, but is
! merely intended to highlight any serious flaws in the coding rather than
! numerical design.
!
program random
   use spral_random
   implicit none

   integer, parameter :: long = selected_int_kind(18)
   integer, parameter :: wp = kind(0d0)

   integer, parameter :: nsamples = 50000
   integer, parameter :: nbins = 100
   real, parameter :: require_confidence = 0.99

   integer :: errors

   errors = 0

   call test_real_dist()
   call test_integer32_dist()
   call test_integer64_dist()
   call test_logical_dist()

   if(errors.eq.0) then
      write(*, "(/a)") "==================="
      write(*, "(a)") "All tests suceeeded"
      write(*, "(a)") "==================="
   else
      write(*, "(/a)") "==================="
      write(*, "(a, i4)") "Failed ", errors
      write(*, "(a)") "==================="
   endif

contains

subroutine test_real_dist
   type(random_state) :: state

   integer :: i, j
   integer :: bin(nbins)
   real(wp) :: sample, chisq

   write(*, "(/a)") "====================================="
   write(*, "(a)") "Testing random_real()"
   write(*, "(a)") "====================================="


   !
   ! Test (-1,1) distribution
   !
   write(*, "(a)", advance="no") "Sampling Unif(-1,1)...... "
   ! Acquire sample
   bin(:) = 0
   do i = 1, nsamples
      sample = random_real(state)
      j = int(nbins * ((sample+1.0_wp)/2.0_wp)) + 1
      bin(j) = bin(j) + 1
   end do
   chisq = 1.0_wp/nsamples * sum( (bin(:)**2.0_wp) * nbins ) - nsamples
   if(chisq < chisq_pval(nbins-1, require_confidence)) then
      write(*, "(a)") "pass"
   else
      write(*, "(a)") "fail"
      write(*, "(a,es12.4)") "chisq statistic = ", chisq
      write(*, "(a,es12.4)") "chisq required  < ", &
         chisq_pval(nbins-1, require_confidence)
      errors = errors + 1
   endif

   !
   ! Test (0,1) distribution
   !
   write(*, "(a)", advance="no") "Sampling Unif(0,1)....... "
   ! Acquire sample
   bin(:) = 0
   do i = 1, nsamples
      sample = random_real(state, positive=.true.)
      j = int(nbins * sample) + 1
      bin(j) = bin(j) + 1
   end do
   chisq = 1.0_wp/nsamples * sum( (bin(:)**2.0_wp) * nbins ) - nsamples
   if(chisq < chisq_pval(nbins-1, require_confidence)) then
      write(*, "(a)") "pass"
   else
      write(*, "(a)") "fail"
      write(*, "(a,es12.4)") "chisq statistic = ", chisq
      write(*, "(a,es12.4)") "chisq required  < ", &
         chisq_pval(nbins-1, require_confidence)
      errors = errors + 1
   endif

end subroutine test_real_dist

subroutine test_integer32_dist
   type(random_state) :: state

   integer :: i
   integer :: bin(nbins)
   real(wp) :: chisq
   integer :: sample

   write(*, "(/a)") "====================================="
   write(*, "(a)") "Testing random_integer() 32-bit"
   write(*, "(a)") "====================================="

   !
   ! Test (1,...,n) distribution
   !
   write(*, "(a)", advance="no") "Sampling Unif(1,...,n)... "
   ! Acquire sample
   bin(:) = 0
   do i = 1, nsamples
      sample = random_integer(state,nbins)
      bin(sample) = bin(sample) + 1
   end do
   chisq = 1.0_wp/nsamples * sum( (bin(:)**2.0_wp) * nbins ) - nsamples
   if(chisq < chisq_pval(nbins-1, require_confidence)) then
      write(*, "(a)") "pass"
   else
      write(*, "(a)") "fail"
      write(*, "(a,es12.4)") "chisq statistic = ", chisq
      write(*, "(a,es12.4)") "chisq required  < ", &
         chisq_pval(nbins-1, require_confidence)
      errors = errors + 1
   endif
end subroutine test_integer32_dist

subroutine test_integer64_dist
   type(random_state) :: state

   integer :: i
   integer :: bin(nbins)
   real(wp) :: chisq
   integer(long) :: sample

   write(*, "(/a)") "====================================="
   write(*, "(a)") "Testing random_integer() 64-bit"
   write(*, "(a)") "====================================="

   !
   ! Test (1,...,n) distribution
   !
   write(*, "(a)", advance="no") "Sampling Unif(1,...,n)... "
   ! Acquire sample
   bin(:) = 0
   do i = 1, nsamples
      sample = random_integer(state,int(nbins,long))
      bin(sample) = bin(sample) + 1
   end do
   chisq = 1.0_wp/nsamples * sum( (bin(:)**2.0_wp) * nbins ) - nsamples
   if(chisq < chisq_pval(nbins-1, require_confidence)) then
      write(*, "(a)") "pass"
   else
      write(*, "(a)") "fail"
      write(*, "(a,es12.4)") "chisq statistic = ", chisq
      write(*, "(a,es12.4)") "chisq required  < ", &
         chisq_pval(nbins-1, require_confidence)
      errors = errors + 1
   endif
end subroutine test_integer64_dist

subroutine test_logical_dist
   type(random_state) :: state

   integer :: i, j
   integer :: bin(2)
   real(wp) :: chisq
   logical :: sample

   write(*, "(/a)") "====================================="
   write(*, "(a)") "Testing random_logical()"
   write(*, "(a)") "====================================="


   !
   ! Test (1,...,n) distribution
   !
   write(*, "(a)", advance="no") "Sampling B(1,0.5)........ "
   ! Acquire sample
   bin(:) = 0
   do i = 1, nsamples
      sample = random_logical(state)
      if(sample) then
         j = 1
      else
         j = 2
      endif
      bin(j) = bin(j) + 1
   end do
   chisq = 1.0_wp/nsamples * sum( (bin(:)**2.0_wp) * 2 ) - nsamples
   if(chisq < chisq_pval(1, require_confidence)) then
      write(*, "(a)") "pass"
   else
      write(*, "(a)") "fail"
      write(*, "(a,es12.4)") "chisq statistic = ", chisq
      write(*, "(a,es12.4)") "chisq required  < ", &
         chisq_pval(1, require_confidence)
      errors = errors + 1
   endif

end subroutine test_logical_dist

real(wp) function chisq_pval(dof, p)
   integer, intent(in) :: dof
   real, intent(in) :: p

   real(wp) :: xp

   if(p.eq.0.99) then
      xp = 2.33
   else
      write(*, "(a)") "Uncoded pval for chisq_pval"
      stop
   endif

   chisq_pval = dof + sqrt(2.0_wp*dof)*xp + (2*xp**2)/3 - 2/3.0_wp
end function chisq_pval

end program random
