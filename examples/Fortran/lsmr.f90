! examples/Fortran/lsmr.f90 - Example code for SPRAL_LSMR package
program spral_lsmr_example
    use  spral_lsmr
    implicit none

    integer, parameter :: wp = kind( 1.0d+0 )

    type ( lsmr_keep )    :: keep
    type ( lsmr_options ) :: options
    type ( lsmr_inform )  :: inform

    integer :: ptr(4), row(9)
    real(wp) :: b(5), u(5), v(3), val(9), x(3)

    integer :: action, m, n, stat

    ! Data for matrix:
    ! ( 1.0     -1.0 )
    ! (     2.0      )
    ! ( 2.0      2.0 )
    ! ( 5.0 3.0 -2.0 )
    ! (          6.0 )
    m = 5; n = 3
    ptr(1:n+1)        = (/ 1,               4,         6,                 10 /)
    row(1:ptr(n+1)-1) = (/ 1,     3,   4,   2,   4,    1,   3,    4,   5 /)
    val(1:ptr(n+1)-1) = (/ 1.0, 2.0, 5.0, 2.0, 3.0, -1.0, 2.0, -2.0, 6.0 /)
    ! Data for rhs b
    b(1:m) = (/ 1.0, 1.0, 1.0, 1.0, 1.0 /)

    ! prepare for LSMR calls (using no preconditioning and default
    ! settings, except switch off diagnostic printing)
    options%unit_diagnostics = -1
    action = 0
    u(1:m) = b(1:m)

    do

       call lsmr_solve(action, m, n, u, v, x, keep, options, inform)

       if (action.eq.0) then
          ! we are done.
          write (*,'(a,i3,a,i3)') ' Exit LSMR with inform%flag = ',inform%flag,&
            ' and inform%itn = ',inform%itn
          write (*,'(a)') ' LS solution is:'
          write (*,'(10f10.2)') x(1:n)
          exit

      else if (action.eq.1) then

            ! Compute v = v + A'*u without altering u
             call matrix_mult_trans(m,n,ptr,row,val,u,v)

      else if (action.eq.2) then

            ! Compute u = u + A*v  without altering v
            call matrix_mult(m,n,ptr,row,val,v,u)

      end if

   end do

   call lsmr_free(keep,stat)

  contains
!**************************************************************
      ! Takes b and computes u = u + A * v (A in CSC format)

      subroutine matrix_mult(m,n,ptr,row,val,v,u)

      integer,  intent(in) :: m,n
      integer,  intent(in) :: ptr(n+1),row(:)
      real(wp), intent(in) :: val(:),v(n)
      real(wp), intent(inout) :: u(m)

      integer:: i,j,k
      real(wp) :: temp

         do j = 1,n
            temp = v(j)
            do k = ptr(j),ptr(j+1)-1
                i = row(k)
                u(i) = u(i) + val(k)*temp
            end do
         end do

      end subroutine matrix_mult

!**************************************************************
      ! Takes b and computes v = v + A^T * u (A in CSC format)

      subroutine matrix_mult_trans(m,n,ptr,row,val,u,v)

      integer,  intent(in) :: m,n
      integer,  intent(in) :: ptr(n+1),row(:)
      real(wp), intent(in) :: val(:),u(m)
      real(wp), intent(inout) :: v(n)

      integer:: i,j,k
      real(wp) :: sum

         do j = 1,n
            sum = 0.0_wp
            do k = ptr(j),ptr(j+1)-1
               i = row(k)
               sum = sum + val(k)*u(i)
            end do
            v(j) = v(j) + sum
         end do

      end subroutine matrix_mult_trans

end program spral_lsmr_example
