module laplace2d_double

  interface apply_laplacian
    module procedure apply_laplacian_double, apply_lap_double
  end interface

  interface set_laplacian_matrix
    module procedure set_laplacian_matrix_double
  end interface

  interface apply_gauss_seidel_step
    module procedure apply_gauss_seidel_step_double, apply_gs_double
  end interface

contains

  subroutine apply_laplacian_double( mx, my, nx, x, Ax)

    implicit none

    integer, intent(in) :: mx, my, nx
    double precision, intent(in) :: x(mx*my, nx)
    double precision, intent(out) :: Ax(mx*my, nx)

    call apply_lap_double( mx, my, nx, x, Ax)

  end subroutine apply_laplacian_double

  subroutine apply_lap_double( mx, my, nx, x, Ax)
! Multiply x by Laplacian matrix for a square of size mx by my

    implicit none

    integer, intent(in) :: mx, my, nx
    double precision, intent(in) :: x(mx, my, nx)
    double precision, intent(out) :: Ax(mx, my, nx)
    integer :: i, j, k
    double precision :: z

    do k = 1, nx
      do i = 1, mx
        do j = 1, my
          z = 4*x(i, j, k)
          if ( i > 1 ) z = z - x(i - 1, j, k)
          if ( j > 1 ) z = z - x(i, j - 1, k)
          if ( i < mx ) z = z - x(i + 1, j, k)
          if ( j < my ) z = z - x(i, j + 1, k)
          Ax(i, j, k) = z
        end do
      end do
    end do

  end subroutine apply_lap_double

  subroutine set_laplacian_matrix_double( nx, ny, a, lda )
! Laplacian matrix for a rectangle of size nx x ny
    implicit none
    integer, intent(in) :: nx, ny, lda
    double precision, intent(out) :: a(lda, nx*ny)
    integer :: i, ix, iy
    a = 0
    do ix = 1, nx
      do iy = 1, ny
        i = ix + (iy - 1)*nx
        a(i, i) = 4
        if ( ix >  1 ) a(i, i -  1) = -1
        if ( ix < nx ) a(i, i +  1) = -1
        if ( iy > 1  ) a(i, i - nx) = -1
        if ( iy < ny ) a(i, i + nx) = -1
      end do
    end do
  end subroutine set_laplacian_matrix_double

  subroutine apply_gauss_seidel_step_double( mx, my, nx, x, Tx )
! Apply one Gauss-Seidel step for preconditioning

    implicit none

    integer, intent(in) :: mx, my, nx
    double precision, intent(in) :: x(mx*my, nx)
    double precision, intent(out) :: Tx(mx*my, nx)

    call apply_gs_double( mx, my, nx, x, Tx )

  end subroutine apply_gauss_seidel_step_double

  subroutine apply_gs_double( mx, my, nx, x, Tx )
! Apply one Gauss-Seidel step for preconditioning

    implicit none

    integer, intent(in) :: mx, my, nx
    double precision, intent(in) :: x(mx, my, nx)
    double precision, intent(out) :: Tx(mx, my, nx)
    integer :: i, j, k
    double precision :: z
    Tx = x/4
    do k = 1, nx
      do i = 1, mx ! forward update
        do j = 1, my
          z = 0
          if ( i > 1 ) z = z + Tx(i - 1, j, k)
          if ( j > 1 ) z = z + Tx(i, j - 1, k)
          if ( i < mx ) z = z + Tx(i + 1, j, k)
          if ( j < my ) z = z + Tx(i, j + 1, k)
          Tx(i, j, k) = (Tx(i, j, k) + z/4)/4
        end do
      end do
      do i = mx, 1, -1 ! backward update
        do j = my, 1, -1
          z = 0
          if ( i > 1 ) z = z + Tx(i - 1, j, k)
          if ( j > 1 ) z = z + Tx(i, j - 1, k)
          if ( i < mx ) z = z + Tx(i + 1, j, k)
          if ( j < my ) z = z + Tx(i, j + 1, k)
          Tx(i, j, k) = (Tx(i, j, k) + z/4)/4
        end do
      end do
    end do
  end subroutine apply_gs_double

end module laplace2d_double

module laplace2d

  use laplace2d_double

end module laplace2d


