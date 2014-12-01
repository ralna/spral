! Define BLAS API in Fortran module
module spral_blas_iface
  implicit none

  private
  public :: daxpy, dcopy, ddot, dnrm2, dscal
  public :: dgemm, dtrsm

  ! Level 1 BLAS
  interface
    subroutine daxpy( n, a, x, incx, y, incy )
      implicit none
      integer, intent(in) :: n, incx, incy
      double precision, intent(in) :: a
      double precision, intent(in   ), dimension(*) :: x
      double precision, intent(inout), dimension(*) :: y
    end subroutine daxpy
    subroutine dcopy( n, x, incx, y, incy )
      implicit none
      integer, intent(in) :: n, incx, incy
      double precision, intent(in ), dimension(*) :: x
      double precision, intent(out), dimension(*) :: y
    end subroutine dcopy
    double precision function ddot( n, x, incx, y, incy )
      implicit none
      integer, intent(in) :: n, incx, incy
      double precision, intent(in), dimension(*) :: x
      double precision, intent(in), dimension(*) :: y
    end function ddot
    double precision function dnrm2( n, x, incx )
      implicit none
      integer, intent(in) :: n, incx
      double precision, intent(in), dimension(*) :: x
    end function dnrm2
    subroutine dscal( n, a, x, incx )
      implicit none
      integer, intent(in) :: n, incx
      double precision, intent(in) :: a
      double precision, intent(in), dimension(*) :: x
    end subroutine dscal
  end interface

  ! Level 3 BLAS
  interface
    subroutine dgemm( ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
      implicit none
      character, intent(in) :: ta, tb
      integer, intent(in) :: m, n, k
      integer, intent(in) :: lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in   ), dimension(lda, *) :: a
      double precision, intent(in   ), dimension(ldb, *) :: b
      double precision, intent(inout), dimension(ldc, *) :: c
    end subroutine dgemm
    subroutine dtrsm( side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb )
      implicit none
      character, intent(in) :: side, uplo, trans, diag
      integer, intent(in) :: m, n, lda, ldb
      double precision, intent(in   ) :: alpha
      double precision, intent(in   ) :: a(lda, *)
      double precision, intent(inout) :: b(ldb, n)
    end subroutine dtrsm
  end interface

end module spral_blas_iface
