! Definition of LAPACK API in module
module spral_lapack_iface
  implicit none

  private
  public :: dpotrf, dlacpy

  interface
    subroutine dpotrf( uplo, n, a, lda, info )
      implicit none
      character, intent(in) :: uplo
      integer, intent(in) :: n, lda
      double precision, intent(inout) :: a(lda, n)
      integer, intent(out) :: info
    end subroutine dpotrf
    subroutine dlacpy( uplo, m, n, a, lda, b, ldb )
      implicit none
      character, intent(in) :: uplo
      integer, intent(in) :: m, n, lda, ldb
      double precision, intent(in ) :: a(lda, n)
      double precision, intent(out) :: b(ldb, n)
    end subroutine dlacpy
  end interface

end module spral_lapack_iface
