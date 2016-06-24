module ldltf_double

  interface num_neg_D
    module procedure num_neg_D_double
  end interface

contains

  integer function num_neg_D_double( n, LDLT, ld, ipiv )
! Counts negative eigenvalues of factor D

    implicit none

    double precision, parameter :: ZERO = 0.0D0
    integer, intent(in) :: n, ld
    double precision, intent(in) :: LDLT(ld, n)
    integer, intent(in) :: ipiv(n)
    integer :: nneg
    integer :: i
    double precision :: r, s, t

    i = 1
    nneg = 0
    do while ( i <= n )
      s = LDLT(i, i)
      if ( ipiv(i) < 0 ) then
        t = LDLT(i + 1, i)
        r = LDLT(i + 1, i + 1)
        if ( s*r - t*t < ZERO ) then
          nneg = nneg + 1
        else if ( s*r - t*t > ZERO .and. s + r < ZERO ) then
          nneg = nneg + 2
        end if
        i = i + 2
      else
        if ( s < ZERO ) nneg = nneg + 1
        i = i + 1
      end if
    end do
    num_neg_D_double = nneg

  end function num_neg_D_double

end module ldltf_double

module ldltf
  use ldltf_double
end module ldltf
