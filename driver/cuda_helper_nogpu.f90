module cuda_helper
  use, intrinsic :: iso_c_binding
  implicit none

contains

  subroutine cuda_init(cnt)
    implicit none
    integer(C_INT), intent(out) :: cnt
    cnt = 0_C_INT
  end subroutine cuda_init
end module cuda_helper
