module cuda_helper
  use spral_cuda
  use, intrinsic :: iso_c_binding
  implicit none

contains

  ! Return number of CUDA devices available, or < 0 for error
  subroutine cuda_init(cnt)
    implicit none

    integer(C_INT), intent(out) :: cnt

    integer(C_INT) :: cuda_error

    cuda_error = cudaGetDeviceCount(cnt)
    select case (cuda_error)
    case (cudaSuccess)
       continue
    case (cudaErrorInsufficientDriver)
       cnt = -1_C_INT
    case (cudaErrorNoDevice)
       cnt = 0_C_INT
    case default
       cnt = -2_C_INT
    end select
  end subroutine cuda_init
end module cuda_helper
