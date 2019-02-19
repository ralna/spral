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
    type(C_PTR) :: ptr, cublas

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

    ! Do some pointless CUDA operation to force an initialization
    cuda_error = cudaMalloc(ptr, 1_C_SIZE_T)
    if (cuda_error .ne. cudaSuccess) then
       write(*, "(a)")"[error][cuda_init] Failed to allocate memory on the CUDA device"
       cnt = -2_C_INT
       return
    end if
    cuda_error = cudaFree(ptr)

    cuda_error = cublasCreate(cublas)
    cuda_error = cublasDestroy(cublas)

    cuda_error = cudaDeviceSynchronize()

  end subroutine cuda_init
end module cuda_helper
