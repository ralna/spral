module cuda_helper
   use spral_cuda
   use, intrinsic :: iso_c_binding
   implicit none

contains

   ! Do some pointless CUDA operation to force an initialization
   subroutine cuda_init()
      type(C_PTR) :: ptr, cublas
      integer :: cuda_error

      cuda_error = cudaMalloc(ptr, 1_C_SIZE_T)
      cuda_error = cudaFree(ptr)

      cuda_error = cublasCreate(cublas)
      cuda_error = cublasDestroy(cublas)

      cuda_error = cudaDeviceSynchronize()
   end subroutine cuda_init
end module cuda_helper
