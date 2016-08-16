! Provides limited interface definitions for CUDA functions in the case
! we are not compiled against CUDA libraries
module spral_cuda
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: cudaGetErrorString
   public :: detect_gpu

contains
! Convert a CUDA error code to a Fortran character string
character(len=200) function cudaGetErrorString(error)
   integer(C_INT) :: error

   cudaGetErrorString = "Not compiled with CUDA support"
end function cudaGetErrorString


! Return true if a GPU is present and code is compiled with CUDA support
logical function detect_gpu()
   detect_gpu = .false.
end function detect_gpu

end module spral_cuda
