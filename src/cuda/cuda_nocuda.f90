! Provides limited interface definitions for CUDA functions in the case
! we are not compiled against CUDA libraries
module spral_cuda
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: detect_gpu

contains

! Return true if a GPU is present and code is compiled with CUDA support
logical function detect_gpu()
   detect_gpu = .false.
end function detect_gpu

end module spral_cuda
