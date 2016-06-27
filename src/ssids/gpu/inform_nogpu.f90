! Version of module providing minimum functionality in case of no GPU support
module spral_ssids_gpu_inform
   use spral_ssids_inform, only : ssids_inform_base
   implicit none

   private
   public :: ssids_inform_gpu

   type, extends(ssids_inform_base) :: ssids_inform_gpu
      integer :: cuda_error = 0 ! cuda error value
      integer :: cublas_error = 0 ! cublas error value
   end type ssids_inform_gpu
end module spral_ssids_gpu_inform
