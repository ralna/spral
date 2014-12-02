module spral_ssids_inform_gpu
   use spral_cuda, only : cudaGetErrorString
   use spral_ssids_datatypes, only : SSIDS_ERROR_CUBLAS_UNKNOWN, &
                                     SSIDS_ERROR_CUDA_UNKNOWN
   use spral_ssids_inform, only : ssids_inform_base
   implicit none

   private
   public :: ssids_inform_gpu

   type, extends(ssids_inform_base) :: ssids_inform_gpu
      integer :: cuda_error = 0 ! cuda error value
      integer :: cublas_error = 0 ! cublas error value
   contains
      procedure, pass(this) :: flagToCharacter => flagToCharacter_gpu
      procedure, pass(this) :: set_cuda_error
   end type ssids_inform_gpu

contains

!
! Returns a string representation (gpu extension)
! Member function inform%flagToCharacter
!
function flagToCharacter_gpu(this) result(msg)
   class(ssids_inform_gpu), intent(in) :: this
   character(len=200) :: msg ! return value

   select case(this%flag)
   !
   ! Errors
   !
   case(SSIDS_ERROR_CUDA_UNKNOWN)
      write(msg,'(2a)') ' Unhandled CUDA error: ', &
         cudaGetErrorString(this%cuda_error)
   case(SSIDS_ERROR_CUBLAS_UNKNOWN)
      msg = 'Unhandled CUBLAS error:'
      ! FIXME?
   !
   ! Fall back to base case
   !
   case default
      msg = this%ssids_inform_base%flagToCharacter()
   end select

end function flagToCharacter_gpu

subroutine set_cuda_error(this, cuda_error)
   class(ssids_inform_gpu), intent(inout) :: this
   integer, intent(in) :: cuda_error

   this%flag = SSIDS_ERROR_CUDA_UNKNOWN
   this%cuda_error = cuda_error
end subroutine


end module spral_ssids_inform_gpu
