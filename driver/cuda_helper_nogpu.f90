module cuda_helper
   use, intrinsic :: iso_c_binding
   implicit none

contains

   ! Do some pointless CUDA operation to force an initialization
   subroutine cuda_init()
      ! Noop
   end subroutine cuda_init
end module cuda_helper
