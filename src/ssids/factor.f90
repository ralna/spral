!
! This module defines ssids_fkeep type and associated procedures (CPU version)
!
module spral_ssids_factor
   use spral_ssids_datatypes, only : wp, long, smalloc_type, node_type
   use spral_ssids_cuda_datatypes, only : gpu_type
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: ssids_fkeep

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for data generated in factorise phase
   !
   type ssids_fkeep_base
      integer :: flag ! copy of error flag.
      real(wp), dimension(:), allocatable :: scaling ! Stores scaling for
         ! each entry (in original matrix order)
      type(node_type), dimension(:), allocatable :: nodes ! Stores pointers
         ! to information about nodes
      type(smalloc_type), pointer :: alloc=>null() ! Linked list of memory pages
         ! pointed to by nodes variable
      logical :: pos_def ! set to true if user indicates matrix pos. definite

      ! Info components to copy on solve
      integer :: matrix_rank
      integer :: maxfront
      integer :: num_delay
      integer(long) :: num_factor
      integer(long) :: num_flops
      integer :: num_neg
      integer :: num_two
   end type ssids_fkeep_base

   !
   ! Data type for data generated in factorise phase (gpu version)
   !
   type, extends(ssids_fkeep_base) :: ssids_fkeep
      type(C_PTR), dimension(:), allocatable :: stream_handle
      type(gpu_type), dimension(:), allocatable :: stream_data
      type(gpu_type) :: top_data
      type(C_PTR) :: gpu_rlist_with_delays = C_NULL_PTR
      type(C_PTR) :: gpu_rlist_direct_with_delays = C_NULL_PTR
      type(C_PTR) :: gpu_clists = C_NULL_PTR
      type(C_PTR) :: gpu_clists_direct = C_NULL_PTR
      type(C_PTR) :: gpu_clen = C_NULL_PTR
      logical :: host_factors = .false.
   end type ssids_fkeep

end module spral_ssids_factor
