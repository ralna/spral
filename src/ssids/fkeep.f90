!
! This module defines ssids_fkeep type and associated procedures (CPU version)
!
module spral_ssids_fkeep
   use spral_ssids_alloc, only : smfreeall
   use spral_ssids_datatypes, only : long, node_type, smalloc_type, &
                                     ssids_akeep, ssids_options, ssids_inform, &
                                     wp, SSIDS_ERROR_ALLOCATION
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: ssids_fkeep

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for data generated in factorise phase
   !
   type ssids_fkeep
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

   contains
      procedure, pass(fkeep) :: inner_factor => inner_factor_cpu ! Do actual factorization
      procedure, pass(fkeep) :: free => free_fkeep ! Frees memory
   end type ssids_fkeep

contains

subroutine inner_factor_cpu(fkeep, akeep, val, options, inform)
   class(ssids_akeep), intent(in) :: akeep
   class(ssids_fkeep), intent(inout) :: fkeep
   real(wp), dimension(*), target, intent(in) :: val
   class(ssids_options), intent(in) :: options
   class(ssids_inform), intent(inout) :: inform

   print *, "CPU Factor not implemented yet!"
   stop
end subroutine inner_factor_cpu

subroutine free_fkeep(fkeep, flag)
   class(ssids_fkeep), intent(inout) :: fkeep
   integer, intent(out) :: flag ! not actually used for cpu version, set to 0

   integer :: st, i

   flag = 0 ! Not used for basic SSIDS, just zet to zero

   ! Skip if nothing intialized
   if (.not.allocated(fkeep%nodes)) return

   ! Free memory
   call smfreeall(fkeep%alloc)
   deallocate(fkeep%alloc)
   nullify(fkeep%alloc)

   deallocate(fkeep%nodes, stat=st)
   deallocate(fkeep%scaling, stat=st)

end subroutine free_fkeep

end module spral_ssids_fkeep
