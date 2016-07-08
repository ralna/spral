module spral_ssids_akeep
   use spral_ssids_datatypes, only : long, wp, SSIDS_ERROR_CUDA_UNKNOWN, &
                                     ssids_options
   use spral_ssids_inform, only : ssids_inform_base
   use spral_ssids_subtree, only : symbolic_subtree_base
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: ssids_akeep_base

   type symbolic_subtree_ptr
      class(symbolic_subtree_base), pointer :: ptr
   end type symbolic_subtree_ptr

   !
   ! Data type for information generated in analyse phase
   !
   type ssids_akeep_base
      logical :: check ! copy of check as input to analyse phase
      integer :: flag ! copy of error flag.
      integer :: maxmn ! maximum value of blkm or blkn
      integer :: n ! Dimension of matrix
      integer :: ne ! Set to number of entries input by user.
      integer(long) :: nfactor 
      integer :: nnodes = -1 ! Number of nodes in assembly tree
      integer :: num_two ! This is set to 0 as we ignore any negative signs
         ! that indicate 2x2 pivot (analyse does not exploit them)

      ! Subtree partition
      integer :: nparts
      integer, dimension(:), allocatable :: part
      type(symbolic_subtree_ptr), dimension(:), allocatable :: subtree

      ! child_list(child_ptr(node):child_ptr(node+1)-1) is list of children
      ! of node. Used to ensure we always sum contributions from children
      ! in the same order
      integer, dimension(:), allocatable :: child_ptr 
      integer, dimension(:), allocatable :: child_list
      integer(C_INT), dimension(:), allocatable :: invp ! inverse of pivot order
         ! that is passed to factorize phase
      integer, dimension(:), allocatable :: level ! level(i) of the assembly
         ! tree at which node i sits. root is level 1.
      integer, dimension(:,:), allocatable :: nlist ! map from A to factors
         ! For nodes i, the entries nlist(1:2, nptr(i):nptr(i+1)-1) define
         ! a relationship:
         ! nodes(node)%lcol( nlist(2,j) ) = val( nlist(1,j) )
     integer, dimension(:), allocatable :: nptr ! Entries into nlist for
         ! nodes of the assembly tree. Has length nnodes+1
      integer, dimension(:), allocatable :: rlist ! rlist(rptr(i):rptr(i+1)-1)
         ! contains the row indices for node i of the assembly tree. 
         ! At each node, the list
         ! is in elimination order. Allocated within mc78_analyse.
      integer(long), dimension(:), allocatable :: rptr ! Pointers into rlist
         ! for nodes of assembly tree. Has length nnodes+1. 
         ! Allocated within mc78_analyse.
      integer, dimension(:), allocatable :: sparent ! sparent(i) is parent
         ! of node i in assembly tree. sparent(i)=nnodes+1 if i is a root.
         ! The parent is always numbered higher than each of its children.
         ! Allocated within mc78_analyse.
      integer, dimension(:), allocatable :: sptr ! (super)node pointers.
         ! Supernode i consists of sptr(i) through sptr(i+1)-1.
         ! Allocated within mc78_analyse.
      integer(long), dimension(:), allocatable :: subtree_work ! For each node,
         ! the number of flops involved in the subtree rooted at that node.

      integer, dimension(:), allocatable :: rlist_direct ! List of rows
         ! for direct addressing of children (assuming no delays!)

      ! Following components are for cleaned up matrix data.
      ! LOWER triangle only. We have to retain these for factorize phase
      ! as used if the user wants to do scaling.
      ! These components are NOT used if check is set to .false.
      ! on call to ssids_analyse.
      integer, allocatable :: ptr(:) ! column pointers
      integer, allocatable :: row(:)! row indices
      integer :: lmap ! used by hsl_mc69
      integer, allocatable :: map(:) ! used by hsl_mc69

      ! Following components are cached members of inform
      integer :: matrix_dup
      integer :: matrix_outrange
      integer :: matrix_missing_diag
      integer :: maxdepth
      integer(long) :: num_flops ! not copied to inform in factor, but used to
         ! determine if parallelism should be used
      integer :: num_sup

      ! Scaling from matching-based ordering
      real(wp), dimension(:), allocatable :: scaling
   contains
      procedure, pass(akeep) :: free => free_akeep_base
      procedure, pass(akeep) :: move_data
   end type ssids_akeep_base

contains

subroutine move_data(akeep, options, inform)
   class(ssids_akeep_base), intent(inout) :: akeep
   type(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform

   ! No-op
end subroutine move_data

subroutine free_akeep_base(akeep, flag)
   class(ssids_akeep_base), intent(inout) :: akeep
   integer, intent(out) :: flag

   integer :: st

   flag = 0

   deallocate(akeep%part, stat=st)
   deallocate(akeep%subtree, stat=st)
   deallocate(akeep%child_ptr, stat=st)
   deallocate(akeep%child_list, stat=st)
   deallocate(akeep%invp, stat=st)
   deallocate(akeep%level, stat=st)
   deallocate(akeep%nlist, stat=st)
   deallocate(akeep%nptr, stat=st)
   deallocate(akeep%rlist, stat=st)
   deallocate(akeep%rlist_direct, stat=st)
   deallocate(akeep%rptr, stat=st)
   deallocate(akeep%sparent, stat=st)
   deallocate(akeep%sptr, stat=st)
   deallocate(akeep%subtree_work, stat=st)
   deallocate(akeep%ptr, stat=st)
   deallocate(akeep%row, stat=st)
   deallocate(akeep%map, stat=st)
end subroutine free_akeep_base

end module spral_ssids_akeep
