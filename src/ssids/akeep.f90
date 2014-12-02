module spral_ssids_akeep
   use spral_ssids_datatypes, only : long, wp
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: ssids_akeep

   !
   ! Data type for information generated in analyse phase
   !
   type ssids_akeep
      logical :: check ! copy of check as input to analyse phase
      integer :: flag ! copy of error flag.
      integer :: maxmn ! maximum value of blkm or blkn
      integer :: n ! Dimension of matrix
      integer :: ne ! Set to number of entries input by user.
      integer(long) :: nfactor 
      integer :: nnodes = -1 ! Number of nodes in assembly tree
      integer :: num_two ! This is set to 0 as we ignore any negative signs
         ! that indicate 2x2 pivot (analyse does not exploit them)

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

      ! GPU pointers (copies of the CPU arrays of same name)
      type(C_PTR) :: gpu_nlist = C_NULL_PTR
      type(C_PTR) :: gpu_rlist = C_NULL_PTR
      type(C_PTR) :: gpu_rlist_direct = C_NULL_PTR
   end type ssids_akeep

end module spral_ssids_akeep
