module spral_ssids_cpu_factor_iface
   use, intrinsic :: iso_c_binding
   use spral_ssids_akeep, only : ssids_akeep_base
   use spral_ssids_datatypes, only : ssids_options, node_type, long, &
                                     DEBUG_PRINT_LEVEL
   use spral_ssids_inform, only : ssids_inform_base
   implicit none

   private
   public :: cpu_node_data, cpu_factor_options, cpu_factor_stats
   public :: factor_cpu, setup_cpu_data, extract_cpu_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! See comments in C++ definition in factor_gpu.cxx for detail
   type, bind(C) :: cpu_node_data
      ! Fixed data from analyse
      integer(C_INT) :: nrow_expected
      integer(C_INT) :: ncol_expected
      type(C_PTR) :: first_child
      type(C_PTR) :: next_child
      type(C_PTR) :: rlist
      logical(C_BOOL) :: even

      ! Data about A
      integer(C_INT) :: num_a
      type(C_PTR) :: amap

      ! Data that changes during factorize
      integer(C_INT) :: ndelay_in
      integer(C_INT) :: ndelay_out
      integer(C_INT) :: nelim
      type(C_PTR) :: lcol
      type(C_PTR) :: perm
      type(C_PTR) :: contrib
   end type cpu_node_data

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type, bind(C) :: cpu_factor_options
      real(C_DOUBLE) :: small
      real(C_DOUBLE) :: u
      integer(C_INT) :: print_level
      integer(C_INT) :: cpu_task_block_size
      integer(C_INT) :: cpu_small_subtree_threshold
   end type cpu_factor_options

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type, bind(C) :: cpu_factor_stats
      integer(C_INT) :: flag
      integer(C_INT) :: num_delay
      integer(C_INT) :: num_neg
      integer(C_INT) :: num_two
      integer(C_INT) :: num_zero
      integer(C_INT) :: maxfront
      integer(C_INT) :: elim_at_pass(5)
      integer(C_INT) :: elim_at_itr(5)
   end type cpu_factor_stats

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface
      subroutine factor_cpu(pos_def, n, nnodes, nodes, aval, scaling, &
            alloc, options, stats) &
            bind(C, name="spral_ssids_factor_cpu_dbl")
         use, intrinsic :: iso_c_binding
         import :: cpu_node_data, cpu_factor_options, cpu_factor_stats
         implicit none
         logical(C_BOOL), value :: pos_def
         integer(C_INT), value :: n
         integer(C_INT), value :: nnodes
         type(cpu_node_data), dimension(nnodes), intent(inout) :: nodes
         real(C_DOUBLE), dimension(*), intent(in) :: aval
         type(C_PTR), value :: scaling
         type(C_PTR), value :: alloc
         type(cpu_factor_options), intent(in) :: options
         type(cpu_factor_stats), intent(out) :: stats
      end subroutine factor_cpu
   end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setup_cpu_data(akeep, cnodes, foptions, coptions)
   class(ssids_akeep_base), target, intent(in) :: akeep
   type(cpu_node_data), dimension(akeep%nnodes+1), target, intent(out) :: cnodes
   type(ssids_options), intent(in) :: foptions
   type(cpu_factor_options), intent(out) :: coptions

   integer :: node, parent

   !
   ! Setup node data
   !
   ! Basic data and initialize linked lists
   do node = 1, akeep%nnodes
      ! Data about factors
      cnodes(node)%nrow_expected = int(akeep%rptr(node+1) - akeep%rptr(node))
      cnodes(node)%ncol_expected = akeep%sptr(node+1) - akeep%sptr(node)
      cnodes(node)%first_child = C_NULL_PTR
      cnodes(node)%next_child = C_NULL_PTR
      cnodes(node)%rlist = C_LOC(akeep%rlist(akeep%rptr(node)))

      ! Print out debug info
      if(foptions%print_level > DEBUG_PRINT_LEVEL) then
         print *, "node ", node, " parent ", akeep%sparent(node), &
            int(akeep%rptr(node+1) - akeep%rptr(node)), "x", &
            akeep%sptr(node+1) - akeep%sptr(node)
      endif

      ! Data about A
      cnodes(node)%num_a = akeep%nptr(node+1) - akeep%nptr(node)
      cnodes(node)%amap = C_LOC(akeep%nlist(1,akeep%nptr(node)))
   end do
   cnodes(akeep%nnodes+1)%first_child = C_NULL_PTR
   ! Build linked lists of children
   do node = 1, akeep%nnodes
      parent = akeep%sparent(node)
      cnodes(node)%next_child = cnodes(parent)%first_child
      cnodes(parent)%first_child = C_LOC( cnodes(node) )
   end do
   ! Setup odd/even distance from root
   do node = akeep%nnodes, 1, -1
      parent = akeep%sparent(node)
      if(parent > akeep%nnodes) then
         cnodes(node)%even = .true.
      else
         cnodes(node)%even = .not.cnodes(parent)%even
      endif
   end do

   !
   ! Setup options
   !
   coptions%small       = foptions%small
   coptions%u           = foptions%u
   coptions%print_level = foptions%print_level
   coptions%cpu_small_subtree_threshold = foptions%cpu_small_subtree_threshold
   coptions%cpu_task_block_size         = foptions%cpu_task_block_size
end subroutine setup_cpu_data

subroutine extract_cpu_data(nnodes, cnodes, fnodes, cstats, finform)
   integer, intent(in) :: nnodes
   type(cpu_node_data), dimension(nnodes+1), intent(in) :: cnodes
   type(node_type), dimension(nnodes+1), intent(out) :: fnodes
   type(cpu_factor_stats), intent(in) :: cstats
   class(ssids_inform_base), intent(inout) :: finform

   integer :: node, nrow, ncol
   integer :: rank

   ! Copy factors (scalars and pointers)
   rank = 0
   do node = 1, nnodes
      ! Components we copy from C version
      rank = rank + cnodes(node)%ncol_expected
      nrow = cnodes(node)%nrow_expected + cnodes(node)%ndelay_in
      ncol = cnodes(node)%ncol_expected + cnodes(node)%ndelay_in
      fnodes(node)%nelim = cnodes(node)%nelim
      fnodes(node)%ndelay = cnodes(node)%ndelay_in
      call C_F_POINTER(cnodes(node)%lcol, fnodes(node)%lcol, &
         shape=(/ (2_long+nrow)*ncol /) )
      call C_F_POINTER(cnodes(node)%perm, fnodes(node)%perm, &
         shape=(/ ncol /) )

      ! Components we abdicate
      fnodes(node)%rdptr = -1
      fnodes(node)%ncpdb = -1
      fnodes(node)%gpu_lcol = C_NULL_PTR
      nullify(fnodes(node)%rsmptr)
      nullify(fnodes(node)%ismptr)
      fnodes(node)%rsmsa = -1
      fnodes(node)%ismsa = -1
   end do

   ! Copy stats
   finform%flag         = cstats%flag
   finform%num_delay    = cstats%num_delay
   finform%num_neg      = cstats%num_neg
   finform%num_two      = cstats%num_two
   finform%matrix_rank  = rank - cstats%num_zero
   finform%maxfront     = cstats%maxfront
   !print *, "Elim at (pass) = ", cstats%elim_at_pass(:)
   !print *, "Elim at (itr) = ", cstats%elim_at_itr(:)
end subroutine extract_cpu_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module spral_ssids_cpu_factor_iface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Provide a way to alloc memory using smalloc (double version)
type(C_PTR) function spral_ssids_smalloc_dbl(calloc, len) bind(C)
   use, intrinsic :: iso_c_binding
   use spral_ssids_datatypes, only : long, smalloc_type
   use spral_ssids_alloc, only : smalloc
   implicit none
   type(C_PTR), value :: calloc
   integer(C_SIZE_T), value :: len

   type(smalloc_type), pointer :: falloc, srcptr
   real(C_DOUBLE), dimension(:), pointer :: ptr
   integer(long) :: srchead
   integer :: st

   call c_f_pointer(calloc, falloc)
   call smalloc(falloc, ptr, len, srcptr, srchead, st)
   if(st.ne.0) then
      spral_ssids_smalloc_dbl = C_NULL_PTR
   else
      spral_ssids_smalloc_dbl = C_LOC(srcptr%rmem(srchead))
   endif
end function spral_ssids_smalloc_dbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Provide a way to alloc memory using smalloc (int version)
type(C_PTR) function spral_ssids_smalloc_int(calloc, len) bind(C)
   use, intrinsic :: iso_c_binding
   use spral_ssids_datatypes, only : long, smalloc_type
   use spral_ssids_alloc, only : smalloc
   implicit none
   type(C_PTR), value :: calloc
   integer(C_SIZE_T), value :: len

   type(smalloc_type), pointer :: falloc, srcptr
   integer(C_INT), dimension(:), pointer :: ptr
   integer(long) :: srchead
   integer :: st

   if(len.lt.0) then
      spral_ssids_smalloc_int = C_NULL_PTR
      return
   endif

   call c_f_pointer(calloc, falloc)
   call smalloc(falloc, ptr, len, srcptr, srchead, st)
   if(st.ne.0) then
      spral_ssids_smalloc_int = C_NULL_PTR
   else
      spral_ssids_smalloc_int = C_LOC(srcptr%imem(srchead))
   endif
end function spral_ssids_smalloc_int
