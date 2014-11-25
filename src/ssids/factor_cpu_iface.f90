module factor_cpu_iface
   use, intrinsic :: iso_c_binding
   use spral_ssids_datatypes, only : node_type
   implicit none

   ! See comments in C++ definition in factor_gpu.cxx for detail
   type, bind(C) :: cpu_node_data_double
      ! Fixed data from analyse
      integer(C_INT) :: nrow_expected
      integer(C_INT) :: ncol_expected
      type(C_PTR) :: first_child
      type(C_PTR) :: next_child
      type(C_PTR) :: rlist

      ! Data about A
      integer(C_INT) :: num_a
      type(C_PTR) :: amap
      type(C_PTR) :: aval

      ! Data that changes during factorize
      integer(C_INT) :: ndelay_in
      integer(C_INT) :: ndelay_out
      integer(C_INT) :: nelim
      type(C_PTR) :: lcol
      type(C_PTR) :: perm
      type(C_PTR) :: contrib
   end type cpu_node_data_double

   type, bind(C) :: cpu_factor_stats
      integer(C_INT) :: flag
   end type cpu_factor_stats

   interface
      subroutine factor_cpu_double(pos_def, nnodes, nodes, val, scaling, &
            alloc, stats) &
            bind(C, name="spral_ssids_factor_cpu_dbl")
         use, intrinsic :: iso_c_binding
         import :: cpu_node_data_double, cpu_factor_stats
         implicit none
         logical(C_BOOL), value :: pos_def
         integer(C_INT), value :: nnodes
         type(cpu_node_data_double), dimension(nnodes), intent(inout) :: nodes
         real(C_DOUBLE), dimension(*), intent(in) :: val
         real(C_DOUBLE), dimension(*), intent(in) :: scaling
         type(C_PTR), value :: alloc
         type(cpu_factor_stats), intent(out) :: stats
      end subroutine factor_cpu_double
   end interface

contains

subroutine setup_cpu_data(nnodes, fnodes, cnodes)
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(in) :: fnodes
   type(cpu_node_data_double), dimension(nnodes), intent(out) :: cnodes

end subroutine setup_cpu_data

end module factor_cpu_iface

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

   call c_f_pointer(calloc, falloc)
   call smalloc(falloc, ptr, len, srcptr, srchead, st)
   if(st.ne.0) then
      spral_ssids_smalloc_int = C_NULL_PTR
   else
      spral_ssids_smalloc_int = C_LOC(srcptr%imem(srchead))
   endif
end function spral_ssids_smalloc_int
