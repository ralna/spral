module spral_ssids_cpu_subtree
   use, intrinsic :: iso_c_binding
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_cpu_iface ! fixme only
   use spral_ssids_cpu_solve ! fixme only
   use spral_ssids_datatypes, only : long, wp, ssids_options, &
      SSIDS_ERROR_ALLOCATION
   use spral_ssids_subtree, only : symbolic_subtree_base, numeric_subtree_base
   use spral_ssids_inform, only : ssids_inform_base
   implicit none

   private
   public :: cpu_symbolic_subtree, construct_cpu_symbolic_subtree
   public :: cpu_numeric_subtree

   type, extends(symbolic_subtree_base) :: cpu_symbolic_subtree
      integer :: n
      integer :: nnodes
      integer, dimension(:), pointer :: nptr
      integer, dimension(:,:), pointer :: nlist
      type(C_PTR) :: csubtree
   contains
      procedure :: factor
      procedure :: factor2 ! fixme remove
      final :: symbolic_final
   end type cpu_symbolic_subtree

   type, extends(numeric_subtree_base) :: cpu_numeric_subtree
      type(cpu_symbolic_subtree), pointer :: symbolic
      logical(C_BOOL) :: posdef
      type(cpu_node_data), dimension(:), allocatable :: cnodes
      type(cpu_factor_options) :: coptions ! FIXME: doesn't belong here
      type(contrib_type), pointer :: contrib
      type(C_PTR) :: csubtree
   contains
      procedure :: get_contrib
      procedure :: solve_fwd
      procedure :: solve_diag
      procedure :: solve_bwd
      final :: numeric_final
   end type cpu_numeric_subtree

   interface
      type(C_PTR) function c_create_symbolic_subtree(nnodes, sptr, sparent, &
            rptr, rlist)&
            bind(C, name="spral_ssids_cpu_create_symbolic_subtree")
         use, intrinsic :: iso_c_binding
         implicit none
         integer(C_INT), value :: nnodes
         integer(C_INT), dimension(nnodes+1), intent(in) :: sptr
         integer(C_INT), dimension(nnodes+1), intent(in) :: sparent
         integer(C_LONG), dimension(nnodes+1), intent(in) :: rptr
         integer(C_INT), dimension(*), intent(in) :: rlist
      end function c_create_symbolic_subtree

      subroutine c_destroy_symbolic_subtree(subtree) &
            bind(C, name="spral_ssids_cpu_destroy_symbolic_subtree")
         use, intrinsic :: iso_c_binding
         implicit none
         type(C_PTR), value :: subtree
      end subroutine c_destroy_symbolic_subtree
      subroutine c_factor_cpu(pos_def, subtree, n, nnodes, nodes, aval, &
            scaling, alloc, options, stats) &
            bind(C, name="spral_ssids_factor_cpu_dbl")
         use, intrinsic :: iso_c_binding
         import :: cpu_node_data, cpu_factor_options, cpu_factor_stats
         implicit none
         logical(C_BOOL), value :: pos_def
         type(C_PTR), value :: subtree
         integer(C_INT), value :: n
         integer(C_INT), value :: nnodes
         type(cpu_node_data), dimension(nnodes), intent(inout) :: nodes
         real(C_DOUBLE), dimension(*), intent(in) :: aval
         type(C_PTR), value :: scaling
         type(C_PTR), value :: alloc
         type(cpu_factor_options), intent(in) :: options
         type(cpu_factor_stats), intent(out) :: stats
      end subroutine c_factor_cpu

      type(C_PTR) function c_create_cpu_subtree(posdef, symbolic_subtree, &
            nnodes, nodes) &
            bind(C, name="spral_ssids_create_cpu_subtree_dbl")
         use, intrinsic :: iso_c_binding
         import :: cpu_node_data
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: symbolic_subtree
         integer(C_INT), value :: nnodes
         type(cpu_node_data), dimension(nnodes), intent(inout) :: nodes
      end function c_create_cpu_subtree

      subroutine c_destroy_cpu_subtree(posdef, subtree) &
            bind(C, name="spral_ssids_destroy_cpu_subtree_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
      end subroutine c_destroy_cpu_subtree
   end interface

contains

function construct_cpu_symbolic_subtree(n, nnodes, sptr, sparent, rptr, &
      rlist, nptr, nlist) result(this)
   class(cpu_symbolic_subtree), pointer :: this
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(nnodes+1), target, intent(in) :: nptr
   integer, dimension(2,nptr(nnodes+1)-1), target, intent(in) :: nlist

   integer :: st

   allocate(this, stat=st)
   if(st.ne.0) then
      nullify(this)
      return
   endif
   this%n = n
   this%nnodes = nnodes
   this%nptr => nptr
   this%nlist => nlist
   this%csubtree = &
      c_create_symbolic_subtree(nnodes, sptr, sparent, rptr, rlist)
end function construct_cpu_symbolic_subtree

subroutine symbolic_final(this)
   type(cpu_symbolic_subtree) :: this

   call c_destroy_symbolic_subtree(this%csubtree)
end subroutine symbolic_final

function factor2(this, posdef, aval, options, inform, alloc, cstats, scaling)
   type(cpu_numeric_subtree), pointer :: factor2
   class(cpu_symbolic_subtree), target, intent(inout) :: this
   logical(C_BOOL), intent(in) :: posdef
   real(wp), dimension(*), intent(in) :: aval
   class(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform
   type(C_PTR), intent(in) :: alloc
   type(cpu_factor_stats), intent(out) :: cstats
   real(wp), dimension(*), target, optional, intent(in) :: scaling

   type(cpu_numeric_subtree), pointer :: cpu_factor
   type(C_PTR) :: cscaling
   integer :: st

   ! Setup output
   allocate(cpu_factor, stat=st)
   if(st.ne.0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      inform%stat = st
      nullify(factor2)
      return
   endif
   factor2 => cpu_factor
   cpu_factor%symbolic => this

   ! Allocate cnodes and setup for main call
   allocate(cpu_factor%cnodes(this%nnodes+1), stat=st)
   if(st.ne.0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      inform%stat = st
      deallocate(cpu_factor, stat=st)
      nullify(factor2)
      return
   endif
   call setup_cpu_data(this%nnodes, cpu_factor%cnodes, this%nptr, this%nlist)
   call cpu_copy_options_in(options, cpu_factor%coptions)

   cpu_factor%posdef = posdef
   cpu_factor%csubtree = &
      c_create_cpu_subtree(posdef, this%csubtree, this%nnodes, cpu_factor%cnodes)

   cscaling = C_NULL_PTR
   if(present(scaling)) cscaling = C_LOC(scaling)

   call c_factor_cpu(posdef, cpu_factor%csubtree, this%n, &
      this%nnodes, cpu_factor%cnodes, aval, cscaling, alloc, &
      cpu_factor%coptions, cstats)
end function factor2

function factor(this, posdef, aval, options, inform, scaling)
   class(numeric_subtree_base), pointer :: factor
   class(cpu_symbolic_subtree), target, intent(inout) :: this
   logical(C_BOOL), intent(in) :: posdef
   real(wp), dimension(*), intent(in) :: aval
   class(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform
   real(wp), dimension(*), optional, intent(in) :: scaling

   type(cpu_numeric_subtree), pointer :: cpu_factor
   integer :: st

   ! Setup output
   allocate(cpu_factor, stat=st)
   if(st.ne.0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      inform%stat = st
      nullify(factor)
      return
   endif
   factor => cpu_factor
   cpu_factor%symbolic => this

   ! Allocate cnodes and setup for main call
   allocate(cpu_factor%cnodes(this%nnodes+1), stat=st)
   if(st.ne.0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      inform%stat = st
      deallocate(cpu_factor, stat=st)
      nullify(factor)
      return
   endif
   call setup_cpu_data(this%nnodes, cpu_factor%cnodes, this%nptr, this%nlist)
   call cpu_copy_options_in(options, cpu_factor%coptions)

   ! Allocate subtree
   cpu_factor%csubtree = &
      c_create_cpu_subtree(posdef, this%csubtree, this%nnodes, cpu_factor%cnodes)

   ! Perform actual factorization
   !call cpu_factor%factor3(this%n, this%nnodes, cpu_factor%cnodes, aval, alloc, coptions, cstats, scaling)

   ! Gather information back to Fortran
   !call extract_cpu_data(this%nnodes, cnodes, fkeep%nodes, cstats, inform)
end function factor

subroutine numeric_final(this)
   type(cpu_numeric_subtree) :: this

   call c_destroy_cpu_subtree(this%posdef, this%csubtree)
end subroutine numeric_final

function get_contrib(this)
   class(contrib_type), pointer :: get_contrib
   class(cpu_numeric_subtree), intent(in) :: this

   get_contrib => this%contrib
end function get_contrib

subroutine solve_fwd(this, nrhs, x, ldx)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx

   ! FIXME
end subroutine solve_fwd

subroutine solve_diag(this, nrhs, x, ldx)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx

   ! FIXME
end subroutine solve_diag

subroutine solve_bwd(this, nrhs, x, ldx)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx

   ! FIXME
end subroutine solve_bwd

end module spral_ssids_cpu_subtree
