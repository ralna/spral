module spral_ssids_cpu_subtree
   use, intrinsic :: iso_c_binding
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_cpu_iface
   use spral_ssids_datatypes, only : long, wp
   use spral_ssids_subtree, only : symbolic_subtree_base, numeric_subtree_base
   implicit none

   private
   public :: cpu_symbolic_subtree, construct_cpu_symbolic_subtree
   public :: cpu_numeric_subtree

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

   type, extends(symbolic_subtree_base) :: cpu_symbolic_subtree
      type(C_PTR) :: csubtree
   contains
      procedure :: factor
      procedure :: factor2 ! fixme remove
      final :: symbolic_final
   end type cpu_symbolic_subtree

   type, extends(numeric_subtree_base) :: cpu_numeric_subtree
      logical(C_BOOL) :: posdef
      type(C_PTR) :: csubtree
      type(contrib_type), pointer :: contrib
   contains
      procedure :: get_contrib
      procedure :: solve_fwd
      procedure :: solve_diag
      procedure :: solve_bwd
      procedure :: factor3 ! fixme remove
      final :: numeric_final
   end type cpu_numeric_subtree

contains

function construct_cpu_symbolic_subtree(nnodes, sptr, sparent, rptr, &
      rlist) result(this)
   class(cpu_symbolic_subtree), pointer :: this
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist

   allocate(this)
   this%csubtree = &
      c_create_symbolic_subtree(nnodes, sptr, sparent, rptr, rlist)
end function construct_cpu_symbolic_subtree

subroutine symbolic_final(this)
   type(cpu_symbolic_subtree) :: this

   call c_destroy_symbolic_subtree(this%csubtree)
end subroutine symbolic_final

function factor2(this, posdef, nnodes, nodes)
   type(cpu_numeric_subtree), pointer :: factor2
   class(cpu_symbolic_subtree), intent(inout) :: this
   logical(C_BOOL), intent(in) :: posdef
   integer, intent(in) :: nnodes
   type(cpu_node_data), dimension(nnodes), intent(inout) :: nodes

   type(cpu_numeric_subtree), pointer :: cpu_factor

   ! Setup output
   allocate(cpu_factor)
   factor2 => cpu_factor

   cpu_factor%posdef = posdef
   cpu_factor%csubtree = &
      c_create_cpu_subtree(posdef, this%csubtree, nnodes, nodes)
end function factor2

subroutine factor3(this, n, nnodes, nodes, aval, alloc, options, &
      stats, scaling)
   class(cpu_numeric_subtree) :: this
   integer(C_INT), intent(in) :: n
   integer(C_INT), intent(in) :: nnodes
   type(cpu_node_data), dimension(nnodes), intent(inout) :: nodes
   real(C_DOUBLE), dimension(*), intent(in) :: aval
   type(C_PTR), intent(in) :: alloc
   type(cpu_factor_options), intent(in) :: options
   type(cpu_factor_stats), intent(out) :: stats
   real(C_DOUBLE), dimension(*), target, optional, intent(in) :: scaling

   type(C_PTR) :: cscaling

   cscaling = C_NULL_PTR
   if(present(scaling)) cscaling = C_LOC(scaling)

   call c_factor_cpu(this%posdef, this%csubtree, n, nnodes, nodes, aval, &
      cscaling, alloc, options, stats)
end subroutine factor3

function factor(this, pos_def, aval, scaling)
   class(numeric_subtree_base), pointer :: factor
   class(cpu_symbolic_subtree), intent(inout) :: this
   logical, intent(in) :: pos_def
   real(wp), dimension(*), intent(in) :: aval
   real(wp), dimension(*), optional, intent(in) :: scaling

   type(cpu_numeric_subtree), pointer :: cpu_factor

   ! Setup output
   allocate(cpu_factor)
   factor => cpu_factor

   ! FIXME: remove redundancy
   !cpu_factor%csubtree = &
   !   c_create_cpu_subtree(posdef, this%csubtree, this%nnodes, this%nodes)
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
