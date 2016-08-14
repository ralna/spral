module spral_ssids_subtree
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_datatypes, only : long, wp, ssids_options
   use spral_ssids_inform, only : ssids_inform_base
   implicit none

   private
   public :: symbolic_subtree_base, numeric_subtree_base

   type, abstract :: symbolic_subtree_base
   contains
      procedure(factor_iface), deferred :: factor
   end type symbolic_subtree_base

   ! Type that represent a subtree solver
   type, abstract :: numeric_subtree_base
   contains
      procedure(get_contrib_iface), deferred :: get_contrib
      procedure(solve_proc_iface), deferred :: solve_fwd
      procedure(solve_proc_iface), deferred :: solve_diag
      procedure(solve_proc_iface), deferred :: solve_diag_bwd
      procedure(solve_proc_iface), deferred :: solve_bwd
   end type numeric_subtree_base

   abstract interface
      function factor_iface(this, posdef, aval, child_contrib, options, &
            inform, scaling)
         import symbolic_subtree_base, numeric_subtree_base, wp
         import ssids_inform_base, ssids_options
         import contrib_type
         implicit none
         class(numeric_subtree_base), pointer :: factor_iface
         class(symbolic_subtree_base), target, intent(inout) :: this
         logical, intent(in) :: posdef
         real(wp), dimension(*), target, intent(in) :: aval
         type(contrib_type), dimension(:), target, intent(inout) :: child_contrib
         class(ssids_options), intent(in) :: options
         class(ssids_inform_base), intent(inout) :: inform
         real(wp), dimension(*), target, optional, intent(in) :: scaling
      end function factor_iface
      function get_contrib_iface(this)
         import contrib_type, numeric_subtree_base
         implicit none
         type(contrib_type) :: get_contrib_iface
         class(numeric_subtree_base), intent(in) :: this
      end function get_contrib_iface
      subroutine solve_proc_iface(this, nrhs, x, ldx, inform)
         import numeric_subtree_base, ssids_inform_base, wp
         implicit none
         class(numeric_subtree_base), intent(inout) :: this
         integer, intent(in) :: nrhs
         real(wp), dimension(*), intent(inout) :: x
         integer, intent(in) :: ldx
         class(ssids_inform_base), intent(inout) :: inform
      end subroutine solve_proc_iface
   end interface
end module spral_ssids_subtree
