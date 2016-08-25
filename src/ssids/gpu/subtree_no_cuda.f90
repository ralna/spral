! This is dummy file compiled when there is no CUDA support
module spral_ssids_gpu_subtree
   use, intrinsic :: iso_c_binding
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_datatypes
   use spral_ssids_inform, only : ssids_inform
   use spral_ssids_subtree, only : symbolic_subtree_base, numeric_subtree_base
   implicit none

   private
   public :: gpu_symbolic_subtree, construct_gpu_symbolic_subtree
   public :: gpu_numeric_subtree, gpu_free_contrib

   type, extends(symbolic_subtree_base) :: gpu_symbolic_subtree
   contains
      procedure :: factor
      final :: symbolic_final
   end type gpu_symbolic_subtree

   type, extends(numeric_subtree_base) :: gpu_numeric_subtree
   contains
      procedure :: get_contrib
      procedure :: solve_fwd
      procedure :: solve_diag
      procedure :: solve_diag_bwd
      procedure :: solve_bwd
      procedure :: enquire_posdef
      procedure :: enquire_indef
      procedure :: alter
      final :: numeric_final
   end type gpu_numeric_subtree

contains

function construct_gpu_symbolic_subtree(device, n, sa, en, sptr, sparent, &
      rptr, rlist, nptr, nlist, options) result(this)
   class(gpu_symbolic_subtree), pointer :: this
   integer, intent(in) :: device
   integer, intent(in) :: n
   integer, intent(in) :: sa
   integer, intent(in) :: en
   integer, dimension(*), target, intent(in) :: sptr
   integer, dimension(*), target, intent(in) :: sparent
   integer(long), dimension(*), target, intent(in) :: rptr
   integer, dimension(*), target, intent(in) :: rlist
   integer, dimension(*), target, intent(in) :: nptr
   integer, dimension(2,*), target, intent(in) :: nlist
   class(ssids_options), intent(in) :: options

   print *, "construct_gpu_symbolic_subtree() called without GPU support."
   print *, "This should never happen."
   stop -1
end function construct_gpu_symbolic_subtree

subroutine symbolic_final(this)
   type(gpu_symbolic_subtree) :: this
end subroutine symbolic_final

function factor(this, posdef, aval, child_contrib, options, inform, scaling)
   class(numeric_subtree_base), pointer :: factor
   class(gpu_symbolic_subtree), target, intent(inout) :: this
   logical, intent(in) :: posdef
   real(wp), dimension(*), target, intent(in) :: aval
   type(contrib_type), dimension(:), target, intent(inout) :: child_contrib
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(inout) :: inform
   real(wp), dimension(*), target, optional, intent(in) :: scaling
end function factor

subroutine numeric_final(this)
   type(gpu_numeric_subtree) :: this
end subroutine numeric_final

function get_contrib(this)
   type(contrib_type) :: get_contrib
   class(gpu_numeric_subtree), intent(in) :: this
end function get_contrib

subroutine gpu_free_contrib(contrib)
   type(contrib_type), intent(inout) :: contrib
end subroutine gpu_free_contrib

subroutine solve_fwd(this, nrhs, x, ldx, inform)
   class(gpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   type(ssids_inform), intent(inout) :: inform
end subroutine solve_fwd

subroutine solve_diag(this, nrhs, x, ldx, inform)
   class(gpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   type(ssids_inform), intent(inout) :: inform
end subroutine solve_diag

subroutine solve_diag_bwd(this, nrhs, x, ldx, inform)
   class(gpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   type(ssids_inform), intent(inout) :: inform
end subroutine solve_diag_bwd

subroutine solve_bwd(this, nrhs, x, ldx, inform)
   class(gpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   type(ssids_inform), intent(inout) :: inform
end subroutine solve_bwd

subroutine enquire_posdef(this, d)
   class(gpu_numeric_subtree), target, intent(in) :: this
   real(wp), dimension(*), target, intent(out) :: d
end subroutine enquire_posdef

subroutine enquire_indef(this, piv_order, d)
   class(gpu_numeric_subtree), target, intent(in) :: this
   integer, dimension(*), target, optional, intent(out) :: piv_order
   real(wp), dimension(2,*), target, optional, intent(out) :: d
end subroutine enquire_indef

subroutine alter(this, d)
   class(gpu_numeric_subtree), target, intent(inout) :: this
   real(wp), dimension(2,*), intent(in) :: d
end subroutine alter

end module spral_ssids_gpu_subtree
