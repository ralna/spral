module spral_ssids_cpu_subtree
   use, intrinsic :: iso_c_binding
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_cpu_iface ! fixme only
   use spral_ssids_datatypes
   use spral_ssids_inform, only : ssids_inform_base
   use spral_ssids_subtree, only : symbolic_subtree_base, numeric_subtree_base
   implicit none

   private
   public :: cpu_symbolic_subtree, construct_cpu_symbolic_subtree
   public :: cpu_numeric_subtree

   type, extends(symbolic_subtree_base) :: cpu_symbolic_subtree
      integer :: n
      type(C_PTR) :: csubtree
   contains
      procedure :: factor
      final :: symbolic_final
   end type cpu_symbolic_subtree

   type, extends(numeric_subtree_base) :: cpu_numeric_subtree
      logical(C_BOOL) :: posdef
      type(cpu_symbolic_subtree), pointer :: symbolic
      type(smalloc_type), pointer :: alloc=>null() ! Linked list of memory pages
         ! pointed to by nodes variable
      type(contrib_type), pointer :: contrib
      type(C_PTR) :: csubtree
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
   end type cpu_numeric_subtree

   interface
      type(C_PTR) function c_create_symbolic_subtree(n, nnodes, sptr, sparent, &
            rptr, rlist, nptr, nlist, options) &
            bind(C, name="spral_ssids_cpu_create_symbolic_subtree")
         use, intrinsic :: iso_c_binding
         import :: cpu_factor_options
         implicit none
         integer(C_INT), value :: n
         integer(C_INT), value :: nnodes
         integer(C_INT), dimension(nnodes+1), intent(in) :: sptr
         integer(C_INT), dimension(nnodes+1), intent(in) :: sparent
         integer(C_LONG), dimension(nnodes+1), intent(in) :: rptr
         integer(C_INT), dimension(*), intent(in) :: rlist
         integer(C_INT), dimension(*), intent(in) :: nptr
         integer(C_INT), dimension(*), intent(in) :: nlist
         type(cpu_factor_options), intent(in) :: options
      end function c_create_symbolic_subtree

      subroutine c_destroy_symbolic_subtree(subtree) &
            bind(C, name="spral_ssids_cpu_destroy_symbolic_subtree")
         use, intrinsic :: iso_c_binding
         implicit none
         type(C_PTR), value :: subtree
      end subroutine c_destroy_symbolic_subtree

      type(C_PTR) function c_create_numeric_subtree(posdef, symbolic_subtree, &
            aval, scaling, options, stats) &
            bind(C, name="spral_ssids_cpu_create_num_subtree_dbl")
         use, intrinsic :: iso_c_binding
         import :: cpu_factor_options, cpu_factor_stats
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: symbolic_subtree
         real(C_DOUBLE), dimension(*), intent(in) :: aval
         type(C_PTR), value :: scaling
         type(cpu_factor_options), intent(in) :: options
         type(cpu_factor_stats), intent(out) :: stats
      end function c_create_numeric_subtree

      subroutine c_destroy_numeric_subtree(posdef, subtree) &
            bind(C, name="spral_ssids_cpu_destroy_num_subtree_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
      end subroutine c_destroy_numeric_subtree

      subroutine c_subtree_solve_fwd(posdef, subtree, nrhs, x, ldx) &
            bind(C, name="spral_ssids_cpu_subtree_solve_fwd_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
         integer(C_INT), value :: nrhs
         real(C_DOUBLE), dimension(*), intent(inout) :: x
         integer(C_INT), value :: ldx
      end subroutine c_subtree_solve_fwd

      subroutine c_subtree_solve_diag(posdef, subtree, nrhs, x, ldx) &
            bind(C, name="spral_ssids_cpu_subtree_solve_diag_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
         integer(C_INT), value :: nrhs
         real(C_DOUBLE), dimension(*), intent(inout) :: x
         integer(C_INT), value :: ldx
      end subroutine c_subtree_solve_diag

      subroutine c_subtree_solve_diag_bwd(posdef, subtree, nrhs, x, ldx) &
            bind(C, name="spral_ssids_cpu_subtree_solve_diag_bwd_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
         integer(C_INT), value :: nrhs
         real(C_DOUBLE), dimension(*), intent(inout) :: x
         integer(C_INT), value :: ldx
      end subroutine c_subtree_solve_diag_bwd

      subroutine c_subtree_solve_bwd(posdef, subtree, nrhs, x, ldx) &
            bind(C, name="spral_ssids_cpu_subtree_solve_bwd_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
         integer(C_INT), value :: nrhs
         real(C_DOUBLE), dimension(*), intent(inout) :: x
         integer(C_INT), value :: ldx
      end subroutine c_subtree_solve_bwd

      subroutine c_subtree_enquire(posdef, subtree, piv_order, d) &
            bind(C, name="spral_ssids_cpu_subtree_enquire_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
         type(C_PTR), value :: piv_order
         type(C_PTR), value :: d
      end subroutine c_subtree_enquire

      subroutine c_subtree_alter(posdef, subtree, d) &
            bind(C, name="spral_ssids_cpu_subtree_alter_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
         real(C_DOUBLE), dimension(*), intent(in) :: d
      end subroutine c_subtree_alter
   end interface

contains

function construct_cpu_symbolic_subtree(n, nnodes, sptr, sparent, rptr, &
      rlist, nptr, nlist, options) result(this)
   class(cpu_symbolic_subtree), pointer :: this
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), target, intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), target, intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), target, intent(in) :: rlist
   integer, dimension(nnodes+1), target, intent(in) :: nptr
   integer, dimension(2,nptr(nnodes+1)-1), target, intent(in) :: nlist
   class(ssids_options), intent(in) :: options

   integer :: node
   integer(long) :: blkm, blkn
   integer :: st
   type(cpu_factor_options) :: coptions

   ! Allocate output
   allocate(this, stat=st)
   if(st.ne.0) then
      nullify(this)
      return
   endif

   ! Store basic details
   this%n = n

   ! Call C++ subtree analyse
   call cpu_copy_options_in(options, coptions)
   this%csubtree = &
      c_create_symbolic_subtree(n, nnodes, sptr, sparent, rptr, rlist, nptr, &
         nlist, coptions)

end function construct_cpu_symbolic_subtree

subroutine symbolic_final(this)
   type(cpu_symbolic_subtree) :: this

   call c_destroy_symbolic_subtree(this%csubtree)
end subroutine symbolic_final

function factor(this, posdef, aval, options, inform, scaling)
   class(numeric_subtree_base), pointer :: factor
   class(cpu_symbolic_subtree), target, intent(inout) :: this
   logical(C_BOOL), intent(in) :: posdef
   real(wp), dimension(*), intent(in) :: aval
   class(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform
   real(wp), dimension(*), target, optional, intent(in) :: scaling

   type(cpu_numeric_subtree), pointer :: cpu_factor
   type(cpu_factor_options) :: coptions
   type(cpu_factor_stats) :: cstats
   type(C_PTR) :: cscaling
   integer :: st

   ! Leave output as null until successful exit
   nullify(factor)

   ! Allocate cpu_factor for output
   allocate(cpu_factor, stat=st)
   if(st.ne.0) goto 10
   cpu_factor%symbolic => this

   ! Call C++ factor routine
   cpu_factor%posdef = posdef
   cscaling = C_NULL_PTR
   if(present(scaling)) cscaling = C_LOC(scaling)
   call cpu_copy_options_in(options, coptions)
   cpu_factor%csubtree = &
      c_create_numeric_subtree(posdef, this%csubtree, &
         aval, cscaling, coptions, cstats)
   if(cstats%flag.ne.0) then
      inform%flag = cstats%flag
      return
   endif

   ! Extract to Fortran data structures
   call cpu_copy_stats_out(this%n, cstats, inform)

   ! Success, set result and return
   factor => cpu_factor
   return

   ! Allocation error handler
   10 continue
   inform%flag = SSIDS_ERROR_ALLOCATION
   inform%stat = st
   deallocate(cpu_factor, stat=st)
   return
end function factor

subroutine numeric_final(this)
   type(cpu_numeric_subtree) :: this

   call c_destroy_numeric_subtree(this%posdef, this%csubtree)
end subroutine numeric_final

function get_contrib(this)
   class(contrib_type), pointer :: get_contrib
   class(cpu_numeric_subtree), intent(in) :: this

   !FIXME: make this actually be useful
   get_contrib => this%contrib
end function get_contrib

subroutine solve_fwd(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform

   call c_subtree_solve_fwd(this%posdef, this%csubtree, nrhs, x, ldx)
end subroutine solve_fwd

subroutine solve_diag(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform

   call c_subtree_solve_diag(this%posdef, this%csubtree, nrhs, x, ldx)
end subroutine solve_diag

subroutine solve_diag_bwd(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform

   call c_subtree_solve_diag_bwd(this%posdef, this%csubtree, nrhs, x, ldx)
end subroutine solve_diag_bwd

subroutine solve_bwd(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform

   call c_subtree_solve_bwd(this%posdef, this%csubtree, nrhs, x, ldx)
end subroutine solve_bwd

subroutine enquire_posdef(this, d)
   class(cpu_numeric_subtree), intent(in) :: this
   real(wp), dimension(*), target, intent(out) :: d

   call c_subtree_enquire(this%posdef, this%csubtree, C_NULL_PTR, C_LOC(d))
end subroutine enquire_posdef

subroutine enquire_indef(this, invp, piv_order, d)
   class(cpu_numeric_subtree), intent(in) :: this
   integer, dimension(*), intent(in) :: invp
   integer, dimension(*), optional, intent(out) :: piv_order
   real(wp), dimension(2,*), target, optional, intent(out) :: d

   integer :: i
   integer, dimension(:), target, allocatable :: po
   type(C_PTR) :: dptr, poptr

   ! Setup pointers
   poptr = C_NULL_PTR
   if(present(piv_order)) then
      allocate(po(this%symbolic%n))
      poptr = C_LOC(po)
   endif
   dptr = C_NULL_PTR
   if(present(d)) dptr = C_LOC(d)

   ! Call C++ routine
   call c_subtree_enquire(this%posdef, this%csubtree, poptr, dptr)

   ! Apply invp to piv_order
   if(present(piv_order)) then
      do i = 1, this%symbolic%n
         piv_order( invp(i) ) = po(i)
      end do
   endif
end subroutine enquire_indef

subroutine alter(this, d)
   class(cpu_numeric_subtree), target, intent(inout) :: this
   real(wp), dimension(2,*), intent(in) :: d

   call c_subtree_alter(this%posdef, this%csubtree, d)
end subroutine alter

end module spral_ssids_cpu_subtree
