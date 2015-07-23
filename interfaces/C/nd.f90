module spral_nd_ciface
   use, intrinsic :: iso_c_binding
   use spral_nestd, nd_options => nestd_options, &
                    nd_inform => nestd_inform, &
                    nd_order => nestd_order
   implicit none

   type, bind(C) :: spral_nd_options
      integer(C_INT) :: array_base
      integer(C_INT) :: print_level
      integer(C_INT) :: unit_diagnostics
      integer(C_INT) :: unit_error
      integer(C_INT) :: amd_call
      integer(C_INT) :: amd_switch1
      real(C_DOUBLE) :: balance
      integer(C_INT) :: coarse_partition_method
      logical(C_BOOL) :: find_supervariables
      integer(C_INT) :: max_improve_cycles
      integer(C_INT) :: matching
      real(C_DOUBLE) :: max_reduction
      real(C_DOUBLE) :: min_reduction
      integer(C_INT) :: partition_method
      integer(C_INT) :: refinement_band
      logical(C_BOOL) :: remove_dense_rows
      integer(C_INT) :: stop_coarsening1
      integer(C_INT) :: stop_coarsening2
   end type spral_nd_options

   type, bind(C) :: spral_nd_inform
      integer(C_INT) :: flag
      integer(C_INT) :: dense
      integer(C_INT) :: nsuper
      integer(C_INT) :: nzsuper
      integer(C_INT) :: stat
   end type spral_nd_inform

contains
   subroutine copy_options_in(coptions, foptions, cindexed)
      type(spral_nd_options), intent(in) :: coptions
      type(nd_options), intent(inout) :: foptions ! inherit some defaults!
      logical, intent(out) :: cindexed

      cindexed                      = (coptions%array_base.eq.0)
      foptions%print_level          = coptions%print_level
      foptions%unit_diagnostics     = coptions%unit_diagnostics
      foptions%unit_error           = coptions%unit_error
      foptions%amd_call             = coptions%amd_call
      foptions%amd_switch1          = coptions%amd_switch1
      foptions%balance              = coptions%balance
      foptions%coarse_partition_method = coptions%coarse_partition_method
      foptions%find_supervariables  = coptions%find_supervariables
      foptions%max_improve_cycles   = coptions%max_improve_cycles
      foptions%matching             = coptions%matching
      foptions%max_reduction        = coptions%max_reduction
      foptions%min_reduction        = coptions%min_reduction
      foptions%partition_method     = coptions%partition_method
      foptions%refinement_band      = coptions%refinement_band
      foptions%remove_dense_rows    = coptions%remove_dense_rows
      foptions%stop_coarsening1     = coptions%stop_coarsening1
      foptions%stop_coarsening2     = coptions%stop_coarsening2
   end subroutine copy_options_in

   subroutine copy_inform_out(finform, cinform)
      type(nd_inform), intent(in) :: finform
      type(spral_nd_inform), intent(out) :: cinform

      cinform%flag      = finform%flag
      cinform%dense     = finform%dense
      cinform%nsuper    = finform%nsuper
      cinform%nzsuper   = finform%nzsuper
      cinform%stat      = finform%stat
   end subroutine copy_inform_out

end module spral_nd_ciface

subroutine spral_nd_default_options(coptions) bind(C)
   use spral_nd_ciface
   implicit none

   type(spral_nd_options), intent(out) :: coptions

   type(nd_options) :: default_options

   coptions%array_base           = 0
   coptions%print_level          = default_options%print_level
   coptions%unit_diagnostics     = default_options%unit_diagnostics
   coptions%unit_error           = default_options%unit_error
   coptions%amd_call             = default_options%amd_call
   coptions%amd_switch1          = default_options%amd_switch1
   coptions%balance              = default_options%balance
   coptions%coarse_partition_method = default_options%coarse_partition_method
   coptions%find_supervariables  = default_options%find_supervariables
   coptions%max_improve_cycles   = default_options%max_improve_cycles
   coptions%matching             = default_options%matching
   coptions%max_reduction        = default_options%max_reduction
   coptions%min_reduction        = default_options%min_reduction
   coptions%partition_method     = default_options%partition_method
   coptions%refinement_band      = default_options%refinement_band
   coptions%remove_dense_rows    = default_options%remove_dense_rows
   coptions%stop_coarsening1     = default_options%stop_coarsening1
   coptions%stop_coarsening2     = default_options%stop_coarsening2
end subroutine spral_nd_default_options

subroutine spral_nd_order(mtx, n, cptr, crow, perm, coptions, cinform) bind(C)
   use spral_nd_ciface
   implicit none

   integer(C_INT), value :: mtx
   integer(C_INT), value :: n
   integer(C_INT), dimension(*), target, intent(in) :: cptr
   integer(C_INT), dimension(*), target, intent(in) :: crow
   integer(C_INT), dimension(*), intent(out) :: perm
   type(spral_nd_options), intent(in) :: coptions
   type(spral_nd_inform), intent(out) :: cinform

   integer(C_INT), dimension(:), pointer :: fptr
   integer(C_INT), dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer(C_INT), dimension(:), allocatable, target :: frow_alloc
   type(nd_options) :: foptions
   type(nd_inform) :: finform

   logical :: cindexed

   ! Copy options in to find out whether we're using C or Fortran indexing
   call copy_options_in(coptions, foptions, cindexed)

   ! Translate arguments
   if(cindexed) then
      allocate(fptr_alloc(n+1), frow_alloc(cptr(n)))
      fptr_alloc(1:n+1) = cptr(1:n+1) + 1
      fptr => fptr_alloc
      frow_alloc(1:fptr(n+1)-1) = crow(1:fptr(n+1)-1)
      frow => frow_alloc
   else
      fptr => cptr(1:n+1)
      frow => crow(1:fptr(n+1)-1)
   endif

   ! Call Fortran procedure
   call nd_order(mtx, n, fptr, frow, perm, foptions, finform)

   ! Copy arguments out
   if(cindexed) perm(1:n) = perm(1:n) - 1
   call copy_inform_out(finform, cinform)
end subroutine spral_nd_order
