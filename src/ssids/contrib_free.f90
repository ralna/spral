! This module contains routines for freeing contrib_type.
! As it depends on routines defined by module that use the type, it needs
! to be a seperate module.
module spral_ssids_contrib_free
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_cpu_subtree, only : cpu_free_contrib
   implicit none

contains
   subroutine contrib_free(contrib)
      type(contrib_type), intent(inout) :: contrib

      call cpu_free_contrib(contrib%posdef, contrib%owner_ptr)
   end subroutine contrib_free
end module spral_ssids_contrib_free

! The C prototype for the following routine is in contrib.h
subroutine spral_ssids_contrib_free_dbl(ccontrib) bind(C)
   use, intrinsic :: iso_c_binding
   use spral_ssids_contrib_free
   implicit none

   type(C_PTR), value :: ccontrib

   type(contrib_type), pointer :: fcontrib

   call c_f_pointer(ccontrib, fcontrib)
   call contrib_free(fcontrib)
end subroutine spral_ssids_contrib_free_dbl
