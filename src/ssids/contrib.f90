module spral_ssids_contrib
   use spral_ssids_datatypes, only : wp
   implicit none

   ! This class represents a contribution block being passed between two
   ! subtrees. It exists in CPU memory, but provides a cleanup routine as
   ! memory management may differ between two subtrees being passed
   type :: contrib_type
      integer :: n ! size of block
      real(wp), dimension(:), pointer :: val ! n x n lower triangular matrix
      integer, dimension(:), pointer :: rlist ! row list
   contains
      procedure :: cleanup => default_contrib_cleanup
   end type contrib_type

contains
   ! Default cleanup is just to Fortran deallocate memory
   subroutine default_contrib_cleanup(this)
      class(contrib_type), intent(inout) :: this

      deallocate(this%val); nullify(this%val)
      deallocate(this%rlist); nullify(this%rlist)
   end subroutine default_contrib_cleanup
end module spral_ssids_contrib
