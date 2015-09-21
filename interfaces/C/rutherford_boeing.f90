module spral_rutherford_boeing_ciface
   use, intrinsic :: iso_c_binding
   use spral_rutherford_boeing
   implicit none

   type handle_type
      integer(C_INT), dimension(:), allocatable :: ptr32
      integer(C_LONG), dimension(:), allocatable :: ptr64
      integer(C_INT), dimension(:), allocatable :: row
      integer(C_INT), dimension(:), allocatable :: col
      real(C_DOUBLE), dimension(:), allocatable :: val
   end type handle_type

   type, bind(C) :: spral_rboptions
      integer(C_INT) :: array_base
      logical(C_BOOL) :: add_diagonal
      real(C_FLOAT) :: extra_space
      integer(C_INT) :: format
      integer(C_INT) :: lwr_upr_full
      integer(C_INT) :: values
   end type spral_rboptions

   interface
      !integer(C_SIZE_T) pure function strlen(string) bind(C)
      !   use :: iso_c_binding
      !   character(C_CHAR), dimension(*), intent(in) :: string
      !end function strlen
      integer(C_SIZE_T) pure function strlen(string) bind(C)
         use :: iso_c_binding
         type(C_PTR), value :: string
      end function strlen
   end interface

contains
   subroutine copy_options_in(coptions, foptions, cindexed)
      type(spral_rboptions), intent(in) :: coptions
      type(rb_reader_options), intent(out) :: foptions
      logical, intent(out) :: cindexed

      cindexed                = (coptions%array_base.eq.0)
      foptions%add_diagonal   = coptions%add_diagonal
      foptions%extra_space    = coptions%extra_space
      foptions%format         = coptions%format
      foptions%lwr_upr_full   = coptions%lwr_upr_full
      foptions%values         = coptions%values
   end subroutine copy_options_in
end module spral_rutherford_boeing_ciface

subroutine spral_rb_default_options(coptions) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   type(spral_rboptions), intent(out) :: coptions

   type(rb_reader_options) :: foptions

   coptions%array_base     = 0
   coptions%add_diagonal   = foptions%add_diagonal
   coptions%extra_space    = foptions%extra_space
   coptions%format         = foptions%format
   coptions%lwr_upr_full   = foptions%lwr_upr_full
   coptions%values         = foptions%values
end subroutine spral_rb_default_options

integer(C_INT) function spral_rb_read_i32d(filename, handle, m, n, ptr, row, &
      col, val, options, type_code, title, identifier) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   character(C_CHAR), dimension(*), target, intent(in) :: filename
   type(C_PTR), intent(out) :: handle
   integer(C_INT), intent(out) :: m
   integer(C_INT), intent(out) :: n
   type(C_PTR), intent(out) :: ptr
   type(C_PTR), intent(out) :: row
   type(C_PTR), intent(out) :: col
   type(C_PTR), intent(out) :: val
   type(spral_rboptions), intent(in) :: options
   type(C_PTR), value :: type_code
   type(C_PTR), value :: title
   type(C_PTR), value :: identifier

   integer :: info
   type(handle_type), pointer :: matrix
   character(len=:), allocatable :: ffilename
   character(len=3) :: ftype_code
   character(len=72) :: ftitle
   character(len=8) :: fidentifier
   type(rb_reader_options) :: foptions
   logical :: cindexed

   integer :: i
   character(C_CHAR), dimension(:), pointer :: string

   ! Handle filename
   allocate(character(len=strlen(C_LOC(filename))) :: ffilename)
   do i = 1, int(strlen(C_LOC(filename)))
      ffilename(i:i) = filename(i)
   end do
   ! Create object to store data in
   allocate(matrix)
   handle = C_LOC(matrix)
   ! Copy options
   call copy_options_in(options, foptions, cindexed)

   ! Main function call
   ! FIXME: Add support for random state
   call rb_read(ffilename, m, n, matrix%ptr32, matrix%row, matrix%col, &
      matrix%val, foptions, info, type_code=ftype_code, title=ftitle, &
      identifier=fidentifier)

   ! Convert to C indexing (if required)
   if(cindexed .and. allocated(matrix%ptr32)) matrix%ptr32(:) = matrix%ptr32(:) - 1
   if(cindexed .and. allocated(matrix%row)) matrix%row(:) = matrix%row(:) - 1
   if(cindexed .and. allocated(matrix%col)) matrix%col(:) = matrix%col(:) - 1
   ! Determine pointers for ptr, row, col and val
   if(allocated(matrix%ptr32)) ptr = C_LOC(matrix%ptr32)
   if(allocated(matrix%row)) row = C_LOC(matrix%row)
   if(allocated(matrix%col)) col = C_LOC(matrix%col)
   if(allocated(matrix%val)) val = C_LOC(matrix%val)
   ! Handle optional strings
   if(C_ASSOCIATED(type_code)) then
      call c_f_pointer(type_code, string, shape = (/ 3 /))
      do i = 1, len(ftype_code)
         string(i) = ftype_code(i:i)
      end do
   endif
   if(C_ASSOCIATED(title)) then
      call c_f_pointer(title, string, shape = (/ 72 /))
      do i = 1, len(ftitle)
         string(i) = ftitle(i:i)
      end do
   endif
   if(C_ASSOCIATED(identifier)) then
      call c_f_pointer(identifier, string, shape = (/ 8 /))
      do i = 1, len(fidentifier)
         string(i) = fidentifier(i:i)
      end do
   endif

   ! Set return code
   spral_rb_read_i32d = info
end function spral_rb_read_i32d

integer(C_INT) function spral_rb_read_i64d(filename, handle, m, n, ptr, row, &
      col, val, options, type_code, title, identifier) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   character(C_CHAR), dimension(*), target, intent(in) :: filename
   type(C_PTR), intent(out) :: handle
   integer(C_INT), intent(out) :: m
   integer(C_INT), intent(out) :: n
   type(C_PTR), intent(out) :: ptr
   type(C_PTR), intent(out) :: row
   type(C_PTR), intent(out) :: col
   type(C_PTR), intent(out) :: val
   type(spral_rboptions), intent(in) :: options
   type(C_PTR), value :: type_code
   type(C_PTR), value :: title
   type(C_PTR), value :: identifier

   integer :: info
   type(handle_type), pointer :: matrix
   character(len=:), allocatable :: ffilename
   character(len=3) :: ftype_code
   character(len=72) :: ftitle
   character(len=8) :: fidentifier
   type(rb_reader_options) :: foptions
   logical :: cindexed

   integer :: i
   character(C_CHAR), dimension(:), pointer :: string

   ! Handle filename
   allocate(character(len=strlen(C_LOC(filename))) :: ffilename)
   do i = 1, int(strlen(C_LOC(filename)))
      ffilename(i:i) = filename(i)
   end do
   ! Create object to store data in
   allocate(matrix)
   handle = C_LOC(matrix)
   ! Copy options
   call copy_options_in(options, foptions, cindexed)

   ! Main function call
   ! FIXME: Add support for random state
   call rb_read(ffilename, m, n, matrix%ptr64, matrix%row, matrix%col, matrix%val, foptions, &
      info, type_code=ftype_code, title=ftitle, identifier=fidentifier)

   ! Convert to C indexing (if required)
   if(cindexed .and. allocated(matrix%ptr64)) matrix%ptr64(:) = matrix%ptr64(:) - 1
   if(cindexed .and. allocated(matrix%row)) matrix%row(:) = matrix%row(:) - 1
   if(cindexed .and. allocated(matrix%col)) matrix%col(:) = matrix%col(:) - 1
   ! Determine pointers for ptr, row, col and val
   if(allocated(matrix%ptr64)) ptr = C_LOC(matrix%ptr64)
   if(allocated(matrix%row)) row = C_LOC(matrix%row)
   if(allocated(matrix%col)) col = C_LOC(matrix%col)
   if(allocated(matrix%val)) val = C_LOC(matrix%val)
   ! Handle optional strings
   if(C_ASSOCIATED(type_code)) then
      call c_f_pointer(type_code, string, shape = (/ 3 /))
      do i = 1, len(ftype_code)
         string(i) = ftype_code(i:i)
      end do
   endif
   if(C_ASSOCIATED(title)) then
      call c_f_pointer(title, string, shape = (/ 72 /))
      do i = 1, len(ftitle)
         string(i) = ftitle(i:i)
      end do
   endif
   if(C_ASSOCIATED(identifier)) then
      call c_f_pointer(identifier, string, shape = (/ 8 /))
      do i = 1, len(fidentifier)
         string(i) = fidentifier(i:i)
      end do
   endif

   ! Set return code
   spral_rb_read_i64d = info
end function spral_rb_read_i64d

subroutine spral_rb_free_handle(handle) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   type(C_PTR), intent(inout) :: handle

   type(handle_type), pointer :: fhandle

   if(.not.C_ASSOCIATED(handle)) return
   call c_f_pointer(handle, fhandle)
   deallocate(fhandle)
   handle = C_NULL_PTR
end subroutine spral_rb_free_handle
