module spral_rutherford_boeing_ciface
   use, intrinsic :: iso_c_binding
   use spral_rutherford_boeing
   implicit none

   integer, parameter :: long = selected_int_kind(18)

   type handle_type
      integer(C_INT), dimension(:), allocatable :: ptr32
      integer(C_LONG), dimension(:), allocatable :: ptr64
      integer(C_INT), dimension(:), allocatable :: row
      real(C_DOUBLE), dimension(:), allocatable :: val
   end type handle_type

   type, bind(C) :: spral_rb_read_options
      integer(C_INT) :: array_base
      logical(C_BOOL) :: add_diagonal
      real(C_FLOAT) :: extra_space
      integer(C_INT) :: lwr_upr_full
      integer(C_INT) :: values
   end type spral_rb_read_options

   interface
      !integer(C_SIZE_T) pure function strlen(string) bind(C)
      !   use :: iso_c_binding
      !   character(C_CHAR), dimension(*), intent(in) :: string
      !end function strlen
      integer(C_SIZE_T) pure function strlen(string) bind(C)
         use :: iso_c_binding
         type(C_PTR), value, intent(in) :: string
      end function strlen
   end interface

contains
   subroutine copy_options_in(coptions, foptions, cindexed)
      type(spral_rb_read_options), intent(in) :: coptions
      type(rb_read_options), intent(out) :: foptions
      logical, intent(out) :: cindexed

      cindexed                = (coptions%array_base.eq.0)
      foptions%add_diagonal   = coptions%add_diagonal
      foptions%extra_space    = coptions%extra_space
      foptions%lwr_upr_full   = coptions%lwr_upr_full
      foptions%values         = coptions%values
   end subroutine copy_options_in

   subroutine convert_string_c2f(cstr, fstr)
      type(C_PTR), intent(in) :: cstr
      character(len=:), allocatable, intent(out) :: fstr

      integer :: i
      character(C_CHAR), dimension(:), pointer :: cstrptr

      allocate(character(len=strlen(cstr)) :: fstr)
      call c_f_pointer(cstr, cstrptr, shape = (/ strlen(cstr)+1 /))
      do i = 1, size(cstrptr)
         fstr(i:i) = cstrptr(i)
      end do
   end subroutine convert_string_c2f

   ! WARNING: Assumes cstr points to a sufficiently large buffer.
   subroutine convert_string_f2c(fstr, cstr)
      character(len=*), intent(in) :: fstr
      type(C_PTR), intent(inout) :: cstr

      integer :: i
      character(C_CHAR), dimension(:), pointer :: cstrptr

      call c_f_pointer(cstr, cstrptr, shape = (/ len(fstr)+1 /))
      do i = 1, len(fstr)
         cstrptr(i) = fstr(i:i)
      end do
      cstrptr(len(fstr)+1) = C_NULL_CHAR
   end subroutine convert_string_f2c
end module spral_rutherford_boeing_ciface

integer(C_INT) function spral_rb_peek(filename, m, n, nelt, nvar, nval, &
      matrix_type, type_code, title, identifier) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   type(C_PTR), value :: filename
   type(C_PTR), value :: m
   type(C_PTR), value :: n
   type(C_PTR), value :: nelt
   type(C_PTR), value :: nvar
   type(C_PTR), value :: nval
   type(C_PTR), value :: matrix_type
   type(C_PTR), value :: type_code
   type(C_PTR), value :: title
   type(C_PTR), value :: identifier

   integer :: fm, fn, fmatrix_type
   integer(long) :: fnelt, fnvar, fnval
   character(len=:), allocatable :: ffilename
   character(len=3) :: ftype_code
   character(len=72) :: ftitle
   character(len=8) :: fidentifier

   integer(C_INT), pointer :: temp_int
   integer(C_LONG), pointer :: temp_long

   ! Convert filename to Fortran string
   call convert_string_c2f(filename, ffilename)

   ! Call Fortran routine
   call rb_peek(ffilename, spral_rb_peek, m=fm, n=fn, nelt=fnelt, nvar=fnvar, &
      nval=fnval, matrix_type=fmatrix_type, type_code=ftype_code, &
      title=ftitle, identifier=fidentifier)

   ! Copy results out as approriate
   if(c_associated(m)) then
      call c_f_pointer(m, temp_int)
      temp_int = fm
   endif
   if(c_associated(n)) then
      call c_f_pointer(n, temp_int)
      temp_int = fn
   endif
   if(c_associated(nelt)) then
      call c_f_pointer(nelt, temp_long)
      temp_long = fnelt
   endif
   if(c_associated(nvar)) then
      call c_f_pointer(nvar, temp_long)
      temp_long = fnvar
   endif
   if(c_associated(nval)) then
      call c_f_pointer(nval, temp_long)
      temp_long = fnval
   endif
   if(c_associated(matrix_type)) then
      call c_f_pointer(matrix_type, temp_int)
      temp_int = fmatrix_type
   endif
   if(c_associated(type_code)) &
      call convert_string_f2c(ftype_code, type_code)
   if(c_associated(title)) &
      call convert_string_f2c(ftitle, title)
   if(c_associated(identifier)) &
      call convert_string_f2c(fidentifier, identifier)
end function spral_rb_peek

subroutine spral_rb_default_read_options(coptions) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   type(spral_rb_read_options), intent(out) :: coptions

   type(rb_read_options) :: foptions

   coptions%array_base     = 0
   coptions%add_diagonal   = foptions%add_diagonal
   coptions%extra_space    = foptions%extra_space
   coptions%lwr_upr_full   = foptions%lwr_upr_full
   coptions%values         = foptions%values
end subroutine spral_rb_default_read_options

integer(C_INT) function spral_rb_read_i32d(filename, handle, m, n, ptr, row, &
      val, options, type_code, title, identifier) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   type(C_PTR), value :: filename
   type(C_PTR), intent(out) :: handle
   integer(C_INT), intent(out) :: m
   integer(C_INT), intent(out) :: n
   type(C_PTR), intent(out) :: ptr
   type(C_PTR), intent(out) :: row
   type(C_PTR), intent(out) :: val
   type(spral_rb_read_options), intent(in) :: options
   type(C_PTR), value :: type_code
   type(C_PTR), value :: title
   type(C_PTR), value :: identifier

   integer :: info
   type(handle_type), pointer :: matrix
   character(C_CHAR), dimension(:), pointer :: cfilename
   character(len=:), allocatable :: ffilename
   character(len=3) :: ftype_code
   character(len=72) :: ftitle
   character(len=8) :: fidentifier
   type(rb_read_options) :: foptions
   logical :: cindexed

   integer :: i
   character(C_CHAR), dimension(:), pointer :: string

   ! Handle filename
   allocate(character(len=strlen(filename)) :: ffilename)
   call c_f_pointer(filename, cfilename, shape = (/ strlen(filename)+1 /))
   do i = 1, int(strlen(filename))
      ffilename(i:i) = cfilename(i)
   end do
   ! Create object to store data in
   allocate(matrix)
   handle = C_LOC(matrix)
   ! Copy options
   call copy_options_in(options, foptions, cindexed)

   ! Main function call
   ! FIXME: Add support for random state
   call rb_read(ffilename, m, n, matrix%ptr32, matrix%row, matrix%val, &
      foptions, info, type_code=ftype_code, title=ftitle, &
      identifier=fidentifier)

   ! Convert to C indexing (if required)
   if(cindexed .and. allocated(matrix%ptr32)) matrix%ptr32(:) = matrix%ptr32(:) - 1
   if(cindexed .and. allocated(matrix%row)) matrix%row(:) = matrix%row(:) - 1
   ! Determine pointers for ptr, row, and val
   if(allocated(matrix%ptr32)) ptr = C_LOC(matrix%ptr32)
   if(allocated(matrix%row)) row = C_LOC(matrix%row)
   if(allocated(matrix%val)) val = C_LOC(matrix%val)
   ! Handle optional strings
   if(C_ASSOCIATED(type_code)) then
      call c_f_pointer(type_code, string, shape = (/ 3+1 /))
      do i = 1, len(ftype_code)
         string(i) = ftype_code(i:i)
      end do
      string(len(ftype_code)+1) = C_NULL_CHAR
   endif
   if(C_ASSOCIATED(title)) then
      call c_f_pointer(title, string, shape = (/ 72+1 /))
      do i = 1, len(ftitle)
         string(i) = ftitle(i:i)
      end do
      string(len(ftitle)+1) = C_NULL_CHAR
   endif
   if(C_ASSOCIATED(identifier)) then
      call c_f_pointer(identifier, string, shape = (/ 8+1 /))
      do i = 1, len(fidentifier)
         string(i) = fidentifier(i:i)
      end do
      string(len(fidentifier)+1) = C_NULL_CHAR
   endif

   ! Set return code
   spral_rb_read_i32d = info
end function spral_rb_read_i32d

integer(C_INT) function spral_rb_read_i64d(filename, handle, m, n, ptr, row, &
      val, options, type_code, title, identifier) bind(C)
   use spral_rutherford_boeing_ciface
   implicit none

   type(C_PTR), value :: filename
   type(C_PTR), intent(out) :: handle
   integer(C_INT), intent(out) :: m
   integer(C_INT), intent(out) :: n
   type(C_PTR), intent(out) :: ptr
   type(C_PTR), intent(out) :: row
   type(C_PTR), intent(out) :: val
   type(spral_rb_read_options), intent(in) :: options
   type(C_PTR), value :: type_code
   type(C_PTR), value :: title
   type(C_PTR), value :: identifier

   integer :: info
   type(handle_type), pointer :: matrix
   character(C_CHAR), dimension(:), pointer :: cfilename
   character(len=:), allocatable :: ffilename
   character(len=3) :: ftype_code
   character(len=72) :: ftitle
   character(len=8) :: fidentifier
   type(rb_read_options) :: foptions
   logical :: cindexed

   integer :: i
   character(C_CHAR), dimension(:), pointer :: string

   ! Handle filename
   allocate(character(len=strlen(filename)) :: ffilename)
   call c_f_pointer(filename, cfilename, shape = (/ strlen(filename)+1 /))
   do i = 1, int(strlen(filename))
      ffilename(i:i) = cfilename(i)
   end do
   ! Create object to store data in
   allocate(matrix)
   handle = C_LOC(matrix)
   ! Copy options
   call copy_options_in(options, foptions, cindexed)

   ! Main function call
   ! FIXME: Add support for random state
   call rb_read(ffilename, m, n, matrix%ptr64, matrix%row, matrix%val, &
      foptions, info, type_code=ftype_code, title=ftitle, &
      identifier=fidentifier)

   ! Convert to C indexing (if required)
   if(cindexed .and. allocated(matrix%ptr64)) matrix%ptr64(:) = matrix%ptr64(:) - 1
   if(cindexed .and. allocated(matrix%row)) matrix%row(:) = matrix%row(:) - 1
   ! Determine pointers for ptr, row, and val
   if(allocated(matrix%ptr64)) ptr = C_LOC(matrix%ptr64)
   if(allocated(matrix%row)) row = C_LOC(matrix%row)
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
