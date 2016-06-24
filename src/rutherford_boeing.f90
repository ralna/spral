! COPYRIGHT (c) 2000,2010,2013 Science and Technology Facilities Council
! Authors: Jonathan Hogg and Iain Duff
!
! Based on modified versions of MC56 and HSL_MC56.
module spral_rutherford_boeing
   use spral_matrix_util, only : half_to_full
   use spral_random, only : random_state, random_real
   implicit none

   integer, parameter :: wp = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)
   real(wp), parameter :: zero = 0.0_wp

   private
   public :: rb_peek, &        ! Peeks at the header of a RB file
             rb_read           ! Reads a RB file
   public :: rb_reader_options ! Options that control what rb_read returns

   ! Possible values of control%format
   integer, parameter :: FORMAT_CSC = 1 ! Compressed Sparse Column
   integer, parameter :: FORMAT_CSR = 2 ! Compressed Sparse Row
   integer, parameter :: FORMAT_COO = 3 ! Coordinate

   ! Possible values control%lwr_upr_full
   integer, parameter :: TRI_LWR  = 1 ! Lower triangle
   integer, parameter :: TRI_UPR  = 2 ! Upper triangle
   integer, parameter :: TRI_FULL = 3 ! Both lower and upper triangles

   ! Possible values of control%values
   integer, parameter :: VALUES_FILE       = 0 ! As per file
   integer, parameter :: VALUES_PATTERN    = 1 ! Pattern only
   integer, parameter :: VALUES_SYM        = 2 ! Random values, symmetric
   integer, parameter :: VALUES_DIAG_DOM   = 3 ! Random vals, diag dominant
   integer, parameter :: VALUES_UNSYM      = 4 ! Random values, unsymmetric

   ! Possible error returns
   integer, parameter :: SUCCESS           =  0 ! No errors
   integer, parameter :: ERROR_UNIT        = -1 ! Failed to find a unit
   integer, parameter :: ERROR_BAD_FILE    = -2 ! Failed to open file
   integer, parameter :: ERROR_NOT_RB      = -3 ! Header not valid for RB
   integer, parameter :: ERROR_IO          = -4 ! Error return from io
   integer, parameter :: ERROR_TYPE        = -5 ! Tried to read bad type
   integer, parameter :: ERROR_ELT_ASM     = -6 ! Read elt as asm or v/v
   integer, parameter :: ERROR_EXTRA_SPACE = -10 ! control%extra_space<1.0
   integer, parameter :: ERROR_LWR_UPR_FULL= -11 ! control%lwr_up_full oor
   integer, parameter :: ERROR_FORMAT      = -12 ! control%format oor
   integer, parameter :: ERROR_VALUES      = -13 ! control%values oor
   integer, parameter :: ERROR_ALLOC       = -20 ! failed on allocate

   ! Possible warnings
   integer, parameter :: WARN_AUX_FILE     = 1 ! values in auxiliary file

   type rb_reader_options
      logical  :: add_diagonal = .false.        ! Add missing diagonal entries
      real     :: extra_space = 1.0             ! Array sizes are mult by this
      integer  :: format = FORMAT_CSC      ! Format to manipulate to
      integer  :: lwr_upr_full = TRI_LWR   ! Ensure entries in lwr/upr tri
      integer  :: values = VALUES_FILE     ! As per file
   end type rb_reader_options

   interface rb_peek
      module procedure rb_peek_file, rb_peek_unit
   end interface rb_peek

   interface rb_read
      module procedure rb_read_double_int32, rb_read_double_int64
   end interface rb_read
contains
   !
   ! This subroutine reads the header information for a file.
   !
   subroutine rb_peek_file(filename, info, m, n, nelt, nvar, nval, &
         type_code, title, identifier)
      character(len=*), intent(in) :: filename  ! File to peek at
      integer, intent(out) :: info              ! Return code
      integer, optional, intent(out) :: m       ! number of rows
      integer, optional, intent(out) :: n       ! number of columns
      integer, optional, intent(out) :: nelt    ! number of elements (0 if asm)
      integer, optional, intent(out) :: nvar    ! number of indices in file
      integer, optional, intent(out) :: nval    ! number of values in file
      character(len=3), optional, intent(out) :: type_code ! eg "rsa"
      character(len=72), optional, intent(out) :: title ! title field of file
      character(len=8), optional, intent(out) :: identifier ! id field of file

      integer :: iunit ! unit file is open on
      integer :: iost ! stat parameter for io calls

      info = SUCCESS

      ! Find a free unit and open file on it
      open(newunit=iunit, file=filename, status="old", action="read", &
         iostat=iost)
      if(iost.ne.0) then
         info = ERROR_BAD_FILE
         return
      endif

      ! Call unit version to do hard work, no need to rewind as we will close
      ! file immediately
      call rb_peek_unit(iunit, info, m, n, nelt, nvar, nval, type_code, &
         title, identifier, rewind=.false.)

      ! Close file
      close(iunit, iostat=iost)
      if(iost.ne.0 .and. info.eq.SUCCESS) then
         ! Note: we ignore close errors if info indicates a previous error
         info = ERROR_IO
         return
      endif
   end subroutine rb_peek_file

   subroutine rb_peek_unit(iunit, info, m, n, nelt, nvar, nval, &
         type_code, title, identifier, rewind)
      integer, intent(in) :: iunit           ! unit file is open on
      integer, intent(out) :: info           ! return code
      integer, optional, intent(out) :: m    ! number of rows
      integer, optional, intent(out) :: n    ! number of columns
      integer, optional, intent(out) :: nelt ! number of elements (0 if asm)
      integer, optional, intent(out) :: nvar ! number of indices in file
      integer, optional, intent(out) :: nval ! number of values in file
      character(len=3), optional, intent(out) :: type_code ! eg "rsa"
      character(len=72), optional, intent(out) :: title ! title field of file
      character(len=8), optional, intent(out) :: identifier ! id field of file
      logical, optional, intent(in) :: rewind ! If present and false, don't
         ! backspace unit to start

      ! "shadow" versions of file data - can't rely on arguments being present
      ! so data is read into these and copied to arguments if required
      integer :: r_m
      integer :: r_n
      integer :: r_nelt
      integer :: r_nvar
      integer :: r_nval
      character(len=3) :: r_type_code
      character(len=72) :: r_title
      character(len=8) :: r_identifier
      logical :: r_rewind

      ! Other local variables
      character(len=80) :: buffer1, buffer2 ! Buffers for reading char data
      integer :: t1, t2, t3, t4 ! Temporary variables for reading integer data
      integer :: iost ! stat parameter for io ops

      info = SUCCESS

      r_rewind = .true.
      if(present(rewind)) r_rewind = rewind

      ! Nibble top of file to find desired information, then return to original
      ! position if required
      read (iunit, '(a72,a8/a80/a80)', iostat=iost) &
         r_title, r_identifier, buffer1, buffer2
      if(iost.ne.0) then
         info = ERROR_IO
         return
      endif
      if(r_rewind) then
         backspace(iunit); backspace(iunit); backspace(iunit)
      endif

      read(buffer2, '(a3,11x,4(1x,i13))') r_type_code, t1, t2, t3, t4

      !
      ! Validate type_code code, remap data depending on value of type_code(3:3)
      !
      select case (r_type_code(1:1))
      case("r", "c", "i", "p", "q")
         ! Good, do nothing
      case default
         ! Not a matrix in RB format
         info = ERROR_NOT_RB
         return
      end select

      select case (r_type_code(2:2))
      case("s", "u", "h", "z", "r")
         ! Good, do nothing
      case default
         ! Not a matrix in RB format
         info = ERROR_NOT_RB
         return
      end select

      select case (r_type_code(3:3))
      case("a")
         ! Assembled format
         r_m = t1
         r_n = t2
         r_nvar = t3
         if(t4.ne.0) then
            ! RB format requires t4 to be an explicit zero
            info = ERROR_NOT_RB
            return
         endif
         r_nval = r_nvar ! one-to-one correspondence between integers and reals
         r_nelt = 0 ! no elemental matrices
      case("e")
         ! Elemental format
         r_m = t1
         r_n = r_m ! Elemental matrices are square
         r_nelt = t2
         r_nvar = t3
         r_nval = t4
      case default
         ! Not a valid RB letter code
         info = ERROR_NOT_RB
         return
      end select
      
      !
      ! Copy out data if requested
      !
      if(present(m)) m = r_m
      if(present(n)) n = r_n
      if(present(nelt)) nelt = r_nelt
      if(present(nvar)) nvar = r_nvar
      if(present(nval)) nval = r_nval
      if(present(type_code)) type_code = r_type_code
      if(present(title)) title = r_title
      if(present(identifier)) identifier = r_identifier
   end subroutine rb_peek_unit

   !
   ! This subroutine reads an assembled matrix
   !
   ! FIXME: Should only need two of ptr, row, col, but any change will require
   ! extensive testing.
   subroutine rb_read_double_int32(filename, m, n, ptr, row, col, val, &
         control, info, type_code, title, identifier, state)
      character(len=*), intent(in) :: filename ! File to read
      integer, intent(out) :: m
      integer, intent(out) :: n
      integer, dimension(:), allocatable, intent(out) :: ptr
      integer, dimension(:), allocatable, target, intent(out) :: row
      integer, dimension(:), allocatable, target, intent(out) :: col
      real(wp), dimension(:), allocatable, target, intent(out) :: val
      type(rb_reader_options), intent(in) :: control ! control variables
      integer, intent(out) :: info ! return code
      character(len=3), optional, intent(out) :: type_code ! file data type
      character(len=72), optional, intent(out) :: title ! file title
      character(len=8), optional, intent(out) :: identifier ! file identifier
      type(random_state), optional, intent(inout) :: state ! state to use for
         ! random number generation

      integer(long), dimension(:), allocatable :: ptr64
      integer :: st

      call rb_read_double_int64(filename, m, n, ptr64, row, col, val, &
         control, info, type_code=type_code, title=title, &
         identifier=identifier, state=state)

      ! FIXME: Add an error code if ne > maxint
      if(allocated(ptr64)) then
         allocate(ptr(n+1), stat=st)
         if(st.ne.0) then
            info = ERROR_ALLOC
            return
         endif
         ptr(1:n+1) = int(ptr64(1:n+1)) ! Forced conversion, FIXME: add guard
      endif
   end subroutine rb_read_double_int32
   subroutine rb_read_double_int64(filename, m, n, ptr, row, col, val, &
         control, info, type_code, title, identifier, state)
      character(len=*), intent(in) :: filename ! File to read
      integer, intent(out) :: m
      integer, intent(out) :: n
      integer(long), dimension(:), allocatable, intent(out) :: ptr
      integer, dimension(:), allocatable, target, intent(out) :: row
      integer, dimension(:), allocatable, target, intent(out) :: col
      real(wp), dimension(:), allocatable, target, intent(out) :: val
      type(rb_reader_options), intent(in) :: control ! control variables
      integer, intent(out) :: info ! return code
      character(len=3), optional, intent(out) :: type_code ! file data type
      character(len=72), optional, intent(out) :: title ! file title
      character(len=8), optional, intent(out) :: identifier ! file identifier
      type(random_state), optional, intent(inout) :: state ! state to use for
         ! random number generation

      ! Below variables are required for calling f77 MC56
      integer, dimension(10) :: icntl
      integer, dimension(5) :: f77info
      integer, dimension(:), allocatable :: ival

      ! Shadow variables for type_code, title and identifier (as arguments
      !  are optional we need a real copy)
      character(len=3) :: r_type_code
      character(len=72) :: r_title
      character(len=8) :: r_identifier

      ! Pointers to simplify which array we are reading in to.
      integer, pointer, dimension(:) :: rcptr => null()
      real(wp), pointer, dimension(:) :: vptr => null()

      real(wp), target :: temp(1) ! place holder array
      integer :: i, k ! loop indices
      integer(long) :: j ! loop indices
      integer :: r, c ! loop indices
      integer :: nnz ! number of non-zeroes
      integer :: nelt ! number of elements in file, should be 0
      integer :: len, len2 ! length of arrays to allocate
      integer :: iunit ! unit we open the file on
      integer :: st, iost ! error codes from allocate and file operations
      logical :: symmetric ! .true. if file claims to be (skew) symmetric or H
      logical :: skew ! .true. if file claims to be skew symmetric
      logical :: pattern ! .true. if we are only storing the pattern
      logical :: expanded ! .true. if pattern has been expanded
      type(random_state) :: state2 ! random number state used if state not present
      integer, dimension(:), allocatable :: iw34 ! work array used by mc34

      info = SUCCESS

      ! Initialize variables to avoid compiler warnings
      symmetric = .false.
      skew = .false.

      ! Validate control paramters
      if(control%extra_space < 1.0) then
         info = ERROR_EXTRA_SPACE
         return
      endif
      if(control%lwr_upr_full.lt.1 .or. control%lwr_upr_full.gt.3) then
         info = ERROR_LWR_UPR_FULL
         return
      endif
      if(control%format.lt.FORMAT_CSC .or. &
            control%format.gt.FORMAT_COO) then
         info = ERROR_FORMAT
         return
      endif
      if(control%values.eq.-1 .or. abs(control%values).gt.4) then
         info = ERROR_VALUES
         return
      endif

      ! Find a free unit and open file on it
      open(newunit=iunit, file=filename, status="old", action="read", &
         iostat=iost)
      if(iost.ne.0) then
         info = ERROR_BAD_FILE
         return
      endif

      ! Read top of file (and rewind) to determine space required
      call rb_peek_unit(iunit, info, m=m, n=n, nelt=nelt, &
         nval=nnz, type_code=r_type_code, title=title, &
         identifier=identifier)
      if(info.ne.0) goto 100

      if(nelt.ne.0) then
         ! Attempting to read element file as assembled
         info = ERROR_ELT_ASM
         goto 100
      endif
      
      !
      ! Allocate space for matrix
      !

      ! ptr
      len = n + 1
      len = max(len, int(len * control%extra_space))
      allocate(ptr(len), stat=st)
      if(st.ne.0) goto 200

      ! row and/or col
      len = nnz
      select case (r_type_code(2:2))
      case("s", "h", "z")
         symmetric = .true.
         skew = (r_type_code(2:2) .eq. "z")
         ! Do we need to allow for expansion?
         ! (a) to get both upper and lower triangles
         if(control%lwr_upr_full .eq. TRI_FULL) len = len * 2
         ! (b) to add additional diagonal entries
         if(control%add_diagonal .or. &
               control%values.eq.-VALUES_DIAG_DOM .or. &
               ( control%values.eq.VALUES_DIAG_DOM .and. &
               (r_type_code(1:1).eq."p" .or. r_type_code(1:1).eq."q")) ) then
            len = len + n
         endif
      case("u", "r")
         symmetric = .false.
         ! Unsymmetric or rectangular, no need to worry about upr/lwr, but
         ! may need to add diagonal.
         if(control%add_diagonal) len = len + n
      end select
      len2 = len
      len = max(len, int(len * control%extra_space))
      select case (control%format)
      case(FORMAT_CSC)
         allocate(row(len), stat=st)
         if(st.ne.0) goto 200
         rcptr => row
         if(symmetric .and. control%lwr_upr_full.eq.TRI_UPR) then
            ! We need to read into %col then copy into %row as we flip
            ! from lwr to upr
            allocate(col(len2), stat=st)
            rcptr => col
         endif
      case(FORMAT_CSR)
         allocate(col(len), stat=st)
         if(st.ne.0) goto 200
         rcptr => col
         if(symmetric .and. control%lwr_upr_full.eq.TRI_LWR) then
            ! We need to read into %row then copy into %col as we flip
            ! from upr to lwr
            allocate(row(len2), stat=st)
            rcptr => row
         endif
      case(FORMAT_COO)
         allocate(col(len), stat=st)
         if(st.ne.0) goto 200
         allocate(row(len), stat=st)
         if(control%lwr_upr_full.eq.TRI_UPR) then
            ! Store as CSR to avoid need to flip from lwr to upr
            rcptr => col
         else
            ! Store as CSC as this is natural format
            rcptr => row
         endif
      end select
      if(st.ne.0) goto 200

      ! Allocate val if required
      if(abs(control%values).ge.VALUES_SYM .or. &
            (control%values.eq.0 .and. r_type_code(1:1).ne."p" .and. &
            r_type_code(1:1).ne."q")) then
         ! We are actually going to store some values
         allocate(val(len), stat=st)
         if(st.ne.0) goto 200
         vptr => val
      else
         ! Use a place holder in call to mc56
         vptr => temp
      endif

      !
      ! Read matrix in its native format (real/integer)
      !

      icntl(1) = iunit
      ! Determine whether we need to read values from file or not
      icntl(2) = 0
      if(control%values.lt.0 .or. control%values.eq.VALUES_PATTERN) &
         icntl(2) = 1
      pattern = icntl(2).eq.1 ! .true. if we are only storing the pattern

      select case(r_type_code(1:1))
      case ("r") ! Real
         call read_data_real(iunit, r_title, r_identifier, &
            r_type_code, m, n, nnz, ptr, rcptr, f77info(1), &
            val=vptr)
      case ("c") ! Complex
         info = ERROR_TYPE
         goto 100
      case ("i") ! Integer
         if(icntl(2).ne.1) then
            allocate(ival(nnz), stat=st)
         else
            allocate(ival(1), stat=st)
         endif
         if(st.ne.0) goto 200
         call read_data_integer(iunit, r_title, r_identifier, &
            r_type_code, m, n, nnz, ptr, rcptr, f77info(1), val=ival)
         if(icntl(2).ne.1) val(1:nnz) = real(ival)
      case ("p", "q") ! Pattern only
         ! Note: if "q", then values are in an auxilary file we cannot read,
         ! so warn the user.
         if(r_type_code(1:1).eq."q") then
            icntl(2) = 1 ! Work around MC56 not recognising a 'q' file
            info = WARN_AUX_FILE
         endif
         call read_data_real(iunit, r_title, r_identifier, &
            r_type_code, m, n, nnz, ptr, rcptr, f77info(1))
         pattern = .true.
      end select
      if(f77info(1).ne.0) then
         ! Unexpected error from MC56 call
         info = -99
         goto 100
      endif

      !
      ! Add any missing diagonal entries
      !
      if(control%add_diagonal .or. &
            (symmetric .and. pattern .and. abs(control%values).eq.3)) then
         if(pattern) then
            call add_missing_diag(m, n, ptr, &
               rcptr)
         else
            call add_missing_diag(m, n, ptr, &
               rcptr, val=val)
         endif
      endif

      !
      ! Expand pattern if we need to generate unsymmetric values for it
      !
      if( ( (pattern .and. abs(control%values).eq.VALUES_UNSYM) ) &
            .and. symmetric .and. control%lwr_upr_full.eq.TRI_FULL) then
         allocate(iw34(n),stat=st)
         if(st.ne.0) goto 200
         call half_to_full(n, rcptr, ptr, iw34)
         expanded = .true.
      else
         expanded = .false.
      endif

      !
      ! Generate values if required
      !
      if(pattern .and. (control%values.lt.0 .or. control%values.ge.2)) then
         do c = 1, n
            k = int( ptr(c+1) - ptr(c) )
            if(present(state)) then
               do j = ptr(c), ptr(c+1)-1
                  val(j) = random_real(state, .false.)
                  r = rcptr(j)
                  if(abs(control%values).eq.3 .and. r.eq.c .and. symmetric)&
                     val(j) = max(100, 10*k)
               end do
            else
               do j = ptr(c), ptr(c+1)-1
                  val(j) = random_real(state2, .false.)
                  r = rcptr(j)
                  if(abs(control%values).eq.3 .and. r.eq.c .and. symmetric)&
                     val(j) = max(100, 10*k)
               end do
            endif
         end do
         pattern = .false.
      end if

      !
      ! Expand to full storage or flip lwr/upr as required
      !
      if(symmetric) then
         select case (control%lwr_upr_full)
         case(TRI_LWR)
            if(control%format.eq.FORMAT_CSR) then
               ! Only need to flip from upr to lwr if want to end up as CSR
               if(associated(vptr, val)) then
                  call flip_lwr_upr(n, ptr, row, &
                     col, st, val)
               else
                  call flip_lwr_upr(n, ptr, row, &
                     col, st)
               endif
               if(st.ne.0) goto 200
               deallocate(row)
            endif
         case(TRI_UPR)
            if(control%format.eq.FORMAT_CSC) then
               ! Only need to flip from upr to lwr if want to end up as CSC
               if(pattern) then
                  call flip_lwr_upr(n, ptr, col, &
                     row, st)
               else
                  call flip_lwr_upr(n, ptr, col, &
                     row, st, val)
               endif
               if(st.ne.0) goto 200
            endif
            if(skew .and. associated(vptr, val)) then
               call sym_to_skew(n, ptr, row, col, val)
            endif
         case(TRI_FULL)
            if(.not. allocated(iw34)) allocate(iw34(n),stat=st)
            if(st.ne.0) goto 200
            if(.not. expanded) then
               if(pattern) then
                  call half_to_full(n, rcptr, ptr, iw34)
               else
                  call half_to_full(n, rcptr, ptr, iw34, &
                     a=val)
               endif
               expanded = .true.
               if(skew .and. .not.pattern) then
                  ! HSL_MC34 doesn't cope with skew symmetry, need to flip
                  ! -ive all entries in upper triangle.
                  call sym_to_skew(n, ptr, row, col, val)
               endif
            endif
         end select
      endif

      !
      ! Convert to coordinate storage if desired
      ! Note: To avoid going to/from lower and upper we have read into
      ! %row or %col depending on the value of control%lwr_upr_full
      !
      if(control%format.eq.FORMAT_COO) then
         if(control%lwr_upr_full.eq.TRI_UPR) then
            ! Matrix is currently in CSR
            ! Fill in the row indices and deallocate ptr
            do i = 1, n
               row(ptr(i) : ptr(i+1)-1) = i
            end do
         else
            ! Matrix is currently in CSC
            ! Fill in the column indices and deallocate ptr
            do i = 1, n
               col(ptr(i) : ptr(i+1)-1) = i
            end do
         endif
         deallocate(ptr)
      endif

      100 continue

      if(present(type_code)) type_code = r_type_code
      if(present(title)) title = r_title
      if(present(identifier)) identifier = r_identifier

      close(iunit, iostat=iost)
      if(iost.ne.0 .and. info.eq.SUCCESS) then
         ! Note: we ignore close errors if info indicates a previous error
         info = ERROR_IO
         return
      endif

      !!!!!!!!!!!!!!!!
      return
      !!!!!!!!!!!!!!!!

      !
      ! Error handlers
      !
      200 continue 
         ! Allocation error
         info = ERROR_ALLOC
         goto 100
   end subroutine rb_read_double_int64

   !
   ! This subroutine takes a matrix in CSC full format that is symmetric
   ! and returns the skew symmetric interpreation obtained by setting all
   ! entries in the upper triangle to minus their original value
   !
   subroutine sym_to_skew(n, ptr, row, col, val)
      integer, intent(inout) :: n
      integer(long), dimension(n+1), intent(inout) :: ptr
      integer, dimension(:), allocatable, intent(inout) :: row
      integer, dimension(:), allocatable, intent(inout) :: col
      real(wp), dimension(ptr(n+1)-1), intent(inout) :: val

      integer :: i
      integer(long) :: j

      if(allocated(row)) then
         ! CSC format
         do i = 1, n
            do j = ptr(i), ptr(i+1)-1
               if(row(j).ge.i) cycle ! in lower triangle
               val(j) = -val(j)
            end do
         end do
      else
         ! CSR format
         do i = 1, n
            do j = ptr(i), ptr(i+1)-1
               if(i.ge.col(j)) cycle ! in lower triangle
               val(j) = -val(j)
            end do
         end do
      endif
   end subroutine sym_to_skew

   !
   ! This subroutine will transpose a matrix. To reduce copying we supply
   ! the destination integer matrix distinct from the source. The destination
   ! val and ptr arrays is the same as the source (if required).
   ! The matrix must be symmetric.
   !
   subroutine flip_lwr_upr(n, ptr, row, col, st, val)
      integer, intent(in) :: n ! Number of rows.columns in matrix (is symmetric)
      integer(long), dimension(n+1), intent(inout) :: ptr ! ptrs into row/col
      integer, dimension(ptr(n+1)-1), intent(in) :: row ! source index array
      integer, dimension(ptr(n+1)-1), intent(out) :: col ! destination index a.
      integer, intent(out) :: st ! stat parameter for allocates
      real(wp), dimension(ptr(n+1)-1), optional, intent(inout) :: val ! numeric
         ! values can be flipped as well, if required (indiciated by presence)

      integer(long) :: i ! loop indices
      integer :: r, c ! loop indices
      integer, dimension(:), allocatable :: wptr ! working copy of ptr
      real(wp), dimension(:), allocatable :: wval ! working copy of val

      ! Allocate memory
      allocate(wptr(n+2), stat=st)
      if(st.ne.0) return
      if(present(val)) allocate(wval(ptr(n+1)-1), stat=st)
      if(st.ne.0) return

      ! Count number of entries in row r as wptr(r+2)
      wptr(2:n+2) = 0
      do c = 1, n
         do i = ptr(c), ptr(c+1)-1
            r = row(i)
            wptr(r+2) = wptr(r+2) + 1
         end do
      end do

      ! Determine insert point for row r as wptr(r+1)
      wptr(1:2) = 1
      do r = 1, n
         wptr(r+2) = wptr(r+1) + wptr(r+2)
      end do

      ! Now loop over matrix inserting entries at correct points
      if(present(val)) then
         do c = 1, n
            do i = ptr(c), ptr(c+1)-1
               r = row(i)
               col(wptr(r+1)) = c
               wval(wptr(r+1)) = val(i)
               wptr(r+1) = wptr(r+1) + 1
            end do
         end do
      else
         do c = 1, n
            do i = ptr(c), ptr(c+1)-1
               r = row(i)
               col(wptr(r+1)) = c
               wptr(r+1) = wptr(r+1) + 1
            end do
         end do
      endif

      ! Finally copy data back to where it needs to be
      ptr(1:n+1) = wptr(1:n+1)
      if(present(val)) val(1:ptr(n+1)-1) = wval(1:ptr(n+1)-1)
   end subroutine flip_lwr_upr
   
   !
   ! Add any missing values to matrix (assumed to be in CSC)
   !
   subroutine add_missing_diag(m, n, ptr, row, val)
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer(long), dimension(n+1), intent(inout) :: ptr
      integer, dimension(:), intent(inout) :: row
      real(wp), dimension(*), optional, intent(inout) :: val

      integer :: col 
      integer(long) :: i
      integer :: ndiag
      logical :: found

      ! Count number of missing diagonal entries
      ndiag = 0
      do col = 1, min(m,n)
         do i = ptr(col), ptr(col+1)-1
            if(row(i).eq.col) ndiag = ndiag + 1
         end do
      end do

      ndiag = min(m,n) - ndiag ! Determine number missing

      ! Process matrix, adding diagonal entries as first entry in column if
      ! not otherwise present
      do col = n, 1, -1
         if(ndiag.eq.0) return
         found = .false.
         if(present(val)) then
            do i = ptr(col+1)-1, ptr(col), -1
               found = ( found .or. (row(i).eq.col) )
               row(i+ndiag) = row(i)
               val(i+ndiag) = val(i)
            end do
         else
            do i = ptr(col+1)-1, ptr(col), -1
               found = ( found .or. (row(i).eq.col) )
               row(i+ndiag) = row(i)
            end do
         endif
         ptr(col+1) = ptr(col+1) + ndiag
         if(.not.found .and. col.le.m) then
            ! Note: only adding diagonal if we're in the square submatrix!
            ndiag = ndiag - 1
            i = ptr(col) + ndiag
            row(i) = col
            if(present(val)) val(i) = zero
         endif
      end do
   end subroutine add_missing_diag


   !  ==================================================
   !  Code for reading files in Rutherford-Boeing format:
   !  (originally based on MC56)
   !  ==================================================

   subroutine read_data_real(lunit, title, key, dattyp, m, nvec, ne, ip, &
         ind, flag, val)
      integer, intent(in) :: lunit ! unit from which to read data
      character(len=72), intent(out) :: title   ! Title read from file
      character(len=8), intent(out) :: key      ! Key read from file
      character(len=3), intent(out) :: dattyp   ! Indicates type of data read:
         ! For matrix data this takes the form 'xyz', where
         ! x ... r, c, i, p, or x.
         ! y ... s, u, h, z, or r.
         ! z ... a or e.
      integer, intent(out) :: m ! Number of rows or the largest index used
         ! depending on whether matrix is assembled or unassembled.
      integer, intent(out) :: nvec ! Number of columns or the number of elements
         ! depending on whether matrix is assembled or unassembled.
      integer, intent(out) :: ne ! Number of entries in matrix
      integer(long), dimension(*), intent(out) :: ip ! Column/Element pointers
      integer, dimension(*), intent(out) :: ind ! Row/Element indices
      integer, intent(out) :: flag ! Return code
      real(wp), dimension(*), optional, intent(out) :: val ! If present,
         ! and DATTYP is not equal to ord, ipt, or icv, returns the numerical
         ! data.

      character(len=80) :: buffer1, buffer2
      integer :: i
      integer :: neltvl, np1, nreal
      character(len=16) :: ptrfmt, indfmt
      character(len=20) :: valfmt

      flag = 0

      ! Read in header block
      read (lunit,'(a72,a8/a80/a80)') title, key, buffer1, buffer2

      ! Check we have matrix data
      if (buffer2(3:3).ne.'e' .and. buffer2(3:3).ne.'a') then
         ! Not matrix data
         flag = ERROR_TYPE
         return
      endif

      ! Read matrix header information
      read(buffer2,'(a3,11x,4(1x,i13))') dattyp, m, nvec, ne, neltvl
      read(lunit,'(2a16,a20)') ptrfmt, indfmt, valfmt

      ! Read ip array
      np1 = nvec+1
      if (dattyp(3:3).eq.'e' .and. dattyp(2:2).eq.'r') np1=2*nvec+1
      read(lunit,ptrfmt) (ip(i),i=1,np1)

      ! Read ind array
      read(lunit,indfmt) (ind(i),i=1,ne)

      if (dattyp(1:1).eq.'p' .or. dattyp(1:1).eq.'x') return ! pattern only

      if(present(val)) then
         ! read values
         nreal = ne
         if (neltvl.gt.0) nreal = neltvl
         read(lunit,valfmt) (val(i),i=1,nreal)
      endif

   end subroutine read_data_real

   subroutine read_data_integer(lunit, title, key, dattyp, m, nvec, ne, ip, &
         ind, flag, val)
      integer, intent(in) :: lunit ! unit from which to read data
      character(len=72), intent(out) :: title   ! Title read from file
      character(len=8), intent(out) :: key      ! Key read from file
      character(len=3), intent(out) :: dattyp   ! Indicates type of data read:
         ! For matrix data this takes the form 'xyz', where
         ! x ... r, c, i, p, or x.
         ! y ... s, u, h, z, or r.
         ! z ... a or e.
      integer, intent(out) :: m ! Number of rows or the largest index used
         ! depending on whether matrix is assembled or unassembled.
      integer, intent(out) :: nvec ! Number of columns or the number of elements
         ! depending on whether matrix is assembled or unassembled.
      integer, intent(out) :: ne ! Number of entries in matrix
      integer(long), dimension(*), intent(out) :: ip ! Column/Element pointers
      integer, dimension(*), intent(out) :: ind ! Row/Element indices
      integer, intent(out) :: flag ! Return code
      integer, dimension(*), optional, intent(out) :: val ! If present,
         ! and DATTYP is not equal to ord, ipt, or icv, returns the numerical
         ! data.

      character(len=80) :: buffer1, buffer2
      integer :: i
      integer :: neltvl, np1, nreal
      character(len=16) :: ptrfmt, indfmt
      character(len=20) :: valfmt

      flag = 0

      ! Read in header block
      read (lunit,'(a72,a8/a80/a80)') title, key, buffer1, buffer2

      ! Check we have matrix data
      if (buffer2(3:3).ne.'e' .and. buffer2(3:3).ne.'a') then
         ! Not matrix data
         flag = ERROR_TYPE
         return
      endif

      ! Read matrix header information
      read(buffer2,'(a3,11x,4(1x,i13))') dattyp, m, nvec, ne, neltvl
      read(lunit,'(2a16,a20)') ptrfmt, indfmt, valfmt

      ! Read ip array
      np1 = nvec+1
      if (dattyp(3:3).eq.'e' .and. dattyp(2:2).eq.'r') np1=2*nvec+1
      read(lunit,ptrfmt) (ip(i),i=1,np1)

      ! Read ind array
      read(lunit,indfmt) (ind(i),i=1,ne)

      if (dattyp(1:1).eq.'p' .or. dattyp(1:1).eq.'x') return ! pattern only

      if(present(val)) then
         ! read values
         nreal = ne
         if (neltvl.gt.0) nreal = neltvl
         read(lunit,valfmt) (val(i),i=1,nreal)
      endif

   end subroutine read_data_integer

end module spral_rutherford_boeing
