module spral_ssids_alloc
  use iso_c_binding
!$ use omp_lib
  use spral_cuda, only : cudaMalloc, cudaFree, aligned_size, c_ptr_plus
  use spral_ssids_datatypes, only : smalloc_type, stack_mem_type, long, wp
  implicit none

  private
  ! CUDA stack allocator
  public :: cuda_stack_alloc_type, & ! Data type
            custack_init,          & ! Initialize stack allocator
            custack_alloc,         & ! Allocate from top of stack
            custack_free,          & ! Free from top of stack
            custack_finalize         ! Free memory associated with stack
  ! Small memory allocator (page based, assumes only deallocated all at once)
  public :: smalloc,      & ! Allocate a small amount of memory from pool
            smfreeall,    & ! Free all data associated with smallocator
            smalloc_setup   ! Initialize smallocator

  type cuda_stack_alloc_type
     private
     type(C_PTR) :: stack = C_NULL_PTR ! GPU pointer to memory
     integer(C_SIZE_T) :: stack_sz = 0 ! Size of stack
     integer(C_SIZE_T) :: top = 0 ! Current top of stack
  end type cuda_stack_alloc_type

  ! Make smalloc generic
  interface smalloc
     module procedure real_alloc, int_alloc
  end interface smalloc

  ! Assorted fixed parameters (could be made part of control)
  integer(long), parameter, public :: PAGE_SZ = 2**20 ! 8+4MB pages (real+int)

contains

  subroutine custack_init(stack, bytes, cuda_error)
    implicit none
    type(cuda_stack_alloc_type), intent(inout) :: stack
    integer(C_SIZE_T), intent(in) :: bytes
    integer, intent(out) :: cuda_error

    cuda_error = 0 ! All is good

    ! Check stack not still in use
    if (stack%top .ne. 0) then
       ! We should never reach this point
       print *, "Attempting to resize non-empty stack!"
       stop
    end if

    ! If stack is already large enough, do nothing
    if (bytes .le. stack%stack_sz) return

    ! Free any preexisting stack
    if (C_ASSOCIATED(stack%stack)) then
       cuda_error = cudaFree(stack%stack)
       if (cuda_error .ne. 0) return
    end if

    ! Always align!
    stack%stack_sz = aligned_size(bytes)

    ! Allocate stack to new size
    cuda_error = cudaMalloc(stack%stack, stack%stack_sz)
    if (cuda_error .ne. 0) return
  end subroutine custack_init

  subroutine custack_finalize(stack, cuda_error)
    implicit none
    type(cuda_stack_alloc_type), intent(inout) :: stack
    integer, intent(out) :: cuda_error

    cuda_error = 0 ! All is good

    ! Reset information
    stack%top = 0
    stack%stack_sz = 0

    ! Don't bother trying to free memory if none allocated
    if (.not. C_ASSOCIATED(stack%stack)) return ! Not initialized

    ! Free memory
    cuda_error = cudaFree(stack%stack)
    if (cuda_error .ne. 0) return

    ! Nullify pointer
    stack%stack = C_NULL_PTR
  end subroutine custack_finalize

  type(C_PTR) function custack_alloc(stack, bytes)
    implicit none
    type(cuda_stack_alloc_type), intent(inout) :: stack
    integer(C_SIZE_T), intent(in) :: bytes

    integer(C_SIZE_T) :: bytes_aligned ! Alloc size has to round up to nearest
      ! multiple of 256 bytes

    bytes_aligned = aligned_size(bytes)
    if ((stack%top + bytes_aligned) .gt. stack%stack_sz) then
       ! Insufficient space
       custack_alloc = C_NULL_PTR
       return
    end if

    custack_alloc = c_ptr_plus(stack%stack, stack%top)
    stack%top = stack%top + bytes_aligned
  end function custack_alloc

  subroutine custack_free(stack, bytes)
    implicit none
    type(cuda_stack_alloc_type), intent(inout) :: stack
    integer(C_SIZE_T), intent(in) :: bytes

    integer(C_SIZE_T) :: bytes_aligned ! Alloc size has to round up to nearest
      ! multiple of 256 bytes

    bytes_aligned = aligned_size(bytes)
    stack%top = stack%top - bytes_aligned
  end subroutine custack_free

!*************************************************
!
! Initialize memory storage of current unallocated page with specified sizes
! (minimum of PAGE_SZ)
! MUST be called prior to smalloc() to set up alloc%lock
!
  subroutine smalloc_setup(alloc, int_hint, real_hint, st)
    implicit none
    type(smalloc_type), intent(inout) :: alloc
    integer(long), intent(in) :: int_hint
    integer(long), intent(in) :: real_hint
    integer, intent(out) :: st

!$  call omp_init_lock(alloc%lock)

    st = 0
    if (.not. allocated(alloc%imem)) then
       alloc%imem_size = max(PAGE_SZ+0_long,int_hint)
       allocate(alloc%imem(alloc%imem_size), stat=st)
    end if
    if (.not. allocated(alloc%rmem)) then
       alloc%rmem_size = max(PAGE_SZ+0_long,real_hint)
       allocate(alloc%rmem(alloc%rmem_size), stat=st)
    end if
  end subroutine smalloc_setup

!*************************************************
!
! Grab some real memory of the specified size and return a pointer to it
! Preference given to same page as last allocation.
! If this is not possible then check if we fit on any existing pages
! If not make a new page
! Note: as may be accessed by multiple threads, we have a lock
!
  subroutine real_alloc(alloc_in, ptr, len, srcptr, srchead, st)
    implicit none
    type(smalloc_type), target, intent(inout) :: alloc_in
    real(wp), dimension(:), pointer, intent(out) :: ptr
    integer(long), intent(in) :: len
    type(smalloc_type), pointer, intent(out) :: srcptr
    integer(long), intent(out) :: srchead
    integer, intent(out) :: st

    type(smalloc_type), pointer :: alloc

    st = 0
    if (len .lt. 0) return

!$  call omp_set_lock(alloc_in%lock)

    ! First try same page as last alloc
    if (associated(alloc_in%top_real)) then
       alloc => alloc_in%top_real
       if (alloc%rhead+len .le. alloc%rmem_size) then
          ! Sufficient space, allocate in this block
          ptr => alloc%rmem(alloc%rhead+1:alloc%rhead+len)
          srcptr => alloc
          srchead = alloc%rhead+1
          alloc%rhead = alloc%rhead + len
!$        call omp_unset_lock(alloc_in%lock)
          return
       end if
    end if

    ! Else check all pages, if we reach the end create a new one of sufficient
    ! size to allocate the pointer
    alloc => alloc_in
    do
       if (.not. allocated(alloc%rmem)) then
          alloc%rmem_size = max(len,PAGE_SZ)
          allocate(alloc%rmem(alloc%rmem_size), stat=st)
          if (st .ne. 0) then
!$           call omp_unset_lock(alloc_in%lock)
             return
          end if
       end if
       if (alloc%rhead+len .le. alloc%rmem_size) then
          ! Sufficient space, allocate in this block
          ptr => alloc%rmem(alloc%rhead+1:alloc%rhead+len)
          srcptr => alloc
          srchead = alloc%rhead+1
          alloc%rhead = alloc%rhead + len
          alloc_in%top_real => alloc
          exit
       end if
       ! Insufficent space, move to next block
       if (.not. associated(alloc%next_alloc)) then
          allocate(alloc%next_alloc, stat=st)
          if (st .ne. 0) then
!$           call omp_unset_lock(alloc_in%lock)
             return
          end if
       end if
       alloc => alloc%next_alloc
    end do

!$  call omp_unset_lock(alloc_in%lock)
  end subroutine real_alloc

!*************************************************
!
! Grab some integer memory of the specified size and return a pointer to it
! Preference given to same page as last allocation.
! If this is not possible then check if we fit on any existing pages
! If not make a new page
! Note: as may be accessed by multiple threads, we have a lock
!
  subroutine int_alloc(alloc_in, ptr, len, srcptr, srchead, st)
    implicit none
    type(smalloc_type), target, intent(inout) :: alloc_in
    integer, dimension(:), pointer, intent(out) :: ptr
    integer(long), intent(in) :: len
    type(smalloc_type), pointer, intent(out) :: srcptr
    integer(long), intent(out) :: srchead
    integer, intent(out) :: st

    type(smalloc_type), pointer :: alloc

    st = 0
    if (len .lt. 0) return

!$  call omp_set_lock(alloc_in%lock)

    ! First try same page as last alloc
    if (associated(alloc_in%top_int)) then
       alloc => alloc_in%top_int
       if (alloc%ihead+len .le. alloc%imem_size) then
          ! Sufficient space, allocate in this block
          ptr => alloc%imem(alloc%ihead+1:alloc%ihead+len)
          srcptr => alloc
          srchead = alloc%ihead+1
          alloc%ihead = alloc%ihead + len
!$        call omp_unset_lock(alloc_in%lock)
          return
       end if
    end if

    ! Else check all pages, if we reach the end create a new one of sufficient
    ! size to allocate the pointer
    alloc => alloc_in
    do
       if (.not. allocated(alloc%imem)) then
          alloc%imem_size = max(len,PAGE_SZ)
          allocate(alloc%imem(alloc%imem_size),stat=st)
          if (st .ne. 0) then
!$           call omp_unset_lock(alloc_in%lock)
             return
          end if
       end if
       if (alloc%ihead+len .le. size(alloc%imem)) then
          ! Sufficient space, allocate in this block
          ptr => alloc%imem(alloc%ihead+1:alloc%ihead+len)
          srcptr => alloc
          srchead = alloc%ihead+1
          alloc%ihead = alloc%ihead + len
          alloc_in%top_int => alloc
          exit
       end if
       ! Insufficent space, move to next block
       if (.not. associated(alloc%next_alloc)) then
          allocate(alloc%next_alloc, stat=st)
          if (st .ne. 0) then
!$           call omp_unset_lock(alloc_in%lock)
             return
          end if
       end if
       alloc => alloc%next_alloc
    end do

!$  call omp_unset_lock(alloc_in%lock)
  end subroutine int_alloc

!*************************************************
!
! Free all memory associated with the specified alloc linked list
! (both real and integer)
!
  subroutine smfreeall(alloc_in)
    implicit none
    type(smalloc_type), intent(inout) :: alloc_in

    type(smalloc_type), pointer :: alloc, alloc2
    integer :: st

!$  call omp_destroy_lock(alloc_in%lock)

    deallocate(alloc_in%rmem, stat=st)
    alloc_in%rhead = 0
    deallocate(alloc_in%imem, stat=st)
    alloc_in%ihead = 0
    alloc => alloc_in%next_alloc
    nullify(alloc_in%next_alloc, alloc_in%top_real, alloc_in%top_int)

    do while (associated(alloc))
       deallocate(alloc%rmem, stat=st)
       deallocate(alloc%imem, stat=st)
       alloc2 => alloc%next_alloc
       deallocate(alloc)
       alloc => alloc2
    end do
  end subroutine smfreeall

end module spral_ssids_alloc
