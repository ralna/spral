!
! This module defines ssids_fkeep type and associated procedures (CPU version)
!
module spral_ssids_fkeep_gpu
   use spral_cuda, only : cudaMemcpy_d2h, cudaMemcpy_h2d, cudaMalloc, cudaFree,&
                          c_ptr_plus, cudaStreamCreate, cudaStreamDestroy, &
                          cudaMemcpy2d, cudaMemcpyHostToDevice, &
                          cudaMemcpyDeviceToHost
   use spral_ssids_cuda_datatypes, only : gpu_type, free_gpu_type
   use spral_ssids_datatypes, only : long, ssids_akeep, ssids_options, &
                                     ssids_inform, thread_stats, wp, &
                                     SSIDS_ERROR_ALLOCATION, &
                                     SSIDS_ERROR_CUDA_UNKNOWN
   use spral_ssids_factor_gpu, only : parfactor
   use spral_ssids_fkeep, only : ssids_fkeep
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: ssids_fkeep_gpu

   !
   ! Data type for data generated in factorise phase (gpu version)
   !
   type, extends(ssids_fkeep) :: ssids_fkeep_gpu
      type(C_PTR), dimension(:), allocatable :: stream_handle
      type(gpu_type), dimension(:), allocatable :: stream_data
      type(gpu_type) :: top_data
      type(C_PTR) :: gpu_rlist_with_delays = C_NULL_PTR
      type(C_PTR) :: gpu_rlist_direct_with_delays = C_NULL_PTR
      type(C_PTR) :: gpu_clists = C_NULL_PTR
      type(C_PTR) :: gpu_clists_direct = C_NULL_PTR
      type(C_PTR) :: gpu_clen = C_NULL_PTR
      logical :: host_factors = .false.

   contains
      procedure, pass(fkeep) :: inner_factor => inner_factor_gpu ! Do actual factorization
      procedure, pass(fkeep) :: free => free_fkeep_gpu ! Frees memory
   end type ssids_fkeep_gpu

contains

subroutine inner_factor_gpu(fkeep, akeep, val, options, inform)
   class(ssids_akeep), intent(in) :: akeep
   class(ssids_fkeep_gpu), intent(inout) :: fkeep
   real(wp), dimension(*), target, intent(in) :: val
   class(ssids_options), intent(in) :: options
   class(ssids_inform), intent(inout) :: inform

   integer :: i
   type(C_PTR) :: gpu_val, gpu_scaling
   type(C_PTR), dimension(:), allocatable :: gpu_contribs
   integer(long) :: sz
   integer, dimension(:,:), allocatable :: map ! work array, one copy per
      ! thread. Size (0:n, num_threads), with 0 index used to track which
      ! node current map refers to.
   type(thread_stats), dimension(:), allocatable :: stats ! one copy
      ! per thread, accumulates per thread statistics that are then summed to
      ! obtain global stats in inform.

   integer :: num_threads
   integer :: cuda_error, st

   num_threads = 1
!$ num_threads = omp_get_max_threads()

   fkeep%host_factors = .false.

   allocate(stats(num_threads), map(0:akeep%n, num_threads), stat=st)
   if (st .ne. 0) go to 10
   map(0, :) = -1 ! initally map unassociated with any node

   ! Setup child contribution array
   ! Note: only non-NULL where we're passing contributions between subtrees
   allocate(gpu_contribs(akeep%nnodes), stat=st)
   if(st.ne.0) goto 10
   gpu_contribs(:) = C_NULL_PTR

   ! Copy A values to GPU
   sz = akeep%nptr(akeep%nnodes+1) - 1
   cuda_error = cudaMalloc(gpu_val, C_SIZEOF(val(1:sz)))
   if(cuda_error.ne.0) goto 200
   cuda_error = cudaMemcpy_h2d(gpu_val, C_LOC(val), C_SIZEOF(val(1:sz)))
   if(cuda_error.ne.0) goto 200
   
   ! Allocate and initialize streams
   if(allocated(fkeep%stream_handle)) then
      if(size(fkeep%stream_handle).lt.options%nstream) then
         do i = 1, size(fkeep%stream_handle)
            if(C_ASSOCIATED(fkeep%stream_handle(i))) then
               cuda_error = cudaStreamDestroy(fkeep%stream_handle(i))
               if(cuda_error.ne.0) goto 200
            endif
         end do
         deallocate(fkeep%stream_handle, stat=st)
         if(st.ne.0) goto 10
      endif
   endif
   if(.not.allocated(fkeep%stream_handle)) then
      allocate(fkeep%stream_handle(options%nstream), stat=st)
      if(st.ne.0) goto 10
      do i = 1, options%nstream
         cuda_error = cudaStreamCreate(fkeep%stream_handle(i))
         if(cuda_error.ne.0) goto 200
      end do
   endif

   ! Cleanup/allocate factor datastructures
   ! FIXME: We should move node<->level assignment to analyze then we can
   ! more easily reuse stream_data
   call free_gpu_type(fkeep%top_data, cuda_error)
   if(allocated(fkeep%stream_data)) then
      do i = 1, size(fkeep%stream_data)
         call free_gpu_type(fkeep%stream_data(i), cuda_error)
         if(cuda_error.ne.0) goto 200
      end do
      deallocate(fkeep%stream_data, stat=st)
      if(st.ne.0) goto 10
   endif
   allocate(fkeep%stream_data(options%nstream), stat=st)
   if (st.ne.0) goto 10

   ! Call main factorization routine
   if (allocated(fkeep%scaling)) then
      ! Copy scaling vector to GPU
      cuda_error = cudaMalloc(gpu_scaling, C_SIZEOF(fkeep%scaling(1:akeep%n)))
      if(cuda_error.ne.0) goto 200
      cuda_error = copy_to_gpu_non_target(gpu_scaling, fkeep%scaling, &
         C_SIZEOF(fkeep%scaling(1:akeep%n)))
      if(cuda_error.ne.0) goto 200

      ! Perform factorization
      call parfactor(fkeep%pos_def, akeep%child_ptr, akeep%child_list, akeep%n,&
         akeep%nptr, akeep%gpu_nlist, gpu_val, akeep%nnodes, fkeep%nodes,      &
         akeep%sptr, akeep%sparent, akeep%rptr, akeep%rlist, akeep%invp,       &
         akeep%rlist_direct, akeep%gpu_rlist, akeep%gpu_rlist_direct,          &
         gpu_contribs, fkeep%stream_handle, fkeep%stream_data,                 &
         fkeep%top_data, fkeep%gpu_rlist_with_delays,                          &
         fkeep%gpu_rlist_direct_with_delays, fkeep%gpu_clists,                 &
         fkeep%gpu_clists_direct, fkeep%gpu_clen, fkeep%alloc, options, stats, &
         ptr_scale=gpu_scaling)
      cuda_error = cudaFree(gpu_scaling)
      if(cuda_error.ne.0) goto 200
   else
      call parfactor(fkeep%pos_def, akeep%child_ptr, akeep%child_list, akeep%n,&
         akeep%nptr, akeep%gpu_nlist, gpu_val, akeep%nnodes, fkeep%nodes,      &
         akeep%sptr, akeep%sparent, akeep%rptr, akeep%rlist, akeep%invp,       &
         akeep%rlist_direct, akeep%gpu_rlist, akeep%gpu_rlist_direct,          &
         gpu_contribs, fkeep%stream_handle, fkeep%stream_data,                 &
         fkeep%top_data, fkeep%gpu_rlist_with_delays,                          &
         fkeep%gpu_rlist_direct_with_delays, fkeep%gpu_clists,                 &
         fkeep%gpu_clists_direct, fkeep%gpu_clen, fkeep%alloc, options, stats)
   end if

   cuda_error = cudaFree(gpu_val)
   if(cuda_error.ne.0) goto 200
   
   ! Do reductions
   i = minval(stats(:)%flag)
   if(i.lt.0) then
      inform%flag = i
      inform%stat = maxval(stats(:)%st)
      if(inform%stat.eq.0) inform%stat = minval(stats(:)%st)
      ! Note: cuda_error and cublas_error are actually C enums, so are +ive
      if(inform%cuda_error.eq.0) inform%cuda_error = maxval(stats(:)%cuda_error)
      if(inform%cublas_error.eq.0) &
         inform%cublas_error = maxval(stats(:)%cublas_error)
      st = inform%stat
   end if
   i = max(inform%flag, maxval(stats(:)%flag))
   inform%maxfront = maxval(stats(:)%maxfront)
   inform%num_factor = sum(stats(:)%num_factor)
   inform%num_flops = sum(stats(:)%num_flops)
   inform%num_delay = sum(stats(:)%num_delay)
   inform%num_neg = sum(stats(:)%num_neg)
   inform%num_two = sum(stats(:)%num_two)
   inform%matrix_rank = akeep%sptr(akeep%nnodes+1)-1 - sum(stats(:)%num_zero)
   fkeep%flag = inform%flag
   return
   !!!!!!!!!!!!!!!!!!!!

   !
   ! Error handling
   !
   10 continue
   inform%flag = SSIDS_ERROR_ALLOCATION
   inform%stat = st
   fkeep%flag = inform%flag
   return

   200 continue
   inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
   inform%cuda_error = cuda_error
   fkeep%flag = inform%flag
   return
end subroutine inner_factor_gpu

!****************************************************************************

! Following function used to get around target requirement of C_LOC()
integer(C_INT) function copy_to_gpu_non_target(gpu_ptr, src, sz)
   type(C_PTR) :: gpu_ptr
   real(wp), dimension(*), target, intent(in) :: src
   integer(C_SIZE_T), intent(in) :: sz

   copy_to_gpu_non_target = cudaMemcpy_h2d(gpu_ptr, C_LOC(src), sz)
end function copy_to_gpu_non_target

!****************************************************************************

subroutine free_fkeep_gpu(fkeep, flag)
   class(ssids_fkeep_gpu), intent(inout) :: fkeep
   integer, intent(out) :: flag  ! Returns any cuda error value

   integer :: st, i

   ! Skip if nothing intialized
   if (.not.allocated(fkeep%nodes)) return

   ! Call superclass free first (sets flag to 0)
   call fkeep%ssids_fkeep%free(flag)

   !
   ! Now cleanup GPU-specific stuff
   !
   call free_gpu_type(fkeep%top_data, flag)
   if(allocated(fkeep%stream_data)) then
      do i = 1, size(fkeep%stream_data)
         call free_gpu_type(fkeep%stream_data(i), flag)
      end do
      deallocate(fkeep%stream_data, stat=st)
   endif

   ! Cleanup top-level presolve info
   if(C_ASSOCIATED(fkeep%gpu_rlist_with_delays)) then
      flag = cudaFree(fkeep%gpu_rlist_with_delays)
      fkeep%gpu_rlist_with_delays = C_NULL_PTR
      if(flag.ne.0) return
   endif
   if(C_ASSOCIATED(fkeep%gpu_clists)) then
      flag = cudaFree(fkeep%gpu_clists)
      fkeep%gpu_clists = C_NULL_PTR
      if(flag.ne.0) return
   endif
   if(C_ASSOCIATED(fkeep%gpu_clists_direct)) then
      flag = cudaFree(fkeep%gpu_clists)
      fkeep%gpu_clists = C_NULL_PTR
      if(flag.ne.0) return
   endif
   if(C_ASSOCIATED(fkeep%gpu_clen)) then
      flag = cudaFree(fkeep%gpu_clen)
      fkeep%gpu_clen = C_NULL_PTR
      if(flag.ne.0) return
   endif

   ! Release streams
   if(allocated(fkeep%stream_handle)) then
      do i = 1, size(fkeep%stream_handle)
         flag = cudaStreamDestroy(fkeep%stream_handle(i))
         if(flag.ne.0) return
      end do
      deallocate(fkeep%stream_handle, stat=st)
   endif
end subroutine free_fkeep_gpu

end module spral_ssids_fkeep_gpu
