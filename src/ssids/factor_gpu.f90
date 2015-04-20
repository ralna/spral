! Copyright (c) 2013 Science and Technology Facilities Council (STFC)
! Authors: Evgueni Ovtchinnikov and Jonathan Hogg
!
! Factorize phase to run on GPU
module spral_ssids_factor_gpu
!$ use omp_lib
   use, intrinsic :: iso_c_binding
   use spral_cuda
   use spral_ssids_alloc
   use spral_ssids_cuda_datatypes
   use spral_ssids_cuda_interfaces
   use spral_ssids_datatypes
   use spral_ssids_dense_factor_gpu, only : &
      node_ldlt, node_llt, multinode_llt, multinode_ldlt
   use spral_ssids_solve_gpu, only : setup_gpu_solve
   implicit none

   private
   public :: parfactor ! Performs factorization phase using multiple streams

   type :: ntype
      integer :: level = 1
      integer :: subtree = 0
      integer(long) :: work_here = 0
      integer(long) :: work_below = 0 ! includes work_here
   end type ntype

   type :: asmtype
      integer :: npassed      ! #rows passed to parent total
      integer :: npassl       ! #cols passed to parent's L part
      integer(long) :: offset ! start of rows to pass up in rlist(:)
         ! i.e. rptr(child)+blkn-1
   end type asmtype

contains

subroutine parfactor(pos_def, child_ptr, child_list, n, nptr, gpu_nlist,      &
      ptr_val, nnodes, nodes, sptr, sparent, rptr, rlist, invp, rlist_direct, &
      gpu_rlist, gpu_rlist_direct, gpu_contribs, stream_handle, stream_data,  &
      top_data, gpu_rlist_with_delays, gpu_rlist_direct_with_delays,          &
      gpu_clists, gpu_clists_direct, gpu_clen, alloc, options, stats, ptr_scale)
   logical, intent(in) :: pos_def ! True if problem is supposedly pos-definite
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: nptr
   type(C_PTR), intent(in) :: gpu_nlist
   type(C_PTR), intent(in) :: ptr_val
   ! Note: gfortran-4.3 bug requires explicit size of nodes array
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer, dimension(*), intent(in) :: sparent
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, dimension(*), intent(in), target :: rlist_direct
   type(C_PTR), intent(in) :: gpu_rlist
   type(C_PTR), intent(in) :: gpu_rlist_direct
   type(C_PTR), dimension(*), intent(inout) :: gpu_contribs
   type(C_PTR), dimension(*), intent(in) :: stream_handle
   type(gpu_type), dimension(:), intent(out) :: stream_data
   type(gpu_type), intent(out) :: top_data
   type(C_PTR), intent(out) :: gpu_rlist_with_delays
   type(C_PTR), intent(out) :: gpu_rlist_direct_with_delays
   type(C_PTR), intent(out) :: gpu_clists
   type(C_PTR), intent(out) :: gpu_clists_direct
   type(C_PTR), intent(out) :: gpu_clen
   type(smalloc_type), target, intent(inout) :: alloc ! Contains actual memory
      ! allocations for L. Everything else (within the subtree) is just a
      ! pointer to this.
   ! explicit size required on buf to avoid gfortran-4.3 bug
   type(ssids_options), intent(in) :: options
   type(thread_stats), dimension(*), intent(inout) :: stats
   type(C_PTR), optional, intent(in) :: ptr_scale

   integer :: stream, root
   integer, dimension(:), allocatable :: stptr, stlist
   type(C_PTR), dimension(:), allocatable :: gpu_LDLT
   type(C_PTR) :: gpu_LDLT_top
   
   integer :: st

   type(cuda_settings_type) :: user_settings
   integer :: this_thread
   logical :: abort

   ! Set GPU device settings as we wish
   call push_ssids_cuda_settings(user_settings, stats(1)%cuda_error)
   if(stats(1)%cuda_error.ne.0) goto 200

   ! Find subtrees to run on each stream
   allocate(stptr(options%nstream+1), stlist(nnodes), stat=stats(1)%st)
   if(stats(1)%st.ne.0) goto 100
   call assign_subtrees(options%nstream, nnodes, child_ptr, child_list, &
      sparent, sptr, rptr, stptr, stlist, options%min_loadbalance, stats(1)%st)
   if(stats(1)%st.ne.0) goto 100

   ! Allocate memory
   allocate(gpu_LDLT(options%nstream), stat=stats(1)%st)
   if(stats(1)%st.ne.0) goto 100

   ! Run subtrees
   abort = .false.
!$OMP PARALLEL DEFAULT(NONE)                                                  &
!$OMP SHARED(abort, alloc, gpu_LDLT, gpu_nlist, gpu_rlist, gpu_rlist_direct,  &
!$OMP    n, nnodes, nodes, options, pos_def, ptr_scale, ptr_val, stats, stptr,&
!$OMP    stlist, stream_data)                                                 &
!$OMP PRIVATE(stream, this_thread)
   this_thread = 1
!$ this_thread = omp_get_thread_num() + 1

!$OMP DO
   do stream = 1, options%nstream
!$OMP FLUSH(abort)
      if(abort) cycle
      !print *, "Running stream ", stream, stptr(stream+1)-stptr(stream)
      allocate(stream_data(stream)%lvllist(nnodes), &
         stream_data(stream)%lvlptr(nnodes + 1), stat=stats(this_thread)%st)
      if(stats(this_thread)%st.ne.0) then
         abort=.true.
!$OMP FLUSH(abort)
         cycle
      endif

      ! FIXME: the below multi-GPU stuff should be made to work
      ! (Maybe just needs sync or something?????)
      !if(omp_get_thread_num() .eq. 0) then
      !   cuda_error = cudaSetDevice(0)
      !   cuda_error = cudaDeviceEnablePeerAccess(2, 0)
      !else
      !   cuda_error = cudaSetDevice(2)
      !   cuda_error = cudaDeviceEnablePeerAccess(0, 0)
      !endif

      call assign_nodes_to_levels_bottom(stptr(stream+1)-stptr(stream), &
         stlist(stptr(stream)), nnodes, child_ptr, child_list, &
         sparent, gpu_contribs, stream_data(stream)%num_levels, &
         stream_data(stream)%lvlptr, stream_data(stream)%lvllist, &
         stats(this_thread)%st)
      if(stats(this_thread)%st.ne.0) then
         abort=.true.
!$OMP FLUSH(abort)
         cycle
      endif
      call subtree_factor_gpu(stream_handle(stream), pos_def, child_ptr, &
         child_list, n, nptr, gpu_nlist, ptr_val, nnodes, nodes, sptr, &
         sparent, rptr, rlist_direct, gpu_rlist, gpu_rlist_direct, &
         gpu_contribs, gpu_LDLT(stream), stream_data(stream), alloc, options, &
         stats(this_thread), ptr_scale)
      if(stats(this_thread)%flag.lt.0) then
         abort=.true.
!$OMP FLUSH(abort)
         cycle
      endif
      if(stats(this_thread)%cuda_error.ne.0) then
         abort=.true.
!$OMP FLUSH(abort)
         cycle
      endif

   end do
!$OMP END DO NOWAIT

   if(stats(this_thread)%st.ne.0) &
      stats(this_thread)%flag = SSIDS_ERROR_ALLOCATION
   if(stats(this_thread)%cuda_error.ne.0) &
      stats(this_thread)%flag = SSIDS_ERROR_CUDA_UNKNOWN
   
!$OMP END PARALLEL ! Implicit barrier
   if(abort) then
      call push_ssids_cuda_settings(user_settings, st)
      return
   endif

   stats(1)%cuda_error = cudaDeviceSynchronize() ! Wait for streams to finish
   if(stats(1)%cuda_error.ne.0) goto 200

   !cuda_error = cudaSetDevice(0)

   ! Run root node
   root = nnodes + 1
   allocate(top_data%lvllist(nnodes), top_data%lvlptr(nnodes + 1), &
      stat=stats(1)%st)
   if(stats(1)%st.ne.0) goto 100
   call assign_nodes_to_levels_top(root, nnodes, child_ptr, child_list, &
      sparent, stptr(options%nstream+1)-1, stlist, &
      top_data%num_levels, top_data%lvlptr, top_data%lvllist, stats(1)%st)
   if(stats(1)%st.ne.0) goto 100
   call subtree_factor_gpu(stream_handle(1), pos_def, child_ptr, child_list, n,&
      nptr, gpu_nlist, ptr_val, nnodes, nodes, sptr, sparent, rptr,            &
      rlist_direct, gpu_rlist, gpu_rlist_direct, gpu_contribs, gpu_LDLT_top,   &
      top_data, alloc, options, stats(1), ptr_scale)
   if(stats(1)%cuda_error.ne.0) goto 200

   ! Free gpu_LDLT as required
   do stream = 1, options%nstream
      if(C_ASSOCIATED(gpu_LDLT(stream))) then
         stats(1)%cuda_error = cudaFree(gpu_LDLT(stream))
         if(stats(1)%cuda_error.ne.0) goto 200
      endif
   end do

   ! Apply any presolve as required
   call perform_presolve(pos_def, child_ptr, child_list, n, nnodes, nodes,    &
      sptr, sparent, rptr, rlist, invp, stream_handle, stream_data, top_data, &
      gpu_rlist_with_delays, gpu_rlist_direct_with_delays, gpu_clists,        &
      gpu_clists_direct, gpu_clen, options, stats(1)%st, stats(1)%cuda_error)
   if(stats(1)%st.ne.0) goto 100
   if(stats(1)%cuda_error.ne.0) goto 200

   ! Restore user GPU device settings
   call push_ssids_cuda_settings(user_settings, stats(1)%cuda_error)
   if(stats(1)%cuda_error.ne.0) goto 200

   return

   100 continue ! Fortran Memory allocation error
   stats(1)%flag = SSIDS_ERROR_ALLOCATION
   call push_ssids_cuda_settings(user_settings, st)
   return

   200 continue ! CUDA Error
   stats(1)%flag = SSIDS_ERROR_CUDA_UNKNOWN
   call push_ssids_cuda_settings(user_settings, st)
   return
end subroutine parfactor

subroutine perform_presolve(pos_def, child_ptr, child_list, n, nnodes, nodes, &
      sptr, sparent, rptr, rlist, invp, stream_handle, stream_data, top_data, &
      gpu_rlist_with_delays, gpu_rlist_direct_with_delays, gpu_clists,        &
      gpu_clists_direct, gpu_clen, options, st, cuda_error)
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer, dimension(*), intent(in) :: sptr
   integer, dimension(*), intent(in) :: sparent
   type(node_type), dimension(nnodes), intent(inout) :: nodes
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   type(C_PTR), dimension(*), intent(in) :: stream_handle
   type(gpu_type), dimension(:), intent(inout) :: stream_data
   type(gpu_type), intent(inout) :: top_data
   type(C_PTR), intent(out) :: gpu_rlist_with_delays
   type(C_PTR), intent(out) :: gpu_rlist_direct_with_delays
   type(C_PTR), intent(out) :: gpu_clists
   type(C_PTR), intent(out) :: gpu_clists_direct
   type(C_PTR), intent(out) :: gpu_clen
   type(ssids_options), intent(in) :: options
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   integer :: stream
   integer(long) :: nrd
   integer, dimension(:), allocatable :: rlist_direct_postfact

   st = 0
   cuda_error = 0

   select case(options%presolve)
   case(0)
      ! Setup data-strctures for normal GPU solve
      call setup_gpu_solve(n, child_ptr, child_list, nnodes, nodes, sparent,  &
         sptr, rptr, rlist, options%nstream, stream_handle, stream_data,      &
         top_data, gpu_rlist_with_delays, gpu_clists, gpu_clists_direct,      &
         gpu_clen, st, cuda_error, gpu_rlist_direct_with_delays)
      if(st.ne.0) return
      if(cuda_error.ne.0) return
   case(1:)
      ! Allocate modified version of rlist_direct
      nrd = rlist_direct_size(nnodes, nodes, rptr)
      allocate(rlist_direct_postfact(2*nrd), stat=st)
      if(st.ne.0) return
      call rebuild_rlist_direct(n, nnodes, nodes, sptr, rptr, rlist, child_ptr,&
         child_list, nrd, rlist_direct_postfact, st)
      if(st.ne.0) return

      ! Setup solves on stream subtrees
      do stream = 1, options%nstream
         call solve_setup(stream_handle(stream), pos_def, sparent, child_ptr, &
            child_list, n, invp, nnodes, nodes, sptr, rptr, rlist,            &
            rlist_direct_postfact, stream_data(stream), st, cuda_error)
         if(st.ne.0) return
         if(cuda_error.ne.0) return
      end do

      ! Wait for streams to finish
      cuda_error = cudaDeviceSynchronize()
      if(cuda_error.ne.0) return

      call solve_setup(stream_handle(1), pos_def, sparent, child_ptr,   &
         child_list, n, invp, nnodes, nodes, sptr, rptr, rlist,         &
         rlist_direct_postfact, top_data, st, cuda_error)
      if(st.ne.0) return
      if(cuda_error.ne.0) return
   end select
end subroutine perform_presolve

!*******************************
!
! This subroutine factorises the subtree(s) that include nodes sa through
! en. Any elements being passed to nodes numbered higher than en are allocated
! using Fortran allocate statemenets rather than stack_alloc.
!
! We maintain the factors seperately from the generated elements to avoid
! copying. Factors are stored in alloc, but pointed to by entries of nodes(:)
! for ease of reference.
!
! Generated elements are stored in a pair of stacks (need two so we can copy
! from child to parent). They are not touched until the factorization has
! been performed on columns that we expect to eliminate.
!
! Entries of A are only added to L just before they are expected to be
! eliminated.
!
subroutine subtree_factor_gpu(stream, pos_def, child_ptr, child_list, n,   &
      nptr, gpu_nlist, ptr_val, nnodes, nodes, sptr, sparent, rptr,        &
      rlist_direct, gpu_rlist, gpu_rlist_direct, gpu_contribs, gpu_LDLT,   &
      gpu, alloc, options, stats, ptr_scale)
   type(C_PTR), intent(in) :: stream ! stream handle to execute on
   logical, intent(in) :: pos_def ! True if problem is supposedly pos-definite
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: nptr
   type(C_PTR), intent(in) :: gpu_nlist
   type(C_PTR), intent(in) :: ptr_val
   ! Note: gfortran-4.3 bug requires explicit size of nodes array
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer, dimension(*), intent(in) :: sparent
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in), target :: rlist_direct
   type(C_PTR), intent(in) :: gpu_rlist
   type(C_PTR), intent(in) :: gpu_rlist_direct
   type(C_PTR), dimension(*), intent(inout) :: gpu_contribs
   type(C_PTR), intent(out) :: gpu_LDLT ! ptr to mem that needs freed
   type(gpu_type), intent(inout) :: gpu
   type(smalloc_type), target, intent(inout) :: alloc ! Contains actual memory
      ! allocations for L. Everything else (within the subtree) is just a
      ! pointer to this.
   ! explicit size required on buf to avoid gfortran-4.3 bug
   type(ssids_options), intent(in) :: options
   type(thread_stats), intent(inout) :: stats
   type(C_PTR), optional, intent(in) :: ptr_scale

   integer :: blkm
   integer :: blkn
   integer :: cblkm
   integer :: cblkn
   integer :: cn
   integer :: cnode
   integer :: i
   integer(long) :: ii
   integer :: j
   integer :: k
   integer :: m
   integer :: ndelay

   integer :: node
   integer, dimension(:), pointer :: lperm

   integer(long) :: level_size
   integer :: level_width, level_height
   integer :: total_nch
   integer :: idata_size, max_idata_size
   integer :: LDLT_size, max_LDLT_size
   integer(long) :: pc_size
   integer :: ncb, max_ncb
   integer :: tile_size
   integer :: maxnelm
   integer :: maxtc
   integer :: ntlev
  
   integer :: nch
   logical :: free_contrib
   integer :: p, q ! general purpose indices
  
   ! per-node assembly info
   type(asmtype), dimension(:), allocatable :: asminf

   integer, allocatable :: iwork(:,:), jwork(:)

   real(wp) :: delta = 0.01, eps = tiny(1.0)
   real(wp), target :: s
   real(wp) :: dummy_real

   ! elimination tree data
   integer :: llist, lev
   integer(long), allocatable :: off_LDLT(:) ! node LDLT contribution offset
      ! in levLDLT

   ! GPU work space (reused for many different things)
   type(cuda_stack_alloc_type) :: gwork
   integer(C_SIZE_T) :: lgpu_work

   ! device pointers
   type(C_PTR), dimension(:), allocatable :: gpu_ldcol
   type(C_PTR) :: ptr_levL, ptr_levLD, ptr_levLDLT
   type(C_PTR) :: ptr_u, ptr_v
   type(C_PTR) :: ptr_cval, ptr_ccval
   type(C_PTR) :: cublas_handle

   ! CUDA-side stats
   type(cuda_stats_type), target :: custats
   type(C_PTR) :: gpu_custats

   gpu_LDLT = C_NULL_PTR

   if(gpu%num_levels.eq.0) return ! Shortcut empty streams (v. small matrices)
  
   gpu%n = n
   gpu%nnodes = nnodes

   ! Initialize CUBLAS handle
   stats%cublas_error = cublasCreate(cublas_handle)
   if(stats%cublas_error.ne.0) goto 300
   stats%cublas_error = cublasSetStream(cublas_handle, stream)
   if(stats%cublas_error.ne.0) goto 300

   ! Initialize CUDA stats
   stats%cuda_error = cudaMalloc(gpu_custats, C_SIZEOF(custats))
   if(stats%cuda_error.ne.0) goto 200
   stats%cuda_error = cudaMemsetAsync(gpu_custats, 0, C_SIZEOF(custats), stream)
   if(stats%cuda_error.ne.0) goto 200

   ! Precalculate level information
   max_LDLT_size = 0
   max_ncb = 2
   max_idata_size = 0
   do lev = 1, gpu%num_levels
      total_nch = 0
      LDLT_size = 0
      idata_size = 0
      p = 0
      q = 0
      do llist = gpu%lvlptr(lev), gpu%lvlptr(lev + 1) - 1
         node = gpu%lvllist(llist)
         blkn = sptr(node+1) - sptr(node)
         blkm = int(rptr(node+1) - rptr(node))
         nch = child_ptr(node + 1) - child_ptr(node)
         total_nch = total_nch + nch
         m = blkm - blkn
         LDLT_size = LDLT_size + m*m
         k = (m - 1)/BLOCK_SIZE + 1
         p = p + (k*(k + 1))/2
         q = q + (blkm - 1)/BLOCK_SIZE + 2
      end do
      ncb = gpu%lvlptr(lev + 1) - gpu%lvlptr(lev)
      idata_size = max(10*p, 9*q, total_nch)
      max_LDLT_size = max(max_LDLT_size, LDLT_size)
      max_ncb = max(max_ncb, ncb)
      max_idata_size = max(max_idata_size, idata_size)
   end do
  
   ii = nptr(nnodes + 1) - 1
   stats%cuda_error = cudaMalloc(gpu_LDLT, 2*max_LDLT_size*C_SIZEOF(dummy_real))
   if(stats%cuda_error.ne.0) goto 200
  
   allocate(gpu%values_L(gpu%num_levels), gpu%off_L(nnodes), &
      off_LDLT(nnodes), asminf(nnodes), stat=stats%st)
   if(stats%st.ne.0) goto 100

   do node = 1, nnodes
      blkm = int(rptr(node + 1) - rptr(node))
      blkn = sptr(node + 1) - sptr(node)
      do cn = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(cn)
         cblkn = sptr(cnode + 1) - sptr(cnode)
         cblkm = int(rptr(cnode + 1) - rptr(cnode))
         m = 0
         do ii = rptr(cnode) + cblkn, rptr(cnode + 1) - 1
            if ( rlist_direct(ii) > blkn ) exit
            m = m + 1
         end do
         asminf(cnode)%npassed = cblkm - cblkn
         asminf(cnode)%npassl = m
         asminf(cnode)%offset = rptr(cnode) + cblkn - 1
      end do
   end do

   !
   ! Loop over levels doing work
   !
   do lev = 1, gpu%num_levels
      lgpu_work = level_gpu_work_size(lev, gpu%lvlptr, gpu%lvllist, child_ptr, &
         child_list, nodes, sptr, rptr, asminf)
      call custack_init(gwork, lgpu_work, stats%cuda_error)
      if(stats%cuda_error.ne.0) goto 200

      ncb = gpu%lvlptr(lev + 1) - gpu%lvlptr(lev)

      !
      ! Initialize level information
      !
      level_size = 0
      level_width = 0
      level_height = 0
      total_nch = 0
      pc_size = 0
      do llist = gpu%lvlptr(lev), gpu%lvlptr(lev + 1) - 1
         node = gpu%lvllist(llist)
         ndelay = nodes(node)%ndelay
         blkn = sptr(node+1) - sptr(node) + ndelay
         blkm = int(rptr(node+1) - rptr(node)) + ndelay
      
         gpu%off_L(node) = level_size
         off_LDLT(node) = pc_size

         nch = child_ptr(node+1) - child_ptr(node)
         level_size = level_size + (blkm + 2_long)*blkn
         level_width = level_width + blkn
         level_height = max(level_height, blkm)
         total_nch = total_nch + nch
         m = blkm - blkn
         pc_size = pc_size + m*(m+0_long)
      end do
      
      !
      ! Generate pointers for this level
      !
      if ( mod(lev, 2) > 0 ) then
         ptr_levLDLT = gpu_LDLT
      else
         ptr_levLDLT = c_ptr_plus(gpu_LDLT, max_LDLT_size*C_SIZEOF(dummy_real))
      end if

      stats%cuda_error = &
         cudaMalloc(gpu%values_L(lev)%ptr_levL, level_size*C_SIZEOF(dummy_real))
      if(stats%cuda_error.ne.0) goto 200
      ptr_levL = gpu%values_L(lev)%ptr_levL
      if(.not. pos_def) then
         stats%cuda_error = &
            cudaMalloc(ptr_levLD, level_size*C_SIZEOF(dummy_real))
         if(stats%cuda_error.ne.0) goto 200

         ! Initialize pointers to LD storage
         if(allocated(gpu_ldcol)) deallocate(gpu_ldcol, stat=stats%st)
         if(stats%st.ne.0) goto 100
         allocate(gpu_ldcol(gpu%lvlptr(lev+1)-gpu%lvlptr(lev)), stat=stats%st)
         if(stats%st.ne.0) goto 100
         do llist = gpu%lvlptr(lev), gpu%lvlptr(lev + 1) - 1
            node = gpu%lvllist(llist)        
            i = llist - gpu%lvlptr(lev) + 1
            gpu_ldcol(i) = &
               c_ptr_plus(ptr_levLD, gpu%off_L(node)*C_SIZEOF(dummy_real))
         end do
      endif

      ! Set up node pointers
      level_size = 0
      do llist = gpu%lvlptr(lev), gpu%lvlptr(lev + 1) - 1
         node = gpu%lvllist(llist)
         ndelay = nodes(node)%ndelay
         blkn = sptr(node+1) - sptr(node) + ndelay
         blkm = int(rptr(node+1) - rptr(node)) + ndelay
      
         nodes(node)%gpu_lcol = c_ptr_plus(gpu%values_L(lev)%ptr_levL, &
            level_size*C_SIZEOF(dummy_real))

         level_size = level_size + (blkm + 2_long)*blkn
      end do

      ! Allocate+Initialize lperm for fronts on this level
      do llist = gpu%lvlptr(lev), gpu%lvlptr(lev + 1) - 1

         node = gpu%lvllist(llist)
         ndelay = nodes(node)%ndelay
         blkn = sptr(node+1) - sptr(node) + ndelay
         blkm = int(rptr(node+1) - rptr(node)) + ndelay
    
         stats%maxfront = max(stats%maxfront, blkn)

         ! Allocate memory for the local permutation (lperm).
         call smalloc(alloc, nodes(node)%perm, blkn+0_long, &
            nodes(node)%ismptr, nodes(node)%ismsa, stats%st)
         if (stats%st .ne. 0) go to 100
         lperm => nodes(node)%perm

         ! Initialise lperm
         j = ndelay + 1
         do i = sptr(node), sptr(node+1)-1
            lperm(j) = i
            j = j + 1
         end do
    
      end do
      
      ! Initialize L to 0, and add A to it.
      call init_L_with_A(stream, lev, gpu%lvlptr, gpu%lvllist, nodes, ncb, &
         level_size, nptr, rptr, gpu_nlist, gpu_rlist, ptr_val, ptr_levL, &
         gwork, stats%st, stats%cuda_error, ptr_scale=ptr_scale)
      if(stats%st.ne.0) goto 100
      if(stats%cuda_error.ne.0) goto 200

      ! The pivot is considered to be zero if less than epsilon(ONE) times
      ! the maximum value of A on this level: must use relative threshold
      ! in case the user opts for no scaling
      ! FIXME: Make optional?
      ! FIXME: Use CUBLAS instead? [Check if zeroed]
      k = int( min(65535_long, (level_size - 1)/256 + 1) )
      ptr_u = custack_alloc(gwork, k*C_SIZEOF(dummy_real))
      ptr_v = custack_alloc(gwork, C_SIZEOF(dummy_real))
      call max_abs( stream, k, level_size, ptr_levL, ptr_u, ptr_v )
      stats%cuda_error = cudaMemcpyAsync_D2H(C_LOC(s), ptr_v, C_SIZEOF(s), &
         stream)
      if(stats%cuda_error.ne.0) goto 200
      stats%cuda_error = cudaStreamSynchronize(stream) ! Wait for ptr_u, ptr_v
      if(stats%cuda_error.ne.0) goto 200
      call custack_free(gwork, C_SIZEOF(dummy_real)) ! ptr_v
      call custack_free(gwork, k*C_SIZEOF(dummy_real)) ! ptr_u
      eps = s*epsilon(ONE)

      !
      ! Assemble fully summed columns
      !
      call assemble_fully_summed(stream, total_nch, lev, gpu%lvlptr,       &
         gpu%lvllist, nodes, ptr_ccval, gpu_contribs, ptr_levL,            &
         gpu_rlist_direct, child_ptr,  child_list, off_LDLT, asminf, rptr, &
         sptr, gwork, stats%st, stats%cuda_error)
      if(stats%st.ne.0) goto 100
      if(stats%cuda_error.ne.0) goto 200
    
      !
      ! Perform factorization (of fully summed columns)
      !
      if(pos_def) then
         call factor_posdef(stream, lev, gpu%lvlptr, nodes, gpu%lvllist, &
            sptr, rptr, ptr_levL, cublas_handle, stats, gwork)
      else
         call factor_indef(stream, lev, gpu%lvlptr, nodes, gpu%lvllist, &
            sparent, sptr, rptr, level_height, level_width, delta, eps, &
            gpu_ldcol, gwork, cublas_handle, options, stats, gpu_custats)
      endif
      if(stats%flag.lt.0) goto 20
      if(stats%st.ne.0) goto 100
      if(stats%cuda_error.ne.0) goto 200
      if(stats%cublas_error.ne.0) goto 300

      !
      ! Form contribution block (of non-fully summed columns)
      !
      if(pc_size > 0) then
         if(pos_def) then
            call form_contrib(stream, lev, gpu%lvlptr, nodes, gpu%lvllist, &
               off_LDLT, sptr, rptr, ptr_levLDLT, gwork, stats%st,         &
               stats%cuda_error)
         else
            call form_contrib(stream, lev, gpu%lvlptr, nodes, gpu%lvllist, &
               off_LDLT, sptr, rptr, ptr_levLDLT, gwork, stats%st,         &
               stats%cuda_error, gpu_ldcol=gpu_ldcol)
         endif
         if(stats%st.ne.0) goto 100
         if(stats%cuda_error.ne.0) goto 200
      endif
      if(stats%flag < 0) goto 20

      !
      ! Assemble children into contribution block
      !
      call assemble_contrib(stream, total_nch, lev, gpu%lvlptr, gpu%lvllist, &
         child_ptr, child_list, sptr, rptr, asminf, pc_size, &
         off_LDLT, ptr_ccval, gpu_contribs, ptr_levLDLT, gpu_rlist_direct, &
         gwork, stats%st, stats%cuda_error)
      if(stats%st.ne.0) goto 100
      if(stats%cuda_error.ne.0) goto 200

      ! Free allocs specific to this level
      if(.not. pos_def) then
         stats%cuda_error = cudaFree( ptr_levLD )
         if(stats%cuda_error.ne.0) goto 200
      endif

      ! Store pointers for use on next level
      ptr_cval = ptr_levL
      ptr_ccval = ptr_levLDLT
    
   end do ! lev

   ! Free stack memory
   call custack_finalize(gwork, stats%cuda_error)
   if(stats%cuda_error.ne.0) goto 200

   ! Copy GPU stats back to host and free it
   stats%cuda_error = cudaMemcpy_d2h(C_LOC(custats), gpu_custats, &
      C_SIZEOF(custats))
   if(stats%cuda_error.ne.0) goto 200
   stats%cuda_error = cudaFree(gpu_custats)
   if(stats%cuda_error.ne.0) goto 200
   stats%num_zero = custats%num_zero
   stats%num_neg = custats%num_neg
   stats%num_two = custats%num_two
  
   if ( options%presolve > 0 ) then
      gpu%presolve = options%presolve
      tile_size = 24
      maxnelm = max_nelim(nnodes, nodes, gpu%num_levels, gpu%lvlptr, gpu%lvllist)
      maxtc = (maxnelm - 1)/tile_size + 1
      m = gpu%lvlptr(gpu%num_levels + 1) - 1
      allocate(iwork(maxtc + 1, 2), jwork(m), stat=stats%st)
      if(stats%st.ne.0) goto 100
      call presolve_lwork( nnodes, nodes, sptr, rptr, gpu%num_levels, &
         gpu%lvlptr, gpu%lvllist, tile_size, maxtc, iwork, ntlev, m )
      stats%cuda_error = cudaMalloc(ptr_u, m*C_SIZEOF(dummy_real))
      if(stats%cuda_error.ne.0) go to 200
      call presolve_first( stream, nnodes, nodes, rptr, &
         gpu%num_levels, gpu%lvlptr, &
         gpu%lvllist, tile_size, ntlev, iwork, maxtc, iwork(1,2), jwork, &
         pos_def, cublas_handle, ptr_u, stats%st, stats%cuda_error, &
         stats%cublas_error )
      if(stats%st.ne.0) goto 100
      if(stats%cuda_error.ne.0) goto 200
      if(stats%cublas_error.ne.0) goto 300
      call presolve_second( stream, nnodes, nodes, sptr, rptr, &
         gpu%num_levels, gpu%values_L, &
         gpu%lvlptr, gpu%lvllist, tile_size, ntlev, iwork, maxtc, iwork(1,2), &
         ptr_u, stats%st, stats%cuda_error)
      if(stats%st.ne.0) goto 100
      if(stats%cuda_error.ne.0) go to 200
      deallocate(iwork, jwork, stat=stats%st)
      if(stats%st.ne.0) goto 100
      stats%cuda_error = cudaFree( ptr_u )
      if(stats%cuda_error.ne.0) go to 200
   end if

   20 continue ! start of cleanup

   free_contrib = .true.
   do llist = gpu%lvlptr(gpu%num_levels), gpu%lvlptr(gpu%num_levels+1)-1
      node = gpu%lvllist(llist)
      if(sparent(node) .le. nnodes) then
         gpu_contribs(node) = &
            c_ptr_plus(ptr_ccval, off_LDLT(node)*C_SIZEOF(dummy_real))
         free_contrib = .false.
      endif
   end do
   if(free_contrib) then
      stats%cuda_error = cudaFree(gpu_LDLT)
      if(stats%cuda_error.ne.0) goto 200
      gpu_LDLT = C_NULL_PTR
   endif
  
   deallocate(off_LDLT, asminf, stat=stats%st)
   if(stats%st.ne.0) goto 100
  
   ! Destroy CUBLAS handle
   stats%cublas_error = cublasDestroy(cublas_handle)
   if(stats%cuda_error.ne.0) goto 300

   return ! Normal return

   100 continue ! Fortran Memory allocation error
   stats%flag = SSIDS_ERROR_ALLOCATION
   return

   200 continue ! CUDA failure
   stats%flag = SSIDS_ERROR_CUDA_UNKNOWN
   return

   300 continue ! CUBLAS failure
   stats%flag = SSIDS_ERROR_CUBLAS_UNKNOWN
   return
   
end subroutine subtree_factor_gpu

! Following routine calculates size of gwork required for a given level.
! gwork size is made up of maximum of space required for each individual
! work routine
integer(C_SIZE_T) function level_gpu_work_size(lev, lvlptr, lvllist, &
      child_ptr, child_list, nodes, sptr, rptr, asminf) result(lgpu)
   integer, intent(in) :: lev
   integer, dimension(*), intent(in) :: lvlptr
   integer, dimension(*), intent(in) :: lvllist
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(asmtype), dimension(:), intent(in) :: asminf

   integer(C_SIZE_T) :: sz, sz2, sz3, sz4, sz5, sz6, sz7, sz8
   integer :: li, ci
   integer :: ndelay, blkm, blkn
   integer(long) :: bx, by
   integer :: node, cnode

   ! Dummy datatypes to get sizes
   type(load_nodes_type) :: lnt_dummy
   type(assemble_cp_type) :: acpt_dummy
   type(assemble_blk_type) :: abt_dummy
   type(assemble_delay_type) :: adt_dummy
   type(multisymm_type) :: mst_dummy
   type(multiswap_type) :: mswt_dummy
   type(multisyrk_type) :: msyrt_dymmy
   type(multinode_fact_type) :: mnft_dummy
   type(multiblock_fact_type) :: mbft_dummy
   type(multireorder_data) :: mr_dummy
   type(multielm_data) :: me_dummy
   type(cstat_data_type) :: cdt_dummy
   integer(C_INT) :: int_dummy
   real(C_DOUBLE) :: real_dummy

   ! Initialize space required to 0
   lgpu = 0

   ! Space for lndata in init_L_with_A()
   sz = lvlptr(lev+1)-lvlptr(lev)
   lgpu = max(lgpu, sz*C_SIZEOF(lnt_dummy))

   ! Space for ptr_u, ptr_v with cuda_max_abs()
   sz = 0
   do li = lvlptr(lev), lvlptr(lev+1)-1
      node = lvllist(li)
      ndelay = nodes(node)%ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      sz = sz + (blkm+2_long)*blkn ! calculate 'level_size'
   end do
   sz = int( min(65535_long, (sz - 1)/256 + 1) )
   lgpu = max(lgpu, &
      sz*C_SIZEOF(real_dummy) + & ! ptr_u
      C_SIZEOF(real_dummy))       ! ptr_v

   ! Space for cpdata, blkdata, ddata and sync in assemble_fully_summed()
   sz = 0
   sz2 = 0
   do li = lvlptr(lev), lvlptr(lev+1)-1
      node = lvllist(li)
      sz = sz + child_ptr(node+1)-child_ptr(node)
      do ci = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(ci)
         bx = (asminf(cnode)%npassl-1) / HOGG_ASSEMBLE_TX + 1
         by = (asminf(cnode)%npassed-1) / HOGG_ASSEMBLE_TY + 1
         sz2 = sz2 + bx*by
      end do
   end do
   lgpu = max(lgpu, sz*C_SIZEOF(acpt_dummy) + sz2*C_SIZEOF(abt_dummy) + &
      sz*C_SIZEOF(adt_dummy) + (1+sz)*C_SIZEOF(int_dummy))

   ! Space for msdata in factor_posdef()
   sz = lvlptr(lev+1)-lvlptr(lev) ! number required for msdata
   lgpu = max(lgpu, sz*C_SIZEOF(mst_dummy))

   ! Space for swapdata in factor_indef()
   lgpu = max(lgpu, sz*C_SIZEOF(mswt_dummy))

   ! Space for msdata, ptr_ind and ptr_B in factor_indef()/factor_posdef()
   ! Also gpu_aux, gpu_perm, gpu_rdata, gpu_mdata, gpu_mnfdata, gpu_mbfdata
   !    in multinode_ldlt() / multinode_llt()
   ! Also gpu_r in node_ldlt()
   ! Also gpu_csdata in collect_stats_indef()
   sz = 0 ! number required for ptr_B
   sz2 = 0 ! number required for ptr_ind
   sz4 = 0 ! number required for gpu_perm
   sz5 = 0 ! number required for gpu_rdata
   sz6 = 0 ! number required for gpu_mdata
   sz7 = 0 ! number required for gpu_mbfdata
   sz8 = 0 ! number required for gpu_r
   do li = lvlptr(lev), lvlptr(lev+1)-1
      node = lvllist(li)
      ndelay = nodes(node)%ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      sz = sz + blkm
      sz2 = sz2 + blkn
      sz4 = sz4 + blkn
      sz5 = sz5 + (blkm-1)/(32*BLOCK_SIZE) + 2
      sz6 = sz6 + ((blkm-1)/32+1) * ((blkn-1)/32+1)
      sz7 = sz7 + (blkm - 1)/(BLOCK_SIZE*(MNF_BLOCKS - 1)) + 1
      if(blkm.eq.blkn .or. lvlptr(lev+1)-lvlptr(lev).eq.1) &
         sz8 = max(sz8, blkm*(blkm+0_C_SIZE_T))
   end do
   sz2 = max(int(sz2), (lvlptr(lev+1)-lvlptr(lev))*BLOCK_SIZE)
   sz = sz*BLOCK_SIZE
   sz3 = lvlptr(lev+1)-lvlptr(lev) ! number of blocks for msdata et al
   lgpu = max(lgpu, &
              max( sz*C_SIZEOF(real_dummy) + & ! ptr_B
                      sz2*C_SIZEOF(int_dummy) + & ! ptr_ind
                      8*C_SIZEOF(int_dummy) + & ! gpu_aux
                      sz4*C_SIZEOF(int_dummy) + & ! gpu_perm
                      sz5*4*C_SIZEOF(int_dummy) + & ! gpu_rdata
                      sz5*C_SIZEOF(mr_dummy) + & ! gpu_rdata
                      sz6*2*C_SIZEOF(int_dummy) + & ! gpu_mdata
                      sz6*C_SIZEOF(me_dummy) + & ! gpu_mdata
                      sz3*C_SIZEOF(mnft_dummy) + & ! gpu_mnfdata
                      sz3*C_SIZEOF(int_dummy) + & ! gpu_stat
                      sz7*C_SIZEOF(mbft_dummy) + & ! gpu_mbfdata
                      sz8*C_SIZEOF(real_dummy), & ! gpu_r
                   sz3*C_SIZEOF(cdt_dummy) & ! gpu_csdata
                   ) &
          )

   ! Space for gpu_msdata in form_contrib()
   sz = 0 ! spac for gpu_msdata
   do li = lvlptr(lev), lvlptr(lev+1)-1
      node = lvllist(li)
      ndelay = nodes(node)%ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      bx = (blkm-blkn-1)/32 + 1
      sz = sz + (bx*(bx+1))/2
   end do
   lgpu = max(lgpu, sz*C_SIZEOF(msyrt_dymmy))

   ! Space for cpdata, blkdata and sync in assemble_contrib()
   sz = 0
   sz2 = 0
   do li = lvlptr(lev), lvlptr(lev+1)-1
      node = lvllist(li)
      sz = sz + child_ptr(node+1)-child_ptr(node)
      do ci = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(ci)
         if(asminf(cnode)%npassed-asminf(cnode)%npassl .le. 0) cycle
         bx = (asminf(cnode)%npassed-asminf(cnode)%npassl-1) / &
            HOGG_ASSEMBLE_TX + 1
         by = (asminf(cnode)%npassed-asminf(cnode)%npassl-1) / &
            HOGG_ASSEMBLE_TY + 1
         sz2 = sz2 + bx*by
      end do
   end do
   lgpu = max(lgpu, sz*C_SIZEOF(acpt_dummy) + sz2*C_SIZEOF(abt_dummy) + &
      (1+sz)*C_SIZEOF(int_dummy))

   ! Allow for alignment retification of up to 10 pointers beyond the first
   lgpu = lgpu + 10*256
end function level_gpu_work_size

! Assigns nodes to level lists, with level 1 being the closest to the leaves
! This version handles the top of the tree, stopping at the level set L_th
! and not proceeding below it.
subroutine assign_nodes_to_levels_top(root, nnodes, child_ptr, child_list, &
      sparent, nL_th, L_th, num_levels, lvlptr, lvllist, st)
   integer, intent(in) :: root
   integer, intent(in) :: nnodes
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, dimension(*), intent(in) :: sparent
   integer, intent(in) :: nL_th
   integer, dimension(*), intent(in) :: L_th
   integer, intent(out) :: num_levels
   integer, dimension(*), intent(out) :: lvlptr
   integer, dimension(*), intent(out) :: lvllist
   integer, intent(out) :: st

   integer :: sa
   integer :: node, lvl, ri
   integer, dimension(:), allocatable :: level ! level of node
   integer, dimension(:), allocatable :: lvlcount

   logical, dimension(:), allocatable :: dead

   ! Find first child of subtree
   sa = root
   do while(child_ptr(sa+1)-child_ptr(sa).gt.0)
      sa = child_list(child_ptr(sa))
   end do

   allocate(level(sa:root), lvlcount(root-sa+1), dead(sa:root), stat=st)
   if(st.ne.0) return
   dead(:) = .false.

   ! Find level of each node, with level 1 being a root
   lvlcount(:) = 0
   if(root.le.nnodes) then
      ! Count root
      level(root) = 1
      lvlcount(1) = 1
      num_levels = 1
   else
      ! Don't count root
      level(root) = 0
      num_levels = 0
   endif
   ! Mark as dead all nodes in L_th
   do ri = 1, nL_th
      node = L_th(ri)
      dead(node) = .true.
   end do
   do node = root-1, sa, -1
      ! Mark as dead any nodes in subtrees below L_th
      if(dead(sparent(node))) dead(node) = .true.
      if(dead(node)) cycle
      ! Record non-dead nodes
      lvl = level(sparent(node)) + 1
      level(node) = lvl
      lvlcount(lvl) = lvlcount(lvl) + 1
      num_levels = max(num_levels, lvl)
   end do

   ! Carefully remove any virtual root we have
   dead(nnodes+1:root) = .true.

   ! Setup pointers, note that the final level we want for each node is
   ! num_levels-level(node)+1 as we number from the bottom, not the top!
   ! We use lvlptr(i+1) as the insert position for level i
   lvlptr(1:2) = 1
   do lvl = 2, num_levels
      lvlptr(lvl+1) = lvlptr(lvl) + lvlcount(num_levels-(lvl-1)+1)
   end do

   ! Finally assign nodes to levels
   do node = sa, root
      if(dead(node)) cycle
      lvl = num_levels - level(node) + 1
      lvllist(lvlptr(lvl+1)) = node
      lvlptr(lvl+1) = lvlptr(lvl+1) + 1
   end do
end subroutine assign_nodes_to_levels_top

! Assigns nodes to level lists, with level 1 being the closest to the leaves
! This version handles multiple roots below the level set L_th.
subroutine assign_nodes_to_levels_bottom(nroot, roots, nnodes, child_ptr, &
      child_list, sparent, gpu_contribs, num_levels, lvlptr, lvllist, st)
   integer, intent(in) :: nroot
   integer, dimension(nroot), intent(in) :: roots
   integer, intent(in) :: nnodes
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, dimension(*), intent(in) :: sparent
   type(C_PTR), dimension(*), intent(in) :: gpu_contribs
   integer, intent(out) :: num_levels
   integer, dimension(*), intent(out) :: lvlptr
   integer, dimension(*), intent(out) :: lvllist
   integer, intent(out) :: st

   integer, dimension(:), allocatable :: sa
   integer :: node, lvl, ri, root
   integer, dimension(:), allocatable :: level ! level of node
   integer, dimension(:), allocatable :: lvlcount

   logical, dimension(:), allocatable :: dead

   if(nroot.eq.0) then
      num_levels = 0
      return
   endif

   allocate(level(nnodes+1), lvlcount(nnodes+1), dead(nnodes+1), stat=st)
   if(st.ne.0) return

   ! Find first children of subtrees
   allocate(sa(nroot), stat=st)
   if(st.ne.0) return
   do ri = 1, nroot
      sa(ri) = roots(ri)
      do while(child_ptr(sa(ri)+1)-child_ptr(sa(ri)).gt.0)
         sa(ri) = child_list(child_ptr(sa(ri)))
      end do
   end do

   ! Find level of each node, with level 1 being a root
   num_levels = 0
   dead(:) = .false.
   lvlcount(:) = 0
   do ri = 1, nroot
      root = roots(ri)
      if(root.le.nnodes) then
         ! Count root
         level(root) = 1
         lvlcount(1) = lvlcount(1) + 1
         num_levels = max(num_levels, 1)
      else
         ! Don't count root
         level(root) = 0
      endif
      do node = root-1, sa(ri), -1
         ! Mark as dead nodes in subtrees rooted at nodes with defined contrib
         if(C_ASSOCIATED(gpu_contribs(node))) dead(node) = .true.
         if(dead(sparent(node))) dead(node) = .true.
         if(dead(node)) cycle
         ! Record non-dead nodes
         lvl = level(sparent(node)) + 1
         level(node) = lvl
         lvlcount(lvl) = lvlcount(lvl) + 1
         num_levels = max(num_levels, lvl)
      end do
   end do

   ! Remove any virtual root we have
   dead(nnodes+1) = .true.

   ! Setup pointers, note that the final level we want for each node is
   ! num_levels-level(node)+1 as we number from the bottom, not the top!
   ! We use lvlptr(i+1) as the insert position for level i
   lvlptr(1:2) = 1
   do lvl = 2, num_levels
      lvlptr(lvl+1) = lvlptr(lvl) + lvlcount(num_levels-(lvl-1)+1)
   end do

   ! Finally assign nodes to levels
   do ri = 1, nroot
      root = roots(ri)
      do node = sa(ri), root
         if(dead(node)) cycle
         lvl = num_levels - level(node) + 1
         lvllist(lvlptr(lvl+1)) = node
         lvlptr(lvl+1) = lvlptr(lvl+1) + 1
      end do
   end do
end subroutine assign_nodes_to_levels_bottom

subroutine assign_subtrees(nstream, nnodes, child_ptr, child_list, sparent, &
      sptr, rptr, stptr, stlist, min_loadbalance, st)
   integer, intent(in) :: nstream
   integer, intent(in) :: nnodes
   integer, intent(in) :: child_ptr(*) ! Pointers into child_list for node
   integer, intent(in) :: child_list(*) ! Children of node node are given by:
      ! child_list(child_ptr(node):child_ptr(node+1)-1)
   integer, dimension(*), intent(in) :: sparent
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(out) :: stptr
   integer, dimension(*), intent(out) :: stlist
   real :: min_loadbalance
   integer, intent(out) :: st

   integer :: node
   integer, dimension(:), allocatable :: level ! level of node
   integer, dimension(:), allocatable :: lvlcount

   type(ntype), dimension(:), allocatable :: nwork
   integer, dimension(:), allocatable :: sthead, stnext
   integer(long), dimension(:), allocatable :: stwork
   integer :: n_lth, nbelow, nchild, stream
   integer, dimension(:), allocatable :: L_th
   real :: loadbalance

   integer :: i, j
   integer(long) :: jj, blkm, blkn, w

   allocate(level(nnodes+1), lvlcount(nnodes), nwork(nnodes+1), stat=st)
   if(st.ne.0) return

   ! Calculate work at and below each node, exploiting postorder
   do node = 1, nnodes
      blkm = rptr(node+1) - rptr(node)
      blkn = sptr(node+1) - sptr(node)
      w = 0
      do jj = blkm, blkm-blkn+1, -1
         w = w + jj**2
      end do
      nwork(node)%work_here = w
      nwork(node)%work_below = nwork(node)%work_below + w
      j = sparent(node)
      nwork(j)%work_below = nwork(j)%work_below + nwork(node)%work_below
      nwork(j)%level = max(nwork(j)%level, nwork(node)%level + 1)
   end do

   ! Find L_th
   allocate(L_th(nnodes), sthead(nstream), stnext(nnodes), stwork(nstream), &
      stat=st)
   if(st.ne.0) return
   n_lth = 1
   nbelow = nnodes+1
   L_th(1) = nnodes+1
   loadbalance = 0.0
   do while(loadbalance < min_loadbalance .and. nbelow>0)
      ! Split largest node into its children
      nbelow = nbelow - 1
      node = L_th(1)
      nchild = child_ptr(node+1)-child_ptr(node)
      if(nchild.eq.0) then
         do i = 2, n_lth
            L_th(i-1) = L_th(i)
         end do
         n_lth = n_lth - 1
      else
         L_th(1) = child_list(child_ptr(node))
         do i = child_ptr(node)+1, child_ptr(node+1)-1
            n_lth = n_lth + 1
            L_th(n_lth) = child_list(i)
         end do
      endif
      ! Sort nodes by work
      call sort_nodes(n_lth, L_th, nwork)
      ! Allocate nodes to streams
      stwork(:) = 0
      sthead(:) = -1
      do i = 1, n_lth
         node = L_th(i)
         ! Find stream with least work
         stream = 1
         do j = 2, nstream
            if(stwork(j) < stwork(stream)) stream = j
         end do
         stnext(node) = sthead(stream)
         sthead(stream) = node
         stwork(stream) = stwork(stream) + nwork(node)%work_below
      end do
      loadbalance = real(minval(stwork(:))) / real(maxval(stwork(:)))
   end do

   ! Store stream's work
   stptr(1) = 1
   do stream = 1, nstream
      stptr(stream+1) = stptr(stream)
      node = sthead(stream)
      do while(node.ne.-1)
         stlist(stptr(stream+1)) = node
         stptr(stream+1) = stptr(stream+1) + 1
         node = stnext(node)
      end do
   end do

end subroutine assign_subtrees

subroutine sort_nodes(n_lth, L_th, nwork)
   integer, intent(in) :: n_lth
   integer, dimension(*), intent(inout) :: L_th
   type(ntype), dimension(*), intent(in) :: nwork

   integer :: i, nxchg, n1, n2

   nxchg = 1 ! Force a first pass
   do while(nxchg > 0)
      nxchg = 0
      do i = 1, n_lth-1
         n1 = L_th(i)
         n2 = L_th(i+1)
         if(nwork(n1)%work_below < nwork(n2)%work_below) then
            ! swap n1 and n2
            L_th(i+1) = n1
            L_th(i) = n2
            nxchg = nxchg + 1
         endif
      end do
   end do
end subroutine sort_nodes

subroutine init_L_with_A(stream, lev, lvlptr, lvllist, nodes, ncb, level_size, &
      nptr, rptr, gpu_nlist, gpu_rlist, ptr_val, ptr_levL, &
      gwork, st, cuda_error, ptr_scale)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: lev
   integer, dimension(*), intent(in) :: lvlptr
   integer, dimension(*), intent(in) :: lvllist
   type(node_type), dimension(*), intent(in) :: nodes
   integer, intent(in) :: ncb
   integer(long), intent(in) :: level_size
   integer, dimension(*), intent(in) :: nptr
   integer(long), dimension(*), intent(in) :: rptr
   type(C_PTR), intent(in) :: gpu_nlist
   type(C_PTR), intent(in) :: gpu_rlist
   type(C_PTR), intent(in) :: ptr_val
   type(C_PTR), intent(in) :: ptr_levL ! target data is altered
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error
   type(C_PTR), optional, intent(in) :: ptr_scale

   integer :: llist, node, i
   type(load_nodes_type), dimension(:), allocatable, target :: lndata
   type(C_PTR) :: gpu_lndata
   real(wp) :: dummy_real

   st = 0
   cuda_error = 0

   ! Initialize data for cuda_load_nodes and copy to GPU
   allocate(lndata(ncb), stat=st)
   if(st.ne.0) return
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)
      i = llist - lvlptr(lev) + 1
      lndata(i)%offn = nptr(node) - 1
      lndata(i)%nnz = nptr(node + 1) - nptr(node)
      lndata(i)%lda = int(rptr(node + 1) - rptr(node))
      lndata(i)%offr = rptr(node) - 1
      lndata(i)%ldl = int(rptr(node + 1) - rptr(node)) + nodes(node)%ndelay
      lndata(i)%lcol = c_ptr_plus( nodes(node)%gpu_lcol, &
         nodes(node)%ndelay * (1+lndata(i)%ldl) * C_SIZEOF(dummy_real) )
   end do
   gpu_lndata = custack_alloc(gwork, ncb*C_SIZEOF(lndata(1)))
   cuda_error = cudaMemcpyAsync_H2D(gpu_lndata, C_LOC(lndata), &
      ncb*C_SIZEOF(lndata(1)), stream)
   if(cuda_error.ne.0) return

   ! Initialize frontal matrices to 0
   cuda_error = &
      cudaMemsetAsync(ptr_levL, 0, level_size*C_SIZEOF(dummy_real), stream)
   if(cuda_error.ne.0) return
 
   ! Store values of A into front for fully summed variables
   if ( present(ptr_scale) ) then
      call load_nodes_sc( stream, ncb, gpu_lndata, gpu_nlist, gpu_rlist,&
         ptr_scale, ptr_val )
   else
      call load_nodes( stream, ncb, gpu_lndata, gpu_nlist, ptr_val )
   end if
   call custack_free(gwork, ncb*C_SIZEOF(lndata(1)))
end subroutine init_L_with_A

subroutine form_contrib(stream, lev, lvlptr, nodes, lvllist, off_LDLT,&
      sptr, rptr, ptr_levLDLT, gwork, st, cuda_error, gpu_ldcol)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: lev
   integer, dimension(*), intent(in) :: lvlptr
   type(node_type), dimension(*), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: lvllist
   integer(long), dimension(*), intent(in) :: off_LDLT
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   type(C_PTR), intent(in) :: ptr_levLDLT
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error
   type(C_PTR), dimension(*), optional, intent(in) :: gpu_ldcol

   type(multisyrk_type), dimension(:), allocatable, target :: msdata
   type(C_PTR) :: gpu_msdata
   
   integer :: i, j, k, m
   integer :: llist, ncb, nn
   integer :: ndelay, nelim, blkm, blkn, node
   real(wp) :: dummy_real

   cuda_error = 0
   st = 0

   ncb = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)
      ndelay = nodes(node)%ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      m = blkm - blkn
      k = (m - 1)/32 + 1
      ncb = ncb + (k*(k + 1))/2
   end do
   allocate(msdata(ncb), stat=st)
   if(st.ne.0) return

   ncb = 0
   nn = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      i = llist - lvlptr(lev) + 1
      node = lvllist(llist)
      ndelay = nodes(node)%ndelay
      nelim = nodes(node)%nelim
      blkn = sptr(node + 1) - sptr(node) + ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      m = blkm - blkn
      nn = nn + 1
      k = (m - 1)/32 + 1
      k = (k*(k + 1))/2
      do j = 1, k
         msdata(ncb+j)%first = ncb
         msdata(ncb+j)%lval = c_ptr_plus(nodes(node)%gpu_lcol, &
            blkn*C_SIZEOF(dummy_real))
         if(present(gpu_ldcol)) then
            msdata(ncb+j)%ldval = c_ptr_plus(gpu_ldcol(i), &
               blkn*C_SIZEOF(dummy_real)) ! LD in indef case
         else
            msdata(ncb+j)%ldval = c_ptr_plus(nodes(node)%gpu_lcol, &
               blkn*C_SIZEOF(dummy_real)) ! L in posdef case
         endif
         msdata(ncb+j)%offc = off_LDLT(node)
         msdata(ncb+j)%n = m
         msdata(ncb+j)%k = nelim
         msdata(ncb+j)%lda = blkm
         msdata(ncb+j)%ldb = blkm
      end do
      ncb = ncb + k
   end do

   if ( ncb > 0 ) then
      gpu_msdata = custack_alloc(gwork, ncb*C_SIZEOF(msdata(1)))
      cuda_error = cudaMemcpyAsync_H2D(gpu_msdata, C_LOC(msdata), &
         ncb*C_SIZEOF(msdata(1)), stream)
      if(cuda_error.ne.0) return
      call cuda_multidsyrk_low_col( stream, ncb, gpu_msdata, ptr_levLDLT )
      call custack_free(gwork, ncb*C_SIZEOF(msdata(1))) ! gpu_msdata
   end if

end subroutine form_contrib

subroutine factor_posdef(stream, lev, lvlptr, nodes, lvllist, sptr, rptr, &
      ptr_levL, cublas_handle, stats, gwork)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: lev
   integer, dimension(*), intent(in) :: lvlptr
   type(node_type), dimension(*), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: lvllist
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(C_PTR), intent(in) :: ptr_levL
   type(C_PTR), intent(in) :: cublas_handle
   type(thread_stats), intent(inout) :: stats
   type(cuda_stack_alloc_type), intent(inout) :: gwork

   type(multisymm_type), dimension(:), allocatable, target :: msdata
   type(C_PTR) :: gpu_msdata

   type(C_PTR) :: ptr_B
   integer :: i, j
   integer :: ncb, llist, rmax, rtot
   integer :: blkm, blkn, node

   type(C_PTR) :: ptr_L

   integer, dimension(:), allocatable :: node_m, node_n
   type(C_PTR), dimension(:), allocatable :: node_lcol

   integer, allocatable, target :: nelm(:) ! work arrays
   real(wp) :: dummy_real

   ncb = lvlptr(lev + 1) - lvlptr(lev)
   allocate(nelm(ncb), msdata(ncb), node_m(ncb), node_n(ncb), node_lcol(ncb), &
      stat=stats%st)
   if(stats%st.ne.0) return

   !
   ! Copy lower triangle into upper triangle so we can use access (i,j) or
   ! (j,i) to get the same number while pivoting.
   !
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)        
      i = llist - lvlptr(lev) + 1
      msdata(i)%lcol = nodes(node)%gpu_lcol
      msdata(i)%ncols = sptr(node + 1) - sptr(node)
      msdata(i)%nrows = int(rptr(node + 1) - rptr(node))
   end do
   gpu_msdata = custack_alloc(gwork, ncb*C_SIZEOF(msdata(1)))
   stats%cuda_error = &
      cudaMemcpyAsync_h2d(gpu_msdata, C_LOC(msdata), ncb*C_SIZEOF(msdata(1)), stream)
   if(stats%cuda_error.ne.0) return
   call multisymm(stream, ncb, gpu_msdata)
   call custack_free(gwork, ncb*C_SIZEOF(msdata(1))) ! gpu_msdata

   !
   ! Factor several nodes simultaneously
   !
   ! Setup
   ncb = 0
   rmax = 0
   rtot = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)        
      i = llist - lvlptr(lev) + 1
      blkn = sptr(node + 1) - sptr(node)
      blkm = int(rptr(node + 1) - rptr(node))
      if ( blkn < blkm ) then
         ncb = ncb + 1
         node_m(ncb) = blkm
         node_n(ncb) = blkn
         node_lcol(ncb) = nodes(node)%gpu_lcol
         rtot = rtot + blkm
      else if ( blkn == blkm ) then
         rmax = max(rmax, blkm)
      end if
   end do
   rmax = max(rmax, rtot)
   ptr_B = custack_alloc(gwork, rmax*BLOCK_SIZE*C_SIZEOF(dummy_real))
   if ( ncb > 0 ) then
      ! Perform simultaneous LLT factorization
      call multinode_llt(stream, ncb, node_m, node_n, node_lcol,              &
         cublas_handle, ptr_levL, ptr_B, BLOCK_SIZE, nelm, stats%flag, gwork, &
         stats%st, stats%cuda_error, stats%cublas_error)
      if(stats%st.ne.0) return
      if(stats%cuda_error.ne.0.or.stats%cublas_error.ne.0) return
      if(stats%flag.lt.0) return
      ! Store outcome of factorization
      ncb = 0
      do llist = lvlptr(lev), lvlptr(lev + 1) - 1
         node = lvllist(llist)
         blkn = sptr(node + 1) - sptr(node)
         blkm = int(rptr(node + 1) - rptr(node))
         if ( blkn < blkm ) then
            ncb = ncb + 1
            nodes(node)%nelim = nelm(ncb)
            do j = blkm, blkm-nelm(ncb)+1, -1
               stats%num_factor = stats%num_factor + j
               stats%num_flops = stats%num_flops + j**2_long
            end do
         end if
      end do
   end if
   
   !
   ! Factor root nodes
   !
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)
      ptr_L = nodes(node)%gpu_lcol
      node = lvllist(llist)
      blkn = sptr(node + 1) - sptr(node)
      blkm = int(rptr(node + 1) - rptr(node))
      ! Positive-definite, LL^T factorization, no pivoting
      if ( blkm == blkn ) then
         call node_llt(stream, blkm, blkn, ptr_L, blkm, ptr_B, BLOCK_SIZE, &
            cublas_handle, stats%flag, gwork, stats%cuda_error,            &
            stats%cublas_error)
         if(stats%cuda_error.ne.0.or.stats%cublas_error.ne.0) return
         if(stats%flag.lt.0) return
         nodes(node)%nelim = blkn
      end if
   end do

   call custack_free(gwork, rmax*BLOCK_SIZE*C_SIZEOF(dummy_real)) ! ptr_B

end subroutine factor_posdef

subroutine collect_stats_indef(stream, lev, lvlptr, lvllist, nodes, &
      sptr, rptr, stats, gwork, gpu_custats)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: lev
   integer, dimension(*), intent(in) :: lvlptr
   integer, dimension(*), intent(in) :: lvllist
   type(node_type), dimension(*), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(thread_stats), intent(inout) :: stats
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   type(C_PTR), intent(in) :: gpu_custats

   type(cstat_data_type), dimension(:), allocatable, target :: csdata
   type(C_PTR) :: gpu_csdata

   integer :: llist, llvlptr
   integer :: node, blkm, blkn, ndelay
   real(wp) :: dummy_real

   llvlptr = lvlptr(lev+1)-lvlptr(lev)
   allocate(csdata(llvlptr), stat=stats%st)
   if(stats%st.ne.0) return
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)
      ndelay = nodes(node)%ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      csdata(llist-lvlptr(lev)+1)%nelim = nodes(node)%nelim
      csdata(llist-lvlptr(lev)+1)%dval = c_ptr_plus(nodes(node)%gpu_lcol, &
         blkm*(blkn+0_long)*C_SIZEOF(dummy_real))
   end do

   gpu_csdata = custack_alloc(gwork, llvlptr*C_SIZEOF(csdata(1)))
   stats%cuda_error = cudaMemcpyAsync_H2D(gpu_csdata, C_LOC(csdata), &
      llvlptr*C_SIZEOF(csdata(1)), stream)
   if(stats%cuda_error.ne.0) return

   call cuda_collect_stats(stream, size(csdata), gpu_csdata, gpu_custats)

   call custack_free(gwork, llvlptr*C_SIZEOF(csdata(1))) ! gpu_csdata
end subroutine collect_stats_indef

! Factorize a nodal matrix (not contrib block)
subroutine factor_indef( stream, lev, lvlptr, nodes, lvllist, sparent, sptr, &
      rptr, level_height, level_width, delta, eps, gpu_ldcol, gwork, &
      cublas_handle, options, stats, gpu_custats )
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: lev
   integer, dimension(*), intent(in) :: lvlptr
   type(node_type), dimension(*), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: lvllist
   integer, dimension(*), intent(in) :: sparent
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer, intent(in) :: level_height
   integer, intent(in) :: level_width
   real(wp), intent(inout) :: delta, eps
   type(C_PTR), dimension(:), intent(in) :: gpu_ldcol
   type(C_PTR), intent(in) :: cublas_handle
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   type(ssids_options), intent(in) :: options
   type(thread_stats), intent(inout) :: stats
   type(C_PTR), intent(in) :: gpu_custats

   type(multisymm_type), dimension(:), allocatable, target :: msdata
   type(C_PTR) :: gpu_msdata

   type(multiswap_type), dimension(:), allocatable, target :: swapdata
   type(C_PTR) :: gpu_swapdata

   type(C_PTR) :: ptr_ind, ptr_B
   integer :: ind_len, B_len
   integer :: i, j, k, p
   integer :: ncb, last_ln, llist, maxr
   integer :: ndelay, blkm, blkn, nelim, parent, node
   integer, dimension(:), pointer :: lperm

   type(C_PTR) :: ptr_D, ptr_L, ptr_LD

   integer, dimension(:), allocatable :: node_m, node_n, node_skip
   type(C_PTR), dimension(:), allocatable :: node_lcol, node_ldcol

   integer, allocatable, target :: iwork(:), perm(:), nelm(:) ! work arrays
   real(wp) :: dummy_real
   integer(C_INT) :: dummy_int

   ncb = lvlptr(lev + 1) - lvlptr(lev)
   allocate(perm(level_width), nelm(ncb), msdata(ncb), &
      node_m(ncb), node_n(ncb), node_lcol(ncb), node_ldcol(ncb), &
      node_skip(ncb), swapdata(ncb), stat=stats%st )
   if(stats%st.ne.0) return

   ! Initialize variables to avoid warnings
   last_ln = 0

   !
   ! Copy lower triangle into upper triangle so we can use access (i,j) or
   ! (j,i) to get the same number while pivoting.
   !
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)        
      i = llist - lvlptr(lev) + 1
      ndelay = nodes(node)%ndelay
      msdata(i)%lcol = nodes(node)%gpu_lcol
      msdata(i)%ncols = sptr(node + 1) - sptr(node) + ndelay
      msdata(i)%nrows = int(rptr(node + 1) - rptr(node)) + ndelay
   end do
   gpu_msdata = custack_alloc(gwork, ncb*C_SIZEOF(msdata(1)))
   stats%cuda_error = &
      cudaMemcpyAsync_H2D(gpu_msdata, C_LOC(msdata), ncb*C_SIZEOF(msdata(1)), &
         stream)
   if(stats%cuda_error.ne.0) return
   call multisymm( stream, ncb, gpu_msdata )
   call custack_free(gwork, ncb*C_SIZEOF(msdata(1)))
   ! Note: Done with gpu_msdata

   !
   ! Swap delays to end
   !
   ncb = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)
      i = llist - lvlptr(lev) + 1
      ndelay = nodes(node)%ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      if(ndelay.gt.0 .and. blkn.gt.1) then
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         k = min(ndelay, blkn - ndelay)
         ncb = ncb + 1
         swapdata(ncb)%nrows = blkm
         swapdata(ncb)%ncols = blkn
         swapdata(ncb)%k = k
         swapdata(ncb)%lcol = nodes(node)%gpu_lcol
         swapdata(ncb)%lda = blkm
         swapdata(ncb)%off = blkn - k
         lperm => nodes(node)%perm
         do i = 1, k
            j = blkn - k + i
            p = lperm(i)
            lperm(i) = lperm(j)
            lperm(j) = p
         end do
      end if
   end do
   if(ncb.gt.0) then
      gpu_swapdata = custack_alloc(gwork, ncb*C_SIZEOF(swapdata(1)))
      stats%cuda_error = cudaMemcpyAsync_H2D(gpu_swapdata, C_LOC(swapdata), &
         ncb*C_SIZEOF(swapdata(1)), stream)
      if(stats%cuda_error.ne.0) return
      call swap_ni2Dm( stream, ncb, gpu_swapdata )
      call custack_free(gwork, ncb*C_SIZEOF(swapdata(1)))
   end if

   !
   ! Factor several nodes simultaneously
   !
   ! Setup
   p = 0
   k = 0
   ncb = 0
   maxr = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)        
      i = llist - lvlptr(lev) + 1
      ndelay = nodes(node)%ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      if(blkn.lt.blkm) then
         ncb = ncb + 1
         lperm => nodes(node)%perm
         node_m(ncb) = blkm
         node_n(ncb) = blkn
         node_lcol(ncb) = nodes(node)%gpu_lcol
         node_ldcol(ncb) = gpu_ldcol(i)
         node_skip(ncb) = k
         do j = 1, blkn
            perm(k + j) = lperm(j)
         end do
         k = k + blkn
         maxr = max(maxr, blkm)
         p = p + blkm
      end if
   end do
  
   ind_len = max(level_width, ncb*BLOCK_SIZE)
   B_len = max(p, level_height)*BLOCK_SIZE

   ptr_ind = custack_alloc(gwork, ind_len*C_SIZEOF(dummy_int))
   ptr_B = custack_alloc(gwork, B_len*C_SIZEOF(dummy_real))

   last_ln = ncb

   if(ncb.gt.1) then

      ! Perform simultaneous factorization of several nodes
      call multinode_ldlt(stream, ncb, node_m, node_n, node_lcol, node_ldcol,  &
         node_skip, ptr_B, ptr_ind, delta, eps, BLOCK_SIZE, perm, nelm, gwork, &
         cublas_handle, stats%st, stats%cuda_error, stats%cublas_error)
      if(stats%st.ne.0 .or. stats%cuda_error.ne.0 .or. &
         stats%cublas_error.ne.0) return

      ! Store outcome of factorization
      ncb = 0
      k = 0
      do llist = lvlptr(lev), lvlptr(lev + 1) - 1
         node = lvllist(llist)
         ndelay = nodes(node)%ndelay
         blkn = sptr(node + 1) - sptr(node) + ndelay
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         parent = sparent(node)
         lperm => nodes(node)%perm
         if(blkn.lt.blkm) then
            ncb = ncb + 1
            nelim = nelm(ncb)
            do j = 1, blkn
               lperm(j) = perm(k + j)
            end do
            k = k + blkn
            nodes(node)%nelim = nelim
            nodes(parent)%ndelay = nodes(parent)%ndelay + blkn - nelim
            stats%num_delay &
               = stats%num_delay + blkn - nelim
            do j = blkm, blkm-nelim+1, -1
               stats%num_factor = stats%num_factor + j
               stats%num_flops = stats%num_flops + j**2_long
            end do
         end if
      end do
   end if
   
   allocate(iwork(2*level_width), stat=stats%st)
   if(stats%st.ne.0) return
   
   !
   ! Factor remaining nodes one by one
   !
   ncb = 0
   k = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      i = llist - lvlptr(lev) + 1
      node = lvllist(llist)
      ptr_L = nodes(node)%gpu_lcol
      ptr_LD = gpu_ldcol(i)

      node = lvllist(llist)
      ndelay = nodes(node)%ndelay
      blkn = sptr(node + 1) - sptr(node) + ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      parent = sparent(node)
      lperm => nodes(node)%perm
   
      ptr_D = c_ptr_plus( ptr_L, blkm*blkn*C_SIZEOF(dummy_real) )

      ! Indefinite, LDL^T factorization, with pivoting
      if(blkm.eq.blkn .or. last_ln.eq.1 ) then

         call node_ldlt(stream, blkm, blkn, ptr_L, ptr_LD, blkm, ptr_D, ptr_B, &
            ptr_ind, delta, eps, BLOCK_SIZE, lperm, iwork, nelim, gwork, &
            cublas_handle, stats%cuda_error, stats%cublas_error)
         if(stats%cuda_error.ne.0 .or. stats%cublas_error.ne.0) return

         if(blkm.eq.blkn .and. nelim.lt.blkn) then
            if(options%action) then
               stats%flag = SSIDS_WARNING_FACT_SINGULAR
            else
               stats%flag = SSIDS_ERROR_SINGULAR
               return
            endif
         end if

         ! Record delays
         nodes(node)%nelim = nelim
         if(blkn.lt.blkm) then
!$OMP ATOMIC
            nodes(parent)%ndelay = nodes(parent)%ndelay + (blkn - nelim)
         endif
         stats%num_delay = stats%num_delay + blkn - nelim
         do j = blkm, blkm-nelim+1, -1
            stats%num_factor = stats%num_factor + j
            stats%num_flops = stats%num_flops + j**2_long
         end do

      end if

   end do

   call custack_free(gwork, B_len*C_SIZEOF(dummy_real)) ! ptr_B
   call custack_free(gwork, ind_len*C_SIZEOF(dummy_int)) ! ptr_ind

   call collect_stats_indef(stream, lev, lvlptr, lvllist, nodes, &
      sptr, rptr, stats, gwork, gpu_custats)

end subroutine factor_indef

subroutine setup_assemble_contrib(stream, lev, lvlptr, lvllist, child_ptr, &
      child_list, sptr, rptr, asminf, gpu_ccval, gpu_contribs, ptr_levLDLT, &
      off_LDLT, gpu_rlist_direct, gwork, ncp, gpu_cpdata, nblk, &
      gpu_blkdata, gpu_sync, st, cuda_error)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: lev
   integer, intent(in) :: lvlptr(*) ! Pointers into lvllist for level
   integer, intent(in) :: lvllist(*) ! Nodes at level lev are given by:
      ! lvllist(lvlptr(lev):lvlptr(lev+1)-1)
   integer, intent(in) :: child_ptr(*) ! Pointers into child_list for node
   integer, intent(in) :: child_list(*) ! Children of node node are given by:
      ! child_list(child_ptr(node):child_ptr(node+1)-1)
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(asmtype), dimension(:), intent(in) :: asminf ! Assembly info
   type(C_PTR), intent(in) :: gpu_ccval ! GPU (*gpu_ccval) points to previous
      ! level's contribution blocks
   type(C_PTR), dimension(*), intent(in) :: gpu_contribs ! For each node, is
      ! either NULL or points to a contribution block arriving from a subtree
   type(C_PTR), intent(in) :: ptr_levLDLT
   integer(long), intent(in) :: off_LDLT(*) ! Offsets for children
   type(C_PTR), intent(in) :: gpu_rlist_direct ! GPU pointer to rlist_direct
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   integer, intent(out) :: ncp ! Number of child-parent pairs
   type(C_PTR), intent(out) :: gpu_cpdata ! Ouput child-parent info
   integer, intent(out) :: nblk ! Number of blocks
   type(C_PTR), intent(out) :: gpu_blkdata ! Output block-by-block info
   type(C_PTR), intent(out) :: gpu_sync
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   type(assemble_cp_type), dimension(:), allocatable, target :: cpdata
   type(assemble_blk_type), dimension(:), allocatable, target :: blkdata

   integer :: child, maxchild
   integer :: cpi, bi, blki, blkj
   integer :: j, k, m, npassl, blkm, blkn
   integer :: llist, node, cnode, npassed, bx, by
   integer :: blk
   real(wp) :: dummy_real
   integer(C_INT) :: dummy_int

   ! Ensure all return values initialized (mainly to prevent warnings)
   ncp = 0
   gpu_cpdata = C_NULL_PTR
   nblk = 0
   gpu_blkdata = C_NULL_PTR
   gpu_sync = C_NULL_PTR
   st = 0
   cuda_error = 0

   ! Count number of children with work to do
   ncp = 0
   nblk = 0
   maxchild = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      node = lvllist(llist)
      blkm = int(rptr(node + 1) - rptr(node))
      blkn = sptr(node + 1) - sptr(node)
      m = blkm - blkn
      if(m.gt.0) then
         do child = child_ptr(node), child_ptr(node+1)-1
            cnode = child_list(child)
            npassed = asminf(cnode)%npassed
            npassl = asminf(cnode)%npassl
            if(npassed-npassl.gt.0) then
               ncp = ncp + 1
               bx = (npassed-npassl-1) / HOGG_ASSEMBLE_TX + 1
               by = (npassed-npassl-1) / HOGG_ASSEMBLE_TY + 1
               nblk = nblk + &
                  calc_blks_lwr(bx, by, HOGG_ASSEMBLE_TX, HOGG_ASSEMBLE_TY)
            endif
         end do
      end if
      maxchild = max(maxchild, child_ptr(node+1)-child_ptr(node))
   end do

   ! Fill in child-parent information
   allocate(cpdata(ncp), stat=st)
   if(st.ne.0) return
   cpi = 1
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
      k = llist - lvlptr(lev) + 1
      node = lvllist(llist)
      blkm = int(rptr(node + 1) - rptr(node))
      blkn = sptr(node + 1) - sptr(node)
      blk = 0
      do child = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(child)
         npassed = asminf(cnode)%npassed
         npassl = asminf(cnode)%npassl
         ! We are only doing the contribution block of the parent
         if(npassed-npassl.le.0) cycle
         cpdata(cpi)%cm = npassed - npassl
         cpdata(cpi)%cn = npassed - npassl
         cpdata(cpi)%ldp = blkm - blkn
         ! Note: row rlist(i) of parent is row rlist(i)-blkn of contribution blk
         ! so we alter pval to account for this
         cpdata(cpi)%pval = c_ptr_plus(ptr_levLDLT, &
            (off_LDLT(node) - blkn*(1+cpdata(cpi)%ldp))*C_SIZEOF(dummy_real) )
         cpdata(cpi)%ldc = asminf(cnode)%npassed
            cpdata(cpi)%cvoffset = off_LDLT(cnode) + &
               npassl * (1+cpdata(cpi)%ldc)
         if(C_ASSOCIATED(gpu_contribs(cnode))) then
            cpdata(cpi)%cv = c_ptr_plus(gpu_contribs(cnode), &
                  (npassl * (1+cpdata(cpi)%ldc)) * C_SIZEOF(dummy_real))
         else
            cpdata(cpi)%cv = c_ptr_plus(gpu_ccval, &
                  (off_LDLT(cnode) + npassl * (1+cpdata(cpi)%ldc)) * &
                  C_SIZEOF(dummy_real))
         endif
         cpdata(cpi)%rlist_direct = c_ptr_plus(gpu_rlist_direct, &
               (asminf(cnode)%offset + npassl)*C_SIZEOF(dummy_int) &
            )
         cpdata(cpi)%sync_offset = cpi-1 - 1 ! 0-indexed
         cpdata(cpi)%sync_wait_for = blk
         ! Calulate how many blocks next iteration needs to wait for
         bx = (cpdata(cpi)%cm-1) / HOGG_ASSEMBLE_TX + 1
         by = (cpdata(cpi)%cn-1) / HOGG_ASSEMBLE_TY + 1
         blk = calc_blks_lwr(bx, by, HOGG_ASSEMBLE_TX, HOGG_ASSEMBLE_TY)
         cpi = cpi + 1
      end do
   end do

   ! Setup block information: do all first children, then all second, etc.
   ! This ensures maximum parallelism can be exploited
   allocate(blkdata(nblk), stat=st)
   if(st.ne.0) return
   bi = 1
   do child = 1, maxchild
      cpi = 1
      do llist = lvlptr(lev), lvlptr(lev + 1) - 1
         node = lvllist(llist)
         do j = child_ptr(node), child_ptr(node+1)-1
            cnode = child_list(j)
            ! We are only doing the contribution block of the parent
            if(asminf(cnode)%npassed-asminf(cnode)%npassl .le. 0) cycle
            if(j-child_ptr(node)+1 .eq. child) then
               bx = (cpdata(cpi)%cm-1) / HOGG_ASSEMBLE_TX + 1
               by = (cpdata(cpi)%cn-1) / HOGG_ASSEMBLE_TY + 1
               do blkj = 0, by-1
                  do blki = 0, bx-1
                     if((blki+1)*HOGG_ASSEMBLE_TX<(blkj+1)*HOGG_ASSEMBLE_TY) &
                        cycle ! Entirely in upper triangle
                     blkdata(bi)%cp = cpi-1
                     blkdata(bi)%blk = blkj*bx + blki
                     bi = bi + 1
                  end do
               end do
            endif
            cpi = cpi + 1
         end do
      end do
   end do

   gpu_cpdata = custack_alloc(gwork, ncp*C_SIZEOF(cpdata(1)))
   gpu_blkdata = custack_alloc(gwork, nblk*C_SIZEOF(blkdata(1)))
   gpu_sync = custack_alloc(gwork, (1+ncp)*C_SIZEOF(dummy_int))

   cuda_error = cudaMemcpyAsync_H2D(gpu_cpdata, C_LOC(cpdata), &
      ncp*C_SIZEOF(cpdata(1)), stream)
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpyAsync_H2D(gpu_blkdata, C_LOC(blkdata), &
      nblk*C_SIZEOF(blkdata(1)), stream)
   if(cuda_error.ne.0) return
end subroutine setup_assemble_contrib

! Return number of blocks not entirely in upper triangle
! For nx x ny block matrix with block size nbx x nby
integer function calc_blks_lwr(nx, ny, nbx, nby)
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   integer, intent(in) :: nbx
   integer, intent(in) :: nby

   integer :: i

   calc_blks_lwr = 0
   do i = 1, nx
      calc_blks_lwr = calc_blks_lwr + min(ny,(i*nbx)/nby)
   end do

end function calc_blks_lwr

subroutine assemble_contrib(stream, total_nch, lev, lvlptr, lvllist,    &
      child_ptr, child_list, sptr, rptr, asminf, pc_size, off_LDLT,     &
      gpu_ccval, gpu_contribs, ptr_levLDLT, gpu_rlist_direct, gwork,    &
      st, cuda_error)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: total_nch
   integer, intent(in) :: lev
   integer, intent(in) :: lvlptr(*) ! Pointers into lvllist for level
   integer, intent(in) :: lvllist(*) ! Nodes at level lev are given by:
      ! lvllist(lvlptr(lev):lvlptr(lev+1)-1)
   integer, intent(in) :: child_ptr(*) ! Pointers into child_list for node
   integer, intent(in) :: child_list(*) ! Children of node node are given by:
      ! child_list(child_ptr(node):child_ptr(node+1)-1)
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(asmtype), dimension(:), intent(in) :: asminf ! Assembly info
   integer(long), intent(in) :: pc_size
   integer(long), dimension(*), intent(in) :: off_LDLT ! Offsets for children
   type(C_PTR), intent(in) :: gpu_ccval ! GPU (*gpu_ccval) points to previous
      ! level's contribution blocks
   type(C_PTR), dimension(*), intent(in) :: gpu_contribs ! For each node, is
      ! either NULL or points to a contribution block arriving from a subtree
   type(C_PTR), intent(in) :: ptr_levLDLT ! GPU pointer to contribution blocks
   type(C_PTR), intent(in) :: gpu_rlist_direct! GPU pointer to rlist_direct copy
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   integer :: ncp, nblk
   type(C_PTR) :: gpu_cpdata, gpu_blkdata, gpu_sync
   type(assemble_cp_type) :: act_dummy
   type(assemble_blk_type) :: abt_dummy
   integer(C_INT) :: dummy_int

   if(total_nch.eq.0 .or. pc_size.eq.0 ) return ! Nothing to do

   call setup_assemble_contrib(stream, lev, lvlptr, lvllist, child_ptr, &
      child_list, sptr, rptr, asminf, gpu_ccval, gpu_contribs, ptr_levLDLT, &
      off_LDLT, gpu_rlist_direct, gwork, ncp, gpu_cpdata, nblk, &
      gpu_blkdata, gpu_sync, st, cuda_error)
   if(st.ne.0 .or. cuda_error.ne.0) return

   call assemble(stream, nblk, 0, gpu_blkdata, ncp, gpu_cpdata, &
      gpu_ccval, ptr_levLDLT, gpu_sync)

   ! Free in reverse alloc order
   call custack_free(gwork, (1+ncp)*C_SIZEOF(dummy_int)) ! gpu_sync
   call custack_free(gwork, nblk*C_SIZEOF(abt_dummy)) ! gpu_blkdata
   call custack_free(gwork, ncp*C_SIZEOF(act_dummy)) ! gpu_cpdata
end subroutine assemble_contrib

subroutine setup_assemble_fully_summed(stream, total_nch, lev, lvlptr, lvllist,&
      child_ptr, child_list, nodes, asminf, sptr, rptr, gpu_ccval,             &
      gpu_contribs, off_LDLT, gpu_rlist_direct, gwork, ncp, gpu_cpdata, nblk,  &
      gpu_blkdata, ndblk, gpu_ddata, gpu_sync, st, cuda_error)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: total_nch ! Total number of children for nodes in
      ! this level
   integer, intent(in) :: lev ! Current level
   integer, intent(in) :: lvlptr(*) ! Pointers into lvllist for level
   integer, intent(in) :: lvllist(*) ! Nodes at level lev are given by:
      ! lvllist(lvlptr(lev):lvlptr(lev+1)-1)
   integer, intent(in) :: child_ptr(*) ! Pointers into child_list for node
   integer, intent(in) :: child_list(*) ! Children of node node are given by:
      ! child_list(child_ptr(node):child_ptr(node+1)-1)
   type(node_type), dimension(*), intent(in) :: nodes ! node data
   type(asmtype), dimension(*), intent(in) :: asminf
   integer, intent(in) :: sptr(*)
   integer(long), intent(in) :: rptr(*)
   type(C_PTR), intent(in) :: gpu_ccval ! GPU (*ptr_ccval) points to previous
      ! level's contribution blocks
   type(C_PTR), dimension(*), intent(in) :: gpu_contribs ! For each node, is
      ! either NULL or points to a contribution block arriving from a subtree
   integer(long), intent(in) :: off_LDLT(*) ! Offsets for children
   type(C_PTR), intent(in) :: gpu_rlist_direct ! GPU pointer to rlist_direct
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   integer, intent(out) :: ncp ! Number of child-parent pairs
   type(C_PTR), intent(out) :: gpu_cpdata ! Ouput child-parent info
   integer, intent(out) :: nblk ! Number of blocks
   type(C_PTR), intent(out) :: gpu_blkdata ! Output block-by-block info
   integer, intent(out) :: ndblk
   type(C_PTR), intent(out) :: gpu_ddata
   type(C_PTR), intent(out) :: gpu_sync
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   integer :: i, bi, ni, blk, bx, by, child, ci, cnode, cpi, node, blki, blkj
   integer :: ndelay, ldp, maxchild, blkm
   integer :: cndelay, cnelim, cblkm, cblkn, llist, nd
   type(C_PTR) :: pval
   type(assemble_cp_type), dimension(:), allocatable, target :: cpdata
      ! child-parent data to be copied to GPU
   type(assemble_blk_type), dimension(:), allocatable, target :: blkdata
      ! block-level data to be copied to GPU
   type(assemble_delay_type), dimension(:), allocatable, target :: ddata
      ! delay data to be copied to GPU
   integer(C_INT) :: dummy_int
   real(wp) :: dummy_real
   
   cuda_error = 0
   st = 0

   ! Ensure all output data is initialized
   ncp = 0
   gpu_cpdata = C_NULL_PTR
   nblk = 0
   gpu_blkdata = C_NULL_PTR
   ndblk = 0
   gpu_ddata = C_NULL_PTR
   gpu_sync = C_NULL_PTR
   cuda_error = 0

   ! Check for trivial return
   if(total_nch.le.0) return ! No children to worry about

   ! Count maximum number of children
   maxchild = 0
   do ni = lvlptr(lev), lvlptr(lev+1)-1
      node = lvllist(ni)
      maxchild = max(maxchild, child_ptr(node+1)-child_ptr(node))
   end do

   ! Initialize child-parent data, count number of blocks at each level
   allocate(cpdata(total_nch), stat=st)
   if(st.ne.0) return
   cpi = 1
   nblk = 0
   do ni = lvlptr(lev), lvlptr(lev+1)-1
      node = lvllist(ni)
      i = ni - lvlptr(lev) + 1
      ndelay = nodes(node)%ndelay
      ldp = int(rptr(node+1)-rptr(node)) + ndelay ! adjusted by #delays
      pval = c_ptr_plus(nodes(node)%gpu_lcol, &
         ndelay*(1+ldp)*C_SIZEOF(dummy_real)) ! adjusted past delays
      blk = 0
      do ci = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(ci)
         cpdata(cpi)%pval = pval
         cpdata(cpi)%ldp = ldp
         cpdata(cpi)%cm = asminf(cnode)%npassed
         cpdata(cpi)%cn = asminf(cnode)%npassl
         cpdata(cpi)%ldc = asminf(cnode)%npassed
         cpdata(cpi)%cvoffset = off_LDLT(cnode)
         if(C_ASSOCIATED(gpu_contribs(cnode))) then
            cpdata(cpi)%cv = gpu_contribs(cnode)
         else
            cpdata(cpi)%cv = &
               c_ptr_plus(gpu_ccval, off_LDLT(cnode)*C_SIZEOF(dummy_real))
         endif
         cpdata(cpi)%rlist_direct = c_ptr_plus(gpu_rlist_direct, &
            asminf(cnode)%offset*C_SIZEOF(dummy_int))
         cpdata(cpi)%sync_offset = max(0, cpi-1 - 1)
         cpdata(cpi)%sync_wait_for = blk
         bx = (cpdata(cpi)%cm-1) / HOGG_ASSEMBLE_TX + 1
         by = (cpdata(cpi)%cn-1) / HOGG_ASSEMBLE_TY + 1
         i = ci-child_ptr(node)+1
         blk = calc_blks_lwr(bx, by, HOGG_ASSEMBLE_TX, HOGG_ASSEMBLE_TY)
         nblk = nblk + blk
         cpi = cpi + 1
      end do
   end do
   ncp = size(cpdata)

   ! Initialize blkdata
   allocate(blkdata(nblk), stat=st)
   if(st.ne.0) return
   bi = 1
   do child = 1, maxchild
      cpi = 1
      do ni = lvlptr(lev), lvlptr(lev+1)-1
         node = lvllist(ni)
         if(child_ptr(node)+child-1 .lt. child_ptr(node+1)) then
            bx = (cpdata(cpi+child-1)%cm-1) / HOGG_ASSEMBLE_TX + 1
            by = (cpdata(cpi+child-1)%cn-1) / HOGG_ASSEMBLE_TY + 1
            do blkj = 0, by-1
               do blki = 0, bx-1
                  if((blki+1)*HOGG_ASSEMBLE_TX<(blkj+1)*HOGG_ASSEMBLE_TY) &
                     cycle ! Entirely in upper triangle
                  blkdata(bi)%cp = cpi + (child-1) - 1 ! 0 indexed
                  blkdata(bi)%blk = blkj*bx + blki
                  bi = bi + 1
               end do
            end do
         endif
         cpi = cpi + child_ptr(node+1)-child_ptr(node)
      end do
   end do

   ! Initialize ddata (for copying in any delays)
   allocate(ddata(total_nch), stat=st)
   if(st.ne.0) return
   ndblk = 0
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1
     node = lvllist(llist)
     ndelay = nodes(node)%ndelay
     blkm = int(rptr(node + 1) - rptr(node)) + ndelay
     i = llist - lvlptr(lev) + 1
     nd = 0
     do child = child_ptr(node), child_ptr(node+1)-1
       cnode = child_list(child)
       cblkm = int(rptr(cnode + 1) - rptr(cnode))
       cblkn = sptr(cnode + 1) - sptr(cnode)
       cndelay = nodes(cnode)%ndelay
       cnelim = nodes(cnode)%nelim
       if(cblkn+cndelay .le. cnelim) cycle ! No delays from this child
       ndblk = ndblk + 1
       ddata(ndblk)%ldd = blkm
       ddata(ndblk)%ndelay = ndelay - nd
       ddata(ndblk)%m = cblkm + cndelay - cnelim
       ddata(ndblk)%n = cblkn + cndelay - cnelim
       ddata(ndblk)%lds = cblkm + cndelay
       ddata(ndblk)%dval = c_ptr_plus(nodes(node)%gpu_lcol, &
         nd*(1_C_SIZE_T+blkm)*C_SIZEOF(dummy_real))
       ddata(ndblk)%sval = c_ptr_plus(nodes(cnode)%gpu_lcol, &
         cnelim*(1_C_SIZE_T+ddata(ndblk)%lds)*C_SIZEOF(dummy_real))
       ddata(ndblk)%roffset = rptr(cnode) + cblkn - 1
       nd = nd + ddata(ndblk)%n
     end do
   end do

   ! Copy data to GPU
   gpu_cpdata = custack_alloc(gwork, ncp*C_SIZEOF(cpdata(1)))
   gpu_blkdata = custack_alloc(gwork, nblk*C_SIZEOF(blkdata(1)))
   gpu_ddata = custack_alloc(gwork, ndblk*C_SIZEOF(ddata(1)))
   gpu_sync = custack_alloc(gwork, (ncp + 1)*C_SIZEOF(dummy_int))

   cuda_error = cudaMemcpyAsync_H2D(gpu_cpdata, C_LOC(cpdata), &
      ncp*C_SIZEOF(cpdata(1)), stream)
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpyAsync_H2D(gpu_blkdata, C_LOC(blkdata), &
      nblk*C_SIZEOF(blkdata(1)), stream)
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpyAsync_H2D(gpu_ddata, C_LOC(ddata), &
      ndblk*C_SIZEOF(ddata(1)), stream)
   if(cuda_error.ne.0) return

end subroutine setup_assemble_fully_summed

!
! Perform assembly for fully summed columns
!
! At a given level, we launch one kernel for each set of ith children:
!    1) Kernel that does assembly for all 1st children
!    2) Kernel that does assembly for all 2nd children
!    3) ...
! Information about a particular child-parent assembly is stored in cpdata
! For each block that gets launched we store an offset into cpdata and a
!   the subblock of that assembly this block is to perform.
subroutine assemble_fully_summed(stream, total_nch, lev, lvlptr, lvllist, &
      nodes, gpu_ccval, gpu_contribs, ptr_levL, gpu_rlist_direct, &
      child_ptr, child_list, off_LDLT, asminf, rptr, sptr, gwork, &
      st, cuda_error)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: total_nch ! Total number of children for nodes in
      ! this level
   integer, intent(in) :: lev ! Current level
   integer, intent(in) :: lvlptr(*) ! Pointers into lvllist for level
   integer, intent(in) :: lvllist(*) ! Nodes at level lev are given by:
      ! lvllist(lvlptr(lev):lvlptr(lev+1)-1)
   integer, intent(in) :: child_ptr(*) ! Pointers into child_list for node
   integer, intent(in) :: child_list(*) ! Children of node node are given by:
      ! child_list(child_ptr(node):child_ptr(node+1)-1)
   type(asmtype), dimension(*), intent(in) :: asminf ! Assembly info
   type(C_PTR), intent(in) :: gpu_ccval ! GPU (*gpu_ccval) points to previous
      ! level's contribution blocks
   type(C_PTR), dimension(*), intent(in) :: gpu_contribs ! For each node, is
      ! either NULL or points to a contribution block arriving from a subtree
   type(C_PTR), intent(in) :: ptr_levL ! GPU (*ptr_levL) points to L storage
      ! for current level
   type(C_PTR), intent(in) :: gpu_rlist_direct! GPU pointer to rlist_direct copy
   type(node_type), intent(in) :: nodes(*)
   integer(long), intent(in) :: off_LDLT(*)
   integer(long), intent(in) :: rptr(*)
   integer, intent(in) :: sptr(*)
   type(cuda_stack_alloc_type), intent(inout) :: gwork
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   integer, dimension(:), pointer :: lperm
   integer :: cnd
   integer :: k
   integer :: cblkn
   integer :: cn
   integer :: cnode
   integer :: i, j
   integer :: llist
   integer :: nd
   integer :: node

   integer :: ncp, nblk, ndblk
   type(C_PTR) :: gpu_cpdata, gpu_blkdata, gpu_ddata, gpu_sync
   type(assemble_cp_type) :: act_dummy
   type(assemble_blk_type) :: abt_dummy
   type(assemble_delay_type) :: adt_dummy
   integer(C_INT) :: dummy_int

   if ( total_nch .le. 0 ) return

   ! Setup data structures (allocates gpu_cpdata, gpu_blkdata)
   call setup_assemble_fully_summed(stream, total_nch, lev, lvlptr, lvllist, &
      child_ptr, child_list, nodes, asminf, sptr, rptr, &
      gpu_ccval, gpu_contribs, off_LDLT, gpu_rlist_direct, gwork, &
      ncp, gpu_cpdata, nblk, gpu_blkdata, ndblk, gpu_ddata, gpu_sync, &
      st, cuda_error)
   if(st.ne.0 .or. cuda_error.ne.0) return

   ! Perform assembly of child contributions to fully summed columns
   call assemble(stream, nblk, 0, gpu_blkdata, ncp, gpu_cpdata, &
      gpu_ccval, ptr_levL, gpu_sync)

   ! Copy any delayed columns
   call add_delays(stream, ndblk, gpu_ddata, gpu_rlist_direct)

   ! Release memory (in reverse order of alloc)
   call custack_free(gwork, (ncp + 1)*C_SIZEOF(dummy_int)) ! gpu_sync
   call custack_free(gwork, ndblk*C_SIZEOF(adt_dummy)) ! gpu_ddata
   call custack_free(gwork, nblk*C_SIZEOF(abt_dummy)) ! gpu_blkdata
   call custack_free(gwork, ncp*C_SIZEOF(act_dummy)) ! gpu_cpdata

   ! Set lperm for delayed columns from children
   do llist = lvlptr(lev), lvlptr(lev + 1) - 1

     node = lvllist(llist)
     lperm => nodes(node)%perm
 
     nd = 0
     do cn = child_ptr(node), child_ptr(node + 1) - 1
        cnode = child_list(cn)
        cnd = nodes(cnode)%ndelay
        cblkn = sptr(cnode + 1) - sptr(cnode)
        k = nd
        do i = nodes(cnode)%nelim + 1, cblkn + cnd
           j = nodes(cnode)%perm(i)
           k = k + 1
           lperm(k) = j
        end do
        nd = nd + cblkn + cnd - nodes(cnode)%nelim
     end do
   end do
end subroutine assemble_fully_summed

integer function max_nelim(nnodes, nodes, nlev, tree_le, node_list) result ( m )
   integer, intent(in) :: nnodes
   type(node_type), intent(in) :: nodes(nnodes)
   integer, intent(in) :: nlev ! # of elimination tree levels
   integer, intent(in) :: tree_le(nlev + 1) ! tree levels' entry points into the
                                             ! list of nodes below
   integer, intent(in) :: node_list(tree_le(nlev + 1) - 1) ! list of nodes

   integer :: node
   integer :: nelim
   integer :: i

   m = 0
   do i = 1, tree_le(nlev + 1) - 1
      node = node_list(i)
      nelim = nodes(node)%nelim
      m = max(m, nelim)
   end do
end function max_nelim

!
! Pre-solve subroutines below do some post-processing with the L-factor
! in order to accelerate the subsequent solve(s) with it.
! Each node matrix
!
!  |L_d|
!  |L_o|
!
! where L_d is the square dense diagonal block and L_o is compactly stored
! sparse off-diagonal block, is replaced by
!
!  |    L_d**(-1)   |
!  | -L_o L_d**(-1) |
!
! so that the triangular solves are performed by gemm-like CUDA kernels
! rather than trsv-like kernels.
!
! In order to invert L_d's, we split them into square tiles of size <tile_size>.
! All L_d's that have the same number of tile columns are inverted
! simultaneously.
!

!
! This subroutine computes the number of 'tile-width levels' <ntlev>, 
! the groups of nodes that have the same number of tile columns,
! and the workspace size (in doubles) <lwork> needed for pre-solve.
!
subroutine presolve_lwork(nnodes, nodes, sptr, rptr, nlev, tree_le, &
      node_list, tile_size, maxtc, iwork, ntlev, lwork)
   integer, intent(in) :: nnodes ! this and next three same as above
   type(node_type), intent(in) :: nodes(nnodes)
   integer, intent(in) :: sptr(nnodes + 1)
   integer(long), intent(in) :: rptr(nnodes + 1)
   integer, intent(in) :: nlev ! # of elimination tree levels
   integer, intent(in) :: tree_le(nlev + 1) ! tree levels' entry points into the
      ! list of nodes below
   integer, intent(in) :: node_list(tree_le(nlev + 1) - 1) ! list of nodes
   integer, intent(in) :: tile_size
   integer, intent(in) :: maxtc ! max # of tile columns in a node
   integer, intent(out) :: iwork(maxtc,2) ! counters
   integer, intent(out) :: ntlev
   integer, intent(out) :: lwork
  
   integer :: lev
   integer :: node
   integer :: blkm, blkn
   integer :: nelim
   integer :: ndelay
   integer :: ntiles
   integer :: level_size
  
   integer :: i
  
   iwork(:,:) = 0 ! initialize counters

   ! iwork(i,1) counts the number of nodes that have i tile columns
   ntlev = 0
   do i = 1, tree_le(nlev + 1) - 1
      node = node_list(i)
      nelim = nodes(node)%nelim
      if ( nelim < 1 ) cycle
      ntiles = (nelim - 1)/tile_size + 1 ! # of tile columns for this node
      if ( iwork(ntiles,1) == 0 ) ntlev = ntlev + 1 ! count new tile-width level
      iwork(ntiles,1) = iwork(ntiles,1) + 1 ! increment the number of nodes in
         ! this level
   end do

   ! iwork(i,2) is the amount of memory in doubles for inverting all L_d's
   ! that have i tile columns
   do i = 1, tree_le(nlev + 1) - 1
      node = node_list(i)
      nelim = nodes(node)%nelim
      if ( nelim < 1 ) cycle
      blkm = int(rptr(node + 1) - rptr(node)) + nodes(node)%ndelay
      ntiles = (nelim - 1)/tile_size + 1

      ! the workspace for inverting L_d
      iwork(ntiles,2) = iwork(ntiles,2) + nelim**2

      ! one-node levels are treated differently and require additional workspace
      if ( iwork(ntiles,1) == 1 ) iwork(ntiles,2) = iwork(ntiles,2) + blkm*nelim
   end do
   lwork = maxval(iwork(:,2)) ! the workspace size for the first pre-solve step
  
   ! the workspace for the second step
   do lev = 1, nlev
      level_size = 0
      do i = tree_le(lev), tree_le(lev + 1) - 1
         node = node_list(i)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         blkn = sptr(node + 1) - sptr(node) + ndelay
         level_size = level_size + (blkm + 2)*blkn
      end do
      lwork = max(lwork, level_size)
   end do
end subroutine presolve_lwork

!
! This subroutine performs the first pre-solve step: the inversion of L_d
! for tile-width levels of more than one node + full presolve for one-node
! levels.
!
subroutine presolve_first(stream, nnodes, nodes, rptr, nlev, tree_le,      &
      node_list, tile_size, ntlev, tleve, maxtc, tlev, tlnode, pos_def,    &
      cublas, gpu_work, st, cuda_error, cublas_error)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: nnodes ! this and next 4 lines as in presolve_lwork
   type(node_type), intent(in) :: nodes(nnodes)
   integer(long), intent(in) :: rptr(nnodes + 1)
   integer, intent(in) :: nlev
   integer, intent(in) :: tree_le(nlev + 1)
   integer, intent(in) :: node_list(tree_le(nlev + 1) - 1) ! see presolve_lwork
   integer, intent(in) :: tile_size
   integer, intent(in) :: ntlev ! # of tile-width levels
   integer, intent(out) :: tleve(ntlev + 1) ! tile-width levels' entry points
   integer, intent(in) :: maxtc ! max # of tile columns in a node
   integer, intent(out) :: tlev(maxtc + 1) ! tlev(i) is the tile-width level of
      ! a node that is i tile columns wide
   integer, intent(out) :: tlnode(tree_le(nlev + 1) - 1) ! list of nodes
      ! arranged in tile-width levels
   logical, intent(in) :: pos_def ! true if the matrix is positive definite
   type(C_PTR), intent(in) :: cublas ! CUBLAS handle
   type(C_PTR), intent(in) :: gpu_work ! workspace
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error
   integer, intent(out) :: cublas_error
  
   integer :: ntln
   integer :: node
   integer :: blkm
   integer :: nelim
   integer :: ndelay
   integer :: ntiles
  
   integer :: i, j, k

   integer(C_SIZE_T) :: off

   ! each node's data for presolve, nodes being arranged in tile-width levels
   type(node_data), allocatable, target :: ndata(:)
   type(C_PTR) :: ptr_ndata ! C pointer to a tile-width level in ndata

   type(C_PTR) :: gpu_u, gpu_w
   real(wp) :: dummy_real

   st = 0; cuda_error = 0; cublas_error = 0

   allocate(ndata(tree_le(nlev + 1) - 1), stat=st)
   if(st.ne.0) return

   ! for now, tlev(i) is the number of nodes that are (i - 1) tile columns wide
   tlev = 0
   do k = 1, tree_le(nlev + 1) - 1
      node = node_list(k)
      nelim = nodes(node)%nelim
      if(nelim.lt.1) cycle
      ntiles = (nelim - 1)/tile_size + 1
      tlev(ntiles + 1) = tlev(ntiles + 1) + 1
   end do

   ! move nonzeros in tlev to tleve so that tleve(j) is the number of
   ! nodes on (j - 1)-th tile-width level
   j = 1
   do i = 2, maxtc + 1
      if(tlev(i).ne.0) then
         j = j + 1
         tleve(j) = tlev(i)
         ! a node that is (i - 1) tile columns wide is on level (j - 1)
         tlev(i - 1) = j - 1
      end if
   end do
  
   ! compute tile-width levels' offsets
   tleve(1) = 0  
   do i = 3, ntlev + 1
      tleve(i) = tleve(i) + tleve(i - 1)
   end do

   ! populate ndata and tlnode
   do k = 1, tree_le(nlev + 1) - 1
      node = node_list(k)
      nelim = nodes(node)%nelim
      if(nelim.lt.1) cycle
      blkm = int(rptr(node + 1) - rptr(node)) + nodes(node)%ndelay
      ! find this node's level
      i = tlev((nelim - 1)/TILE_SIZE + 1)
      ! find this node's place in the list
      j = tleve(i) + 1
      tlnode(j) = node
      ! fill in pre-solve data
      ndata(j)%nrows = blkm
      ndata(j)%ncols = nelim
      ndata(j)%ld = blkm
      ndata(j)%ptr_v = nodes(node)%gpu_lcol
      tleve(i) = j ! next node's offset on this level
   end do

   ! compute entry points into tile-width levels
   do i = ntlev + 1, 2, -1
      tleve(i) = tleve(i - 1) + 1
   end do
   tleve(1) = 1

   ! perform the first pre-solve step level-by-level;
   ! note that on each level we are dealing with nodes of the same tile width
   do i = 1, ntlev
      ntln = tleve(i + 1) - tleve(i) ! # of nodes on this level

      ! find this level's offset in ndata
      off = (tleve(i) - 1)*C_SIZEOF(ndata(1))
      ptr_ndata = c_ptr_plus(C_LOC(ndata), off)

      ! we start by inverting all diagonal tiles in each node's L_d
      ! and applying the inverses to tiles below them in L_d
      if(pos_def) then
         ! non-unit diagonal
         cuda_error = multi_Ld_inv_init(stream, ntln, ptr_ndata, tile_size, &
            1, gpu_work)
      else
         ! unit diagonal
         cuda_error = multi_Ld_inv_init(stream, ntln, ptr_ndata, tile_size, &
            0, gpu_work)
      end if
      if(cuda_error.ne.0) return

      if(ntln.gt.1) then
         ! simultaneously invert all L_d's on this level
         cuda_error = multi_Ld_inv(stream, ntln, ptr_ndata, tile_size, gpu_work)
      else
         ! if there is only one node on this level, then do full pre-solve
         j = tleve(i)
         node = tlnode(j)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         gpu_u = nodes(node)%gpu_lcol
         off = (ndata(j)%ncols**2)*C_SIZEOF(dummy_real)
         gpu_w = c_ptr_plus(gpu_work, off)
         cuda_error = &
            cudaMemcpyAsync_D2D(gpu_w, gpu_u, blkm*nelim*C_SIZEOF(dummy_real), &
               stream)
         call presolve_single_node( stream, cublas, tile_size, blkm, &
            nelim, gpu_work, &
            nelim, gpu_w, blkm, gpu_u, blkm, cuda_error, cublas_error )
      end if
      if(cuda_error.ne.0.or.cublas_error.ne.0) return
   end do
end subroutine presolve_first

!
! This subroutine does full pre-solve for a single nrows-by-ncols node
!
subroutine presolve_single_node(stream, cublas, tile_size, nrows, ncols, &
      gpu_Di, ldd, gpu_L, ldl, gpu_Li, ldi, cuda_error, cublas_error)
   type(C_PTR), intent(in) :: stream
   type(C_PTR) :: cublas
   integer, intent(in) :: tile_size
   integer, intent(in) :: nrows
   integer, intent(in) :: ncols
   type(C_PTR) :: gpu_Di ! pointer to the array storing the block-diagonal
      ! matrix with diagonal tiles' inverses on the diagonal
   integer, intent(in) :: ldd ! the array's leading dimension 
   type(C_PTR) :: gpu_L ! pointer to the node matrix
   integer, intent(in) :: ldl ! its leading dimension
   type(C_PTR) :: gpu_Li ! pointer to the pre-solved node matris
   integer, intent(in) :: ldi ! its leading dimension
   integer, intent(out) :: cuda_error ! error flag
   integer, intent(out) :: cublas_error ! error flag
  
   integer :: tiles
   integer :: stride
   integer :: step
   integer :: n, m, k
   integer :: i
   integer(long) :: off
   type(C_PTR) :: gpu_a, gpu_b, gpu_c
   real(wp) :: dummy_real

   cuda_error = 0; cublas_error = 0
  
   off = ncols
   cuda_error = cudaMemcpy2DAsync(gpu_Li, ldi*C_SIZEOF(dummy_real), gpu_Di, &
      ldd*C_SIZEOF(dummy_real), ncols*C_SIZEOF(dummy_real), off, &
      cudaMemcpyDeviceToDevice, stream)
   if(cuda_error.ne.0) return

   if(nrows.gt.ncols) then
      gpu_a = c_ptr_plus(gpu_L, ncols*C_SIZEOF(dummy_real))
      gpu_b = c_ptr_plus(gpu_Li, ncols*C_SIZEOF(dummy_real))
      cublas_error = cublasDgemm(cublas, 'N', 'N', nrows - ncols, ncols, ncols,&
         -ONE, gpu_a, ldl, gpu_Li, ldi, ZERO, gpu_b, ldi)
      if(cublas_error.ne.0) return
   end if
  
   tiles = (ncols - 1)/tile_size + 1
   stride = 6
   i = tiles*tile_size
   do while ( tiles > 0 )
      stride = min(stride, tiles)
      do step = 1, stride - 1
         i = i - tile_size
         n = nrows - i
         m = i - (tiles - stride)*tile_size
         k = min(ncols - i, tile_size)
         off = (i + i*ldi)*C_SIZEOF(dummy_real)
         gpu_a = c_ptr_plus(gpu_Li, off)
         off = (i + (tiles - stride)*tile_size*ldl)*C_SIZEOF(dummy_real)
         gpu_b = c_ptr_plus(gpu_L, off)
         off = (i + (tiles - stride)*tile_size*ldi)*C_SIZEOF(dummy_real)
         gpu_c = c_ptr_plus(gpu_Li, off)
         cublas_error = cublasDgemm(cublas, 'N', 'N', n, m, k, -ONE, gpu_a, &
            ldi, gpu_b, ldl, ONE, gpu_c, ldi)
         if(cublas_error.ne.0) return
      end do
      i = i - tile_size
      if(i.lt.1) exit
      n = nrows - i
      m = i
      k = min(ncols-i, stride*tile_size)
      off = (i + i*ldl)*C_SIZEOF(dummy_real)
      gpu_a = c_ptr_plus(gpu_Li, off)
      gpu_b = c_ptr_plus(gpu_L, i*C_SIZEOF(dummy_real))
      gpu_c = c_ptr_plus(gpu_Li, i*C_SIZEOF(dummy_real))
      cublas_error = cublasDgemm(cublas, 'N', 'N', n, m, k, -ONE, gpu_a, ldl, &
         gpu_b, ldl, ONE, gpu_c, ldi)
      if(cublas_error.ne.0) return
      tiles = tiles - stride
   end do
end subroutine presolve_single_node

!
! This subroutine finalizes presolve by computing -L_o L_d**(-1)
!
subroutine presolve_second(stream, nnodes, nodes, sptr, rptr, nlev, tree_lev, &
      tree_le, node_list, tile_size, ntlev, tleve, maxtc, tlev, gpu_work, st, &
      cuda_error)
   type(C_PTR), intent(in) :: stream
   integer, intent(in) :: nnodes
   type(node_type), intent(in) :: nodes(nnodes)
   integer, intent(in) :: sptr(nnodes + 1)
   integer(long), intent(in) :: rptr(nnodes + 1)
   integer, intent(in) :: nlev ! # of elimination tree levels
   type(eltree_level), intent(in) :: tree_lev(nlev) ! levels' data
   integer, intent(in) :: tree_le(nlev + 1) ! levels' entry points into the
      ! list of nodes below
   integer, intent(in) :: node_list(tree_le(nlev + 1) - 1) ! list of nodes
   integer, intent(in) :: tile_size
   integer, intent(in) :: ntlev ! # of tile-width levels
   integer, intent(out) :: tleve(ntlev + 1) ! tile-width levels' entry points
   integer, intent(in) :: maxtc ! max # of tile columns in a node
   integer, intent(out) :: tlev(maxtc + 1) ! tlev(i) is the tile-width level of
      ! a node that is i tile columns wide
   type(C_PTR), intent(in) :: gpu_work ! workspace
   integer, intent(out) :: st ! allocation status
   integer, intent(out) :: cuda_error
  
   integer :: lev
   integer :: ncb
   integer :: node
   integer :: blkm
   integer :: blkn
   integer :: nelim
   integer :: ndelay
   integer :: ntiles

   integer(long) :: level_size
   integer(long) :: offp
  
   integer :: i, j, k

   integer(C_SIZE_T) :: sz

   type(node_solve_data), allocatable, target :: sdata(:)

   type(C_PTR) :: gpu_levL
   type(C_PTR) :: gpu_solve_data
  
   type(C_PTR) :: gpu_u, gpu_v, gpu_w
   real(wp) :: dummy_real

   cuda_error = 0

   do lev = 1, nlev
      ! compute the number of CUDA blocks for simultaneous computation
      ! of the product -L_o L_d**(-1) for all nodes on this tree level
      ncb = 0
      level_size = 0
      do j = tree_le(lev), tree_le(lev + 1) - 1
         node = node_list(j)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         blkn = sptr(node + 1) - sptr(node) + ndelay
         level_size = level_size + (blkm + 2)*blkn

         if(nelim.lt.1) cycle
         ntiles = (nelim - 1)/tile_size + 1
         i = tlev(ntiles)
         ! one-node tile-width levels has already been processed
         if ( tleve(i + 1) - tleve(i) < 2 ) cycle

         ncb = ncb + ((blkm - nelim - 1)/32 + 1)*((nelim - 1)/32 + 1)
      end do
      if(ncb.lt.1) cycle

      allocate(sdata(ncb), stat=st)
      if(st.ne.0) return

      gpu_levL = tree_lev(lev)%ptr_levL
      sz = level_size*C_SIZEOF(dummy_real)
      cuda_error = cudaMemcpyAsync_D2D(gpu_work, gpu_levL, sz, stream)
      if(cuda_error.ne.0) return

      ! set up matrix multiplication data
      k = 0
      offp = 0
      do j = tree_le(lev), tree_le(lev + 1) - 1
         node = node_list(j)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         blkn = sptr(node + 1) - sptr(node) + ndelay

         if(nelim.lt.1) then
            offp = offp + (blkm + 2)*blkn
            cycle
         end if
         ntiles = (nelim - 1)/tile_size + 1
         i = tlev(ntiles)
         if(tleve(i + 1)-tleve(i) .lt. 2) then
            offp = offp + (blkm + 2)*blkn
            cycle
         end if

         sz = offp*C_SIZEOF(dummy_real)
         gpu_u = c_ptr_plus(gpu_levL, sz)
         offp = offp + nelim
         sz = offp*C_SIZEOF(dummy_real)
         gpu_v = c_ptr_plus(gpu_levL, sz)
         gpu_w = c_ptr_plus(gpu_work, sz)
         k = multinode_dgemm_setup(blkm - nelim, nelim, nelim, gpu_w, blkm, &
            gpu_u, blkm, gpu_v, blkm, gpu_v, 1, sdata, k)
         offp = offp + (blkm + 2)*blkn - nelim
      end do
      
      cuda_error = cudaMalloc(gpu_solve_data, ncb*C_SIZEOF(sdata(1)))
      if(cuda_error.ne.0) return
      cuda_error = cudaMemcpyAsync_H2D(gpu_solve_data, C_LOC(sdata), &
         ncb*C_SIZEOF(sdata(1)), stream)
      if(cuda_error.ne.0) return

      ! perform simultaneous multiplication
      call multinode_dgemm_n( stream, ncb, gpu_solve_data, -ONE )

      deallocate(sdata,stat=st)
      if(st.ne.0) return
      cuda_error = cudaFree( gpu_solve_data )
      if(cuda_error.ne.0) return
   end do
end subroutine presolve_second

!****************************************************************************
!
! Build a post-factorization direct mapping between a node's rlist and that
! of its parent such that updates during forward solve can be done as:
! parent_data(rlist_direct(i,1)) += child_contribution(i),
! where parent_data is either its portion of the rhs or its contribution.
!

!
! This function computes the map size.
!
integer(long) function rlist_direct_size(nnodes, nodes, rptr) result ( nrd )
   integer, intent(in) :: nnodes
   type(node_type), intent(in) :: nodes(nnodes)
   integer(long), dimension(*), intent(in) :: rptr
   
   integer :: node
   
   nrd = 0
   do node = 1, nnodes
      nrd = nrd + rptr(node + 1) - rptr(node) + nodes(node)%ndelay
   end do
end function rlist_direct_size

!
! This subroutine computes rlist_direct.
!
subroutine rebuild_rlist_direct(n, nnodes, nodes, sptr, rptr, rlist, &
      child_ptr, child_list, nrd, rlist_direct, st)
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   type(node_type), intent(inout) :: nodes(nnodes)
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer(long), intent(in) :: nrd
   integer, intent(out) :: rlist_direct(nrd,2)
   integer, intent(out) :: st

   integer(long) :: nr
   integer :: ci
   integer :: node
   integer :: nelim
   integer :: ndelay
   integer :: blkm
   integer :: blkn
   integer :: cnode
   integer :: cnelim
   integer :: cndelay
   integer :: cblkm
   integer :: cblkn
   
   integer :: i, j
   integer(long) :: ii

   integer, dimension(:), pointer :: lperm
   integer, allocatable :: map(:)
   
   allocate(map(n), stat=st)
   if(st.ne.0) return
   map = 0

   nr = 1
   do node = 1, nnodes
      ndelay = nodes(node)%ndelay
      blkm = int(rptr(node + 1) - rptr(node)) + ndelay
      nodes(node)%rdptr = nr
      nr = nr + blkm
   end do
   do node = 1, nnodes
      if ( child_ptr(node + 1) - child_ptr(node) < 1 ) cycle
        nelim = nodes(node)%nelim
        ndelay = nodes(node)%ndelay
        blkn = sptr(node+1) - sptr(node) + ndelay
        blkm = int(rptr(node+1) - rptr(node)) + ndelay
        lperm => nodes(node)%perm
        do i = 1, blkn
          map(lperm(i)) = i
        end do
        j = blkn
        do ii = rptr(node) + blkn - ndelay, rptr(node + 1) - 1
          j = j + 1
          map(rlist(ii)) = j
        end do
        do ci = child_ptr(node), child_ptr(node + 1) - 1
          cnode = child_list(ci)
          cnelim = nodes(cnode)%nelim
          cndelay = nodes(cnode)%ndelay
          cblkn = sptr(cnode + 1) - sptr(cnode)
          cblkm = int(rptr(cnode + 1) - rptr(cnode)) + cndelay
          lperm => nodes(cnode)%perm
          nr = nodes(cnode)%rdptr
          do j = 1, cblkn + cndelay
            rlist_direct(nr + j - 1, 1) = map(lperm(j))
          end do
          j = cblkn + cndelay
          do ii = rptr(cnode) + cblkn, rptr(cnode + 1) - 1
             rlist_direct(nr + j, 1) = map(rlist(ii))
             j = j + 1
          end do
          !
          ! rlist_direct(:,2) is used to determine which row of the cild's
          ! contribution goes where during forward solve
          j = 0
          do i = 1, cblkm
            if ( rlist_direct(nr + i - 1, 1) <= nelim ) then
              j = j + 1
              rlist_direct(nr + j - 1, 2) = i
            end if
          end do
          !
          ! this many rows from those just listed in rlist_direct(:,2) 
          ! will update the parent's rows of the rhs
          nodes(cnode)%ncpdb = j - cnelim
          !
          ! the rest will be added to parent's contribution
          do i = 1, cblkm
            if ( rlist_direct(nr + i - 1, 1) > nelim ) then
              j = j + 1
              rlist_direct(nr + j - 1, 2) = i
            end if
          end do
        end do
        lperm => nodes(node)%perm
        do i = 1, blkn
          map(lperm(i)) = 0
        end do
        j = blkn
        do ii = rptr(node) + blkn - ndelay, rptr(node + 1) - 1
          j = j + 1
          map(rlist(ii)) = 0
        end do
   end do
end subroutine rebuild_rlist_direct

!
! This subroutine sets up the data for the solve phase
!
subroutine solve_setup(stream, pos_def, sparent, child_ptr, child_list, n, &
      invp, nnodes, nodes, sptr, rptr, rlist, rlist_direct, fact_data, st, &
      cuda_error)
   type(C_PTR), intent(in) :: stream
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: sparent
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: invp
   ! Note: gfortran-4.3 bug requires explicit size of nodes array
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   integer(C_INT), dimension(*), intent(in), target :: rlist_direct
   type(gpu_type), intent(inout) :: fact_data
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error
   
   integer :: node
   integer :: blkm
   integer :: blkn
   integer :: ndelay
   integer :: nelim
   integer :: cnode
   integer :: cblkm
   integer :: cblkn
   integer :: cnelim
   integer :: cn
   integer :: cnd
   integer :: lev, li
   integer :: lx_size ! # solution components computed on a (el. tree) level
   integer :: max_lx_size ! max of the above over all levels
   integer :: lc_size ! # contributions computed on a level
   integer :: max_lc_size ! max of the above
   integer :: lcc_size ! # of children contributions from previous level
   integer :: ln_size ! sum of nodes' heights on a level
   integer :: offx, offc ! offsets into level's x and contributions
   integer :: ncb ! # of CUDA blocks for multisolve
   integer :: nr, nc, nd
   integer :: i, j, k, m

   integer, dimension(:), pointer :: lperm   

   integer(long) :: rdptr ! node's entry point into rlist_direct
   integer(long) :: rd_size ! # rows in rlist_direct
   integer(long) :: sz ! SIZE_T integer
   integer(long) :: off ! an offset
   integer(long) :: ii

   integer, allocatable :: mask(:) ! this subtree mask
   integer(C_INT), allocatable, target :: col_ind(:) ! cols order after fact.
   integer(C_INT), allocatable, target :: row_ind(:,:) ! source and destination
      ! indices for cuda_gather_dx kernel
   integer(C_LONG), allocatable, target :: lwork(:) ! work array

   integer :: nch, total_nch, max_nch ! # node's children, level's total and max
   integer :: ncc, nccd ! child's contrib: total, to the parent's diag block
   integer :: ncp, max_ncp ! # child-parent pairs, its max
   integer :: cpi ! a child-parent pair index
   type(assemble_cp_type), dimension(:), allocatable, target :: cpdata
      ! child-parent data to be copied to GPU
   type(C_PTR) :: gpu_cpdata ! child-parent info

   integer :: nblk ! # CUDA blocks for assembly
   integer :: blk, bx, by
   type(assemble_blk_type), dimension(:), allocatable, target :: blkdata
      ! block-level data to be copied to GPU
   type(C_PTR) :: gpu_blkdata ! Output block-by-block info

   type(node_solve_data), allocatable, target :: sdata(:)
   type(C_PTR) :: gpu_solve_data
  
   type(C_PTR) :: gpu_u, gpu_v ! tmp pointers
   real(wp) :: dummy_real
   integer(C_INT) :: dummy_int
   
   ! post-factorization rlist_direct has two columns: see rebuild_rlist_direct

   ! post-factorization rlist_direct height == 
   ! total number of rows in nodal matrices
   rd_size = rlist_direct_size(nnodes, nodes, rptr)
   fact_data%rd_size = rd_size

   ! copy rlist_direct to GPU
   cuda_error = cudaMalloc(fact_data%gpu_rlist_direct, &
      2*rd_size*C_SIZEOF(rlist_direct(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpyAsync_H2D(fact_data%gpu_rlist_direct, &
      C_LOC(rlist_direct), 2*rd_size*C_SIZEOF(rlist_direct(1)), stream)
   if(cuda_error.ne.0) return

   allocate(col_ind(n), row_ind(rd_size,2), fact_data%off_lx(nnodes), &
      fact_data%off_lc(nnodes), fact_data%off_ln(nnodes), mask(nnodes), stat=st)
   if(st.ne.0) return
   
   ! mask the nodes of processed subtree
   mask = 0
   do lev = 1, fact_data%num_levels
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         mask(node) = 1
      end do
   end do

   ! find the children from child subtree
   do lev = 1, fact_data%num_levels
      m = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         do cn = child_ptr(node), child_ptr(node + 1) - 1
            cnode = child_list(cn)
            if(mask(cnode).eq.0) m = m + 1
         end do
      end do
      fact_data%values_L(lev)%nimp = m
      if(m.eq.0) cycle
      allocate(fact_data%values_L(lev)%import(m), stat=st)
      if(st.ne.0) return
      m = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         do cn = child_ptr(node), child_ptr(node + 1) - 1
            cnode = child_list(cn)
            if(mask(cnode).eq.0) then
               m = m + 1
               fact_data%values_L(lev)%import(m) = cnode
            end if
         end do
      end do
   end do

   ! find the nodes that are top subtree's nodes' children
   do lev = 1, fact_data%num_levels
      m = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         if(sparent(node).gt.nnodes ) cycle
         if(mask(sparent(node)).eq.0) m = m + 1
      end do
      fact_data%values_L(lev)%nexp = m
      if(m.eq.0) cycle
      allocate(fact_data%values_L(lev)%export(m), stat=st)
      if(st.ne.0) return
      m = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         if(sparent(node).gt.nnodes) cycle
         if(mask(sparent(node)).eq.0) then
            m = m + 1
            fact_data%values_L(lev)%export(m) = node
         end if
      end do
   end do

   max_lx_size = 0
   max_lc_size = 0
   max_ncp = 0
   nc = 0 ! column number
   nr = 0 ! nodes' row number

   do lev = 1, fact_data%num_levels

      lx_size = 0
      lc_size = 0
      total_nch = 0
      max_nch = 0
      fact_data%values_L(lev)%off_col_ind = nc
      fact_data%values_L(lev)%off_row_ind = nr
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkn = sptr(node + 1) - sptr(node)
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         lperm => nodes(node)%perm
         fact_data%off_ln(node) = nr - fact_data%values_L(lev)%off_row_ind
         do i = 1, nelim
            nc = nc + 1
            j = invp(lperm(i))
            col_ind(nc) = j
            nr = nr + 1
            row_ind(nr,1) = j ! source index for cuda_gather_dx
            row_ind(nr,2) = nc ! destination index
         end do
         do i = nelim + 1, blkn + ndelay
            nr = nr + 1
            row_ind(nr,1) = invp(lperm(i))
            row_ind(nr,2) = 0 ! multiplication by D flag is off
         end do
         do ii = rptr(node) + blkn, rptr(node + 1) - 1
            nr = nr + 1
            row_ind(nr,1) = invp(rlist(ii))
            row_ind(nr,2) = 0
         end do
         nch = child_ptr(node + 1) - child_ptr(node)
         total_nch = total_nch + nch
         max_nch = max(max_nch, nch)
         fact_data%off_lx(node) = lx_size
         fact_data%off_lc(node) = lc_size
         lx_size = lx_size + nelim
         lc_size = lc_size + blkm - nelim
      end do
      max_lx_size = max(max_lx_size, lx_size)
      max_lc_size = max(max_lc_size, lc_size)
      fact_data%values_L(lev)%lx_size = lx_size
      fact_data%values_L(lev)%lc_size = lc_size
      fact_data%values_L(lev)%ln_size = lx_size + lc_size
      fact_data%values_L(lev)%total_nch = total_nch

      lcc_size = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         do cn = child_ptr(node), child_ptr(node + 1) - 1
            cnode = child_list(cn)
            if ( mask(cnode) == 0 ) cycle
            cnelim = nodes(cnode)%nelim
            cnd = nodes(cnode)%ndelay
            cblkm = int(rptr(cnode + 1) - rptr(cnode)) + cnd
            fact_data%off_lc(cnode) = lcc_size
            lcc_size = lcc_size + cblkm - cnelim
         end do
      end do
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         do cn = child_ptr(node), child_ptr(node + 1) - 1
            cnode = child_list(cn)
            if ( mask(cnode) /= 0 ) cycle
            cnelim = nodes(cnode)%nelim
            cnd = nodes(cnode)%ndelay
            cblkm = int(rptr(cnode + 1) - rptr(cnode)) + cnd
            fact_data%off_lc(cnode) = lcc_size
            lcc_size = lcc_size + cblkm - cnelim
         end do
      end do
      fact_data%values_L(lev)%lcc_size = lcc_size
      if ( lev > 1 ) fact_data%values_L(lev - 1)%lc_size = lcc_size
      max_lc_size = max(max_lc_size, lcc_size)
      
   end do
   
   deallocate(mask, stat=st)
   if(st.ne.0) return
    
   do lev = 1, fact_data%num_levels

      max_nch = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         nch = child_ptr(node + 1) - child_ptr(node)
         max_nch = max(max_nch, nch)
      end do

      lx_size = fact_data%values_L(lev)%lx_size
      lc_size = fact_data%values_L(lev)%lc_size
      ln_size = fact_data%values_L(lev)%ln_size
      lcc_size = fact_data%values_L(lev)%lcc_size
      total_nch = fact_data%values_L(lev)%total_nch

      if(total_nch.gt.0) then
         allocate(cpdata(total_nch), stat=st)
         if(st.ne.0) return

         cpi = 1
         nblk = 0
         do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
            node = fact_data%lvllist(li)
            offx = fact_data%off_lx(node)
            ndelay = nodes(node)%ndelay
            blkm = int(rptr(node + 1) - rptr(node)) + ndelay
            blk = 0
            do cn = child_ptr(node), child_ptr(node + 1) - 1
               cnode = child_list(cn)
               cnelim = nodes(cnode)%nelim
               offc = fact_data%off_lc(cnode)
               cnd = nodes(cnode)%ndelay
               cblkn = sptr(cnode + 1) - sptr(cnode)
               cblkm = int(rptr(cnode + 1) - rptr(cnode)) + cnd
               rdptr = nodes(cnode)%rdptr
               cpdata(cpi)%pvoffset = offx
               cpdata(cpi)%ldp = lx_size
               cpdata(cpi)%cm = nodes(cnode)%ncpdb
               cpdata(cpi)%cn = 1
               cpdata(cpi)%ldc = lcc_size
               cpdata(cpi)%cvoffset = offc - cnelim
               off = (rd_size + rdptr + cnelim - 1)*C_SIZEOF(dummy_int)
               cpdata(cpi)%ind = c_ptr_plus(fact_data%gpu_rlist_direct, off)
               off = (rdptr - 1)*C_SIZEOF(dummy_int)
               cpdata(cpi)%rlist_direct = &
                  c_ptr_plus(fact_data%gpu_rlist_direct, off)
               cpdata(cpi)%sync_offset = max(0, cpi-1 - 1)
               cpdata(cpi)%sync_wait_for = blk
               bx = (cpdata(cpi)%cm - 1) / HOGG_ASSEMBLE_TX + 1
               by = 1
               nblk = nblk + bx*by
               blk = bx*by
               cpi = cpi + 1
            end do
         end do
    
         allocate(blkdata(nblk), stat=st)
         if(st.ne.0) return

         i = 1
         do nch = 1, max_nch
            cpi = 1
            do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
               node = fact_data%lvllist(li)
               if( child_ptr(node) + nch - 1 .lt. child_ptr(node + 1) ) then
                  bx = (cpdata(cpi + nch - 1)%cm - 1) / HOGG_ASSEMBLE_TX + 1
                  by = 1
                  do blk = 0, bx*by - 1
                     blkdata(i)%cp = cpi + (nch - 1) - 1 ! 0 indexed
                     blkdata(i)%blk = blk ! 0 indexed
                     i = i + 1
                  end do
               endif
               cpi = cpi + child_ptr(node + 1) - child_ptr(node)
            end do
         end do
    
         ncp = size(cpdata)
         max_ncp = max(max_ncp, ncp)
         sz = ncp*C_SIZEOF(cpdata(1))
         cuda_error = cudaMalloc(gpu_cpdata, sz)
         if(cuda_error.ne.0) return
         cuda_error = cudaMemcpyAsync_H2D(gpu_cpdata, C_LOC(cpdata), sz, stream)
         if(cuda_error.ne.0) return
    
         sz = nblk*C_SIZEOF(blkdata(1))
         cuda_error = cudaMalloc(gpu_blkdata, sz)
         if(cuda_error.ne.0) return
         cuda_error = &
            cudaMemcpyAsync_H2D(gpu_blkdata, C_LOC(blkdata), sz, stream)
         if(cuda_error.ne.0) return

         fact_data%values_L(lev)%ncp_pre = ncp
         fact_data%values_L(lev)%gpu_cpdata_pre = gpu_cpdata
         fact_data%values_L(lev)%ncb_asm_pre = nblk
         fact_data%values_L(lev)%gpu_blkdata_pre = gpu_blkdata

         deallocate(cpdata, blkdata, stat=st)
         if(st.ne.0) return

      else

         fact_data%values_L(lev)%ncp_pre = 0
         fact_data%values_L(lev)%ncb_asm_pre = 0

      end if

      ncb = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         blkn = sptr(node + 1) - sptr(node) + ndelay
         if ( nelim < 1 ) cycle
         ncb = ncb + (blkm - 1)/64 + 1
      end do
      fact_data%values_L(lev)%ncb_slv_n = ncb

      if ( ncb > 0 ) then

         allocate (sdata(ncb), stat=st)
         if(st.ne.0) return

         cuda_error = cudaMalloc(gpu_solve_data, ncb*C_SIZEOF(sdata(1)))
         if(cuda_error.ne.0) return

         k = 0
         do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
            node = fact_data%lvllist(li)
            offx = fact_data%off_lx(node)
            offc = fact_data%off_lc(node)
            ndelay = nodes(node)%ndelay
            nelim = nodes(node)%nelim
            blkm = int(rptr(node + 1) - rptr(node)) + ndelay
            blkn = sptr(node + 1) - sptr(node) + ndelay
            if ( nelim < 1 ) cycle
            do i = 1, (blkm - 1)/64 + 1
               sdata(k + i)%off_a = fact_data%off_L(node)
               sdata(k + i)%off_b = offx
               sdata(k + i)%off_u = offx
               sdata(k + i)%off_v = offc
               sdata(k + i)%lda = blkm
               sdata(k + i)%ldb = lx_size
               sdata(k + i)%ldu = lx_size
               sdata(k + i)%ldv = lc_size
               sdata(k + i)%nrows = blkm
               sdata(k + i)%ncols = nelim
               sdata(k + i)%offb = k
            end do
            k = k + (blkm - 1)/64 + 1
         end do
         cuda_error = cudaMemcpyAsync_H2D(gpu_solve_data, C_LOC(sdata), &
            ncb*C_SIZEOF(sdata(1)), stream)
         if(cuda_error.ne.0) return
         fact_data%values_L(lev)%gpu_solve_n_data = gpu_solve_data

         deallocate(sdata, stat=st)
         if(st.ne.0) return
      end if

      ncb = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         blkn = sptr(node + 1) - sptr(node) + ndelay
         if ( nelim < 1 ) cycle
         ncb = ncb + (nelim - 1)/8 + 1
      end do
      fact_data%values_L(lev)%ncb_slv_t = ncb

      if(ncb.gt.0) then

         allocate(sdata(ncb), stat=st)
         if(st.ne.0) return

         cuda_error = cudaMalloc(gpu_solve_data, ncb*C_SIZEOF(sdata(1)))
         if(cuda_error.ne.0) return

         k = 0
         do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
            node = fact_data%lvllist(li)
            offx = fact_data%off_lx(node)
            offc = fact_data%off_ln(node)
            ndelay = nodes(node)%ndelay
            nelim = nodes(node)%nelim
            blkm = int(rptr(node + 1) - rptr(node)) + ndelay
            blkn = sptr(node + 1) - sptr(node) + ndelay
            if ( nelim < 1 ) cycle
            do i = 1, (nelim - 1)/8 + 1
               sdata(k + i)%off_a = fact_data%off_L(node)
               sdata(k + i)%off_b = offc
               sdata(k + i)%off_u = offx
               sdata(k + i)%off_v = offx
               sdata(k + i)%lda = blkm
               sdata(k + i)%ldb = ln_size
               sdata(k + i)%ldu = lx_size
               sdata(k + i)%ldv = lx_size
               sdata(k + i)%nrows = blkm
               sdata(k + i)%ncols = nelim
               sdata(k + i)%offb = k
            end do
            k = k + (nelim - 1)/8 + 1
         end do
         cuda_error = &
            cudaMemcpyAsync_H2D(gpu_solve_data, C_LOC(sdata), &
               ncb*C_SIZEOF(sdata(1)), stream)
         if(cuda_error.ne.0) return
         fact_data%values_L(lev)%gpu_solve_t_data = gpu_solve_data

         deallocate(sdata, stat=st)
         if(st.ne.0) return

      end if

      ncp = 0
      nblk = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         blkm = int(rptr(node + 1) - rptr(node))
         blkn = sptr(node + 1) - sptr(node)
         m = blkm - blkn
         if ( m > 0 ) then
            do nch = child_ptr(node), child_ptr(node+1)-1
               cnode = child_list(nch)
               cnelim = nodes(cnode)%nelim
               cnd = nodes(cnode)%ndelay
               cblkm = int(rptr(cnode + 1) - rptr(cnode)) + cnd
               rdptr = nodes(cnode)%rdptr
               ncc = cblkm - cnelim
               nccd = nodes(cnode)%ncpdb
               if ( ncc - nccd > 0 ) then
                  ncp = ncp + 1
                  bx = (ncc - nccd - 1) / HOGG_ASSEMBLE_TX + 1
                  by = 1
                  nblk = nblk + bx*by
               endif
            end do
         end if
      end do
      max_ncp = max(max_ncp, ncp)
    
      if(ncp.gt.0) then
         allocate(cpdata(ncp), blkdata(nblk), stat=st)
         if(st.ne.0) return
    
         cpi = 1
         do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
            node = fact_data%lvllist(li)
            offx = fact_data%off_lc(node)
            nelim = nodes(node)%nelim
            blkm = int(rptr(node + 1) - rptr(node))
            blkn = sptr(node + 1) - sptr(node)
            blk = 0
            do nch = child_ptr(node), child_ptr(node+1)-1
               cnode = child_list(nch)
               cnelim = nodes(cnode)%nelim
               offc = fact_data%off_lc(cnode)
               cnd = nodes(cnode)%ndelay
               cblkn = sptr(cnode + 1) - sptr(cnode)
               cblkm = int(rptr(cnode + 1) - rptr(cnode)) + cnd
               rdptr = nodes(cnode)%rdptr
               ncc = cblkm - cnelim
               nccd = nodes(cnode)%ncpdb
               if(ncc-nccd.le.0) cycle
               cpdata(cpi)%cm = ncc - nccd
               cpdata(cpi)%cn = 1 !nrhs
               cpdata(cpi)%ldp = lc_size
               cpdata(cpi)%pvoffset = offx - nelim
               cpdata(cpi)%ldc = lcc_size
               cpdata(cpi)%cvoffset = offc - cnelim
               off = (rd_size + rdptr + cnelim + nccd - 1)*C_SIZEOF(dummy_int)
               cpdata(cpi)%ind = c_ptr_plus(fact_data%gpu_rlist_direct, off)
               off = (rdptr - 1)*C_SIZEOF(dummy_int)
               cpdata(cpi)%rlist_direct = &
                  c_ptr_plus(fact_data%gpu_rlist_direct, off)
               cpdata(cpi)%sync_offset = max(0, cpi-1 - 1) ! 0-indexed
               cpdata(cpi)%sync_wait_for = blk
               bx = (cpdata(cpi)%cm - 1) / HOGG_ASSEMBLE_TX + 1
               by = 1
               blk = bx*by ! Store for next iteration
               cpi = cpi + 1
            end do
         end do

         i = 1
         do nch = 1, max_nch
            cpi = 1
            do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
               node = fact_data%lvllist(li)
               do j = child_ptr(node), child_ptr(node + 1) - 1
                  cnode = child_list(j)
                  cnelim = nodes(cnode)%nelim
                  cnd = nodes(cnode)%ndelay
                  cblkm = int(rptr(cnode + 1) - rptr(cnode)) + cnd
                  ncc = cblkm - cnelim
                  nccd = nodes(cnode)%ncpdb
                  if(ncc-nccd.le.0) cycle
                  if(j-child_ptr(node)+1.eq.nch) then
                     bx = (cpdata(cpi)%cm - 1) / HOGG_ASSEMBLE_TX + 1
                     by = 1
                     do blk = 1, bx*by
                        blkdata(i)%cp = cpi - 1
                        blkdata(i)%blk = blk - 1
                        i = i + 1
                     end do
                  endif
                  cpi = cpi + 1
               end do
            end do
         end do

         sz = ncp*C_SIZEOF(cpdata(1))
         cuda_error = cudaMalloc(gpu_cpdata, sz)
         if(cuda_error.ne.0) return
         cuda_error = cudaMemcpyAsync_H2D(gpu_cpdata, C_LOC(cpdata), sz, stream)
         if(cuda_error.ne.0) return
    
         sz = nblk*C_SIZEOF(blkdata(1))
         cuda_error = cudaMalloc(gpu_blkdata, sz)
         if(cuda_error.ne.0) return
         cuda_error = &
            cudaMemcpyAsync_H2D(gpu_blkdata, C_LOC(blkdata), sz, stream)
         if(cuda_error.ne.0) return
    
         fact_data%values_L(lev)%ncp_post = ncp
         fact_data%values_L(lev)%gpu_cpdata_post = gpu_cpdata
         fact_data%values_L(lev)%ncb_asm_post = nblk
         fact_data%values_L(lev)%gpu_blkdata_post = gpu_blkdata

         deallocate(cpdata, blkdata, stat=st)
         if(st.ne.0) return

      else

         fact_data%values_L(lev)%ncp_post = 0
         fact_data%values_L(lev)%ncb_asm_post = 0

      end if

   end do

   fact_data%max_lx_size = max_lx_size
   fact_data%max_lc_size = max_lc_size
  
   cuda_error = cudaMalloc(fact_data%gpu_sync, (max_ncp+1)*C_SIZEOF(dummy_int))
   if(cuda_error.ne.0) return

   ! move col_ind and row_ind to GPU
   cuda_error = cudaMalloc(fact_data%gpu_col_ind, n*C_SIZEOF(dummy_int))
   if(cuda_error.ne.0) return
   cuda_error = cudaMalloc(fact_data%gpu_row_ind, &
      2*rd_size*C_SIZEOF(rlist_direct(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpyAsync_H2D(fact_data%gpu_col_ind, C_LOC(col_ind), &
      n*C_SIZEOF(col_ind(1)), stream)
   if(cuda_error.ne.0) return

   cuda_error = cudaMemcpyAsync_H2D(fact_data%gpu_row_ind, C_LOC(row_ind), &
      2*rd_size*C_SIZEOF(row_ind(1,1)), stream)
   if(cuda_error.ne.0) return
   cuda_error = cudaDeviceSynchronize()
   if(cuda_error.ne.0) return

   deallocate(col_ind, row_ind, stat=st)
   if(st.ne.0) return
   
   ! collect D factor from nodes into an array

   cuda_error = cudaMalloc(fact_data%gpu_diag, (2*n + 1)*C_SIZEOF(dummy_real))
   if(cuda_error.ne.0) return
   if ( pos_def ) return
   
   cuda_error = cudaMalloc(gpu_v, 2*max_lx_size*C_SIZEOF(dummy_int))
   if(cuda_error.ne.0) return

   allocate(lwork(n), stat=st)
   if(st.ne.0) return

   nd = 0
   do lev = 1, fact_data%num_levels
      nc = 0
      do li = fact_data%lvlptr(lev), fact_data%lvlptr(lev + 1) - 1
         node = fact_data%lvllist(li)
         ndelay = nodes(node)%ndelay
         nelim = nodes(node)%nelim
         blkn = sptr(node + 1) - sptr(node) + ndelay
         blkm = int(rptr(node + 1) - rptr(node)) + ndelay
         off = fact_data%off_L(node) + (blkm+0_long)*blkn
         do i = 1, nelim
            nc = nc + 1
            lwork(nc) = off + 2*i - 1
         end do
      end do
      off = (2*nd + 1)*C_SIZEOF(dummy_real)
      gpu_u = c_ptr_plus(fact_data%gpu_diag, off)
      cuda_error = cudaMemcpyAsync_H2D(gpu_v, C_LOC(lwork), &
         nc*C_SIZEOF(lwork(1)), stream)
      if(cuda_error.ne.0) return
      call gather_diag( stream, nc, fact_data%values_L(lev)%ptr_levL, &
         gpu_u, gpu_v )
      nd = nd + nc
   end do

   cuda_error = cudaFree(gpu_v)
   if(cuda_error.ne.0) return

end subroutine solve_setup
   
end module spral_ssids_factor_gpu
