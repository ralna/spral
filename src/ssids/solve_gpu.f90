module spral_ssids_solve_gpu
!$ use omp_lib
   use iso_c_binding
   use spral_cuda
   use spral_ssids_alloc, only : cuda_stack_alloc_type, custack_init, &
      custack_alloc, custack_free, custack_finalize
   use spral_ssids_cuda_datatypes
   use spral_ssids_cuda_interfaces
   use spral_ssids_datatypes
   implicit none

   private
   public :: bwd_solve_gpu,      & ! Backwards solve on GPU (presolve=0)
             bwd_multisolve_gpu, & ! Backwards solve on GPU (presolve=1)
             fwd_solve_gpu,      & ! Forwards solve on GPU (presolve=0)
             fwd_multisolve_gpu, & ! Forwards solve on GPU (presolve=1)
             d_solve_gpu,        & ! D solve on GPU
             setup_gpu_solve,    & ! Setup data strucutres prior to solve
             free_lookup_gpu       ! Cleanup data structures

   interface free_lookup_gpu
      module procedure free_lookup_gpu_fwd, free_lookup_gpu_bwd
   end interface free_lookup_gpu

contains

subroutine bwd_solve_gpu(job, posdef, nnodes, sptr, nstream, stream_handle, &
      stream_data, top_data, invp, x, st, cuda_error)
   integer, intent(in) :: job
   logical, intent(in) :: posdef
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, intent(in) :: nstream
   type(C_PTR), dimension(*), intent(in) :: stream_handle
   type(gpu_type), dimension(nstream), intent(in) :: stream_data
   type(gpu_type), intent(in) :: top_data
   integer, dimension(*), intent(in) :: invp
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(out) :: cuda_error
   integer, intent(out) :: st  ! stat parameter

   integer :: stream, i, n

   real(C_DOUBLE), dimension(:), allocatable, target :: xlocal
   type(C_PTR) :: gpu_x

   st = 0
   cuda_error = 0

   ! Push x on to GPU
   n = sptr(nnodes+1)-1
   allocate(xlocal(n), stat=st)
   if(st.ne.0) return
   do i = 1, n
      xlocal(i) = x(invp(i))
   end do
   cuda_error = cudaMalloc(gpu_x, n*C_SIZEOF(xlocal(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpy_h2d(gpu_x, C_LOC(xlocal), n*C_SIZEOF(xlocal(1)))
   if(cuda_error.ne.0) return

   ! Backwards solve for top part of tree
   call subtree_bwd_solve_gpu(job, posdef, top_data%num_levels,                &
      top_data%bwd_slv_lookup, top_data%bwd_slv_lwork, top_data%bwd_slv_nsync, &
      gpu_x, st, cuda_error, stream_handle(1))
   if(cuda_error.ne.0) return

   ! Backwards solve for lower parts of tree
   ! (As they are independant, order of loop doesn't matter)
   ! FIXME: Should probably be an OpenMP parallel do eventually
   do stream = 1, nstream
      call subtree_bwd_solve_gpu(job, posdef, stream_data(stream)%num_levels,  &
         stream_data(stream)%bwd_slv_lookup, stream_data(stream)%bwd_slv_lwork,&
         stream_data(stream)%bwd_slv_nsync, gpu_x, st, cuda_error,             &
         stream_handle(stream))
      if(cuda_error.ne.0) return
   end do 
   cuda_error = cudaDeviceSynchronize() ! Wait for streams to finish
   if(cuda_error.ne.0) return

   ! Bring x back from GPU
   cuda_error = cudaMemcpy_d2h(C_LOC(xlocal), gpu_x, n*C_SIZEOF(xlocal(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaFree(gpu_x)
   if(cuda_error.ne.0) return
   do i = 1, n
      x(invp(i)) = xlocal(i)
   end do
end subroutine bwd_solve_gpu

!*************************************************************************
!
! This subroutine performs a backwards solve on the chunk of nodes specified
! by sa:en.
!
subroutine subtree_bwd_solve_gpu(job, posdef, num_levels, bwd_slv_lookup, &
      lwork, nsync, gpu_x, st, cuda_error, stream)
   integer, intent(in) :: job
   logical, intent(in) :: posdef
   integer, intent(in) :: num_levels
   type(lookups_gpu_bwd), dimension(:) :: bwd_slv_lookup
   integer, intent(in) :: lwork
   integer, intent(in) :: nsync
   type(C_PTR), intent(inout) :: gpu_x
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error
   type(C_PTR), intent(in) :: stream

   integer(C_INT) :: dummy_int
   real(C_DOUBLE), target :: dummy_real
   type(C_PTR) :: gpu_work, gpu_sync
   integer(long) :: nrect, ndiag
   integer :: lvl

   logical(C_BOOL) :: dsolve, unit_diagonal

   nrect = 0; ndiag = 0
   cuda_error = 0
   st = 0

   if(posdef) then
      dsolve = .false. ! Never solve with D if we have an LL^T factorization
      unit_diagonal = .false.
   else ! indef
      dsolve = (job.ne.SSIDS_SOLVE_JOB_BWD) ! Do we solve L^T or (DL^T)?
      unit_diagonal = .true.
   endif

   ! Allocate workspace
   cuda_error = cudaMalloc(gpu_work, lwork*C_SIZEOF(dummy_real))
   if(cuda_error.ne.0) return
   cuda_error = cudaMalloc(gpu_sync, 2*nsync*C_SIZEOF(dummy_int))
   if(cuda_error.ne.0) return
   
   ! Backwards solve DL^Tx = z or L^Tx = z
   do lvl = num_levels, 1, -1
      call run_bwd_solve_kernels(dsolve, unit_diagonal, gpu_x, gpu_work, &
         nsync, gpu_sync, bwd_slv_lookup(lvl), stream)
   end do

   ! Free workspace
   cuda_error = cudaFree(gpu_work)
   if(cuda_error.ne.0) return
   cuda_error = cudaFree(gpu_sync)
   if(cuda_error.ne.0) return
end subroutine subtree_bwd_solve_gpu

subroutine d_solve_gpu(nnodes, sptr, nstream, stream_handle, &
      stream_data, top_data, n, invp, x, st, cuda_error)
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, intent(in) :: nstream
   type(C_PTR), dimension(*), intent(in) :: stream_handle
   type(gpu_type), dimension(nstream), intent(in) :: stream_data
   type(gpu_type), intent(in) :: top_data
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: invp
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(out) :: cuda_error
   integer, intent(out) :: st  ! stat parameter

   integer :: stream, i

   real(C_DOUBLE), dimension(:), allocatable, target :: xlocal
   type(C_PTR) :: gpu_x, gpu_y

   st = 0
   cuda_error = 0

   allocate(xlocal(sptr(nnodes+1)-1), stat=st)
   if(st.ne.0) return

   ! Allocate workspace on GPU (code doesn't work in place, so need in and out)
   cuda_error = cudaMalloc(gpu_x, 2*aligned_size(n*C_SIZEOF(xlocal(n))))
   if(cuda_error.ne.0) return
   gpu_y = c_ptr_plus_aligned(gpu_x, n*C_SIZEOF(xlocal(n)))

   ! Push x on to GPU
   do i = 1, n
      xlocal(i) = x(invp(i))
   end do
   cuda_error = cudaMemcpy_h2d(gpu_x, C_LOC(xlocal), n*C_SIZEOF(xlocal(1)))
   if(cuda_error.ne.0) return

   ! Backwards solve for top part of tree
   call subtree_d_solve_gpu(top_data%num_levels, top_data%bwd_slv_lookup, &
      gpu_x, gpu_y, stream_handle(1))
   if(cuda_error.ne.0) return

   ! Backwards solve for lower parts of tree
   ! (As they are independant, order of loop doesn't matter)
   ! FIXME: Should probably be an OpenMP parallel do eventually
   do stream = 1, nstream
      call subtree_d_solve_gpu(stream_data(stream)%num_levels,  &
         stream_data(stream)%bwd_slv_lookup, gpu_x, gpu_y, &
         stream_handle(stream))
   end do 
   cuda_error = cudaDeviceSynchronize() ! Wait for streams to finish
   if(cuda_error.ne.0) return

   ! Bring x back from GPU
   cuda_error = cudaMemcpy_d2h(C_LOC(xlocal), gpu_y, n*C_SIZEOF(xlocal(1)))
   if(cuda_error.ne.0) return
   ! Free GPU memory
   cuda_error = cudaFree(gpu_x)
   if(cuda_error.ne.0) return
   ! Undo permtuation
   do i = 1, n
      x(invp(i)) = xlocal(i)
   end do
end subroutine d_solve_gpu

!*************************************************************************
!
! This subroutine performs a D solve on the specified subtree
!
subroutine subtree_d_solve_gpu(num_levels, bwd_slv_lookup, gpu_x, gpu_y, stream)
   integer, intent(in) :: num_levels
   type(lookups_gpu_bwd), dimension(:) :: bwd_slv_lookup
   type(C_PTR), intent(inout) :: gpu_x
   type(C_PTR), intent(inout) :: gpu_y
   type(C_PTR), intent(in) :: stream

   integer :: lvl
   
   ! Diagonal solve Dy = x
   do lvl = num_levels, 1, -1
      call run_d_solve_kernel(gpu_x, gpu_y, bwd_slv_lookup(lvl), stream)
   end do
end subroutine subtree_d_solve_gpu

subroutine free_lookup_gpu_bwd(gpul, cuda_error)
   type(lookups_gpu_bwd), intent(inout) :: gpul
   integer, intent(out) :: cuda_error

   ! Note: only gpul%gemv is a cudaMalloc'd address, all others are just
   ! pointer addition from that location
   if(C_ASSOCIATED(gpul%gemv)) then
      cuda_error = cudaFree(gpul%gemv); gpul%gemv = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif
end subroutine free_lookup_gpu_bwd

subroutine fwd_solve_gpu(posdef, child_ptr, child_list, n, invp, nnodes, nodes,&
      rptr, nstream, stream_handle, stream_data, top_data, x, st, cuda_error)
   logical, intent(in) :: posdef
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(in) :: nodes
   integer(long), dimension(*), intent(in) :: rptr
   integer, intent(in) :: nstream
   type(C_PTR), dimension(*), intent(in) :: stream_handle
   type(gpu_type), dimension(*), intent(in) :: stream_data
   type(gpu_type), intent(in) :: top_data
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   integer :: i

   real(C_DOUBLE), dimension(:), allocatable, target :: xlocal
   type(C_PTR) :: gpu_x

   integer, dimension(:), allocatable :: cvmap
   type(C_PTR) :: gpu_cvalues

   integer :: blkm
   integer :: cn
   integer :: cnode
   integer :: ndelay
   integer :: nelim

   integer :: node

   integer :: stream
   type(C_PTR), dimension(:), allocatable, target :: cvalues

   integer(long) :: stack_ptr
   type(C_PTR) :: gpu_stack

   real(C_DOUBLE) :: dummy_real

   st = 0
   cuda_error = 0

   allocate(xlocal(n), stat=st)
   if(st.ne.0) return

   ! Push x on to GPU
   do i = 1, n
      xlocal(i) = x(invp(i))
   end do
   cuda_error = cudaMalloc(gpu_x, n*C_SIZEOF(xlocal(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpy_h2d(gpu_x, C_LOC(xlocal), n*C_SIZEOF(xlocal(1)))
   if(cuda_error.ne.0) return
   x(1:n) = xlocal(:)

   ! Build map from nodes to position in child lists
   allocate(cvmap(nnodes), stat=st)
   if(st.ne.0) return
   do node = 1, nnodes
      do cn = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(cn)
         cvmap(cnode) = cn-1
      end do
   end do

   ! Determine size of "stack" and allocate it
   ! (Not an actual stack?)
   stack_ptr = 0
   do node = 1, nnodes
      ndelay = nodes(node)%ndelay
      nelim = nodes(node)%nelim
      blkm = int(rptr(node+1) - rptr(node)) + ndelay
      stack_ptr = stack_ptr + blkm-nelim
   end do
   cuda_error = cudaMalloc(gpu_stack, stack_ptr*C_SIZEOF(dummy_real))
   if(cuda_error.ne.0) return

   ! Build index map for use of "stack" and copy to GPU
   allocate(cvalues(nnodes), stat=st)
   if(st.ne.0) return
   cuda_error = cudaMalloc(gpu_cvalues, nnodes*C_SIZEOF(cvalues(1)))
   if(cuda_error.ne.0) return
   stack_ptr = 0
   do node = 1, nnodes
      ndelay = nodes(node)%ndelay
      nelim = nodes(node)%nelim
      blkm = int(rptr(node+1) - rptr(node)) + ndelay
      if(blkm-nelim>0) cvalues(cvmap(node)+1) = &
         c_ptr_plus(gpu_stack, stack_ptr*C_SIZEOF(dummy_real))
      stack_ptr = stack_ptr + blkm-nelim
   end do
   cuda_error = cudaMemcpy_h2d(gpu_cvalues, C_LOC(cvalues), &
      nnodes*C_SIZEOF(cvalues(1)))
   if(cuda_error.ne.0) return
   deallocate(cvalues, stat=st)
   if(st.ne.0) return

   ! Forwards solve for lower parts of tree
   ! (As they are independant, order of loop doesn't matter)
   ! FIXME: Should probably be an OpenMP do eventually
   do stream = 1, nstream
      call subtree_fwd_solve_gpu_lvl(posdef, stream_data(stream)%num_levels,   &
         gpu_x, stream_data(stream)%fwd_slv_lookup,                            &
         stream_data(stream)%fwd_slv_lwork, stream_data(stream)%fwd_slv_nlocal,&
         stream_data(stream)%fwd_slv_nsync, stream_data(stream)%fwd_slv_nasync,&
         gpu_cvalues, cuda_error, stream_handle(stream))
      if(cuda_error.ne.0) return
   end do 

   cuda_error = cudaDeviceSynchronize() ! Wait for streams to finish
   if(cuda_error.ne.0) return

   ! Forwards solve for top part of tree
   call subtree_fwd_solve_gpu_lvl(posdef, top_data%num_levels, gpu_x,          &
      top_data%fwd_slv_lookup, top_data%fwd_slv_lwork, top_data%fwd_slv_nlocal,&
      top_data%fwd_slv_nsync, top_data%fwd_slv_nasync, gpu_cvalues, cuda_error,&
      stream_handle(1))
   if(cuda_error.ne.0) return

   cuda_error = cudaDeviceSynchronize() ! Wait for streams to finish
   if(cuda_error.ne.0) return

   ! Free memory
   cuda_error = cudaFree(gpu_cvalues)
   if(cuda_error.ne.0) return
   cuda_error = cudaFree(gpu_stack)
   if(cuda_error.ne.0) return
   
   ! Bring x back from GPU
   cuda_error = cudaMemcpy_d2h(C_LOC(xlocal), gpu_x, n*C_SIZEOF(dummy_real))
   if(cuda_error.ne.0) return
   cuda_error = cudaFree(gpu_x)
   if(cuda_error.ne.0) return
   do i = 1, n
      x(invp(i)) = xlocal(i)
   end do

end subroutine fwd_solve_gpu

subroutine subtree_fwd_solve_gpu_lvl(posdef, num_levels, gpu_x, fwd_slv_lookup,&
      lwork, nlocal, nsync, nasm_sync, gpu_cvalues, cuda_error, stream)
   logical, intent(in) :: posdef
   integer, intent(in) :: num_levels
   type(C_PTR), intent(inout) :: gpu_x
   type(lookups_gpu_fwd), dimension(*) :: fwd_slv_lookup
   integer, intent(in) :: lwork
   integer, intent(in) :: nlocal
   integer, intent(in) :: nsync
   integer, intent(in) :: nasm_sync
   type(C_PTR), intent(in) :: gpu_cvalues
   integer, intent(out) :: cuda_error
   type(C_PTR), intent(in) :: stream

   integer :: lvl

   type(C_PTR) :: gpu_work, gpu_sync, gpu_asm_sync, gpu_xlocal
   logical(C_BOOL) :: cposdef

   integer(C_INT) :: dummy_int
   real(C_DOUBLE) :: dummy_real

   cposdef = posdef ! convert from Fortran to C logical size

   ! Do actual work
   cuda_error = cudaMalloc(gpu_sync, 2*nsync*C_SIZEOF(dummy_int))
   if(cuda_error.ne.0) return
   cuda_error = cudaMalloc(gpu_asm_sync, (1+nasm_sync)*C_SIZEOF(dummy_int))
   if(cuda_error.ne.0) return
   cuda_error = cudaMalloc(gpu_xlocal, nlocal*C_SIZEOF(dummy_real))
   if(cuda_error.ne.0) return
   cuda_error = cudaMalloc(gpu_work, lwork*C_SIZEOF(dummy_real))
   if(cuda_error.ne.0) return
   do lvl = 1, num_levels
      call run_fwd_solve_kernels(cposdef, fwd_slv_lookup(lvl), gpu_xlocal, &
         gpu_cvalues, gpu_x, gpu_cvalues, gpu_work, nsync, gpu_sync,       &
         nasm_sync, gpu_asm_sync, stream)
   end do
   cuda_error = cudaFree(gpu_work)
   if(cuda_error.ne.0) return
   cuda_error = cudaFree(gpu_xlocal)
   if(cuda_error.ne.0) return
   cuda_error = cudaFree(gpu_asm_sync)
   if(cuda_error.ne.0) return
   cuda_error = cudaFree(gpu_sync)
   if(cuda_error.ne.0) return
end subroutine subtree_fwd_solve_gpu_lvl

subroutine free_lookup_gpu_fwd(gpu, cuda_error)
   type(lookups_gpu_fwd), intent(inout) :: gpu
   integer, intent(out) :: cuda_error

   ! Note: only gpu%assemble is a cudaMalloc'd location, others are all
   ! just pointer addition from that address
   cuda_error = cudaFree(gpu%assemble)
   if(cuda_error.ne.0) return
end subroutine free_lookup_gpu_fwd

subroutine create_gpu_lookup_fwd(nlvl, lvllist, nodes, child_ptr, child_list, &
      cvmap, sptr, rptr, rptr_with_delays, gpu_rlist_with_delays, gpu_clen,   &
      gpu_clists, gpu_clists_direct, gpul, nsync, nlocal, lwork, nasm_sync,   &
      stream, st, cuda_error)
   integer, intent(in) :: nlvl
   integer, dimension(*), intent(in) :: lvllist
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, dimension(*), intent(in) :: cvmap
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer(C_SIZE_T), dimension(*), intent(in) :: rptr_with_delays
   type(C_PTR), intent(in) :: gpu_rlist_with_delays
   type(C_PTR), intent(in) :: gpu_clen
   type(C_PTR), intent(in) :: gpu_clists
   type(C_PTR), intent(in) :: gpu_clists_direct
   type(lookups_gpu_fwd), intent(out) :: gpul
   integer, intent(out) :: nsync
   integer, intent(out) :: nlocal
   integer, intent(out) :: lwork
   integer, intent(out) :: nasm_sync
   type(C_PTR), intent(in) :: stream
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   integer :: i, j, ci, ni
   integer :: node, blkm, blkn, nelim, gldl, nchild, ndelay, syncblk
   integer :: child, cnelim, cblkm, cblkn, cndelay

   integer :: nx, ny
   integer :: nassemble, nassemble2, nasmblk, ntrsv, ngemv, nreduce, nscatter
   type(assemble_lookup_type), dimension(:), allocatable, target :: &
      assemble_lookup
   type(assemble_lookup2_type), dimension(:), allocatable, target :: &
      assemble_lookup2
   type(assemble_blk_type), dimension(:), allocatable, target :: asmblkdata
   type(trsv_lookup_type), dimension(:), allocatable, target :: trsv_lookup
   type(gemv_notrans_lookup), dimension(:), allocatable, target :: gemv_lookup
   type(reduce_notrans_lookup), dimension(:), allocatable, target :: &
      reduce_lookup
   type(scatter_lookup_type), dimension(:), allocatable, target :: &
      scatter_lookup

   integer(C_SIZE_T) :: sz
   type(C_PTR) :: pmem

   type(C_PTR) :: ptdummy_int
   integer(C_INT) :: dummy_int
   real(C_DOUBLE) :: dummy_real

   ! Initialize outputs
   st = 0; cuda_error = 0
   nsync = 0
   nlocal = 0
   lwork = 0

   ! Calculate size of lookups and allocate
   nassemble = 0
   nasm_sync = 0
   nassemble2 = 0
   nasmblk = 0
   ntrsv = 0
   ngemv = 0
   nreduce = 0
   nscatter = 0
   do ni = 1, nlvl
      node = lvllist(ni)
      ! Setup basic data about node
      ndelay = nodes(node)%ndelay
      nelim = nodes(node)%nelim
      blkn = sptr(node+1) - sptr(node) + ndelay
      blkm = int(rptr(node+1) - rptr(node)) + ndelay
      nchild = child_ptr(node+1) - child_ptr(node)

      nassemble = nassemble + (blkm-1) / SLV_ASSEMBLE_NB + 1
      nasm_sync = nasm_sync + 1
      nassemble2 = nassemble2 + nchild
      do ci = child_ptr(node), child_ptr(node+1)-1
         child = child_list(ci)
         cndelay = nodes(child)%ndelay
         cnelim = nodes(child)%nelim
         cblkm = int(rptr(child+1) - rptr(child)) + cndelay

         nasmblk = nasmblk + (cblkm-cnelim-1) / SLV_ASSEMBLE_NB + 1
      end do

      if(nelim.gt.0) &
         ntrsv = ntrsv + (nelim-1)/SLV_TRSV_NB_TASK + 1

      if(blkm-nelim.gt.0 .and. nelim.ne.0) then
         nx = (blkm-nelim-1)/SLV_GEMV_NX + 1
         ny = (nelim-1)/SLV_GEMV_NY + 1
         ngemv = ngemv + nx*ny
         nreduce = nreduce + nx
      endif

      nscatter = nscatter + (nelim-1)/SLV_SCATTER_NB + 1
   end do
   allocate(assemble_lookup(nassemble), trsv_lookup(ntrsv), gemv_lookup(ngemv),&
      reduce_lookup(nreduce), scatter_lookup(nscatter), &
      assemble_lookup2(nassemble2), asmblkdata(nasmblk), stat=st)
   if(st.ne.0) return

   sz = aligned_size(nassemble*C_SIZEOF(assemble_lookup(1))) + &
      aligned_size(ntrsv*C_SIZEOF(trsv_lookup(1))) + &
      aligned_size(ngemv*C_SIZEOF(gemv_lookup(1))) + &
      aligned_size(nreduce*C_SIZEOF(reduce_lookup(1))) + &
      aligned_size(nscatter*C_SIZEOF(scatter_lookup(1))) + &
      aligned_size(nassemble2*C_SIZEOF(assemble_lookup2(1))) + &
      aligned_size(nasmblk*C_SIZEOF(asmblkdata(1)))
   cuda_error = cudaMalloc(pmem, sz)
   if(cuda_error.ne.0) return

   ! Setup lookups
   nsync = 0
   nlocal = 0
   lwork = 0
   nassemble = 0
   nassemble2 = 0
   nasmblk = 0
   ntrsv = 0
   ngemv = 0
   nreduce = 0
   nscatter = 0
   do ni = 1, nlvl
      node = lvllist(ni)
      ! Setup basic data about node
      ndelay = nodes(node)%ndelay
      nelim = nodes(node)%nelim
      blkn = sptr(node+1) - sptr(node) + ndelay
      blkm = int(rptr(node+1) - rptr(node)) + ndelay
      gldl = blkm

      ! Add contributions
      nx = (blkm-1) / SLV_ASSEMBLE_NB + 1
      nchild = child_ptr(node+1) - child_ptr(node)
      do i = 0, nx-1
         nassemble = nassemble + 1
         assemble_lookup(nassemble)%m = &
            min(SLV_ASSEMBLE_NB, blkm-i*SLV_ASSEMBLE_NB)
         assemble_lookup(nassemble)%xend = &
            max(0, min(SLV_ASSEMBLE_NB, nelim-i*SLV_ASSEMBLE_NB))
         assemble_lookup(nassemble)%list = c_ptr_plus( gpu_rlist_with_delays, &
            rptr_with_delays(node) + i*SLV_ASSEMBLE_NB*C_SIZEOF(dummy_int) )
         assemble_lookup(nassemble)%x_offset = nlocal + i*SLV_ASSEMBLE_NB
         assemble_lookup(nassemble)%contrib_idx = cvmap(node)
         assemble_lookup(nassemble)%contrib_offset = i*SLV_ASSEMBLE_NB-nelim
         assemble_lookup(nassemble)%nchild = nchild
         assemble_lookup(nassemble)%clen = &
            c_ptr_plus( gpu_clen, (child_ptr(node)-1)*C_SIZEOF(dummy_int) )
         assemble_lookup(nassemble)%clists = &
            c_ptr_plus( gpu_clists, (child_ptr(node)-1)*C_SIZEOF(ptdummy_int) )
         assemble_lookup(nassemble)%clists_direct = &
            c_ptr_plus( gpu_clists_direct, (child_ptr(node)-1)*C_SIZEOF(ptdummy_int) )
         assemble_lookup(nassemble)%cvalues_offset = child_ptr(node)-1
      end do

      ! Add contributions (new)
      syncblk = 0
      do ci = child_ptr(node), child_ptr(node+1)-1
         child = child_list(ci)
         cndelay = nodes(child)%ndelay
         cnelim = nodes(child)%nelim
         cblkn = sptr(child+1) - sptr(child) + cndelay
         cblkm = int(rptr(child+1) - rptr(child)) + cndelay

         nassemble2 = nassemble2 + 1
         assemble_lookup2(nassemble2)%m = cblkm-cnelim
         assemble_lookup2(nassemble2)%nelim = nelim
         assemble_lookup2(nassemble2)%list = &
            c_ptr_plus( gpu_clists_direct, (ci-1)*C_SIZEOF(ptdummy_int) )
         if(blkm > nelim) then
            assemble_lookup2(nassemble2)%cvparent = cvmap(node)
         else
            assemble_lookup2(nassemble2)%cvparent = 0 ! Avoid OOB error
         endif
         assemble_lookup2(nassemble2)%cvchild = cvmap(child)
         assemble_lookup2(nassemble2)%sync_offset = ni - 1
         assemble_lookup2(nassemble2)%sync_waitfor = syncblk
         assemble_lookup2(nassemble2)%x_offset = nlocal

         syncblk = syncblk + (cblkm-cnelim-1)/SLV_ASSEMBLE_NB + 1
         do i = 1, (cblkm-cnelim-1)/SLV_ASSEMBLE_NB + 1
            nasmblk = nasmblk + 1
            asmblkdata(nasmblk)%cp = nassemble2-1
            asmblkdata(nasmblk)%blk = i-1
         end do
      end do

      ! Solve with diagonal block
      if(nelim.gt.0) then
         nx = (nelim-1)/SLV_TRSV_NB_TASK + 1
         do i = 0, nx-1
            ntrsv = ntrsv + 1
            trsv_lookup(ntrsv)%n = nelim
            trsv_lookup(ntrsv)%a = nodes(node)%gpu_lcol
            trsv_lookup(ntrsv)%lda = gldl
            trsv_lookup(ntrsv)%x_offset = nlocal
            trsv_lookup(ntrsv)%sync_offset = 2*nsync
         end do
      endif

      ! Update with off-diagonal block
      if(blkm-nelim.gt.0 .and. nelim.ne.0) then
         nx = (blkm-nelim-1)/SLV_GEMV_NX + 1
         ny = (nelim-1)/SLV_GEMV_NY + 1
         do j = 0, ny-1
            do i = 0, nx-1
               ngemv = ngemv + 1
               gemv_lookup(ngemv)%m = &
                  min(SLV_GEMV_NX, (blkm-nelim)-i*SLV_GEMV_NX)
               gemv_lookup(ngemv)%n = min(SLV_GEMV_NY, nelim-j*SLV_GEMV_NY)
               gemv_lookup(ngemv)%a = c_ptr_plus( nodes(node)%gpu_lcol, &
                  nelim*C_SIZEOF(dummy_real) + &
                  (i*SLV_GEMV_NX + j*SLV_GEMV_NY*gldl)*C_SIZEOF(dummy_real) )
               gemv_lookup(ngemv)%lda = gldl
               gemv_lookup(ngemv)%x_offset = nlocal+j*SLV_GEMV_NY
               gemv_lookup(ngemv)%y_offset = lwork+(blkm-nelim)*j+SLV_GEMV_NX*i
            end do
         end do
         do i = 0, nx-1
            nreduce = nreduce + 1
            reduce_lookup(nreduce)%m = min(SLV_GEMV_NX,blkm-nelim-i*SLV_GEMV_NX)
            reduce_lookup(nreduce)%n = ny
            reduce_lookup(nreduce)%src_offset = lwork + i*SLV_GEMV_NX
            reduce_lookup(nreduce)%ldsrc = blkm-nelim
            reduce_lookup(nreduce)%dest_idx = cvmap(node)
            reduce_lookup(nreduce)%dest_offset = i*SLV_GEMV_NX
         end do
         lwork = lwork + (blkm-nelim)*ny
      end if

      ! Copy eliminated (and delayed) variables back to x
      nx = (nelim-1)/SLV_SCATTER_NB + 1
      j = nscatter+1
      do i=0, nx-1
         nscatter = nscatter + 1
         scatter_lookup(nscatter)%n = min(SLV_SCATTER_NB,nelim-i*SLV_SCATTER_NB)
         scatter_lookup(nscatter)%src_offset = nlocal+i*SLV_SCATTER_NB
         scatter_lookup(nscatter)%index = c_ptr_plus( gpu_rlist_with_delays, &
            rptr_with_delays(node) + i*SLV_SCATTER_NB*C_SIZEOF(dummy_int) )
         scatter_lookup(nscatter)%dest_offset = 0
      end do

      nsync = nsync + 1
      nlocal = nlocal + nelim
   end do

   gpul%nassemble = nassemble
   sz = gpul%nassemble*C_SIZEOF(assemble_lookup(1))
   gpul%assemble = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%assemble, C_LOC(assemble_lookup), sz, &
      stream)
   if(cuda_error.ne.0) return

   gpul%ntrsv = ntrsv
   sz = gpul%ntrsv*C_SIZEOF(trsv_lookup(1))
   gpul%trsv = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%trsv, C_LOC(trsv_lookup), sz, stream)
   if(cuda_error.ne.0) return

   gpul%ngemv = ngemv
   sz = gpul%ngemv*C_SIZEOF(gemv_lookup(1))
   gpul%gemv = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%gemv, C_LOC(gemv_lookup), sz, stream)
   if(cuda_error.ne.0) return

   gpul%nreduce = nreduce
   sz = gpul%nreduce*C_SIZEOF(reduce_lookup(1))
   gpul%reduce = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%reduce, C_LOC(reduce_lookup), sz, &
      stream)
   if(cuda_error.ne.0) return

   gpul%nscatter = nscatter
   sz = gpul%nscatter*C_SIZEOF(scatter_lookup(1))
   gpul%scatter = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%scatter, C_LOC(scatter_lookup), sz, &
      stream)
   if(cuda_error.ne.0) return

   gpul%nassemble2 = nassemble2
   gpul%nasm_sync = nasm_sync
   sz = gpul%nassemble2*C_SIZEOF(assemble_lookup2(1))
   gpul%assemble2 = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%assemble2, C_LOC(assemble_lookup2), &
      sz, stream)
   if(cuda_error.ne.0) return

   gpul%nasmblk = nasmblk
   sz = gpul%nasmblk*C_SIZEOF(asmblkdata(1))
   gpul%asmblk = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%asmblk, C_LOC(asmblkdata), &
      sz, stream)
   if(cuda_error.ne.0) return
end subroutine create_gpu_lookup_fwd

subroutine create_gpu_lookup_bwd(nlvl, lvllist, nodes, sptr, rptr,         &
      rptr_with_delays, gpu_rlist_with_delays, gpul, lwork, nsync, stream, &
      st, cuda_error)
   integer, intent(in) :: nlvl
   integer, dimension(*), intent(in) :: lvllist
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer(C_SIZE_T), dimension(*), intent(in) :: rptr_with_delays
   type(C_PTR), intent(in) :: gpu_rlist_with_delays
   type(lookups_gpu_bwd), intent(out) :: gpul
   integer, intent(out) :: lwork
   integer, intent(out) :: nsync
   type(C_PTR), intent(in) :: stream
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error

   integer :: i, j, xi, yi, li
   integer :: node, nelim, nd, blkm, blkn
   integer :: nx, ny, nrds, ntrsv, nscatter
   integer :: blk_gemv, blk_rds, blk_trsv, blk_scatter
   integer :: woffset, lupd
   type(gemv_transpose_lookup), dimension(:), allocatable, target :: gemv_lookup
   type(reducing_d_solve_lookup), dimension(:), allocatable, target :: &
      rds_lookup
   type(trsv_lookup_type), dimension(:), allocatable, target :: trsv_lookup
   type(scatter_lookup_type), dimension(:), allocatable, target :: &
      scatter_lookup

   integer(C_INT) :: dummy_int
   real(C_DOUBLE) :: dummy_real

   type(C_PTR) :: pmem
   integer(C_SIZE_T) :: sz

   st = 0; cuda_error = 0

   blk_gemv = 0
   blk_rds = 0
   blk_trsv = 0
   blk_scatter = 0
   lwork = 0
   nsync = 0
   do li = 1, nlvl
      node = lvllist(li)
      nsync = nsync + 1

      nelim = nodes(node)%nelim
      nd = nodes(node)%ndelay
      blkm = int(rptr(node+1) - rptr(node)) + nd

      nx = (blkm-nelim-1)/SLV_TRSM_TR_NBX + 1
      ny = (nelim-1)/SLV_TRSM_TR_NBY + 1
      blk_gemv = blk_gemv + nx*ny
      
      nrds = (nelim-1)/SLV_REDUCING_D_SOLVE_THREADS_PER_BLOCK + 1
      blk_rds = blk_rds + nrds
      lwork = lwork + nx*ny*SLV_TRSM_TR_NBY

      ntrsv = (nelim-1)/SLV_TRSV_NB_TASK + 1
      blk_trsv = blk_trsv + ntrsv

      nscatter = (nelim-1)/SLV_SCATTER_NB + 1
      blk_scatter = blk_scatter + nscatter
   end do
   allocate(gemv_lookup(blk_gemv), rds_lookup(blk_rds), trsv_lookup(blk_trsv), &
      scatter_lookup(blk_scatter), stat=st)
   if(st.ne.0) return

   sz = aligned_size(blk_gemv*C_SIZEOF(gemv_lookup(1))) + &
      aligned_size(blk_rds*C_SIZEOF(rds_lookup(1))) + &
      aligned_size(blk_trsv*C_SIZEOF(trsv_lookup(1))) + &
      aligned_size(blk_scatter*C_SIZEOF(scatter_lookup(1)))
   cuda_error = cudaMalloc(pmem, sz)
   if(cuda_error.ne.0) return

   woffset = 0
   blk_gemv = 1
   blk_rds = 1
   blk_trsv = 1
   blk_scatter = 1
   j = 0 ! Count of nodes at this level used for determing sync
   do li = 1, nlvl
      node = lvllist(li)
      nelim = nodes(node)%nelim
      if(nelim.eq.0) cycle
      nd = nodes(node)%ndelay
      blkn = sptr(node+1) - sptr(node) + nd
      blkm = int(rptr(node+1) - rptr(node)) + nd

      ! Setup tasks for gemv
      nx = (blkm-nelim-1)/SLV_TRSM_TR_NBX + 1
      ny = (nelim-1)/SLV_TRSM_TR_NBY + 1
      lupd = ny * SLV_TRSM_TR_NBY
      do yi = 0, ny-1
         do xi = 0, nx-1
            gemv_lookup(blk_gemv)%m = &
               min(SLV_TRSM_TR_NBX, blkm-nelim-xi*SLV_TRSM_TR_NBX)
            gemv_lookup(blk_gemv)%n = &
               min(SLV_TRSM_TR_NBY, nelim-yi*SLV_TRSM_TR_NBY)
            gemv_lookup(blk_gemv)%lda = blkm
            gemv_lookup(blk_gemv)%a = &
               c_ptr_plus(nodes(node)%gpu_lcol, &
               C_SIZEOF(dummy_real)*( & 
                  gemv_lookup(blk_gemv)%lda*yi*SLV_TRSM_TR_NBY + &
                  nelim+xi*SLV_TRSM_TR_NBX) &
               )
            gemv_lookup(blk_gemv)%rlist = &
               c_ptr_plus( gpu_rlist_with_delays, rptr_with_delays(node) + &
               C_SIZEOF(dummy_int)*(nelim+xi*SLV_TRSM_TR_NBX) )
            gemv_lookup(blk_gemv)%yoffset = &
               woffset + xi*lupd + yi*SLV_TRSM_TR_NBY;
            blk_gemv = blk_gemv + 1
         end do
      end do

      ! Setup tasks for reducing_d_solve()
      nrds = (nelim-1)/SLV_REDUCING_D_SOLVE_THREADS_PER_BLOCK + 1
      do xi = 0, nrds-1
         rds_lookup(blk_rds)%first_idx = &
            xi*SLV_REDUCING_D_SOLVE_THREADS_PER_BLOCK
         rds_lookup(blk_rds)%m = min(SLV_REDUCING_D_SOLVE_THREADS_PER_BLOCK, &
            nelim-xi*SLV_REDUCING_D_SOLVE_THREADS_PER_BLOCK)
         rds_lookup(blk_rds)%n = nx
         rds_lookup(blk_rds)%updoffset = woffset
         rds_lookup(blk_rds)%d = c_ptr_plus(nodes(node)%gpu_lcol, &
            blkm*blkn*C_SIZEOF(dummy_real) )
         rds_lookup(blk_rds)%perm = &
            c_ptr_plus(gpu_rlist_with_delays, rptr_with_delays(node))
         rds_lookup(blk_rds)%ldupd = lupd
         blk_rds = blk_rds+1
      end do

      ! Setup tasks for trsv
      ntrsv = (nelim-1)/SLV_TRSV_NB_TASK + 1
      do i = 1, ntrsv
         trsv_lookup(blk_trsv)%n = nelim
         trsv_lookup(blk_trsv)%a = nodes(node)%gpu_lcol
         trsv_lookup(blk_trsv)%lda = blkm
         trsv_lookup(blk_trsv)%x_offset = woffset
         trsv_lookup(blk_trsv)%sync_offset = 2*j
         blk_trsv = blk_trsv + 1
      end do
      j = j + 1

      nscatter = (nelim-1)/SLV_SCATTER_NB + 1
      do i=0, nscatter-1
         scatter_lookup(blk_scatter)%n = &
            min(SLV_SCATTER_NB, nelim-i*SLV_SCATTER_NB)
         scatter_lookup(blk_scatter)%src_offset = &
            woffset+i*SLV_SCATTER_NB
         scatter_lookup(blk_scatter)%index = c_ptr_plus(gpu_rlist_with_delays, &
            rptr_with_delays(node) + i*SLV_SCATTER_NB*C_SIZEOF(dummy_int))
         scatter_lookup(blk_scatter)%dest_offset = 0
         blk_scatter = blk_scatter + 1
      end do

      woffset = woffset + nx*ny*SLV_TRSM_TR_NBY
   end do

   gpul%ngemv = blk_gemv-1
   gpul%nrds = blk_rds-1
   gpul%ntrsv = blk_trsv-1
   gpul%nscatter = blk_scatter-1

   sz = gpul%ngemv*C_SIZEOF(gemv_lookup(1))
   gpul%gemv = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%gemv, C_LOC(gemv_lookup), sz, stream)
   if(cuda_error.ne.0) return

   sz = gpul%nrds*C_SIZEOF(rds_lookup(1))
   gpul%rds = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%rds, C_LOC(rds_lookup), sz, stream)
   if(cuda_error.ne.0) return

   sz = gpul%ntrsv*C_SIZEOF(trsv_lookup(1))
   gpul%trsv = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%trsv, C_LOC(trsv_lookup), sz, stream)
   if(cuda_error.ne.0) return

   sz = gpul%nscatter*C_SIZEOF(scatter_lookup(1))
   gpul%scatter = pmem
   pmem = c_ptr_plus_aligned(pmem, sz)
   cuda_error = cudaMemcpyAsync_h2d(gpul%scatter, C_LOC(scatter_lookup), sz, &
      stream)
   if(cuda_error.ne.0) return
end subroutine create_gpu_lookup_bwd

subroutine setup_gpu_solve(n, child_ptr, child_list, nnodes, nodes, sparent,  &
      sptr, rptr, rlist, nstream, stream_handle, stream_data, top_data,       &
      gpu_rlist_with_delays, gpu_clists, gpu_clists_direct, gpu_clen, st,     &
      cuda_error, gpu_rlist_direct_with_delays)
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(nnodes), intent(in) :: sparent
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), target, intent(in) :: rlist
   integer, intent(in) :: nstream
   type(C_PTR), dimension(*), intent(in) :: stream_handle
   type(gpu_type), dimension(*), intent(inout) :: stream_data
   type(gpu_type), intent(inout) :: top_data
   type(C_PTR), intent(out) :: gpu_rlist_with_delays
   type(C_PTR), intent(out) :: gpu_clists
   type(C_PTR), intent(out) :: gpu_clists_direct
   type(C_PTR), intent(out) :: gpu_clen
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error
   type(C_PTR), intent(out) :: gpu_rlist_direct_with_delays

   integer :: node, nd, nelim
   integer(long) :: i, k
   integer(long) :: blkm, blkn, ip

   integer(C_INT) :: dummy_int ! used for size of element of rlist
   integer, dimension(:), allocatable, target :: rlist2, rlist_direct2
   integer, dimension(:), pointer :: lperm

   integer :: stream
   integer(long) :: int_data

   integer :: cn, cnelim, cnode, parent
   integer, dimension(:), allocatable :: cvmap, rmap
   integer(C_INT), dimension(:), allocatable, target :: clen
   type(C_PTR), dimension(:), allocatable, target :: clists
   type(C_PTR), dimension(:), allocatable, target :: clists_direct
   integer(C_SIZE_T), dimension(:), allocatable :: rptr_with_delays

   st = 0; cuda_error = 0

   !
   ! Copy rlist, but apply invp and delays while doing so
   !
   ! Count space
   allocate(rptr_with_delays(nnodes+1), stat=st)
   if(st.ne.0) return
   rptr_with_delays(1) = 0
   do node = 1, nnodes
      nd = nodes(node)%ndelay
      blkm = int(rptr(node+1) - rptr(node)) + nd
      rptr_with_delays(node+1) = rptr_with_delays(node) + &
         C_SIZEOF(dummy_int)*blkm
   end do
   cuda_error = cudaMalloc(gpu_rlist_with_delays, rptr_with_delays(nnodes+1))
   if(cuda_error.ne.0) return
   int_data = rptr_with_delays(nnodes+1)
   allocate(rlist2(rptr_with_delays(nnodes+1)/C_SIZEOF(dummy_int)), stat=st)
   if(st.ne.0) return
   ip = 0
   do node = 1, nnodes
      nelim = nodes(node)%nelim
      nd = nodes(node)%ndelay
      blkn = sptr(node+1) - sptr(node) + nd
      blkm = int(rptr(node+1) - rptr(node)) + nd
      lperm => nodes(node)%perm
      do i = 1, blkn
         rlist2(ip+i) = lperm(i) - 1
      end do
      k = rptr(node)+blkn-nd
      do i = blkn+1, blkm
         rlist2(ip+i) = rlist(k) - 1
         k = k + 1
      end do
      ip = ip + blkm
   end do
   cuda_error = cudaMemcpy_h2d(gpu_rlist_with_delays, C_LOC(rlist2), &
      rptr_with_delays(nnodes+1))
   if(cuda_error.ne.0) return

   !
   ! Setup rlist_direct_with_delays
   !
   allocate(rlist_direct2(rptr_with_delays(nnodes+1)/C_SIZEOF(dummy_int)), &
      rmap(0:n-1), stat=st)
   if(st.ne.0) return
   cuda_error = cudaMalloc(gpu_rlist_direct_with_delays, &
      rptr_with_delays(nnodes+1))
   if(cuda_error.ne.0) return
   do node = 1, nnodes
      nelim = nodes(node)%nelim
      parent = sparent(node)
      if(parent<0 .or. parent>nnodes) cycle ! root
      ! drop parent locs into map
      do i = rptr_with_delays(parent)/C_SIZEOF(dummy_int), rptr_with_delays(parent+1)/C_SIZEOF(dummy_int) -1
         rmap(rlist2(i+1)) = int(i - rptr_with_delays(parent)/C_SIZEOF(dummy_int))
      end do
      ! build rlist_direct2
      do i = rptr_with_delays(node)/C_SIZEOF(dummy_int)+nelim, rptr_with_delays(node+1)/C_SIZEOF(dummy_int) -1
         rlist_direct2(i+1) = rmap(rlist2(i+1))
      end do
   end do
   cuda_error = cudaMemcpy_h2d(gpu_rlist_direct_with_delays, &
      C_LOC(rlist_direct2), rptr_with_delays(nnodes+1))
   if(cuda_error.ne.0) return

   !
   ! Setup solve info
   !
   allocate(clists(nnodes), clists_direct(nnodes), cvmap(nnodes), &
      clen(nnodes), stat=st)
   if(st.ne.0) return
   cuda_error = cudaMalloc(gpu_clists, nnodes*C_SIZEOF(clists(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMalloc(gpu_clists_direct, nnodes*C_SIZEOF(clists_direct(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMalloc(gpu_clen, nnodes*C_SIZEOF(dummy_int))
   if(cuda_error.ne.0) return
   do node = 1, nnodes
      do cn = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(cn)
         cnelim = nodes(cnode)%nelim

         cvmap(cnode) = cn-1
         clists(cn) = c_ptr_plus( gpu_rlist_with_delays, &
            rptr_with_delays(cnode) + cnelim*C_SIZEOF(dummy_int) )
         clists_direct(cn) = c_ptr_plus( gpu_rlist_direct_with_delays, &
            rptr_with_delays(cnode) + cnelim*C_SIZEOF(dummy_int) )
         clen(cn) = int((rptr_with_delays(cnode+1) - rptr_with_delays(cnode)) /&
            C_SIZEOF(dummy_int)) - cnelim
      end do
   end do
   cuda_error = cudaMemcpy_h2d(gpu_clists, C_LOC(clists), &
      nnodes*C_SIZEOF(clists(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpy_h2d(gpu_clists_direct, C_LOC(clists_direct), &
      nnodes*C_SIZEOF(clists_direct(1)))
   if(cuda_error.ne.0) return
   cuda_error = cudaMemcpy_h2d(gpu_clen, C_LOC(clen), nnodes*C_SIZEOF(clen(1)))
   if(cuda_error.ne.0) return

   !
   ! Setup solve lookups
   !
   do stream = 1, nstream
      call setup_bwd_slv(nodes, sptr, rptr, stream_data(stream)%num_levels,    &
         stream_data(stream)%lvlptr, stream_data(stream)%lvllist,              &
         rptr_with_delays, gpu_rlist_with_delays,                              &
         stream_data(stream)%bwd_slv_lookup, stream_data(stream)%bwd_slv_lwork,&
         stream_data(stream)%bwd_slv_nsync, st, cuda_error,                    &
         stream_handle(stream))
      if(st.ne.0) return
      if(cuda_error.ne.0) return
      call setup_fwd_slv(child_ptr, child_list, nodes, sptr, rptr,             &
         stream_data(stream)%num_levels, stream_data(stream)%lvlptr,           &
         stream_data(stream)%lvllist, rptr_with_delays, gpu_rlist_with_delays, &
         gpu_clen, gpu_clists, gpu_clists_direct, cvmap,                       &
         stream_data(stream)%fwd_slv_lookup,                                   &
         stream_data(stream)%fwd_slv_lwork, stream_data(stream)%fwd_slv_nlocal,&
         stream_data(stream)%fwd_slv_nsync, stream_data(stream)%fwd_slv_nasync,&
         st, cuda_error, stream_handle(stream))
      if(st.ne.0) return
      if(cuda_error.ne.0) return
   end do
   call setup_bwd_slv(nodes, sptr, rptr, top_data%num_levels, top_data%lvlptr, &
      top_data%lvllist, rptr_with_delays, gpu_rlist_with_delays,               &
      top_data%bwd_slv_lookup, top_data%bwd_slv_lwork, top_data%bwd_slv_nsync, &
      st, cuda_error, stream_handle(1))
   if(st.ne.0) return
   if(cuda_error.ne.0) return
   call setup_fwd_slv(child_ptr, child_list, nodes, sptr, rptr,                &
      top_data%num_levels, top_data%lvlptr, top_data%lvllist, rptr_with_delays,&
      gpu_rlist_with_delays, gpu_clen, gpu_clists, gpu_clists_direct, cvmap,   &
      top_data%fwd_slv_lookup, top_data%fwd_slv_lwork, top_data%fwd_slv_nlocal,&
      top_data%fwd_slv_nsync, top_data%fwd_slv_nasync, st, cuda_error,         &
      stream_handle(1))
   if(st.ne.0) return
   if(cuda_error.ne.0) return
end subroutine setup_gpu_solve

subroutine setup_fwd_slv(child_ptr, child_list, nodes, sptr, rptr, num_levels,&
      lvlptr, lvllist, rptr_with_delays, gpu_rlist_with_delays, gpu_clen,     &
      gpu_clists, gpu_clists_direct, cvmap, fwd_slv_lookup, lwork, nlocal,    &
      nsync, nasm_sync, st, cuda_error, stream)
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer, intent(in) :: num_levels
   integer, dimension(*), intent(in) :: lvlptr
   integer, dimension(*), intent(in) :: lvllist
   integer(C_SIZE_T), dimension(*), intent(in) :: rptr_with_delays
   type(C_PTR), intent(in) :: gpu_rlist_with_delays
   type(C_PTR), intent(in) :: gpu_clen
   type(C_PTR), intent(in) :: gpu_clists
   type(C_PTR), intent(in) :: gpu_clists_direct
   integer, dimension(*), intent(in) :: cvmap
   type(lookups_gpu_fwd), dimension(:), allocatable, intent(out) :: &
      fwd_slv_lookup
   integer, intent(out) :: lwork
   integer, intent(out) :: nlocal
   integer, intent(out) :: nsync
   integer, intent(out) :: nasm_sync
   integer, intent(out) :: st
   integer, intent(out) :: cuda_error
   type(C_PTR), intent(in) :: stream

   integer :: lvl, i, j, k, p

   st = 0
   cuda_error = 0

   allocate(fwd_slv_lookup(num_levels), stat=st)
   if(st.ne.0) return

   nsync = 0
   nlocal = 0
   lwork = 0
   nasm_sync = 0
   do lvl = 1, num_levels
      call create_gpu_lookup_fwd(lvlptr(lvl+1)-lvlptr(lvl),                   &
         lvllist(lvlptr(lvl):lvlptr(lvl+1)-1), nodes, child_ptr, child_list,  &
         cvmap, sptr, rptr, rptr_with_delays, gpu_rlist_with_delays, gpu_clen,&
         gpu_clists, gpu_clists_direct, fwd_slv_lookup(lvl), i, j, k, p,      &
         stream, st, cuda_error)
      if(st.ne.0) return
      if(cuda_error.ne.0) return
      nsync = max(nsync, i)
      nlocal = max(nlocal, j)
      lwork = max(lwork, k)
      nasm_sync = max(nasm_sync, p)
   end do
end subroutine setup_fwd_slv

subroutine setup_bwd_slv(nodes, sptr, rptr, num_levels, lvlptr, lvllist,      &
      rptr_with_delays, gpu_rlist_with_delays, bwd_slv_lookup, lwork, nsync,  &
      st, cuda_error, stream)
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer, intent(in) :: num_levels
   integer, dimension(*), intent(in) :: lvlptr
   integer, dimension(*), intent(in) :: lvllist
   integer(C_SIZE_T), dimension(*), intent(in) :: rptr_with_delays
   type(C_PTR), intent(in) :: gpu_rlist_with_delays
   type(lookups_gpu_bwd), dimension(:), allocatable, intent(out) :: &
      bwd_slv_lookup
   integer, intent(out) :: lwork
   integer, intent(out) :: nsync
   integer, intent(out) :: cuda_error
   integer, intent(out) :: st
   type(C_PTR), intent(in) :: stream

   integer :: i, j, lvl

   st = 0

   allocate(bwd_slv_lookup(num_levels), stat=st)
   if(st.ne.0) return

   lwork = 0
   nsync = 0
   do lvl = 1, num_levels
      call create_gpu_lookup_bwd(lvlptr(lvl+1)-lvlptr(lvl),             &
         lvllist(lvlptr(lvl):lvlptr(lvl+1)-1), nodes, sptr, rptr,       &
         rptr_with_delays, gpu_rlist_with_delays, bwd_slv_lookup(lvl),  &
         i, j, stream, cuda_error, st)
      if(cuda_error.ne.0) return
      if(st.ne.0) return
      lwork = max(lwork, i)
      nsync = max(nsync, j)
   end do
end subroutine setup_bwd_slv

subroutine fwd_multisolve_gpu( nnodes, nodes, rptr, &
      nstream, stream_handle, stream_data, top_data, &
      nrhs, gpu_x, cuda_error, st )
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(in) :: nodes
   integer(long), dimension(nnodes + 1), intent(in) :: rptr
   integer, intent(in) :: nstream
   type(C_PTR), intent(in) :: stream_handle(*)
   type(gpu_type), dimension(*), intent(in) :: stream_data
   type(gpu_type), intent(in) :: top_data
   integer, intent(in) :: nrhs
   type(C_PTR) :: gpu_x
   integer, intent(out) :: cuda_error
   integer, intent(out) :: st
   
   integer :: stream
   integer :: node
   type(C_PTR), allocatable :: gpu_contrib(:)
   
   allocate ( gpu_contrib(nnodes), stat = st )
   if ( st /= 0 ) return
   gpu_contrib = C_NULL_PTR

   do stream = 1, nstream
      call subtree_fwd_multisolve_gpu &
         ( nnodes, nodes, rptr, &
         stream_handle(stream), stream_data(stream), &
         nrhs, gpu_x, gpu_contrib, &
         cuda_error )
      if ( cuda_error .ne. 0 ) return
   end do 

   cuda_error = cudaDeviceSynchronize() ! Wait for streams to finish
   if(cuda_error.ne.0) return

   call subtree_fwd_multisolve_gpu( nnodes, nodes, rptr, &
      stream_handle(1), top_data, nrhs, gpu_x, gpu_contrib, cuda_error )

   do node = 1, nnodes
      if ( C_ASSOCIATED(gpu_contrib(node)) ) then
         cuda_error = cudaFree(gpu_contrib(node))
      end if
   end do

end subroutine fwd_multisolve_gpu

subroutine subtree_fwd_multisolve_gpu(nnodes, nodes, rptr, stream_handle, &
      fact_data, nrhs, gpu_x, gpu_contrib, cuda_error)
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(in) :: nodes
   integer(long), dimension(nnodes + 1), intent(in) :: rptr
   type(C_PTR), intent(in) :: stream_handle
   type(gpu_type), intent(in) :: fact_data
   integer, intent(in) :: nrhs
   type(C_PTR) :: gpu_x
   type(C_PTR) :: gpu_contrib(nnodes)
   integer, intent(out) :: cuda_error
   
   integer :: n
   integer :: lev
   integer :: node
   integer :: nelim
   integer :: nrows
   integer :: offc
   integer :: lx_size
   integer :: lc_size
   integer :: max_lx_size
   integer :: max_lc_size
   integer :: offi
   integer :: ncp ! Number of child-parent pairs
   integer :: nblk ! Number of blocks
   integer :: ncb

   integer :: i

   integer(long) :: rd_size
   integer(long) :: upitch, vpitch
   integer(long) :: ncols

   type(C_PTR) :: gpu_cpdata ! Ouput child-parent info
   type(C_PTR) :: gpu_blkdata ! Output block-by-block info
  
   integer(C_SIZE_T) :: lgpu_work

   type(C_PTR) :: gpu_lx1, gpu_lx2
   type(C_PTR) :: gpu_lc1, gpu_lc2
   type(C_PTR) :: gpu_lcp
   type(C_PTR) :: gpu_lcc
   type(C_PTR) :: gpu_ind
   type(C_PTR) :: gpu_sync
   type(C_PTR) :: gpu_u
   type(C_PTR) :: gpu_v

   type(cuda_stack_alloc_type) :: gwork
   real(wp) :: dummy_real
   integer(C_INT) :: dummy_int
   
   n = fact_data%n

   rd_size = fact_data%rd_size

   gpu_ind = fact_data%gpu_col_ind

   max_lx_size = fact_data%max_lx_size*nrhs
   max_lc_size = fact_data%max_lc_size*nrhs
   lgpu_work = 2*max_lx_size*C_SIZEOF(dummy_real) + &
               2*max_lc_size*C_SIZEOF(dummy_real) + &
               4*256 ! To allow for alignment
   call custack_init(gwork, lgpu_work, cuda_error)
   if(cuda_error.ne.0) return

   gpu_lx1 = custack_alloc(gwork, max_lx_size*C_SIZEOF(dummy_real))
   gpu_lx2 = custack_alloc(gwork, max_lx_size*C_SIZEOF(dummy_real))
   gpu_lc1 = custack_alloc(gwork, max_lc_size*C_SIZEOF(dummy_real))
   gpu_lc2 = custack_alloc(gwork, max_lc_size*C_SIZEOF(dummy_real))
   gpu_lcp = gpu_lc2

   do lev = 1, fact_data%num_levels

      gpu_lcc = gpu_lcp
      if( mod(lev,2) .gt. 0 ) then
         gpu_lcp = gpu_lc1
      else
         gpu_lcp = gpu_lc2
      end if

      lx_size   = fact_data%values_L(lev)%lx_size
      lc_size   = fact_data%values_L(lev)%lc_size
      offi      = fact_data%values_L(lev)%off_col_ind

      gpu_u = c_ptr_plus(gpu_ind, offi*C_SIZEOF(dummy_int))
      call gather( stream_handle, lx_size, nrhs, gpu_x, n, &
         gpu_lx1, lx_size, gpu_u )

      if ( fact_data%values_L(lev)%total_nch > 0 ) then
         do i = 1, fact_data%values_L(lev)%nimp
            node = fact_data%values_L(lev)%import(i)
            nelim = nodes(node)%nelim
            offc = fact_data%off_lc(node)
            nrows = int(rptr(node + 1) - rptr(node)) + nodes(node)%ndelay
            gpu_u = gpu_contrib(node)
            gpu_v = c_ptr_plus(gpu_lcc, offc*C_SIZEOF(dummy_real))
            upitch = (nrows - nelim)*C_SIZEOF(dummy_real)
            vpitch = fact_data%values_L(lev)%lcc_size*C_SIZEOF(dummy_real)
            ncols = nrhs
            cuda_error = cudaMemcpy2DAsync(gpu_v, vpitch, gpu_u, upitch, &
               upitch, ncols, cudaMemcpyDeviceToDevice, stream_handle)
            if ( cuda_error /= 0 ) return
         end do
         ncp         = fact_data%values_L(lev)%ncp_pre
         gpu_cpdata  = fact_data%values_L(lev)%gpu_cpdata_pre
         nblk        = fact_data%values_L(lev)%ncb_asm_pre
         gpu_blkdata = fact_data%values_L(lev)%gpu_blkdata_pre
         gpu_sync    = fact_data%gpu_sync
         call assemble_solve_phase( stream_handle, nrhs, &
            nblk, gpu_blkdata, ncp, &
            gpu_cpdata, gpu_lx1, gpu_lcc, gpu_sync )
      end if

      cuda_error = cudaMemset(gpu_lcp, 0, lc_size*nrhs*C_SIZEOF(dummy_real))
      if ( cuda_error .ne. 0 ) return
      ncb = fact_data%values_L(lev)%ncb_slv_n
      if ( ncb .gt. 0 ) then
         call multinode_solve_n( stream_handle, ncb, nrhs, &
            fact_data%values_L(lev)%ptr_levL, gpu_lx1, gpu_lx2, gpu_lcp, &
            fact_data%values_L(lev)%gpu_solve_n_data )
         gpu_v = c_ptr_plus(gpu_ind, offi*C_SIZEOF(dummy_int))
         call scatter( stream_handle, lx_size, nrhs, &
            gpu_lx2, lx_size, gpu_x, n, gpu_v )
      end if

      ncp = fact_data%values_L(lev)%ncp_post
      if ( ncp .gt. 0 ) then
         gpu_cpdata  = fact_data%values_L(lev)%gpu_cpdata_post
         nblk        = fact_data%values_L(lev)%ncb_asm_post
         gpu_blkdata = fact_data%values_L(lev)%gpu_blkdata_post
         gpu_sync    = fact_data%gpu_sync
         call assemble_solve_phase( stream_handle, nrhs, &
            nblk, gpu_blkdata, ncp, &
            gpu_cpdata, gpu_lcp, gpu_lcc, gpu_sync )
      end if

      vpitch = lc_size*C_SIZEOF(dummy_real)
      do i = 1, fact_data%values_L(lev)%nexp
         node = fact_data%values_L(lev)%export(i)
         nelim = nodes(node)%nelim
         offc = fact_data%off_lc(node)
         nrows = int(rptr(node + 1) - rptr(node)) + nodes(node)%ndelay
         gpu_v = c_ptr_plus(gpu_lcp, offc*C_SIZEOF(dummy_real))
         upitch = (nrows - nelim)*C_SIZEOF(dummy_real)
         cuda_error = cudaMalloc(gpu_u, upitch*nrhs)
         if ( cuda_error /= 0 ) return
         gpu_contrib(node) = gpu_u
         ncols = nrhs
         cuda_error = cudaMemcpy2DAsync(gpu_u, upitch, gpu_v, vpitch, &
            upitch, ncols, cudaMemcpyDeviceToDevice, stream_handle)
         if ( cuda_error /= 0 ) return
      end do

   end do

   call custack_free(gwork, max_lc_size*C_SIZEOF(dummy_real)) ! gpu_lc2
   call custack_free(gwork, max_lc_size*C_SIZEOF(dummy_real)) ! gpu_lc1
   call custack_free(gwork, max_lx_size*C_SIZEOF(dummy_real)) ! gpu_lx2
   call custack_free(gwork, max_lx_size*C_SIZEOF(dummy_real)) ! gpu_lx1

   call custack_finalize(gwork, cuda_error)
   if(cuda_error.ne.0) return

end subroutine subtree_fwd_multisolve_gpu

subroutine bwd_multisolve_gpu(pos_def, job, nstream, stream_handle, &
      stream_data, top_data, nrhs, gpu_x, cuda_error)
   logical, intent(in) :: pos_def
   integer, intent(in) :: job
   integer, intent(in) :: nstream
   type(C_PTR), dimension(*), intent(in) :: stream_handle
   type(gpu_type), dimension(*), intent(in) :: stream_data
   type(gpu_type), intent(in) :: top_data
   integer, intent(in) :: nrhs
   type(C_PTR) :: gpu_x
   integer, intent(out) :: cuda_error
   
   integer :: stream

   call subtree_bwd_multisolve_gpu( pos_def, job, stream_handle(1), top_data, &
      nrhs, gpu_x, cuda_error)
   if(cuda_error.ne.0) return

   do stream = 1, nstream
      call subtree_bwd_multisolve_gpu(pos_def, job, stream_handle(stream), &
         stream_data(stream), nrhs, gpu_x, cuda_error)
      if(cuda_error.ne.0) return
   end do 

   cuda_error = cudaDeviceSynchronize() ! Wait for streams to finish
   if(cuda_error.ne.0) return

end subroutine bwd_multisolve_gpu

subroutine subtree_bwd_multisolve_gpu(pos_def, job, stream_handle, fact_data, &
      nrhs, gpu_x, cuda_error)
   logical, intent(in) :: pos_def
   integer, intent(in) :: job
   type(C_PTR), intent(in) :: stream_handle
   type(gpu_type), intent(in) :: fact_data
   integer, intent(in) :: nrhs
   type(C_PTR) :: gpu_x   
   integer, intent(out) :: cuda_error
   
   integer :: n
   integer :: lev
   integer :: ncb
   
   integer(long) :: sz
   
   integer :: ln_len, lx_len, ldln, ldlx

   type(C_PTR) :: gpu_ln
   type(C_PTR) :: gpu_ind
   type(C_PTR) :: gpu_indx

   type(C_PTR) :: gpu_u
   type(C_PTR) :: gpu_v

   type(cuda_stack_alloc_type) :: gwork
   integer(C_SIZE_T) :: lgpu_work
   real(wp) :: dummy_real
   integer(C_INT) :: dummy_int
   
   n = fact_data%n
   
   ln_len = nrhs*(fact_data%max_lx_size + fact_data%max_lc_size)
   lx_len = nrhs*fact_data%max_lx_size
   lgpu_work = ln_len*C_SIZEOF(dummy_real) + &
      2*lx_len*C_SIZEOF(dummy_real) + 3*256 ! Allow for alignment on 3 ptrs
   call custack_init(gwork, lgpu_work, cuda_error)
   if(cuda_error.ne.0) return
   gpu_ln = custack_alloc(gwork, ln_len*C_SIZEOF(dummy_real))
   gpu_u = custack_alloc(gwork, lx_len*C_SIZEOF(dummy_real))
   gpu_v = custack_alloc(gwork, lx_len*C_SIZEOF(dummy_real))
   
   ! Backwards solve DL^Tx = z or L^Tx = z or Dx = z
   do lev = fact_data%num_levels, 1, -1

      ldln = fact_data%values_L(lev)%ln_size

      ldlx = fact_data%values_L(lev)%lx_size

      if ( job == SSIDS_SOLVE_JOB_DIAG ) then
         sz = fact_data%values_L(lev)%off_col_ind * C_SIZEOF(dummy_int)
         gpu_indx = c_ptr_plus(fact_data%gpu_col_ind, sz)
         gpu_ind = gpu_indx
      else
         sz = fact_data%values_L(lev)%off_row_ind * C_SIZEOF(dummy_int)
         gpu_indx = c_ptr_plus(fact_data%gpu_row_ind, sz)
         gpu_ind = c_ptr_plus(gpu_indx, fact_data%rd_size*C_SIZEOF(dummy_int))
      end if
     
      if ( pos_def .or. job == SSIDS_SOLVE_JOB_BWD ) then
         call gather( stream_handle, ldln, nrhs, &
            gpu_x, n, gpu_ln, ldln, gpu_indx )
      else if ( job == SSIDS_SOLVE_JOB_DIAG ) then
         sz = (2*fact_data%values_L(lev)%off_col_ind) * C_SIZEOF(dummy_real)
         gpu_v = c_ptr_plus(fact_data%gpu_diag, sz)
         call apply_d( stream_handle, ldlx, nrhs, gpu_v, &
            gpu_x, n, gpu_u, ldlx, gpu_indx )
         call scatter( stream_handle, ldlx, nrhs, &
            gpu_u, ldlx, gpu_x, n, gpu_indx )
         cycle
      else
         call gather_dx( stream_handle, ldln, nrhs, fact_data%gpu_diag, &
            gpu_x, n, gpu_ln, ldln, gpu_ind, gpu_indx )
      end if

      ncb = fact_data%values_L(lev)%ncb_slv_t
      if ( ncb .gt. 0 ) then
         call multinode_solve_t( stream_handle, ncb, nrhs, &
            fact_data%values_L(lev)%ptr_levL, gpu_ln, gpu_u, gpu_v, &
            fact_data%values_L(lev)%gpu_solve_t_data )
      end if

      sz = fact_data%values_L(lev)%off_col_ind * C_SIZEOF(dummy_int)
      gpu_indx = c_ptr_plus(fact_data%gpu_col_ind, sz)
      call scatter_sum( stream_handle, ldlx, nrhs, &
         gpu_u, ldlx, gpu_v, ldlx, &
         gpu_x, n, gpu_indx )

   end do
   
   call custack_free(gwork, lx_len*C_SIZEOF(dummy_real)) ! gpu_v
   call custack_free(gwork, lx_len*C_SIZEOF(dummy_real)) ! gpu_u
   call custack_free(gwork, ln_len*C_SIZEOF(dummy_real)) ! gpu_ln

   call custack_finalize(gwork, cuda_error)
   if(cuda_error.ne.0) return
   
end subroutine subtree_bwd_multisolve_gpu

end module spral_ssids_solve_gpu
