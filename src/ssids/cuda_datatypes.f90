! Copyright (c) 2013 Science and Technology Facilities Council (STFC)
! Authors: Evgueni Ovtchinnikov and Jonathan Hogg
!
! Interoperable datatypes for passing structured data to CUDA
! (Done as separate module from spral_ssids_cuda_interfaces so we can USE it
!  in interface blocks)
module spral_ssids_cuda_datatypes
   use iso_c_binding
   implicit none

   private
   ! Data types
   public :: load_nodes_type, assemble_cp_type, assemble_delay_type, &
      assemble_blk_type, smblk, gemv_transpose_lookup, reducing_d_solve_lookup,&
      trsv_lookup_type, scatter_lookup_type, gemv_notrans_lookup, &
      reduce_notrans_lookup, assemble_lookup_type, lookups_gpu_fwd, &
      lookups_gpu_bwd, multiblock_fact_type, multinode_fact_type, &
      cuda_stats_type, cstat_data_type, multisymm_type, &
      multiswap_type, &
      multisyrk_type, &
      node_data, node_solve_data, multireorder_data, multielm_data, &
      assemble_lookup2_type
   ! Constants
   public :: SLV_ASSEMBLE_NB, SLV_GEMV_NX, SLV_GEMV_NY, SLV_TRSM_TR_NBX, &
      SLV_TRSM_TR_NBY, SLV_REDUCING_D_SOLVE_THREADS_PER_BLOCK, &
      SLV_TRSV_NB_TASK, SLV_SCATTER_NB, GPU_ALIGN


   type, bind(C) :: load_nodes_type
      integer(C_INT) :: nnz   ! Number of entries to map
      integer(C_INT) :: lda   ! Leading dimension of A
      integer(C_INT) :: ldl   ! Leading dimension of L
      type(C_PTR) :: lcol     ! Pointer to non-delay part of L
      integer(C_INT) :: offn  ! Offset into nlist
      integer(C_LONG) :: offr ! Offset into rlist
   end type load_nodes_type

   type, bind(C) :: assemble_cp_type
      integer(C_INT) :: pvoffset
      type(C_PTR) :: pval
      integer(C_INT) :: ldp
      integer(C_INT) :: cm
      integer(C_INT) :: cn
      integer(C_INT) :: ldc
      integer(C_LONG) :: cvoffset
      type(C_PTR) :: cv
      type(C_PTR) :: rlist_direct
      type(C_PTR) :: ind
      integer(C_INT) :: sync_offset
      integer(C_INT) :: sync_wait_for
   end type assemble_cp_type

   type, bind(C) :: assemble_delay_type
      integer(C_INT) :: ndelay
      integer(C_INT) :: m
      integer(C_INT) :: n
      integer(C_INT) :: ldd
      integer(C_INT) :: lds
      type(C_PTR) :: dval
      type(C_PTR) :: sval
      integer(C_LONG) :: roffset
   end type assemble_delay_type

   type, bind(C) :: assemble_blk_type
      integer(C_INT) :: cp
      integer(C_INT) :: blk
   end type assemble_blk_type

   type, bind(C) :: smblk
      integer(C_INT) :: bcol
      integer(C_INT) :: blkm
      integer(C_INT) :: nelim
      type(C_PTR) :: lcol
      integer(C_SIZE_T) :: lpitch
      type(C_PTR) :: d
      type(C_PTR) :: rlist
      integer(C_SIZE_T) :: ndone
      integer(C_SIZE_T) :: upd
   end type

   type, bind(C) :: gemv_transpose_lookup
      integer(C_INT) :: m
      integer(C_INT) :: n
      type(C_PTR) :: a
      integer(C_INT) :: lda
      type(C_PTR) :: rlist
      integer(C_INT) :: yoffset
   end type

   type, bind(C) :: reducing_d_solve_lookup
      integer(C_INT) :: first_idx
      integer(C_INT) :: m
      integer(C_INT) :: n
      integer(C_INT) :: ldupd
      integer(C_INT) :: updoffset
      type(C_PTR) :: d
      type(C_PTR) :: perm
   end type

   type, bind(C) :: trsv_lookup_type
      integer(C_INT) :: n
      type(C_PTR) :: a
      integer(C_INT) :: lda
      integer(C_INT) :: x_offset
      integer(C_INT) :: sync_offset
   end type

   type, bind(C) :: scatter_lookup_type
      integer(C_INT) :: n
      integer(C_INT) :: src_offset
      type(C_PTR) :: index
      integer(C_INT) :: dest_offset
   end type

   type, bind(C) :: gemv_notrans_lookup
      integer(C_INT) :: m
      integer(C_INT) :: n
      type(C_PTR) :: a
      integer(C_INT) :: lda
      integer(C_INT) :: x_offset
      integer(C_INT) :: y_offset
   end type

   type, bind(C) :: reduce_notrans_lookup
      integer(C_INT) :: m
      integer(C_INT) :: n
      integer(C_INT) :: src_offset
      integer(C_INT) :: ldsrc
      integer(C_INT) :: dest_idx
      integer(C_INT) :: dest_offset
   end type

   type, bind(C) :: assemble_lookup_type
      integer(C_INT) :: m
      integer(C_INT) :: xend
      type(C_PTR) :: list
      integer(C_INT) :: x_offset
      integer(C_INT) :: contrib_idx
      integer(C_INT) :: contrib_offset
      integer(C_INT) :: nchild
      type(C_PTR) :: clen
      type(C_PTR) :: clists
      type(C_PTR) :: clists_direct
      integer(C_INT) :: cvalues_offset
      integer(C_INT) :: first
   end type

   type, bind(C) :: assemble_lookup2_type
      integer(C_INT) :: m
      integer(C_INT) :: nelim
      integer(C_INT) :: x_offset
      type(C_PTR) :: list
      integer(C_INT) :: cvparent
      integer(C_INT) :: cvchild
      integer(C_INT) :: sync_offset
      integer(C_INT) :: sync_waitfor
   end type

   type, bind(C) :: lookups_gpu_fwd ! for fwd slv
      integer(C_INT) :: nassemble
      integer(C_INT) :: nasm_sync
      integer(C_INT) :: nassemble2
      integer(C_INT) :: nasmblk
      integer(C_INT) :: ntrsv
      integer(C_INT) :: ngemv
      integer(C_INT) :: nreduce
      integer(C_INT) :: nscatter
      type(C_PTR) :: assemble
      type(C_PTR) :: assemble2
      type(C_PTR) :: asmblk
      type(C_PTR) :: trsv
      type(C_PTR) :: gemv
      type(C_PTR) :: reduce
      type(C_PTR) :: scatter
   end type

   type, bind(C) :: lookups_gpu_bwd ! for bwd slv
      integer(C_INT) :: ngemv
      integer(C_INT) :: nrds
      integer(C_INT) :: ntrsv
      integer(C_INT) :: nscatter
      type(C_PTR) :: gemv
      type(C_PTR) :: rds
      type(C_PTR) :: trsv
      type(C_PTR) :: scatter
      type(C_PTR) :: gemv_times
      type(C_PTR) :: rds_times
      type(C_PTR) :: trsv_times
      type(C_PTR) :: scatter_times
   end type

   type, bind(C) :: multiblock_fact_type
      integer(C_INT) :: nrows
      integer(C_INT) :: ncols
      integer(C_INT) :: ld
      integer(C_INT) :: p
      type(C_PTR) :: aptr
      type(C_PTR) :: ldptr
      integer(C_INT) :: offf
      type(C_PTR) :: dptr
      integer(C_INT) :: node
      integer(C_INT) :: offb
   end type multiblock_fact_type

   type, bind(C) :: multinode_fact_type
      integer(C_INT) :: nrows
      integer(C_INT) :: ncols
      type(C_PTR) :: lval
      type(C_PTR) :: ldval
      type(C_PTR) :: dval
      integer(C_INT) :: offp
      integer(C_INT) :: ib
      integer(C_INT) :: jb
      integer(C_INT) :: done
      integer(C_INT) :: rght
      integer(C_INT) :: lbuf
   end type multinode_fact_type

   type, bind(C) :: cuda_stats_type
      integer(C_INT) :: num_two
      integer(C_INT) :: num_neg
      integer(C_INT) :: num_zero
   end type cuda_stats_type

   type, bind(C) :: cstat_data_type
      integer(C_INT) :: nelim
      type(C_PTR) :: dval
   end type cstat_data_type

   type, bind(C) :: multisymm_type
      type(C_PTR) :: lcol
      integer(C_INT) :: ncols
      integer(C_INT) :: nrows
   end type multisymm_type

   type, bind(C) :: multiswap_type
      integer(C_INT) :: nrows
      integer(C_INT) :: ncols
      integer(C_INT) :: k
      type(C_PTR) :: lcol
      integer(C_INT) :: lda
      integer(C_INT) :: off
   end type multiswap_type

   type, bind(C) :: multisyrk_type
      integer(C_INT) :: first
      type(C_PTR) :: lval
      type(C_PTR) :: ldval
      integer(C_LONG) :: offc
      integer(C_INT) :: n
      integer(C_INT) :: k
      integer(C_INT) :: lda
      integer(C_INT) :: ldb
   end type

   type, bind(C) :: node_data
    type(C_PTR) :: ptr_v
    integer(C_INT) :: ld
    integer(C_INT) :: nrows
    integer(C_INT) :: ncols
   end type node_data

   type, bind(C) :: node_solve_data
    type(C_PTR) :: ptr_a;
    type(C_PTR) :: ptr_b;
    type(C_PTR) :: ptr_u;
    type(C_PTR) :: ptr_v;
    integer(C_INT) :: lda
    integer(C_INT) :: ldb
    integer(C_INT) :: ldu
    integer(C_INT) :: ldv
    integer(C_INT) :: nrows
    integer(C_INT) :: ncols
    integer(C_INT) :: nrhs
    integer(C_INT) :: offb
    integer(C_LONG) :: off_a
    integer(C_INT) :: off_b
    integer(C_INT) :: off_u
    integer(C_INT) :: off_v
   end type node_solve_data

   type, bind(C) :: multireorder_data
    integer(C_INT) :: node
    integer(C_INT) :: block
    integer(C_INT) :: nblocks
   end type multireorder_data

   type, bind(C) :: multielm_data
    integer(C_INT) :: node
    integer(C_INT) :: offb
   end type multielm_data

   ! Preprocessor constants
   integer, parameter :: SLV_ASSEMBLE_NB = 128 ! MUST be same as C #define
   integer, parameter :: SLV_GEMV_NX = 32 ! MUST be same as C #define
   integer, parameter :: SLV_GEMV_NY = 32 ! MUST be same as C #define
   integer, parameter :: SLV_TRSM_TR_NBX = 256 ! MUST be same as C #define
   integer, parameter :: SLV_TRSM_TR_NBY = 32 ! MUST be same as C #define
   integer, parameter :: SLV_REDUCING_D_SOLVE_THREADS_PER_BLOCK = 256 
      ! MUST be same as C #define
   integer, parameter :: SLV_TRSV_NB_TASK = 32 ! MUST be same as C #define
   integer, parameter :: SLV_SCATTER_NB = 256 ! MUST be same as C #define

   integer, parameter :: GPU_ALIGN = 256 ! Align on this byte boundary

end module spral_ssids_cuda_datatypes
