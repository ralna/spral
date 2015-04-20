! COPYRIGHT (c) 2011-4 The Science and Technology Facilities Council (STFC)
! Original date 20 December 2011, Version 1.0.0
!
! Written by: Jonathan Hogg and Jennifer Scott
!
! Originally based on HSL_MA97 v2.2.0
!

module spral_ssids
!$ use omp_lib
   use, intrinsic :: iso_c_binding
   use spral_cuda, only : cudaMemcpy_d2h, cudaMemcpy_h2d, cudaMalloc, cudaFree,&
                          c_ptr_plus, cudaStreamCreate, cudaStreamDestroy, &
                          cudaMemcpy2d, cudaMemcpyHostToDevice, &
                          cudaMemcpyDeviceToHost
   use spral_ssids_cuda_interfaces
   use spral_match_order, only : match_order_metis
   use spral_matrix_util, only : SPRAL_MATRIX_REAL_SYM_INDEF, &
                                 SPRAL_MATRIX_REAL_SYM_PSDEF, &
                                 convert_coord_to_cscl, clean_cscl_oop, &
                                 apply_conversion_map
   use spral_metis_wrapper, only : metis_order
   use spral_scaling, only : auction_scale_sym, equilib_scale_sym, &
                             hungarian_scale_sym, &
                             equilib_options, equilib_inform, &
                             hungarian_options, hungarian_inform
   use spral_ssids_alloc, only : smalloc_setup, smalloc, smfreeall
   use spral_ssids_analyse, only : analyse_phase, check_order, expand_matrix, &
                                   expand_pattern
   use spral_ssids_datatypes
   use spral_ssids_factor_gpu, only : parfactor
   use spral_ssids_solve_cpu, only : solve_calc_chunk, inner_solve, &
                                     subtree_bwd_solve
   use spral_ssids_solve_gpu, only : bwd_solve_gpu, fwd_solve_gpu, &
                                     fwd_multisolve_gpu, bwd_multisolve_gpu, &
                                     d_solve_gpu, free_lookup_gpu
   implicit none

   private
   ! Data types
   public :: ssids_akeep, ssids_fkeep, ssids_options, ssids_inform
   ! User interface routines
   public :: ssids_analyse,         & ! Analyse phase, CSC-lower input
             ssids_analyse_coord,   & ! Analyse phase, Coordinate input
             ssids_factor,          & ! Factorize phase
             ssids_solve,           & ! Solve phase
             ssids_free,            & ! Free akeep and/or fkeep
             ssids_enquire_posdef,  & ! Pivot information in posdef case
             ssids_enquire_indef,   & ! Pivot information in indef case
             ssids_alter              ! Alter diagonal

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Make interfaces generic.
   interface ssids_analyse
      module procedure analyse_double
   end interface ssids_analyse

   interface ssids_analyse_coord
      module procedure ssids_analyse_coord_double
   end interface ssids_analyse_coord

   interface ssids_factor
      module procedure ssids_factor_double
   end interface ssids_factor

   interface ssids_solve
      module procedure ssids_solve_one_double
      module procedure ssids_solve_mult_double
   end interface ssids_solve

   interface ssids_free
      module procedure free_akeep_double
      module procedure free_fkeep_double
      module procedure free_both_double
   end interface ssids_free

   interface ssids_enquire_posdef
      module procedure ssids_enquire_posdef_double
   end interface ssids_enquire_posdef

   interface ssids_enquire_indef
      module procedure ssids_enquire_indef_double
   end interface ssids_enquire_indef

   interface ssids_alter
      module procedure ssids_alter_double
   end interface ssids_alter

contains

!****************************************************************************
!
! Analyse phase.
! Matrix entered in CSC format (lower triangle).
! The user optionally inputs the pivot order. If not, metis called. 
! Structure is then expanded.
! Supervariables are computed
! and then the assembly tree is constructed and the data structures
! required by the factorization are set up.
! There is no checking of the user's data if check = .false.
! Otherwise, matrix_util routines are used to clean data.
!
subroutine analyse_double(check, n, ptr, row, akeep, options, inform, &
      order, val)
   logical, intent(in) :: check ! if set to true, matrix data is checked
     ! and cleaned data stored in akeep
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   type(ssids_akeep), intent(out) :: akeep ! See derived-type declaration
   type(ssids_options), intent(in) :: options ! See derived-type declaration
   type(ssids_inform), intent(out) :: inform  ! See derived-type declaration
   integer, optional, intent(inout) :: order(:)
     ! Must be present and set on entry if options%ordering = 0.
     ! If options%ordering = 0 and i is used to index a variable, |order(i)|
     ! must hold its position in the pivot sequence. If a 1x1 pivot i is
     ! required, the user must set order(i)>0. If a 2x2 pivot involving
     ! variables i and j is required, the user must set
     ! order(i)<0, order(j)<0 and |order(j)| = |order(i)|+1.
     ! If i is not used to index a variable,
     ! order(i) must be set to zero.
     ! On exit, holds the pivot order to be used by factorization.
     ! Note: this input is consistent with our out-of-core solvers.
     !!!!! Note: although we allow 2x2 pivots to be input, we actually ignore 
     ! the signs (we reset signs of order after call to hsl_mc68 or hsl_mc80)
   real(wp), optional, intent(in) :: val(:) ! must be present
     ! if a matching-based elimination ordering is required
     ! (options%ordering 2).
     ! If present,  val(k) must hold the value of the entry in row(k).

   character(50)  :: context      ! Procedure name (used when printing).
   integer :: mu_flag       ! error flag for matrix_util routines
   integer :: mp           ! stream number for diagnostic messages
   integer :: nout         ! stream for errors
   integer :: nout1        ! stream for warnings
   integer :: nz           ! entries in expanded matrix
   integer :: st           ! stat parameter
   integer :: flag         ! error flag for metis

   integer, dimension(:), allocatable :: order2
   integer, dimension(:), allocatable :: ptr2 ! col. pointers for expanded mat
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix

   ! The following are only used for matching-based orderings
   real(wp), dimension(:), allocatable :: val_clean ! cleaned values if
     ! val is present and checking is required. 
   real(wp), dimension(:), allocatable :: val2 ! expanded matrix if
     ! val is present.

   integer :: mo_flag

   ! Initialise
   context = 'ssids_analyse'
   call ssids_free(akeep, inform%cuda_error)
   if(inform%cuda_error.ne.0) then
      inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
      akeep%flag = inform%flag
      call ssids_print_flag(context,nout,inform%flag, &
         cuda_error=inform%cuda_error)
      return
   endif
   inform%flag = 0
   inform%matrix_missing_diag = 0
   inform%matrix_outrange = 0
   inform%matrix_dup = 0
   inform%matrix_rank = n
   inform%maxdepth = 0
   inform%num_sup = 0
   inform%stat = 0

   ! Set stream numbers
   mp = options%unit_diagnostics
   nout = options%unit_error
   if (options%print_level < 0) nout = -1
   nout1 = options%unit_warning
   if (options%print_level < 0) nout1 = -1

   ! Print status on entry
   if (options%print_level>=1 .and. mp>=0) then
     write (mp,'(/a)') ' On entry to ssids_analyse:'
     write (mp,'(a,i15)') ' options%print_level       =  ', &
        options%print_level
     write (mp,'(a,i15)') ' options%unit_diagnostics  =  ',mp
     write (mp,'(a,i15)') ' options%unit_error        =  ',options%unit_error
     write (mp,'(a,i15)') ' options%unit_warning      =  ',options%unit_warning
     write (mp,'(a,i15)') ' options%nemin             =  ',options%nemin
     write (mp,'(a,i15)') ' options%ordering          =  ',options%ordering
     write (mp,'(a,i15)') ' n                         =  ',n
   end if

   akeep%check = check
   akeep%n = n
   akeep%flag = 0

   ! Checking of matrix data
   if (n < 0) then
      inform%flag = SSIDS_ERROR_A_N_OOR
      call ssids_print_flag(context,nout,inform%flag)
      akeep%flag = inform%flag
      return
   end if

   if (n .eq. 0) then
      akeep%nnodes = 0
      allocate(akeep%sptr(0), stat=st) ! used to check if analyse has been run
      if (st .ne. 0) go to 490
      akeep%matrix_dup = 0
      akeep%matrix_missing_diag = 0
      akeep%matrix_outrange = 0
      akeep%maxdepth = 0
      akeep%num_sup = 0
      return
   end if

   ! check options%ordering has a valid value
   if (options%ordering < 0 .or. options%ordering > 2) then
      inform%flag = SSIDS_ERROR_ORDER
      call ssids_print_flag(context,nout,inform%flag)
      akeep%flag = inform%flag
      return
   end if

   ! check val present when expected
   if (options%ordering.eq.2) then
     if (.not.present(val)) then
        inform%flag = SSIDS_ERROR_VAL
        call ssids_print_flag(context,nout,inform%flag)
        akeep%flag = inform%flag
        return
     end if
   end if

   akeep%ne = ptr(n+1)-1

   st = 0
   if (check) then
      deallocate (akeep%ptr,stat=st)
      allocate (akeep%ptr(n+1),stat=st)
      if (st .ne. 0) go to 490

      if (present(val)) then
         call clean_cscl_oop(SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, &
            akeep%ptr, akeep%row, mu_flag, val_in=val, val_out=val_clean,   &
            lmap=akeep%lmap, map=akeep%map,  &
            noor=inform%matrix_outrange, ndup=inform%matrix_dup)
      else
         call clean_cscl_oop(SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, &
            akeep%ptr, akeep%row, mu_flag, lmap=akeep%lmap, map=akeep%map,  &
            noor=inform%matrix_outrange, ndup=inform%matrix_dup)
      end if
      ! Check for errors
      if (mu_flag < 0) then
         if (mu_flag .eq. -1) inform%flag  = SSIDS_ERROR_ALLOCATION
         if (mu_flag .eq. -5) inform%flag  = SSIDS_ERROR_A_PTR
         if (mu_flag .eq. -6) inform%flag  = SSIDS_ERROR_A_PTR
         if (mu_flag .eq. -10) inform%flag = SSIDS_ERROR_A_ALL_OOR
         call ssids_print_flag(context,nout,inform%flag)
         akeep%flag = inform%flag
         return
      end if

      ! Check whether warning needs to be raised
      ! Note: same numbering of positive flags as in matrix_util
      if (mu_flag > 0) then
         inform%flag = mu_flag
         call ssids_print_flag(context,nout1,inform%flag)
      end if
      nz = akeep%ptr(n+1) - 1
   else
      nz = ptr(n+1)-1
   end if

   !
   ! If the pivot order is not supplied, we need to compute an order.
   ! Otherwise, we check the supplied order.
   !

   deallocate (akeep%invp,stat=st)
   allocate (akeep%invp(n),order2(n),ptr2(n+1),row2(2*nz),stat=st)
   if (st .ne. 0) go to 490
   if(options%ordering.eq.2) then
      allocate(akeep%scaling(n), val2(2*nz), stat=st)
      if (st .ne. 0) go to 490
   end if

   select case(options%ordering)
   case(0)
      if (.not.present(order)) then
         ! we have an error since user should have supplied the order
         inform%flag = SSIDS_ERROR_ORDER
         akeep%flag = inform%flag
         call ssids_print_flag(context,nout,inform%flag)
         return
      end if
      call check_order(n,order,akeep%invp,akeep,options,inform)
      if (inform%flag < 0) go to 490
      order2(1:n) = order(1:n)
      if (check) then
         call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)
      else
         call expand_pattern(n, nz, ptr, row, ptr2, row2)
      end if
   case(1)
      ! METIS ordering
      if (check) then
         call metis_order(n, akeep%ptr, akeep%row, order2, akeep%invp, &
            flag, inform%stat)
         call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)
      else
         call metis_order(n, ptr, row, order2, akeep%invp, &
            flag, inform%stat)
         call expand_pattern(n, nz, ptr, row, ptr2, row2)
      end if
      if (flag < 0) go to 490
   case(2)
      ! matching-based ordering required
      ! Expand the matrix as more efficient to do it and then
      ! call match_order_metis() with full matrix supplied

      if (check) then
         call expand_matrix(n, nz, akeep%ptr, akeep%row, val_clean, ptr2, &
            row2, val2)
         deallocate (val_clean,stat=st)
      else
         call expand_matrix(n, nz, ptr, row, val, ptr2, row2, val2)
      end if

      call match_order_metis(n, ptr2, row2, val2, order2, akeep%scaling, &
         mo_flag, inform%stat)

      select case(mo_flag)
      case(0)
         ! Success; do nothing
      case(1)
         ! singularity warning required
         inform%flag = SSIDS_WARNING_ANAL_SINGULAR
         call ssids_print_flag(context,nout1,inform%flag)
      case(-1)
         inform%flag = SSIDS_ERROR_ALLOCATION
         call ssids_print_flag(context,nout,inform%flag)
         akeep%flag = inform%flag
         return
      case default
         inform%flag = SSIDS_ERROR_UNKNOWN
         call ssids_print_flag(context,nout,inform%flag)
         akeep%flag = inform%flag
         return
      end select

      deallocate (val2,stat=st)
   end select


   ! perform rest of analyse
   if (check) then
      call analyse_phase(n, akeep%ptr, akeep%row, ptr2, row2, order2,  &
         akeep%invp, akeep, options, inform)
   else
      call analyse_phase(n, ptr, row, ptr2, row2, order2, akeep%invp, &
         akeep, options, inform)
   end if

   if (present(order)) order(1:n) = abs(order2(1:n))

   490 continue
   inform%stat = st
   if (inform%stat .ne. 0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
   end if
   akeep%flag = inform%flag

end subroutine analyse_double

!****************************************************************************
!
! Analyse phase.
! Matrix entered in coordinate format.
! matrix_util routine is used to convert the data to CSC format.
! The user optionally inputs the pivot order. If not, metis called. 
! Structure is then expanded.
! Supervariables are computed
! and then the assembly tree is constructed and the data structures
! required by the factorization are set up.
!
subroutine ssids_analyse_coord_double(n, ne, row, col, akeep, options, &
      inform, order, val)
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: ne ! entries to be input by user
   integer, intent(in) :: row(:) ! row indices
   integer, intent(in) :: col(:) ! col indices
   type(ssids_akeep), intent(out) :: akeep ! See derived-type declaration
   type(ssids_options), intent(in) :: options ! See derived-type declaration
   type(ssids_inform), intent(out) :: inform      ! See derived-type declaration
   integer, intent(inout), optional  :: order(:)
      ! Must be present and set on entry if options%ordering = 0 
      ! i is used to index a variable, order(i) must
      ! hold its position in the pivot sequence.
      ! If i is not used to index a variable,
      ! order(i) must be set to zero.
      ! On exit, holds the pivot order to be used by factorization.
   real(wp), optional, intent(in) :: val(:) ! must be present
     ! if a matching-based elimination ordering is required 
     ! (options%ordering = 2).
     ! If present, val(k) must hold value of entry in row(k) and col(k).

   integer, dimension(:), allocatable :: ptr2 ! col. pointers for expanded mat
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix
   integer, dimension(:), allocatable :: order2 ! pivot order

   integer :: mo_flag

   real(wp), dimension(:), allocatable :: val_clean ! cleaned values if
     ! val is present.
   real(wp), dimension(:), allocatable :: val2 ! expanded matrix (val present)

   character(50)  :: context      ! Procedure name (used when printing).
   integer :: mu_flag       ! error flag for matrix_util routines
   integer :: mp           ! stream number for diagnostic messages
   integer :: nout         ! stream for errors
   integer :: nout1        ! stream for warnings
   integer :: nz           ! entries in expanded matrix
   integer :: flag         ! error flag for metis
   integer :: st           ! stat parameter

   ! Initialise
   context = 'ssids_analyse_coord'
   call ssids_free(akeep, inform%cuda_error)
   if(inform%cuda_error.ne.0) then
      inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
      akeep%flag = inform%flag
      call ssids_print_flag(context,nout,inform%flag, &
         cuda_error=inform%cuda_error)
      return
   endif
   inform%flag = 0
   inform%matrix_missing_diag = 0
   inform%matrix_outrange = 0
   inform%matrix_dup = 0
   inform%matrix_rank = n
   inform%maxdepth = 0
   inform%num_sup = 0
   inform%stat = 0

   ! Set stream numbers
   mp = options%unit_diagnostics
   nout = options%unit_error
   if (options%print_level < 0) nout = -1
   nout1 = options%unit_warning
   if (options%print_level < 0) nout1 = -1

   ! Output status on entry
   if (options%print_level>=1 .and. mp>=0) then
     write (mp,'(/a)') ' On entry to ssids_analyse_coord:'
     write (mp,'(a,i15)') ' options%print_level       =  ', &
        options%print_level
     write (mp,'(a,i15)') ' options%unit_diagnostics  =  ',mp
     write (mp,'(a,i15)') ' options%unit_error        =  ',options%unit_error
     write (mp,'(a,i15)') ' options%unit_warning      =  ',options%unit_warning
     write (mp,'(a,i15)') ' options%nemin             =  ',options%nemin
     write (mp,'(a,i15)') ' options%ordering          =  ',options%ordering
     write (mp,'(a,i15)') ' n                         =  ',n
     write (mp,'(a,i15)') ' ne                        =  ',ne
   end if

   akeep%check = .true.
   akeep%n = n
   akeep%ne = ne
   akeep%flag = 0

   !
   ! Checking of matrix data
   !
   if (n < 0 .or. ne < 0) then
      inform%flag = SSIDS_ERROR_A_N_OOR
      akeep%flag = inform%flag
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   if (n .eq. 0) then
      akeep%nnodes = 0
      allocate(akeep%sptr(0), stat=st) ! used to check if analyse has been run
      if (st .ne. 0) go to 490
      akeep%matrix_dup = 0
      akeep%matrix_missing_diag = 0
      akeep%matrix_outrange = 0
      akeep%maxdepth = 0
      akeep%num_sup = 0
      return
   end if

   ! check options%ordering has a valid value
   if (options%ordering < 0 .or. options%ordering > 2) then
      inform%flag = SSIDS_ERROR_ORDER
      call ssids_print_flag(context,nout,inform%flag)
      akeep%flag = inform%flag
      return
   end if

   ! check val present when expected
   if (options%ordering.eq.2) then
     if (.not.present(val)) then
        inform%flag = SSIDS_ERROR_VAL
        call ssids_print_flag(context,nout,inform%flag)
        akeep%flag = inform%flag
        return
     end if
   end if

   st = 0
   deallocate (akeep%ptr,stat=st)
   allocate (akeep%ptr(n+1),stat=st)
   if (st .ne. 0) go to 490

   if (present(val)) then
      call convert_coord_to_cscl(SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ne, row, &
         col, akeep%ptr, akeep%row, mu_flag, val_in=val, val_out=val_clean,  &
         lmap=akeep%lmap, map=akeep%map,                                     &
         noor=inform%matrix_outrange,  ndup=inform%matrix_dup)
   else
      call convert_coord_to_cscl(SPRAL_MATRIX_REAL_SYM_INDEF, n, n, ne, row, &
         col, akeep%ptr, akeep%row, mu_flag, lmap=akeep%lmap, map=akeep%map, &
         noor=inform%matrix_outrange,  ndup=inform%matrix_dup)
   end if

   ! Check for errors
   if (mu_flag < 0) then
      if (mu_flag .eq. -1)  inform%flag = SSIDS_ERROR_ALLOCATION
      if (mu_flag .eq. -10) inform%flag = SSIDS_ERROR_A_ALL_OOR
      call ssids_print_flag(context,nout,inform%flag)
      akeep%flag = inform%flag
      return
   end if

   ! Check whether warning needs to be raised
   ! Note: same numbering of positive flags as in matrix_util
   if (mu_flag > 0) then
      inform%flag = mu_flag
      call ssids_print_flag(context,nout1,inform%flag)
      akeep%flag = inform%flag
   end if

   nz = akeep%ptr(n+1) - 1

   ! If the pivot order is not supplied, we need to compute an order
   ! here, before we expand the matrix structure.
   ! Otherwise, we must check the supplied order.

   deallocate(akeep%invp, stat=st)
   allocate (akeep%invp(n),order2(n),ptr2(n+1),row2(2*nz),stat=st)
   if (st .ne. 0) go to 490
   if(options%ordering.eq.2) then
      allocate(val2(2*nz),akeep%scaling(n),stat=st)
      if (st .ne. 0) go to 490
   end if

   select case(options%ordering)
   case(0)
      if (.not.present(order)) then
         ! we have an error since user should have supplied the order
         inform%flag = SSIDS_ERROR_ORDER
         akeep%flag = inform%flag
         call ssids_print_flag(context,nout,inform%flag)
         return
      end if
      call check_order(n,order,akeep%invp,akeep,options,inform)
      if (inform%flag < 0) go to 490
      order2(1:n) = order(1:n)
      call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)

   case(1)
      ! METIS ordering
      call metis_order(n, akeep%ptr, akeep%row, order2, akeep%invp, &
         flag, inform%stat)
      if (flag < 0) go to 490
      call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)

   case(2)
      ! matching-based ordering required

      call expand_matrix(n, nz, akeep%ptr, akeep%row, val_clean, ptr2, row2, &
         val2)
      deallocate (val_clean,stat=st)

      call match_order_metis(n,ptr2,row2,val2,order2,akeep%scaling,mo_flag, &
         inform%stat)

      select case(mo_flag)
      case(0)
         ! Success; do nothing
      case(1)
         ! singularity warning required
         inform%flag = SSIDS_WARNING_ANAL_SINGULAR
         call ssids_print_flag(context,nout1,inform%flag)
      case(-1)
         inform%flag = SSIDS_ERROR_ALLOCATION
         call ssids_print_flag(context,nout,inform%flag)
         akeep%flag = inform%flag
         return
      case default
         inform%flag = SSIDS_ERROR_UNKNOWN
         call ssids_print_flag(context,nout,inform%flag)
         akeep%flag = inform%flag
         return
      end select

      deallocate (val2,stat=st)
   end select


   ! we now have the expanded structure held using ptr2, row2
   ! and we are ready to get on with the analyse phase.
   call analyse_phase(n, akeep%ptr, akeep%row, ptr2, row2, order2,  &
      akeep%invp, akeep, options, inform)
   if (inform%flag < 0) go to 490

   if (present(order)) order(1:n) = abs(order2(1:n))

   490 continue
   inform%stat = st
   if (inform%stat .ne. 0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
   end if
   akeep%flag = inform%flag
    
end subroutine ssids_analyse_coord_double

!****************************************************************************

! Following function used to get around target requirement of C_LOC()
integer(C_INT) function copy_to_gpu_non_target(gpu_ptr, src, sz)
   type(C_PTR) :: gpu_ptr
   real(wp), dimension(*), target, intent(in) :: src
   integer(C_SIZE_T), intent(in) :: sz

   copy_to_gpu_non_target = cudaMemcpy_h2d(gpu_ptr, C_LOC(src), sz)
end function copy_to_gpu_non_target


!****************************************************************************
!
! Factorize phase
!
subroutine ssids_factor_double(posdef, val, akeep, fkeep, options, inform, &
      scale, ptr, row)
   logical, intent(in) :: posdef 
   real(wp), dimension(*), target, intent(in) :: val ! A values (lwr triangle)
   type(ssids_akeep), intent(in) :: akeep
   type(ssids_fkeep), intent(inout) :: fkeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(out) :: inform
   real(wp), dimension(:), optional, intent(inout) :: scale ! used to hold
      ! row and column scaling factors. Must be set on entry if
      ! options%scaling <= 0
      ! Note: Has to be assumed shape, not assumed size or fixed size to work
      ! around funny compiler bug
   integer, dimension(akeep%n+1), optional, intent(in) :: ptr ! must be
      ! present if on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.
   integer, dimension(*), optional, intent(in) :: row ! must be present if
      ! on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.

   integer(long) :: sz
   real(wp), dimension(:), allocatable, target :: val2
   type(C_PTR), dimension(:), allocatable :: gpu_contribs
   type(C_PTR) :: gpu_val, gpu_scaling
   character(len=50) :: context

   integer :: i
   integer :: n, nz
   integer :: nout, nout1
   integer :: mp
   logical :: sing
   integer :: cuda_error, st
   type(thread_stats), dimension(:), allocatable :: stats ! one copy
      ! per thread, accumulates per thread statistics that are then summed to
      ! obtain global stats in inform.
   integer, dimension(:,:), allocatable :: map ! work array, one copy per
      ! thread. Size (0:n, num_threads), with 0 index used to track which
      ! node current map refers to.
   integer :: num_threads
   ! Solve parameters. Tree is broken up into multiple chunks. Parent-child
   ! relations between chunks are stored in fwd_ptr and fwd (see solve routine
   ! comments)
   type(smalloc_type), pointer :: next_alloc
   integer :: matrix_type
   real(wp), dimension(:), allocatable :: scaling

   ! Types related to scaling routines
   type(hungarian_options) :: hsoptions
   type(hungarian_inform) :: hsinform
   type(equilib_options) :: esoptions
   type(equilib_inform) :: esinform
   
   integer :: flag

   ! Setup for any printing we may require
   context = 'ssids_factor'
   mp = options%unit_diagnostics
   nout = options%unit_error
   if (options%print_level < 0) nout = -1
   nout1 = options%unit_warning
   if (options%print_level < 0) nout1 = -1

   ! Perform appropriate printing
   if (options%print_level >= 1 .and. mp >= 0) then
      if (posdef) then
         write (mp,'(//a,i2,a)') &
            ' Entering ssids_factor with posdef = .true. and :'
         write (mp,'(a,5(/a,i12),5(/a,es12.4))') &
            ' options parameters (options%) :', &
            ' print_level         Level of diagnostic printing           = ', &
            options%print_level,      &
            ' unit_diagnostics    Unit for diagnostics                   = ', &
            options%unit_diagnostics, &
            ' unit_error          Unit for errors                        = ', &
            options%unit_error,       &
            ' unit_warning        Unit for warnings                      = ', &
            options%unit_warning,     &
            ' scaling             Scaling control                        = ', &
            options%scaling
      else ! indef
         write (mp,'(//a,i2,a)') &
            ' Entering ssids_factor with posdef = .false. and :'
         write (mp,'(a,5(/a,i12),5(/a,es12.4))') &
            ' options parameters (options%) :', &
            ' print_level         Level of diagnostic printing           = ', &
            options%print_level,      &
            ' unit_diagnostics    Unit for diagnostics                   = ', &
            options%unit_diagnostics, &
            ' unit_error          Unit for errors                        = ', &
            options%unit_error,       &
            ' unit_warning        Unit for warnings                      = ', &
            options%unit_warning,     &
            ' scaling             Scaling control                        = ', &
            options%scaling,          &
            ' small               Small pivot size                       = ', &
            options%small,           &
            ' u                   Initial relative pivot tolerance       = ', &
            options%u,               &
            ' multiplier          Multiplier for increasing array sizes  = ', &
            options%multiplier
      end if
   end if

   if (.not.allocated(akeep%sptr) .or. akeep%flag < 0) then
      ! Analyse cannot have been run
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      fkeep%flag = inform%flag
      return
   end if

   ! Initialize inform output
   inform%flag = SSIDS_SUCCESS
   inform%matrix_dup = akeep%matrix_dup
   inform%matrix_missing_diag = akeep%matrix_missing_diag
   inform%matrix_outrange = akeep%matrix_outrange
   inform%maxdepth = akeep%maxdepth
   inform%num_sup = akeep%num_sup
   inform%maxfront = 0
   inform%num_neg = 0
   inform%num_delay = 0
   inform%num_factor = 0
   inform%num_flops = 0
   inform%num_sup = akeep%nnodes
   inform%num_two = 0
   inform%stat = 0

   fkeep%flag = 0
   fkeep%host_factors = .false.
   st = 0

   n = akeep%n

   if (akeep%nnodes.eq.0) then
      inform%flag = SSIDS_SUCCESS
      inform%matrix_rank = 0
      fkeep%flag = inform%flag
      return
   end if

   fkeep%pos_def = posdef
   if(posdef) then
      matrix_type = SPRAL_MATRIX_REAL_SYM_PSDEF
   else
      matrix_type = SPRAL_MATRIX_REAL_SYM_INDEF
   endif

   ! If matrix has been checked, produce a clean version of val in val2
   if (akeep%check) then
      nz = akeep%ptr(n+1) - 1
      allocate (val2(nz),stat=st)
      if (st .ne. 0) go to 10
      call apply_conversion_map(matrix_type, akeep%lmap, akeep%map, val, &
         nz, val2)
   else
      ! analyse run with no checking so must have ptr and row present
      if (.not.present(ptr)) inform%flag = SSIDS_ERROR_PTR_ROW
      if (.not.present(row)) inform%flag = SSIDS_ERROR_PTR_ROW
      if (inform%flag < 0) then
         call ssids_print_flag(context,nout,inform%flag)
         fkeep%flag = inform%flag
         return
      end if
      nz = akeep%ne
   end if

   ! At this point, either  ptr, row, val   
   !                  or    akeep%ptr, akeep%row, val2
   ! hold the lower triangular part of A

   !
   ! Perform scaling if required
   !
   if (options%scaling.gt.0 .or. present(scale)) then
      if(allocated(fkeep%scaling)) then
         if(size(fkeep%scaling).lt.n) then
            deallocate(fkeep%scaling, stat=st)
            allocate(fkeep%scaling(n), stat=st)
         end if
      else
         allocate(fkeep%scaling(n), stat=st)
      end if
      if (st.ne.0) go to 10
   else
      deallocate(fkeep%scaling, stat=st)
   end if

   if(allocated(akeep%scaling) .and. options%scaling.ne.3) then
      inform%flag = SSIDS_WARNING_MATCH_ORD_NO_SCALE
      call ssids_print_flag(context,nout1,inform%flag)
   end if

   select case (options%scaling)
   case(:0) ! User supplied or none
      if (present(scale)) then
         do i = 1, n
            fkeep%scaling(i) = scale(akeep%invp(i))
         end do
      end if
   case(1) ! Matching-based scaling by Hungarian Algorithm (MC64 algorithm)
      ! Allocate space for scaling
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 10
      ! Run Hungarian algorithm
      hsoptions%scale_if_singular = options%action
      if (akeep%check) then
         call hungarian_scale_sym(n, akeep%ptr, akeep%row, val2, scaling, &
            hsoptions, hsinform)
      else
         call hungarian_scale_sym(n, ptr, row, val, scaling, &
            hsoptions, hsinform)
      end if
      select case(hsinform%flag)
      case(-1)
         ! Allocation error
         st = hsinform%stat
         go to 10
      case(-2)
         ! Structually singular matrix and control%action=.false.
         inform%flag = SSIDS_ERROR_SINGULAR
         call ssids_print_flag(context,nout,inform%flag)
         fkeep%flag = inform%flag
         return
      end select
      ! Permute scaling to correct order
      do i = 1, n
         fkeep%scaling(i) = scaling(akeep%invp(i))
      end do
      ! Copy scaling(:) to user array scale(:) if present
      if (present(scale)) then
         scale(1:n) = scaling(1:n)
      end if
      ! Cleanup memory
      deallocate(scaling, stat=st)

   case(2) ! Matching-based scaling by Auction Algorithm
      ! Allocate space for scaling
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 10
      ! Run auction algorithm
      if (akeep%check) then
         call auction_scale_sym(n, akeep%ptr, akeep%row, val2, scaling, &
            options%auction, inform%auction)
      else
         call auction_scale_sym(n, ptr, row, val, scaling, &
            options%auction, inform%auction)
      end if
      if (inform%auction%flag .ne. 0) then
         ! only possible error is allocation failed
         st = inform%auction%stat
         go to 10 
      endif
      ! Permute scaling to correct order
      do i = 1, n
         fkeep%scaling(i) = scaling(akeep%invp(i))
      end do
      ! Copy scaling(:) to user array scale(:) if present
      if (present(scale)) then
         scale(1:n) = scaling(1:n)
      end if
      ! Cleanup memory
      deallocate(scaling, stat=st)

   case(3) ! Scaling generated during analyse phase for matching-based order
      if (.not.allocated(akeep%scaling)) then
         ! No scaling saved from analyse phase
         inform%flag = SSIDS_ERROR_NO_SAVED_SCALING
         call ssids_print_flag(context,nout,inform%flag)
         fkeep%flag = inform%flag
         return
      end if
      do i = 1, n
         fkeep%scaling(i) = akeep%scaling(akeep%invp(i))
      end do

   case(4:) ! Norm equilibriation algorithm (MC77 algorithm)
      ! Allocate space for scaling
      allocate(scaling(n), stat=st)
      if(st.ne.0) goto 10
      ! Run equilibriation algorithm
      if (akeep%check) then
         call equilib_scale_sym(n, akeep%ptr, akeep%row, val2, scaling, &
            esoptions, esinform)
      else
         call equilib_scale_sym(n, ptr, row, val, scaling, &
            esoptions, esinform)
      end if
      if (esinform%flag .ne. 0) then
         ! Only possible error is memory allocation failure
         st = esinform%stat
         go to 10
      endif
      ! Permute scaling to correct order
      do i = 1, n
         fkeep%scaling(i) = scaling(akeep%invp(i))
      end do
      ! Copy scaling(:) to user array scale(:) if present
      if (present(scale)) then
         do i = 1, n
            scale(akeep%invp(i)) = fkeep%scaling(i)
         end do
      end if
      ! Cleanup memory
      deallocate(scaling, stat=st)
   end select

   !if(allocated(fkeep%scaling)) &
   !   print *, "minscale, maxscale = ", minval(fkeep%scaling), &
   !      maxval(fkeep%scaling)

   ! Setup data storage
   if(allocated(fkeep%nodes)) then
      if(size(fkeep%nodes).lt.akeep%nnodes+1) then
         deallocate(fkeep%nodes,stat=st)
         allocate(fkeep%nodes(akeep%nnodes+1), stat=st)
         if (st .ne. 0) go to 10
      end if
   else
      allocate(fkeep%nodes(akeep%nnodes+1), stat=st)
      if (st .ne. 0) go to 10
   end if
   fkeep%nodes(1:akeep%nnodes+1)%ndelay = 0
         
   ! Do an inital storage allocation:
   ! * options%multiplier * n             integers (for nodes(:)%perm)
   ! * options%multiplier * (nfactor+2*n) reals    (for nodes(:)%lcol)
   ! FIXME: do we really need this memory????
   if(associated(fkeep%alloc)) then
      if(fkeep%alloc%imem_size.lt. &
            max(n+0_long, int(options%multiplier*n,kind=long))) then
         deallocate(fkeep%alloc%imem, stat=st)
         fkeep%alloc%imem_size = &
            max(n+0_long, int(options%multiplier*n,kind=long))
         allocate(fkeep%alloc%imem(fkeep%alloc%imem_size),stat=st)
         if (st .ne. 0) go to 10
      end if
      if(fkeep%alloc%rmem_size.lt. max(akeep%nfactor+2*n, &
            int(options%multiplier*real(akeep%nfactor,wp)+2*n,kind=long))) then
         deallocate(fkeep%alloc%rmem, stat=st)
         fkeep%alloc%rmem_size = max(akeep%nfactor+2*n, &
            int(options%multiplier*real(akeep%nfactor,wp)+2*n,kind=long))
         allocate(fkeep%alloc%rmem(fkeep%alloc%rmem_size), stat=st)
         if (st .ne. 0) go to 10
      end if
      next_alloc => fkeep%alloc
      do while(associated(next_alloc))
         next_alloc%rhead = 0
         next_alloc%ihead = 0
         next_alloc => next_alloc%next_alloc
      end do
      nullify(fkeep%alloc%top_real, fkeep%alloc%top_int)

   else
      allocate(fkeep%alloc, stat=st)
      if (st .ne. 0) go to 10
      call smalloc_setup(fkeep%alloc, &
         max(n+0_long, int(options%multiplier*n,kind=long)), &
         max(akeep%nfactor+2*n, &
            int(options%multiplier*real(akeep%nfactor,wp)+2*n,kind=long)), st)
      if (st .ne. 0) go to 10
   end if

   num_threads = 1
!$ num_threads = omp_get_max_threads()

   allocate(stats(num_threads), map(0:n, num_threads), stat=st)
   if (st .ne. 0) go to 10
   map(0, :) = -1 ! initally map unassociated with any node

   !
   ! GPU-based factorization
   !

   ! Setup child contribution array
   ! Note: only non-NULL where we're passing contributions between subtrees
   allocate(gpu_contribs(akeep%nnodes), stat=st)
   if(st.ne.0) goto 10
   gpu_contribs(:) = C_NULL_PTR

   ! Copy A values to GPU
   sz = akeep%nptr(akeep%nnodes+1) - 1
   if (akeep%check) then
      cuda_error = cudaMalloc(gpu_val, sz*C_SIZEOF(val2(1)))
      if(cuda_error.ne.0) goto 200
      cuda_error = cudaMemcpy_h2d(gpu_val, C_LOC(val2), sz*C_SIZEOF(val2(1)))
   else
      cuda_error = cudaMalloc(gpu_val, sz*C_SIZEOF(val(1)))
      if(cuda_error.ne.0) goto 200
      cuda_error = cudaMemcpy_h2d(gpu_val, C_LOC(val), sz*C_SIZEOF(val(1)))
   endif
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
      cuda_error = cudaMalloc(gpu_scaling, n*C_SIZEOF(fkeep%scaling(1)))
      if(cuda_error.ne.0) goto 200
      cuda_error = copy_to_gpu_non_target(gpu_scaling, fkeep%scaling, &
         n*C_SIZEOF(fkeep%scaling(1)))
      if(cuda_error.ne.0) goto 200

      ! Perform factorization
      call parfactor(fkeep%pos_def, akeep%child_ptr, akeep%child_list, n,  &
         akeep%nptr, akeep%gpu_nlist, gpu_val, akeep%nnodes, fkeep%nodes,  &
         akeep%sptr, akeep%sparent, akeep%rptr, akeep%rlist, akeep%invp,   &
         akeep%rlist_direct, akeep%gpu_rlist, akeep%gpu_rlist_direct,      &
         gpu_contribs, fkeep%stream_handle, fkeep%stream_data,             &
         fkeep%top_data, fkeep%gpu_rlist_with_delays,                      &
         fkeep%gpu_rlist_direct_with_delays, fkeep%gpu_clists,             &
         fkeep%gpu_clists_direct, fkeep%gpu_clen, fkeep%alloc, options, stats, &
         ptr_scale=gpu_scaling)
      cuda_error = cudaFree(gpu_scaling)
      if(cuda_error.ne.0) goto 200
   else
      call parfactor(fkeep%pos_def, akeep%child_ptr, akeep%child_list, n,  &
         akeep%nptr, akeep%gpu_nlist, gpu_val, akeep%nnodes, fkeep%nodes,  &
         akeep%sptr, akeep%sparent, akeep%rptr, akeep%rlist, akeep%invp,   &
         akeep%rlist_direct, akeep%gpu_rlist, akeep%gpu_rlist_direct,      &
         gpu_contribs, fkeep%stream_handle, fkeep%stream_data,             &
         fkeep%top_data, fkeep%gpu_rlist_with_delays,                      &
         fkeep%gpu_rlist_direct_with_delays, fkeep%gpu_clists,             &
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

   if (inform%flag < 0) then
      call ssids_print_flag(context, nout, inform%flag, &
         cuda_error=inform%cuda_error, st=inform%stat)
      fkeep%flag = inform%flag
      return
   end if

   if(akeep%n.ne.inform%matrix_rank) then
      ! Rank deficient
      ! Note: If we reach this point then must be options%action=.true.
      if ( options%action ) then
        inform%flag = SSIDS_WARNING_FACT_SINGULAR
      else
        inform%flag = SSIDS_ERROR_SINGULAR
      end if
      call ssids_print_flag(context, nout1, inform%flag)
   end if

   if (options%print_level >= 1 .and. options%unit_diagnostics >= 0) then
      write (options%unit_diagnostics,'(/a)') &
         ' Completed factorisation with:'
      write (options%unit_diagnostics, &
         '(a,2(/a,i12),2(/a,es12.4),5(/a,i12))') &
         ' information parameters (inform%) :', &
         ' flag                   Error flag                               = ',&
         inform%flag, &
         ' maxfront               Maximum frontsize                        = ',&
         inform%maxfront, &
         ' num_factor             Number of entries in L                   = ',&
         real(inform%num_factor), &
         ' num_flops              Number of flops performed                = ',&
         real(inform%num_flops), &
         ' num_two                Number of 2x2 pivots used                = ',&
         inform%num_two, &
         ' num_delay              Number of delayed eliminations           = ',&
         inform%num_delay, &
         ' rank                   Computed rank                            = ',&
         inform%matrix_rank, &
         ' num_neg                Computed number of negative eigenvalues  = ',&
         inform%num_neg
 
   end if

   fkeep%flag = inform%flag
   fkeep%matrix_rank = inform%matrix_rank
   fkeep%maxfront = inform%maxfront
   fkeep%num_delay = inform%num_delay
   fkeep%num_factor = inform%num_factor
   fkeep%num_flops = inform%num_flops
   fkeep%num_neg = inform%num_neg
   fkeep%num_two = inform%num_two
   return
   !!!!!!!!!!!!!!!!!!!!

   !
   ! Error handling
   !
   10 continue
   inform%flag = SSIDS_ERROR_ALLOCATION
   inform%stat = st
   fkeep%flag = inform%flag
   call ssids_print_flag(context, nout, inform%flag, st=inform%stat)
   return

   200 continue
   inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
   inform%cuda_error = cuda_error
   call ssids_print_flag(context, nout, inform%flag, &
      cuda_error=inform%cuda_error)
   return

end subroutine ssids_factor_double

!
! Copies all gpu data back to host
!
subroutine ssids_move_data_inner(akeep, fkeep, options, nout, context, inform)
   type(ssids_akeep), intent(in) :: akeep
   type(ssids_fkeep), intent(inout) :: fkeep
   type(ssids_options), intent(in) :: options
   integer, intent(in) :: nout
   character(len=*), intent(in) :: context
   type(ssids_inform), intent(inout) :: inform

   !integer :: lev
   integer :: cuda_error
   integer :: st

   ! We assume that the factor has been done on the GPU. Do we need to copy
   ! data back to host?
   if(options%use_gpu_solve) return ! Solve to be done on GPU, no movement
   if(fkeep%host_factors) return ! Data already moved

   ! Copy data as desired
   call copy_back_to_host(fkeep%host_factors, options%nstream, &
      fkeep%stream_data, &
      fkeep%top_data, fkeep%nodes, akeep%sptr, &
      akeep%rptr, fkeep%alloc, &
      cuda_error, st)
   if(st.ne.0) goto 100
   if(cuda_error.ne.0) goto 200

   return ! Normal return

   100 continue ! Fortran allocation error
   inform%flag = SSIDS_ERROR_ALLOCATION
   call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
   return

   200 continue ! CUDA error
   inform%cuda_error = cuda_error
   inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
   call ssids_print_flag(context,nout,inform%flag,cuda_error=inform%cuda_error)
   return
end subroutine ssids_move_data_inner

subroutine copy_back_to_host(host_factors, nstream, stream_data, top_data, &
      nodes, sptr, rptr, alloc, cuda_error, st)
   logical, intent(out) :: host_factors
   integer, intent(in) :: nstream
   type(gpu_type), dimension(:), intent(in) :: stream_data
   type(gpu_type), intent(in) :: top_data
   type(node_type), dimension(*), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(smalloc_type), intent(inout) :: alloc ! Contains actual memory
      ! allocations for L. Everything else (within the subtree) is just a
      ! pointer to this.
   integer, intent(out) :: cuda_error
   integer, intent(out) :: st

   integer :: stream

   host_factors = .true. ! Record that data has been copied to host

   st = 0

   do stream = 1, nstream
      call copy_stream_data_to_host(stream_data(stream), nodes, sptr, &
         rptr, alloc, cuda_error, st)
      if(cuda_error.ne.0 .or. st.ne.0) return
   end do
   call copy_stream_data_to_host(top_data, nodes, sptr, &
      rptr, alloc, cuda_error, st)
   if(cuda_error.ne.0 .or. st.ne.0) return

end subroutine copy_back_to_host

subroutine copy_stream_data_to_host(stream_data, &
      nodes, sptr, rptr, alloc, cuda_error, st)
   type(gpu_type), intent(in) :: stream_data
   type(node_type), dimension(*), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   type(smalloc_type), intent(inout) :: alloc ! Contains actual memory
      ! allocations for L. Everything else (within the subtree) is just a
      ! pointer to this.
   integer, intent(out) :: cuda_error
   integer, intent(out) :: st

   integer :: llist, lev, node, ndelay, blkn, blkm
   integer(long) :: offp
   integer(long) :: level_size
   real(wp), dimension(:), allocatable, target :: work
   real(wp), dimension(:), pointer :: lcol
   type(C_PTR) :: ptr_levL
   
   ! Initialize return values
   cuda_error = 0
   st = 0

   ! Shortcut empty streams (occurs for v. small matrices)
   if(stream_data%num_levels.eq.0) return

   ! Copy one level at a time, then split into nodes
   do lev = 1, stream_data%num_levels
    
      level_size = 0
      do llist = stream_data%lvlptr(lev), stream_data%lvlptr(lev + 1) - 1
         node = stream_data%lvllist(llist)
         ndelay = nodes(node)%ndelay
         blkn = sptr(node+1) - sptr(node) + ndelay
         blkm = int(rptr(node+1) - rptr(node)) + ndelay
         call smalloc(alloc, nodes(node)%lcol, (blkn+0_long)*blkm+2*blkn, &
            nodes(node)%rsmptr, nodes(node)%rsmsa, st)
         if (st .ne. 0) return
         level_size = level_size + (blkm + 2_long)*blkn
      end do
      
      allocate(work(level_size), stat=st)
      if(st.ne.0) return
      
      ptr_levL = c_ptr_plus( stream_data%values_L(lev)%ptr_levL, 0_C_SIZE_T )
      cuda_error = cudaMemcpy_d2h(C_LOC(work), ptr_levL, &
         level_size*C_SIZEOF(work(1)))
      if(cuda_error.ne.0) return

      do llist = stream_data%lvlptr(lev), stream_data%lvlptr(lev + 1) - 1
         node = stream_data%lvllist(llist)
         ndelay = nodes(node)%ndelay
         blkn = sptr(node+1) - sptr(node) + ndelay
         blkm = int(rptr(node+1) - rptr(node)) + ndelay
         lcol => nodes(node)%lcol
         offp = stream_data%off_L(node)
         call dcopy( (blkm + 2)*blkn, work(offp + 1), 1, lcol, 1 )
      end do
      
      deallocate ( work )
      
   end do
end subroutine copy_stream_data_to_host

!*************************************************************************
!
! Solve phase single x. 
!
subroutine ssids_solve_one_double(x1, akeep, fkeep, options, inform, job)
   real(wp), dimension(:), intent(inout) :: x1 ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i) is the corresponding component of the
      ! right-hand side.On exit, if i has been used to index a variable,
      ! x(i) holds solution for variable i
   type(ssids_akeep), intent(in) :: akeep
   type(ssids_fkeep), intent(inout) :: fkeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(out) :: inform
   integer, optional, intent(in) :: job

   integer :: ldx

   ldx = size(x1)
   call ssids_solve_mult_double(1, x1, ldx, akeep, fkeep, options, inform, &
      job=job)
end subroutine ssids_solve_one_double

!*************************************************************************

! Note: GPU solve + nrhs>1 doesn't work without presolve
subroutine ssids_solve_mult_double(nrhs, x, ldx, akeep, fkeep, options, &
      inform, job)
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,nrhs), intent(inout), target :: x
   type(ssids_akeep), intent(in) :: akeep ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i,j) is the corresponding component of the
      ! right-hand side for the jth system (j = 1,2,..., nrhs).
      ! On exit, if i has been used to index a variable,
      ! x(i,j) holds solution for variable i to system j
   ! For details of keep, options, inform : see derived type description
   type(ssids_fkeep), intent(inout) :: fkeep !inout for moving data
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(out) :: inform
   integer, optional, intent(in) :: job ! used to indicate whether
      ! partial solution required
      ! job = 1 : forward eliminations only (PLX = B)
      ! job = 2 : diagonal solve (DX = B) (indefinite case only)
      ! job = 3 : backsubs only ((PL)^TX = B)
      ! job = 4 : diag and backsubs (D(PL)^TX = B) (indefinite case only)
      ! job absent: complete solve performed

   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i, r
   integer :: local_job ! local job parameter
   integer :: n
   integer :: nout
   integer, dimension(:,:), allocatable :: map

   integer :: nchunk, num_threads
   integer, dimension(:), allocatable :: chunk_sa, chunk_en, fwd_ptr, fwd
   integer :: cuda_error

   type(C_PTR) :: gpu_x
   type(C_PTR) :: gpu_scale
   type(C_PTR) :: gpu_invp

   type(cuda_settings_type) :: user_settings ! Stores user values we change
      ! temporarily

   inform%flag = SSIDS_SUCCESS

   ! Perform appropriate printing
   if (options%print_level >= 1 .and. options%unit_diagnostics >= 0) then
      write (options%unit_diagnostics,'(//a)') &
         ' Entering ssids_solve with:'
      write (options%unit_diagnostics,'(a,4(/a,i12),(/a,i12))') &
         ' options parameters (options%) :', &
         ' print_level         Level of diagnostic printing        = ', &
         options%print_level, &
         ' unit_diagnostics    Unit for diagnostics                = ', &
         options%unit_diagnostics, &
         ' unit_error          Unit for errors                     = ', &
         options%unit_error, &
         ' unit_warning        Unit for warnings                   = ', &
         options%unit_warning, &
         ' nrhs                                                    = ', &
         nrhs
      if (nrhs > 1) write (options%unit_diagnostics,'(/a,i12)') &
         ' ldx                                                     = ', &
         ldx
   end if

   nout = options%unit_error
   if (options%print_level < 0) nout = -1
   context = 'ssids_solve'

   if (akeep%nnodes.eq.0) return

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   inform%flag = max(SSIDS_SUCCESS, fkeep%flag) ! Preserve warnings
   ! immediate return if already had an error
   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   n = akeep%n
   if (ldx .lt. n) then
      inform%flag = SSIDS_ERROR_X_SIZE
      call ssids_print_flag(context,nout,inform%flag)
      if (nout .ge. 0) write (nout,'(a,i8,a,i8)') &
         ' Increase ldx from ', ldx, ' to at least ', n
      return
   end if

   if (nrhs .lt. 1) then
      inform%flag = SSIDS_ERROR_X_SIZE
      call ssids_print_flag(context,nout,inform%flag)
      if (nout .ge. 0) write (nout,'(a,i8,a,i8)') &
         ' nrhs must be at least 1. nrhs = ', nrhs
      return
   end if

   call push_ssids_cuda_settings(user_settings, cuda_error)
   if(cuda_error.ne.0) goto 200

   ! Copy previous phases' inform data from akeep and fkeep
   inform%matrix_dup = akeep%matrix_dup
   inform%matrix_missing_diag = akeep%matrix_missing_diag
   inform%matrix_outrange = akeep%matrix_outrange
   inform%maxdepth = akeep%maxdepth
   inform%matrix_rank = fkeep%matrix_rank
   inform%maxfront = fkeep%maxfront
   inform%num_delay = fkeep%num_delay
   inform%num_factor = fkeep%num_factor
   inform%num_flops = fkeep%num_flops
   inform%num_neg = fkeep%num_neg
   inform%num_sup = akeep%num_sup
   inform%num_two = fkeep%num_two

   ! Set local_job
   local_job = 0
   if (present(job)) then
      if (job .lt. SSIDS_SOLVE_JOB_FWD .or. job .gt. SSIDS_SOLVE_JOB_DIAG_BWD) &
         inform%flag = SSIDS_ERROR_JOB_OOR
      if (fkeep%pos_def .and. job.eq.SSIDS_SOLVE_JOB_DIAG) &
         inform%flag = SSIDS_ERROR_JOB_OOR
      if (fkeep%pos_def .and. job.eq.SSIDS_SOLVE_JOB_DIAG_BWD) &
         inform%flag = SSIDS_ERROR_JOB_OOR
      if (inform%flag.eq.SSIDS_ERROR_JOB_OOR) then
         call ssids_print_flag(context,nout,inform%flag)
         return
      end if
      local_job = job
   end if

   if(.not.options%use_gpu_solve .and. options%presolve.ne.0) then
      ! Presolve and CPU solve incompatible
      inform%flag = SSIDS_ERROR_PRESOLVE_INCOMPAT
      call ssids_print_flag(context,nout,inform%flag)
      return
   endif

   if ( options%presolve == 0 ) then

     if (allocated(fkeep%scaling)) then
        if (local_job == SSIDS_SOLVE_JOB_ALL .or. &
              local_job == SSIDS_SOLVE_JOB_FWD) then
           do r = 1, nrhs
              !x(1:n,r) = x(1:n,r) * fkeep%scaling(1:n)
              do i = 1, n
                 x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
              end do
           end do
        end if
     end if

   else

     if (allocated(fkeep%scaling)) then
       cuda_error = cudaMalloc(gpu_scale, n*C_SIZEOF(fkeep%scaling(1)))
       if(cuda_error.ne.0) goto 200
       cuda_error = cudaMemcpy_h2d(gpu_scale, n, fkeep%scaling)
       if(cuda_error.ne.0) goto 200
       cuda_error = cudaMalloc(gpu_invp, n*C_SIZEOF(akeep%invp(1)))
       if(cuda_error.ne.0) goto 200
       cuda_error = cudaMemcpy_h2d(gpu_invp, n, akeep%invp)
       if(cuda_error.ne.0) goto 200
     end if

     cuda_error = cudaMalloc(gpu_x, nrhs*n*C_SIZEOF(x(1,1)))
     if(cuda_error.ne.0) goto 200
     if(n.eq.ldx) then
       cuda_error = cudaMemcpy_h2d(gpu_x, C_LOC(x), nrhs*n*C_SIZEOF(x(1,1)))
       if(cuda_error.ne.0) goto 200
     else
       cuda_error = cudaMemcpy2d(gpu_x, n*C_SIZEOF(x(1,1)), C_LOC(x), &
         ldx*C_SIZEOF(x(1,1)), n*C_SIZEOF(x(1,1)), int(nrhs, C_SIZE_T), &
         cudaMemcpyHostToDevice)
       if(cuda_error.ne.0) goto 200
     end if

     if(allocated(fkeep%scaling) .and. &
         (local_job == SSIDS_SOLVE_JOB_ALL .or. &
          local_job == SSIDS_SOLVE_JOB_FWD) ) then
       call scale( n, nrhs, gpu_x, n, gpu_scale, gpu_invp )
     end if
     
   end if

   ! Copy factor data to/from GPU as approriate (if required!)
   call ssids_move_data_inner(akeep, fkeep, options, nout, context, inform)

   ! We aim to have 4 chunks per thread to hopefully provide sufficient
   ! tree-level parallelism.
   num_threads = 1
!$ num_threads = omp_get_max_threads()
   call solve_calc_chunk(akeep%nnodes, fkeep%nodes, akeep%sparent, akeep%rptr, &
      4*num_threads, nchunk, chunk_sa, chunk_en, fwd_ptr, fwd, inform%stat)
   if(inform%stat.ne.0) goto 100

   if(options%use_gpu_solve .and. ( local_job.eq.SSIDS_SOLVE_JOB_FWD .or. &
         local_job.eq.SSIDS_SOLVE_JOB_ALL)) then
      allocate(map(0:akeep%n, num_threads), &
         stat=inform%stat)
      if(inform%stat.ne.0) goto 100
      if(options%presolve.eq.0) then
        call fwd_solve_gpu(fkeep%pos_def, akeep%child_ptr, akeep%child_list,   &
           akeep%n, akeep%invp, akeep%nnodes, fkeep%nodes, akeep%rptr,         &
           options%nstream, fkeep%stream_handle, fkeep%stream_data,            &
           fkeep%top_data, x, inform%stat, cuda_error)
        if(inform%stat.ne.0) goto 100
        if(cuda_error.ne.0) goto 200
      else
        call fwd_multisolve_gpu(akeep%nnodes, fkeep%nodes, akeep%rptr,     &
           options%nstream, fkeep%stream_handle, fkeep%stream_data,        &
           fkeep%top_data, nrhs, gpu_x, cuda_error, inform%stat)
        if(inform%stat.ne.0) goto 100
        if(cuda_error.ne.0) goto 200
      end if
      ! Fudge local_job if required to perform backwards solve
      if(local_job.eq.SSIDS_SOLVE_JOB_ALL) then
         if(fkeep%pos_def) then
            local_job = SSIDS_SOLVE_JOB_BWD
         else
            local_job = SSIDS_SOLVE_JOB_DIAG_BWD
         end if
      elseif(local_job.eq.SSIDS_SOLVE_JOB_FWD) then
         local_job = -1 ! done
      end if
   endif

   if(options%use_gpu_solve .and. local_job.eq.SSIDS_SOLVE_JOB_DIAG) then
      if(options%presolve.eq.0) then
         call d_solve_gpu(akeep%nnodes, akeep%sptr, options%nstream, &
            fkeep%stream_handle, fkeep%stream_data, fkeep%top_data, akeep%n, &
            akeep%invp, x, inform%stat, cuda_error)
         if(inform%stat.ne.0) goto 100
      else
         call bwd_multisolve_gpu(fkeep%pos_def, local_job, options%nstream, &
            fkeep%stream_handle, fkeep%stream_data, fkeep%top_data,   &
            nrhs, gpu_x, cuda_error)
      end if
      if(cuda_error.ne.0) goto 200
      local_job = -1 ! done
   endif

   ! Perform supernodal forward solve or diagonal solve (both in serial)
   call inner_solve(fkeep%pos_def, local_job, akeep%nnodes, &
      fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, akeep%invp, nrhs, &
      x, ldx, inform%stat)
   if (inform%stat .ne. 0) goto 100

   if( local_job.eq.SSIDS_SOLVE_JOB_DIAG_BWD .or. &
         local_job.eq.SSIDS_SOLVE_JOB_BWD .or. &
         local_job.eq.SSIDS_SOLVE_JOB_ALL ) then
      if(options%use_gpu_solve) then
        if(options%presolve.eq.0) then
           call bwd_solve_gpu(local_job, fkeep%pos_def, akeep%nnodes,      &
              akeep%sptr, options%nstream, fkeep%stream_handle,            &
              fkeep%stream_data, fkeep%top_data, akeep%invp, x,            &
              inform%stat, cuda_error)
           if(cuda_error.ne.0) goto 200
        else
           call bwd_multisolve_gpu(fkeep%pos_def, local_job, options%nstream, &
              fkeep%stream_handle, fkeep%stream_data, fkeep%top_data,   &
              nrhs, gpu_x, cuda_error)
           if(cuda_error.ne.0) goto 200
        end if
      else
         call subtree_bwd_solve(akeep%nnodes, 1, local_job, fkeep%pos_def,  &
            akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, &
            akeep%invp, nrhs, x, ldx, inform%stat)
      endif
   end if
   if (inform%stat .ne. 0) goto 100

   if ( options%presolve == 0 ) then

     if (allocated(fkeep%scaling)) then
        if (local_job == SSIDS_SOLVE_JOB_ALL .or. &
              local_job == SSIDS_SOLVE_JOB_BWD .or. &
              local_job == SSIDS_SOLVE_JOB_DIAG_BWD) then
           do r = 1, nrhs
              !x(1:n,r) = x(1:n,r) * fkeep%scaling(1:n)
              do i = 1, n
                 x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
              end do
           end do
        end if
     end if

   else

      if ( allocated(fkeep%scaling) .and. &
             (local_job == SSIDS_SOLVE_JOB_ALL .or. &
              local_job == SSIDS_SOLVE_JOB_BWD .or. &
              local_job == SSIDS_SOLVE_JOB_DIAG_BWD) ) then
         call scale( n, nrhs, gpu_x, n, gpu_scale, gpu_invp )
      end if

      if (allocated(fkeep%scaling)) then
         cuda_error = cudaFree( gpu_scale )
         if(cuda_error.ne.0) goto 200
         cuda_error = cudaFree( gpu_invp )
         if(cuda_error.ne.0) goto 200
      end if

      if(n.eq.ldx) then
        cuda_error = cudaMemcpy_d2h(C_LOC(x), gpu_x, nrhs*n*C_SIZEOF(x(1,1)))
        if(cuda_error.ne.0) goto 200
      else
        cuda_error = cudaMemcpy2d(C_LOC(x), ldx*C_SIZEOF(x(1,1)), gpu_x, &
          n*C_SIZEOF(x(1,1)), n*C_SIZEOF(x(1,1)), int(nrhs, C_SIZE_T), &
          cudaMemcpyDeviceToHost)
        if(cuda_error.ne.0) goto 200
      end if
      cuda_error = cudaFree(gpu_x)
      if(cuda_error.ne.0) goto 200

   end if

   call pop_ssids_cuda_settings(user_settings, cuda_error)
   if(cuda_error.ne.0) goto 200

   return

   100 continue
   inform%flag = SSIDS_ERROR_ALLOCATION
   call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
   call pop_ssids_cuda_settings(user_settings, cuda_error)
   return

   200 continue ! CUDA error
   inform%cuda_error = cuda_error
   inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
   call ssids_print_flag(context,nout,inform%flag,cuda_error=inform%cuda_error)
   call pop_ssids_cuda_settings(user_settings, cuda_error)
   return
end subroutine ssids_solve_mult_double

!*************************************************************************
!
! Return diagonal entries to user
!
subroutine ssids_enquire_posdef_double(akeep, fkeep, options, inform, d)
   type(ssids_akeep), intent(in) :: akeep
   type(ssids_fkeep), target, intent(in) :: fkeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(out) :: inform
   real(wp), dimension(*), intent(out) :: d

   integer :: blkn, blkm
   character(50)  :: context      ! Procedure name (used when printing).
   integer(long) :: i
   integer :: j
   integer :: n
   integer :: node
   integer :: nout
   integer :: piv

   type(node_type), pointer :: nptr

   context = 'ssids_enquire_posdef' 
   inform%flag = SSIDS_SUCCESS

   nout = options%unit_error
   if (options%print_level < 0) nout = -1

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      ! immediate return if had an error
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   if (.not.fkeep%pos_def) then
      inform%flag = SSIDS_ERROR_NOT_LLT
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   n = akeep%n
   ! ensure d is not returned undefined
   d(1:n) = zero ! ensure do not returned with this undefined
   if ( .not. fkeep%host_factors ) return
   
   piv = 1
   do node = 1, akeep%nnodes
      nptr => fkeep%nodes(node)
      blkn = akeep%sptr(node+1) - akeep%sptr(node)
      blkm = int(akeep%rptr(node+1) - akeep%rptr(node))
      i = 1
      do j = 1, blkn
         d(piv) = nptr%lcol(i)
         i = i + blkm + 1
         piv = piv + 1
      end do
   end do
end subroutine ssids_enquire_posdef_double

!*************************************************************************
! In indefinite case, the pivot sequence used will not necessarily be
!
! the same as that passed to ssids_factor (because of delayed pivots). This
! subroutine allows the user to obtain the pivot sequence that was
! actually used.
! also the entries of D^{-1} are returned using array d.
!
subroutine ssids_enquire_indef_double(akeep, fkeep, options, inform, &
      piv_order, d)
   type(ssids_akeep), intent(in) :: akeep
   type(ssids_fkeep), target, intent(in) :: fkeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(out) :: inform
   integer, dimension(*), optional, intent(out) :: piv_order
      ! If i is used to index a variable, its position in the pivot sequence
      ! will be placed in piv_order(i), with its sign negative if it is
      ! part of a 2 x 2 pivot; otherwise, piv_order(i) will be set to zero.
   real(wp), dimension(2,*), optional, intent(out) :: d ! The diagonal
      ! entries of D^{-1} will be placed in d(1,:i) and the off-diagonal
      ! entries will be placed in d(2,:). The entries are held in pivot order.

   integer :: blkn, blkm
   character(50)  :: context      ! Procedure name (used when printing).
   integer :: i, j, k
   integer :: n
   integer :: nd
   integer :: node
   integer :: nout
   integer(long) :: offset
   integer :: piv
   real(C_DOUBLE), dimension(:), allocatable, target :: d2
   type(C_PTR) :: srcptr
   integer, dimension(:), allocatable :: lvllookup
   integer :: st, cuda_error
   real(wp) :: real_dummy

   type(node_type), pointer :: nptr

   context = 'ssids_enquire_indef' 
   inform%flag = SSIDS_SUCCESS

   nout = options%unit_error
   if (options%print_level < 0) nout = -1

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      ! immediate return if had an error
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   if (fkeep%pos_def) then
      inform%flag = SSIDS_ERROR_NOT_LDLT
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   n = akeep%n
   if(present(d)) then
      ! ensure d is not returned undefined
      d(1:2,1:n) = zero
   end if
   
   if(fkeep%host_factors) then
      piv = 1
      do node = 1, akeep%nnodes
         nptr => fkeep%nodes(node)
         j = 1
         nd = nptr%ndelay
         blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
         blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
         offset = blkm*(blkn+0_long)
         do while(j .le. nptr%nelim)
            if (nptr%lcol(offset+2*j).ne.0) then
               ! 2x2 pivot
               if(present(piv_order))  then
                  k = akeep%invp( nptr%perm(j) )
                  piv_order(k) = -piv
                  k = akeep%invp( nptr%perm(j+1) )
                  piv_order(k) = -(piv+1)
               end if
               if(present(d)) then
                  d(1,piv) = nptr%lcol(offset+2*j-1)
                  d(2,piv) = nptr%lcol(offset+2*j)
                  d(1,piv+1) = nptr%lcol(offset+2*j+1)
                  d(2,piv+1) = 0
               end if
               piv = piv + 2
               j = j + 2
            else
               ! 1x1 pivot
               if(present(piv_order)) then
                  k = akeep%invp( nptr%perm(j) )
                  piv_order(k) = piv
               end if
               if(present(d)) then
                  d(1,piv) = nptr%lcol(offset+2*j-1)
                  d(2,piv) = 0
               end if
               piv = piv + 1
               j = j + 1
            end if
         end do
      end do
   else
      allocate(lvllookup(akeep%nnodes), d2(2*akeep%n), stat=st)
      if(st.ne.0) goto 100
      do i = 1, fkeep%top_data%num_levels
         do j = fkeep%top_data%lvlptr(i), fkeep%top_data%lvlptr(i+1)-1
            node = fkeep%top_data%lvllist(j)
            lvllookup(node) = i
         end do
      end do

      piv = 1
      do node = 1, akeep%nnodes
         nptr => fkeep%nodes(node)
         j = 1
         nd = nptr%ndelay
         blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
         blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
         offset = blkm*(blkn+0_long)
         srcptr = c_ptr_plus(nptr%gpu_lcol, offset*C_SIZEOF(real_dummy))
         cuda_error = cudaMemcpy_d2h(C_LOC(d2), srcptr, &
            2*nptr%nelim*C_SIZEOF(d2(1)))
         if(cuda_error.ne.0) goto 200
         do while(j .le. nptr%nelim)
            if (d2(2*j).ne.0) then
               ! 2x2 pivot
               if(present(piv_order))  then
                  k = akeep%invp( nptr%perm(j) )
                  piv_order(k) = -piv
                  k = akeep%invp( nptr%perm(j+1) )
                  piv_order(k) = -(piv+1)
               end if
               if(present(d)) then
                  d(1,piv) = d2(2*j-1)
                  d(2,piv) = d2(2*j)
                  d(1,piv+1) = d2(2*j+1)
                  d(2,piv+1) = 0
               end if
               piv = piv + 2
               j = j + 2
            else
               ! 1x1 pivot
               if(present(piv_order)) then
                  k = akeep%invp( nptr%perm(j) )
                  piv_order(k) = piv
               end if
               if(present(d)) then
                  d(1,piv) = d2(2*j-1)
                  d(2,piv) = 0
               end if
               piv = piv + 1
               j = j + 1
            end if
         end do
      end do
   endif

   return ! Normal return

   100 continue ! Memory allocation error
   inform%stat = st
   inform%flag = SSIDS_ERROR_ALLOCATION
   call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
   return

   200 continue ! CUDA error
   inform%cuda_error = cuda_error
   inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
   call ssids_print_flag(context,nout,inform%flag,cuda_error=inform%cuda_error)
   return
end subroutine ssids_enquire_indef_double

!*************************************************************************
!
! In indefinite case, the entries of D^{-1} may be changed using this routine.
!
subroutine ssids_alter_double(d, akeep, fkeep, options, inform)
   real(wp), dimension(2,*), intent(in) :: d  ! The required diagonal entries
     ! of D^{-1} must be placed in d(1,i) (i = 1,...n)
     ! and the off-diagonal entries must be placed in d(2,i) (i = 1,...n-1).
   type(ssids_akeep), intent(in) :: akeep
   type(ssids_fkeep), target, intent(inout) :: fkeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(out) :: inform

   integer :: blkm, blkn
   character(50)  :: context      ! Procedure name (used when printing).
   integer(long) :: ip
   integer :: i, j
   integer :: nd
   integer :: node
   integer :: nout
   integer :: piv
   real(wp), dimension(:), allocatable, target :: d2
   integer, dimension(:), allocatable :: lvllookup
   type(C_PTR) :: srcptr
   integer :: st, cuda_error
   real(wp) :: real_dummy

   type(node_type), pointer :: nptr

   inform%flag = SSIDS_SUCCESS

   context = 'ssids_alter'

   nout = options%unit_error
   if (options%print_level < 0) nout = -1

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
    end if

   ! immediate return if already had an error
   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      inform%flag = SSIDS_ERROR_CALL_SEQUENCE
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if

   if (fkeep%pos_def) then
      inform%flag = SSIDS_ERROR_NOT_LDLT
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if 

   if (options%presolve.ne.0) then
      inform%flag = SSIDS_ERROR_PRESOLVE_INCOMPAT
      call ssids_print_flag(context,nout,inform%flag)
      return
   end if
   
   if(fkeep%host_factors) then
      piv = 1
      do node = 1, akeep%nnodes
         nptr => fkeep%nodes(node)
         nd = nptr%ndelay
         blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
         blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
         ip = blkm*(blkn+0_long) + 1
         do j = 1, nptr%nelim
            nptr%lcol(ip)   = d(1,piv)
            nptr%lcol(ip+1) = d(2,piv)
            ip = ip + 2
            piv = piv + 1
         end do
      end do
   endif

   ! FIXME: Move to gpu_factors a la host_factors?
   if(options%use_gpu_solve) then
      allocate(lvllookup(akeep%nnodes), d2(2*akeep%n), stat=st)
      if(st.ne.0) goto 100
      do i = 1, fkeep%top_data%num_levels
         do j = fkeep%top_data%lvlptr(i), fkeep%top_data%lvlptr(i+1)-1
            node = fkeep%top_data%lvllist(j)
            lvllookup(node) = i
         end do
      end do

      piv = 1
      do node = 1, akeep%nnodes
         nptr => fkeep%nodes(node)
         nd = nptr%ndelay
         blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
         blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
         ip = 1
         do j = 1, nptr%nelim
            d2(ip)   = d(1,piv)
            d2(ip+1) = d(2,piv)
            ip = ip + 2
            piv = piv + 1
         end do
         srcptr = c_ptr_plus(fkeep%nodes(node)%gpu_lcol, &
            blkm*(blkn+0_long)*C_SIZEOF(real_dummy))
         cuda_error = cudaMemcpy_h2d(srcptr, C_LOC(d2), &
            2*nptr%nelim*C_SIZEOF(d2(1)))
         if(cuda_error.ne.0) goto 200
      end do
   endif

   return ! Normal return

   100 continue ! Memory allocation error
   inform%stat = st
   inform%flag = SSIDS_ERROR_ALLOCATION
   call ssids_print_flag(context,nout,inform%flag,st=inform%stat)
   return

   200 continue ! CUDA error
   inform%cuda_error = cuda_error
   inform%flag = SSIDS_ERROR_CUDA_UNKNOWN
   call ssids_print_flag(context,nout,inform%flag,cuda_error=inform%cuda_error)
   return

end subroutine ssids_alter_double

!*************************************************************************

subroutine free_akeep_double(akeep, cuda_error)
   type(ssids_akeep), intent(inout) :: akeep
   integer, intent(out) :: cuda_error

   integer :: st

   deallocate(akeep%child_ptr, stat=st)
   deallocate(akeep%child_list, stat=st)
   deallocate(akeep%invp, stat=st)
   deallocate(akeep%level, stat=st)
   deallocate(akeep%nlist, stat=st)
   deallocate(akeep%nptr, stat=st)
   deallocate(akeep%rlist, stat=st)
   deallocate(akeep%rlist_direct, stat=st)
   deallocate(akeep%rptr, stat=st)
   deallocate(akeep%sparent, stat=st)
   deallocate(akeep%sptr, stat=st)
   deallocate(akeep%subtree_work, stat=st)
   deallocate(akeep%ptr, stat=st)
   deallocate(akeep%row, stat=st)
   deallocate(akeep%map, stat=st)

   ! Free GPU arrays if needed
   if(C_ASSOCIATED(akeep%gpu_nlist)) then
      cuda_error = cudaFree(akeep%gpu_nlist)
      akeep%gpu_nlist = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif
   if(C_ASSOCIATED(akeep%gpu_rlist)) then
      cuda_error = cudaFree(akeep%gpu_rlist)
      akeep%gpu_rlist = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif
   if(C_ASSOCIATED(akeep%gpu_rlist_direct)) then
      cuda_error = cudaFree(akeep%gpu_rlist_direct)
      akeep%gpu_rlist_direct = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif
end subroutine free_akeep_double

!*******************************

subroutine free_fkeep_double(fkeep, cuda_error)
   type(ssids_fkeep), intent(inout) :: fkeep
   integer, intent(out) :: cuda_error

   integer :: st, i

   cuda_error = 0

   if (.not.allocated(fkeep%nodes)) return

   call smfreeall(fkeep%alloc)
   deallocate(fkeep%alloc)
   nullify(fkeep%alloc)

   deallocate(fkeep%nodes, stat=st)
   deallocate(fkeep%scaling, stat=st)

   ! And clean up factors
   call free_gpu_type(fkeep%top_data, cuda_error)
   if(allocated(fkeep%stream_data)) then
      do i = 1, size(fkeep%stream_data)
         call free_gpu_type(fkeep%stream_data(i), cuda_error)
      end do
      deallocate(fkeep%stream_data, stat=st)
   endif

   ! Cleanup top-level presolve info
   if(C_ASSOCIATED(fkeep%gpu_rlist_with_delays)) then
      cuda_error = cudaFree(fkeep%gpu_rlist_with_delays)
      fkeep%gpu_rlist_with_delays = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif
   if(C_ASSOCIATED(fkeep%gpu_clists)) then
      cuda_error = cudaFree(fkeep%gpu_clists)
      fkeep%gpu_clists = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif
   if(C_ASSOCIATED(fkeep%gpu_clists_direct)) then
      cuda_error = cudaFree(fkeep%gpu_clists)
      fkeep%gpu_clists = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif
   if(C_ASSOCIATED(fkeep%gpu_clen)) then
      cuda_error = cudaFree(fkeep%gpu_clen)
      fkeep%gpu_clen = C_NULL_PTR
      if(cuda_error.ne.0) return
   endif

   ! Release streams
   if(allocated(fkeep%stream_handle)) then
      do i = 1, size(fkeep%stream_handle)
         cuda_error = cudaStreamDestroy(fkeep%stream_handle(i))
         if(cuda_error.ne.0) return
      end do
      deallocate(fkeep%stream_handle, stat=st)
   endif
end subroutine free_fkeep_double

subroutine free_gpu_type(sdata, cuda_error)
   type(gpu_type), intent(inout) :: sdata
   integer, intent(out) :: cuda_error

   integer :: lev
   integer :: st

   if(allocated(sdata%values_L)) then
      do lev = 1, sdata%num_levels
         cuda_error = cudaFree(sdata%values_L(lev)%ptr_levL)
         if(cuda_error.ne.0) return
         if ( sdata%values_L(lev)%ncp_pre > 0 ) then
            sdata%values_L(lev)%ncp_pre = 0
            cuda_error = cudaFree(sdata%values_L(lev)%gpu_cpdata_pre)
            if(cuda_error.ne.0) return
            cuda_error = cudaFree(sdata%values_L(lev)%gpu_blkdata_pre)
            if(cuda_error.ne.0) return
         end if
         if ( sdata%values_L(lev)%ncp_post > 0 ) then
            sdata%values_L(lev)%ncp_post = 0
            cuda_error = cudaFree(sdata%values_L(lev)%gpu_cpdata_post)
            if(cuda_error.ne.0) return
            cuda_error = cudaFree(sdata%values_L(lev)%gpu_blkdata_post)
            if(cuda_error.ne.0) return
         end if
         if ( sdata%values_L(lev)%ncb_slv_n > 0 ) then
            sdata%values_L(lev)%ncb_slv_n = 0
            cuda_error = cudaFree(sdata%values_L(lev)%gpu_solve_n_data)
            if(cuda_error.ne.0) return
         end if
         if ( sdata%values_L(lev)%ncb_slv_t > 0 ) then
            sdata%values_L(lev)%ncb_slv_t = 0
            cuda_error = cudaFree(sdata%values_L(lev)%gpu_solve_t_data)
            if(cuda_error.ne.0) return
         end if
         if ( sdata%values_L(lev)%nexp > 0 ) then
            sdata%values_L(lev)%nexp = 0
            deallocate( sdata%values_L(lev)%export )
         end if
         if ( sdata%values_L(lev)%nimp > 0 ) then
            sdata%values_L(lev)%nimp = 0
            deallocate( sdata%values_L(lev)%import )
         end if
      end do
      deallocate( sdata%values_L, stat=st )
   endif
   deallocate( sdata%lvlptr, stat=st )
   deallocate( sdata%lvllist, stat=st )
   deallocate( sdata%off_L, stat=st )
   if(allocated(sdata%bwd_slv_lookup)) then
      do lev = 1, sdata%num_levels
         call free_lookup_gpu(sdata%bwd_slv_lookup(lev), cuda_error)
         if(cuda_error.ne.0) return
         call free_lookup_gpu(sdata%fwd_slv_lookup(lev), cuda_error)
         if(cuda_error.ne.0) return
      end do
      deallocate(sdata%bwd_slv_lookup, stat=st)
   endif
   if(sdata%presolve.ne.0) then
      deallocate( sdata%off_lx, stat=st )
      deallocate( sdata%off_lc, stat=st )
      deallocate( sdata%off_ln, stat=st )
      cuda_error = cudaFree(sdata%gpu_rlist_direct)
      sdata%gpu_rlist_direct = C_NULL_PTR
      if(cuda_error.ne.0) return
      cuda_error = cudaFree(sdata%gpu_sync)
      sdata%gpu_sync = C_NULL_PTR
      if(cuda_error.ne.0) return
      cuda_error = cudaFree(sdata%gpu_row_ind)
      sdata%gpu_row_ind = C_NULL_PTR
      if(cuda_error.ne.0) return
      cuda_error = cudaFree(sdata%gpu_col_ind)
      sdata%gpu_col_ind = C_NULL_PTR
      if(cuda_error.ne.0) return
      cuda_error = cudaFree(sdata%gpu_diag)
      sdata%gpu_diag = C_NULL_PTR
      if(cuda_error.ne.0) return
   end if
   
end subroutine free_gpu_type

!*******************************

subroutine free_both_double(akeep, fkeep, cuda_error)
   type(ssids_akeep), intent(inout) :: akeep
   type(ssids_fkeep), intent(inout) :: fkeep
   integer, intent(out) :: cuda_error

   call free_akeep_double(akeep, cuda_error)
   if(cuda_error.ne.0) return
   call free_fkeep_double(fkeep, cuda_error)
   if(cuda_error.ne.0) return
end subroutine free_both_double

end module spral_ssids
