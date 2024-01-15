!
! To convert from double:
! * Change wp
! * Change _double
! * Change err_tol to 5e-5
! * Change err_tol_scale to 5e-3
! * Change fred_small to 1e-6
!
program main
!$ use omp_lib
   use spral_hw_topology, only : numa_region
   use spral_matrix_util, only : SPRAL_MATRIX_REAL_SYM_PSDEF, &
      SPRAL_MATRIX_REAL_SYM_INDEF, print_matrix
   use spral_metis_wrapper, only : metis_order
   use spral_random
   use spral_random_matrix, only : random_matrix_generate
   use spral_scaling, only : hungarian_scale_sym, hungarian_options, &
      hungarian_inform
   use spral_ssids
   implicit none

   integer, parameter :: long = selected_int_kind(18)


   integer, parameter :: wp = kind(0d0)
!   real(wp), parameter :: err_tol = 5e-12
!   real(wp), parameter :: err_tol_scale = 2e-10
   real(wp), parameter :: err_tol = 5e-11
   real(wp), parameter :: err_tol_scale = 1e-08
   real(wp), parameter :: fred_small = 1e-14
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: zero = 0.0_wp

   integer :: errors

   integer, parameter :: we_unit = 11
   character(len=6) :: we_file = "we.out"
   integer, parameter :: dl_unit = 12
   character(len=6) :: dl_file = "dl.out"

   type :: matrix_type
      integer :: n
      integer(long) :: ne
      integer, dimension(:), allocatable :: ptr, row, col
      real(wp), dimension(:), allocatable :: val
   end type matrix_type

   ! Error flags
   integer, parameter :: SSIDS_SUCCESS                   = 0
   integer, parameter :: SSIDS_ERROR_CALL_SEQUENCE       = -1
   integer, parameter :: SSIDS_ERROR_A_N_OOR             = -2
   integer, parameter :: SSIDS_ERROR_A_PTR               = -3
   integer, parameter :: SSIDS_ERROR_A_ALL_OOR           = -4
   integer, parameter :: SSIDS_ERROR_SINGULAR            = -5
   integer, parameter :: SSIDS_ERROR_NOT_POS_DEF         = -6
   integer, parameter :: SSIDS_ERROR_PTR_ROW             = -7
   integer, parameter :: SSIDS_ERROR_ORDER               = -8
   integer, parameter :: SSIDS_ERROR_VAL                 = -9
   integer, parameter :: SSIDS_ERROR_X_SIZE              = -10
   integer, parameter :: SSIDS_ERROR_JOB_OOR             = -11
   integer, parameter :: SSIDS_ERROR_PRESOLVE_INCOMPAT   = -12
   integer, parameter :: SSIDS_ERROR_NOT_LLT             = -13
   integer, parameter :: SSIDS_ERROR_NOT_LDLT            = -14
   integer, parameter :: SSIDS_ERROR_NO_SAVED_SCALING    = -15
   integer, parameter :: SSIDS_ERROR_ALLOCATION          = -50
   integer, parameter :: SSIDS_ERROR_CUDA_UNKNOWN        = -51
   integer, parameter :: SSIDS_ERROR_CUBLAS_UNKNOWN      = -52
   integer, parameter :: SSIDS_ERROR_UNIMPLEMENTED       = -98
   integer, parameter :: SSIDS_ERROR_UNKNOWN             = -99

   ! warning flags
   integer, parameter :: SSIDS_WARNING_IDX_OOR          = 1
   integer, parameter :: SSIDS_WARNING_DUP_IDX          = 2
   integer, parameter :: SSIDS_WARNING_DUP_AND_OOR      = 3
   integer, parameter :: SSIDS_WARNING_MISSING_DIAGONAL = 4
   integer, parameter :: SSIDS_WARNING_MISS_DIAG_OORDUP = 5
   integer, parameter :: SSIDS_WARNING_ANAL_SINGULAR    = 6
   integer, parameter :: SSIDS_WARNING_FACT_SINGULAR    = 7
   integer, parameter :: SSIDS_WARNING_MATCH_ORD_NO_SCALE=8

   ! warning flags

   if(we_unit.gt.6) open(unit=we_unit,file=we_file,status="replace")
   if(dl_unit.gt.6) open(unit=dl_unit,file=dl_file,status="replace")

   errors = 0

   call test_warnings
   call test_errors
   call test_special
   call test_random
   call test_random_scale
   call test_big

   write(*, "(/a)") "=========================="
   write(*, "(a,i4)") "Total number of errors = ", errors

   if(we_unit.gt.6) close(we_unit)
   if(dl_unit.gt.6) close(dl_unit)

   if(errors.ne.0) stop 1 ! ERROR CODE for make check script

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_errors
   type(matrix_type) :: a
   type(ssids_options) :: options
   type(ssids_akeep) :: akeep
   type(ssids_fkeep) :: fkeep
   type(ssids_inform) :: info

   logical :: check
   integer :: i
   integer :: n
   integer(long) :: ne
   logical :: posdef
   integer :: nrhs
   integer :: temp
   integer :: cuda_error
   integer, dimension(:), allocatable :: order
   real(wp), dimension(:,:), allocatable :: d, x
   real(wp), dimension(:), allocatable :: x1, d1
   type(random_state) :: state

   options%unit_error = we_unit
   options%unit_warning = we_unit

   write(*,"(/a)") "======================"
   write(*,"(a)") "Testing errors:"
   write(*,"(a)") "======================"

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ssids_analyse (entry by columns)

   write(*,"(/a)") " * Testing bad arguments ssids_analyse (columns)"

   check = .true.
   options%ordering = 0

   call simple_mat_lower(a)
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do
   write(*,"(a)",advance="no") " * Testing n<0..............................."
   n = -1
   call ssids_analyse(check, n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_A_N_OOR)

   call simple_mat_lower(a)
   write(*,"(a)",advance="no") " * Testing ptr with zero component..........."
   temp = a%ptr(1)
   a%ptr(1) = 0
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_A_PTR)
   a%ptr(1) = temp

   write(*,"(a)",advance="no") " * Testing non-monotonic ptr................."
   temp = a%ptr(2)
   a%ptr(2) = a%ptr(3)
   a%ptr(3) = temp
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_A_PTR)
   a%ptr(3) = a%ptr(2)
   a%ptr(2) = temp

   write(*,"(a)",advance="no") " * Testing all A%row oor....................."
   A%row(1:a%ptr(2)-1) = 0
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_A_ALL_OOR)

   write(*,"(a)",advance="no") " * Testing nemin oor........................."
   call simple_mat_lower(a)
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do
   options%nemin = -1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_SUCCESS)
   call ssids_free(akeep, cuda_error)
   options%nemin = 8

   write(*,"(a)",advance="no") " * Testing order absent......................"
   options%ordering = 0
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing order too short..................."
   if (allocated(order)) deallocate(order)
   allocate(order(a%n-1))
   order(1:a%n-1) = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   call ssids_free(akeep, cuda_error)
   deallocate(order)

   call simple_mat(a)
   write(*,"(a)",advance="no") " * Testing order out of range above.........."
   allocate(order(a%n))
   order(1) = a%n+1
   do i = 2,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   deallocate(order)

   call simple_mat(a)
   write(*,"(a)",advance="no") " * Testing order out of range below.........."
   allocate(order(a%n))
   order(1) = 0
   do i = 2,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)

   write(*,"(a)",advance="no") " * Testing options%ordering out of range....."
   options%ordering = -1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)

   call simple_mat(a)
   write(*,"(a)",advance="no") " * Testing options%ordering oor.............."
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   options%ordering = 3
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   deallocate(order)
   options%ordering = 0

   call simple_mat(a)
   write(*,"(a)",advance="no") " * Testing val absent........................"
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   options%ordering = 2
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_ERROR_VAL)
   deallocate(order)
   options%ordering = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ssids_analyse (entry using coordinate form)

   write(*,"(/a)")&
       " * Testing bad arguments ssids_analyse_coord (coordinate form)"

   call simple_mat_lower(a)
   write(*,"(a)",advance="no") " * Testing order out of range above.........."
   allocate(order(a%n))
   order(1) = a%n+1
   do i = 2,a%n
     order(i) = i
   end do
   options%ordering = 0
   call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, info, &
       order=order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   deallocate(order)

   call simple_mat_lower(a)
   write(*,"(a)",advance="no") " * Testing order out of range below.........."
   allocate(order(a%n))
   order(1) = 0
   do i = 2,a%n
     order(i) = i
   end do
   options%ordering = 0
   call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, info, &
       order=order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   deallocate(order)

   call simple_mat_lower(a)
   write(*,"(a)",advance="no") " * Testing options%ordering oor.............."
   allocate(order(a%n))
   options%ordering = 25
   call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, info, &
       order=order)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   deallocate(order)

   call simple_mat_lower(a)
   allocate(order(a%n))
   write(*,"(a)",advance="no") " * Testing val absent........................"
   options%ordering = 2
   call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, info, &
       order=order)
   call print_result(info%flag, SSIDS_ERROR_VAL)
   deallocate(order)

   call simple_mat_lower(a)
   write(*,"(a)",advance="no") " * Testing order absent......................"
   options%ordering = 0
   call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, info)
   call print_result(info%flag, SSIDS_ERROR_ORDER)
   call ssids_free(akeep, cuda_error)

   call simple_mat_lower(a)
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do
   write(*,"(a)",advance="no") " * Testing n<0..............................."
   n = -1
   call ssids_analyse_coord(n, a%ne, a%row, a%col, akeep, options, info, &
      order=order)
   call print_result(info%flag, SSIDS_ERROR_A_N_OOR)

   write(*,"(a)",advance="no") " * Testing ne < 0............................"
   n = 1
   ne = -1
   call ssids_analyse_coord(n, ne, a%row, a%col, akeep, options, info, &
      order=order)
   call print_result(info%flag, SSIDS_ERROR_A_N_OOR)

   write(*,"(a)",advance="no") " * Testing all oor..........................."
   n = 1
   ne = 1
   a%row(1) = -1
   a%col(1) = -1
   call ssids_analyse_coord(n, ne, a%row, a%col, akeep, options, info, &
      order=order)
   call print_result(info%flag, SSIDS_ERROR_A_ALL_OOR)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ssids_factor

   posdef = .false.
   write(*,"(/a)") " * Testing errors from ssids_factor"

   write(*,"(a)",advance="no") " * Testing after analyse error..............."
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!!!
   call simple_mat(a)
   write(*,"(a)",advance="no") " * Testing not calling analyse..............."
   posdef = .false.
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
!!!!!!!!!!

   call simple_mat(a)
   check = .false.
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing ptr absent........................"
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   if(info%flag.ne.0) then
      write(*, "(a,i4)") &
     "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      call ssids_free(akeep,fkeep,cuda_error)
      return
   endif

   call ssids_factor(posdef,a%val,akeep,fkeep,options,info,row=a%row)
   call print_result(info%flag, SSIDS_ERROR_PTR_ROW)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!!!!

   call simple_mat(a)
   check = .false.
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing row absent........................"
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   if(info%flag.ne.0) then
      write(*, "(a,i4)") &
     "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      call ssids_free(akeep,fkeep,cuda_error)
      return
   endif

   call ssids_factor(posdef,a%val,akeep,fkeep,options,info,ptr=a%ptr)
   call print_result(info%flag, SSIDS_ERROR_PTR_ROW)
   call ssids_free(akeep, cuda_error)
   call ssids_free(fkeep, cuda_error)

!!!!!!!!!

   write(*,"(a)",advance="no") " * Testing factor with singular matrix......."
   call simple_sing_mat(a)
   check = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   if(info%flag < 0) then
      write(*, "(a,i4)") &
     "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif

   options%action = .false.
   posdef = .false.
   call ssids_factor(posdef,a%val, akeep, fkeep, options, info)
   call print_result(info%flag, SSIDS_ERROR_SINGULAR)
   call ssids_free(akeep, cuda_error)
   call ssids_free(fkeep, cuda_error)

!!!!!!!!!

   write(*,"(a)",advance="no") &
      " * Testing factor with singular matrix (MC64 scale)."
   call simple_sing_mat(a)
   check = .true.
   options%ordering = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info)
   if(info%flag < 0) then
      write(*, "(a,i4)") &
     "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif

   options%action = .false.
   options%scaling = 1 ! MC64
   posdef = .false.
   call ssids_factor(posdef,a%val, akeep, fkeep, options, info)
   call print_result(info%flag, SSIDS_ERROR_SINGULAR)
   call ssids_free(akeep, cuda_error)
   call ssids_free(fkeep, cuda_error)

!!!!!!!!!

   write(*,"(a)",advance="no") " * Testing factor psdef with indef..........."
   call simple_sing_mat(a)
   a%val(1) = -a%val(1)
   check = .true.
   options%ordering = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info)
   if(info%flag < 0) then
      write(*, "(a,i4)") &
     "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   options%scaling = 0
   call ssids_factor(.true., a%val, akeep, fkeep, options, info)
   call print_result(info%flag, SSIDS_ERROR_NOT_POS_DEF)
   call ssids_free(akeep, fkeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing factor psdef with indef, large..."
   call gen_bordered_block_diag(.false., (/ 15, 455, 10 /), 20, a%n, a%ptr, &
      a%row, a%val, state)
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, val=a%val)
   call ssids_factor(.true., a%val, akeep, fkeep, options, info)
   call print_result(info%flag, SSIDS_ERROR_NOT_POS_DEF)
   call ssids_free(akeep, fkeep, cuda_error)

!!!!!!!!!
   write(*,"(a)",advance="no") " * Testing u oor............................."
   call simple_mat(a)
   check = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   if(info%flag < 0) then
      write(*, "(a,i4)") &
     "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   posdef = .false.
   options%u = -0.1
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   options%u = 0.01
   call print_result(info%flag, SSIDS_SUCCESS)
   call ssids_free(akeep, fkeep, cuda_error)

!!!!!!!!!
   write(*,"(a)",advance="no") " * Testing options%scaling=3 no matching....."
   call simple_mat(a)
   check = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   if(info%flag < 0) then
      write(*, "(a,i4)") &
     "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   posdef = .false.
   options%scaling = 3 ! MC64 from matching-based order
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   options%scaling = 0
   call print_result(info%flag, SSIDS_ERROR_NO_SAVED_SCALING)
   call ssids_free(akeep, fkeep, cuda_error)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! tests on call to ssids_solve

   write(*,"(/a)") " * Testing bad arguments ssids_solve"

   write(*,"(a)",advance="no") " * Testing solve after factor error.........."
   call simple_sing_mat(a)
   check = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)

   options%action = .false.
   options%scaling = 0
   posdef = .false.
   call ssids_factor(posdef,a%val, akeep, fkeep, options, info)
   if(info%flag /= SSIDS_ERROR_SINGULAR) then
      write(*, "(a,i4)") &
     "Unexpected flag returned by ssids_factor. flag = ", info%flag
      errors = errors + 1
      return
   endif

   if (allocated(x1)) deallocate(x1)
   allocate(x1(a%n))
   x1(1:a%n) = one
   call ssids_solve(x1,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!
   call simple_mat(a)
   deallocate(x1)
   allocate(x1(a%n))
   x1(1:a%n) = one
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing solve out of sequence............."
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)

   call ssids_solve(x1,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!

   call simple_mat(a)
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing job out of range below............"
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   call ssids_solve(x1,akeep,fkeep,options,info,job=-1)
   call print_result(info%flag, SSIDS_ERROR_JOB_OOR)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!

   call simple_mat(a)
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing job out of range above............"
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   call ssids_solve(x1,akeep,fkeep,options,info,job=5)
   call print_result(info%flag, SSIDS_ERROR_JOB_OOR)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!

   call simple_mat(a)
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing error in x (one rhs).............."
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   if (allocated(x1)) deallocate(x1)
   allocate(x1(a%n-1))
   call ssids_solve(x1,akeep,fkeep,options,info,job=5)
   call print_result(info%flag, SSIDS_ERROR_X_SIZE)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!

   call simple_mat(a)
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing error in lx......................."
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   if (allocated(x)) deallocate(x)
   allocate(x(a%n-1,1))
   call ssids_solve(1,x,a%n-1,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_X_SIZE)
   call ssids_free(akeep,fkeep,cuda_error)

!!!!!!

   call simple_mat(a)
   posdef = .false.
   write(*,"(a)",advance="no") " * Testing error in nrhs....................."
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call ssids_factor(posdef,a%val,akeep,fkeep,options,info)
   nrhs = -2
   if (allocated(x)) deallocate(x)
   allocate(x(a%n,1))
   call ssids_solve(nrhs,x,a%n,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_X_SIZE)
   call ssids_free(akeep,fkeep,cuda_error)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! tests on call to ssids_enquire_posdef

   write(*,"(/a)") " * Testing bad arguments ssids_enquire_posdef"

   write(*,"(a)",advance="no") " * Testing call to enquire after error......."
   call simple_sing_mat(a)
   check = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call assert_flag('analyse', info%flag, SSIDS_WARNING_MISSING_DIAGONAL)

   options%action = .false.
   posdef = .true.
   call ssids_factor(posdef, a%val, akeep, fkeep, options, info)
   call assert_flag('factor', info%flag, SSIDS_ERROR_NOT_POS_DEF)

   if (allocated(d1)) deallocate(d1)
   allocate(d1(a%n))
   call ssids_enquire_posdef(akeep,fkeep,options,info,d1)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep, fkeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing call to enquire out of seq........"
   call simple_mat(a)
   posdef = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, &
        akeep, options, info, order=order)

   if (allocated(d1)) deallocate(d1)
   allocate(d1(a%n))
   call ssids_enquire_posdef(akeep,fkeep,options,info,d1)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing call to enquire_posdef with indef."
   call simple_mat(a)
   call ssids_analyse(check, a%n, a%ptr, a%row, &
        akeep, options, info, order=order)
   call ssids_factor(.false., a%val, akeep, fkeep, options, &
      info)
   if (allocated(d1)) deallocate(d1)
   allocate(d1(a%n))
   call ssids_enquire_posdef(akeep,fkeep,options,info,d1)
   call print_result(info%flag, SSIDS_ERROR_NOT_LLT)
   call ssids_free(akeep, fkeep, cuda_error)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! tests on call to ssids_enquire_indef

   write(*,"(/a)") " * Testing bad arguments ssids_enquire_indef"

   write(*,"(a)",advance="no") " * Testing call to enquire after error......."
   call simple_sing_mat(a)
   check = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)

   options%action = .false.
   posdef = .false.
   call ssids_factor(posdef,a%val, akeep, fkeep, options, info)

   if (allocated(d)) deallocate(d)
   allocate(d(2,a%n))
   call ssids_enquire_indef(akeep,fkeep,options,info,piv_order=order,d=d)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep, fkeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing call to enquire out of seq........"
   call simple_mat(a)
   posdef = .false.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, &
        akeep, options, info, order=order)

   if (allocated(d)) deallocate(d)
   allocate(d(2,a%n))
   call ssids_enquire_indef(akeep,fkeep,options,info,piv_order=order,d=d)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing call to enquire_indef with posdef."
   call simple_mat(a)
   posdef = .false.
   options%ordering = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call ssids_factor(.true., a%val, akeep, fkeep, options, info)
   if (allocated(d)) deallocate(d)
   allocate(d(2,a%n))
   call ssids_enquire_indef(akeep,fkeep,options,info,piv_order=order,d=d)
   call print_result(info%flag, SSIDS_ERROR_NOT_LDLT)
   call ssids_free(akeep, fkeep, cuda_error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(*,"(/a)") " * Testing bad arguments ssids_alter"

   ! tests on call to ssids_alter
   write(*,"(a)",advance="no") " * Testing call to alter after error........."
   call simple_sing_mat(a)
   check = .true.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)

   options%action = .false.
   posdef = .false.
   call ssids_factor(posdef,a%val, akeep, fkeep, options, info)

   if (allocated(d)) deallocate(d)
   allocate(d(2,a%n))
   call ssids_alter(d,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep, fkeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing call to alter out of seq.........."
   call simple_mat(a)
   posdef = .false.
   options%ordering = 0
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i = 1,a%n
     order(i) = i
   end do

   call ssids_analyse(check, a%n, a%ptr, a%row, &
        akeep, options, info, order=order)

   if (allocated(d)) deallocate(d)
   allocate(d(2,a%n))
   call ssids_alter(d,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_CALL_SEQUENCE)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing call to alter_indef with posdef..."
   call simple_mat(a)
   posdef = .false.
   call ssids_analyse(check, a%n, a%ptr, a%row, &
      akeep, options, info, order=order)
   call ssids_factor(.true., a%val, akeep, fkeep, options, info)
   if (allocated(d)) deallocate(d)
   allocate(d(2,a%n))
   call ssids_alter(d,akeep,fkeep,options,info)
   call print_result(info%flag, SSIDS_ERROR_NOT_LDLT)
   call ssids_free(akeep, fkeep, cuda_error)

end subroutine test_errors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_warnings
   type(matrix_type) :: a
   type(ssids_options) :: options
   type(ssids_akeep) :: akeep
   type(ssids_inform) :: info

   logical :: check
   integer :: i
   logical :: posdef
   integer(long) :: ne
   integer, dimension(:), allocatable :: order
   real(wp), dimension(:,:), allocatable :: rhs
   real(wp), dimension(:,:), allocatable :: res
   real(wp), dimension(:,:), allocatable :: x
   real(wp), dimension(:), allocatable :: x1
   integer :: cuda_error

   write(*,"(a)")
   write(*,"(a)") "================"
   write(*,"(a)") "Testing warnings"
   write(*,"(a)") "================"

   options%unit_warning = we_unit
   options%unit_diagnostics = dl_unit
   options%print_level = 2

   options%ordering = 0 ! supply the ordering

   check = .true.

   write(*,"(/a/)") " * Testing warnings (columns)"

   call simple_mat(a)
   ne = a%ne
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do

   write(*,"(a)",advance="no") " * Testing out of range above............"
   call simple_mat(a,2)
   ne = a%ptr(a%n+1)-1
   a%ptr(a%n+1) = a%ptr(a%n+1) + 1
   a%row(ne+1) = -1
   a%val(ne+1) = 1.
   a%ne = ne + 1
   posdef = .true.
   call gen_rhs(a, rhs, x1, x, res, 1)
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_WARNING_IDX_OOR)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_IDX_OOR)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing out of range below............"
   call simple_mat_lower(a,1)
   ne = a%ne
   a%ptr(a%n+1) = a%ptr(a%n+1) + 1
   a%row(ne+1) = a%n + 1
   a%val(ne+1) = 1.
   a%ne = ne + 1
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag,SSIDS_WARNING_IDX_OOR)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_IDX_OOR)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing duplicates...................."
   call simple_mat(a,2)
   ne = a%ptr(a%n+1)-1
   a%ne = ne
   a%ptr(a%n+1) = a%ptr(a%n+1) + 1
   a%row(ne+1) = a%n
   a%val(ne+1) = 10.
   a%ne = ne + 1
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag,SSIDS_WARNING_DUP_IDX)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_DUP_IDX)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing out of range and duplicates..."
   call simple_mat_lower(a,2)
   ne = a%ptr(a%n+1)-1
   a%ne = ne
   a%ptr(a%n+1) = a%ptr(a%n+1) + 2
   a%row(ne+1) = a%n + 1
   a%val(ne+1) = 10.
   a%row(ne+2) = a%n
   a%val(ne+2) = 10.
   a%ne = ne + 2
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag,SSIDS_WARNING_DUP_AND_OOR)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_DUP_AND_OOR, fs=.true.)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)", advance="no") " * Testing missing diagonal entry (indef)....."
   a%ptr = (/ 1, 4, 5, 6, 7 /)
   a%row(1:6) = (/ 1, 2, 4,     2,    4,    4 /)
   a%val(1:6) = (/   10.0, 2.0, 3.0, &
                     10.0, &
                     4.0, &
                     10.0 /)
   posdef = .false.
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do
   a%ne = a%ptr(a%n+1) - 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_WARNING_MISSING_DIAGONAL)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_MISSING_DIAGONAL)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)", advance="no") " * Testing missing diagonal and out of range.."
   call simple_mat_lower(a)
   a%ptr = (/ 1, 4, 5, 6, 8 /)
   a%row(1:7) = (/ 1, 2, 4,     2,    4,    4,   -1 /)
   a%val(1:7) = (/   10.0, 2.0, 3.0, &
                     10.0, &
                     4.0, &
                     10.0, 1.0 /)
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do
   a%ne = a%ptr(a%n+1) - 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_WARNING_MISS_DIAG_OORDUP)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_MISS_DIAG_OORDUP)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") " * Testing arrays min size (zero diag)........"
   call simple_mat_zero_diag(a)
   a%ne = a%ptr(a%n+1) - 1
   call gen_rhs(a, rhs, x1, x, res, 1)
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_WARNING_MISSING_DIAGONAL)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_MISSING_DIAGONAL)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)", advance="no") " * Testing missing diagonal and duplicate....."
   call simple_mat(a)
   a%ptr = (/ 1, 5, 6, 7, 8 /)
   a%row(1:7) = (/ 1, 2, 2, 4,     2,    4,    4 /)
   a%val(1:7) = (/   10.0, 2.0, 3.0, 1.0, &
                     10.0, &
                     4.0, &
                     10.0 /)
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do
   a%ne = a%ptr(a%n+1) - 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, order)
   call print_result(info%flag, SSIDS_WARNING_MISS_DIAG_OORDUP)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_MISS_DIAG_OORDUP, fs=.true.)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)", advance="no") " * Testing analyse with structurally singular."
   call simple_sing_mat2(a)
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i=1,a%n
     order(i) = i
   end do
   a%ne = a%ptr(a%n+1) - 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call print_result(info%flag, SSIDS_WARNING_ANAL_SINGULAR)
   options%action = .true.
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_FACT_SINGULAR, fs=.true.)

   write(*,"(a)", advance="no") " * Testing analyse with structurally singular."
   call simple_sing_mat2(a)
   call gen_rhs(a, rhs, x1, x, res, 1)
   if (allocated(order)) deallocate(order)
   allocate(order(a%n))
   do i=1,a%n
     order(i) = i ! not set order(3) = 0 but code should still be ok
                  ! provided order holds permutation
   end do
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call print_result(info%flag, SSIDS_WARNING_ANAL_SINGULAR)
   options%action = .true.
   a%ne = a%ptr(a%n+1) - 1
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_FACT_SINGULAR, fs=.true.)
   call ssids_free(akeep, cuda_error)
!!!!!!!!!

   posdef = .false.
   write(*,"(a)") " * Testing factor with singular matrix......."
   call simple_sing_mat(a)
   if (allocated(order)) deallocate(order)
   allocate (order(1:a%n))
   do i = 1,a%n
     order(i) = i
   end do
   a%ne = a%ptr(a%n+1) - 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   if(info%flag < 0) then
      write(*, "(a,i4)") &
         "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   options%action = .true.
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
        SSIDS_WARNING_FACT_SINGULAR, fs=.true.)
   call ssids_free(akeep, cuda_error)

   posdef = .false.
   write(*,"(a)") " * Testing factor with match ord no scale...."
   call simple_mat(a)
   options%ordering = 2
   a%ne = a%ptr(a%n+1) - 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, val=a%val)
   if(info%flag < 0) then
      write(*, "(a,i4)") &
         "Unexpected error during analyse. flag = ", info%flag
      errors = errors + 1
      return
   endif
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
        SSIDS_WARNING_MATCH_ORD_NO_SCALE, fs=.true.)
   call ssids_free(akeep, cuda_error)
   options%ordering = 0 ! restore

   write(*,"(/a/)") " * Testing warnings (coord)"

   call simple_mat(a)
   ne = a%ne

   write(*,"(a)",advance="no") " * Testing out of range above............"
   call simple_mat(a,2)
   ne = a%ptr(a%n+1)-1
   a%ptr(a%n+1) = a%ptr(a%n+1) + 1
   a%row(ne+1) = -1
   a%col(ne+1) = 1
   a%val(ne+1) = 1.
   a%ne = ne + 1
   posdef = .true.
   call gen_rhs(a, rhs, x1, x, res, 1)
   options%ordering = 1
   call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, info)
   call print_result(info%flag,SSIDS_WARNING_IDX_OOR)
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_IDX_OOR)

   write(*,"(a)", advance="no") &
      " * Testing analyse struct singular and MC80.."
   posdef = .false.
   call simple_sing_mat(a)
   call gen_rhs(a, rhs, x1, x, res, 1)
   options%ordering = 2
   options%scaling = 3 ! scaling from Matching-based ordering
   call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, info, &
      val=a%val)
   call print_result(info%flag, SSIDS_WARNING_ANAL_SINGULAR)
   options%action = .true.
   a%ne = a%ptr(a%n+1) - 1
   call chk_answer(posdef, a, akeep, options, rhs, x, res, &
      SSIDS_WARNING_FACT_SINGULAR, fs=.true.)
   call ssids_free(akeep, cuda_error)
end subroutine test_warnings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_special
   type(matrix_type) :: a
   type(ssids_options) :: options, default_options
   type(ssids_akeep) :: akeep
   type(ssids_fkeep) :: fkeep
   type(ssids_inform) :: info

   integer :: i
   logical :: check
   logical :: posdef
   integer :: st, cuda_error
   integer :: test
   integer, dimension(:), allocatable :: order
   real(wp), dimension(:), allocatable :: scale
   real(wp), dimension(:), allocatable :: x1
   real(wp), dimension(:,:), allocatable :: rhs, x, res
   type(random_state) :: state

   integer :: big_test_n = int(1e5 + 5)

   options%unit_error = we_unit; default_options%unit_error = we_unit
   options%unit_warning = we_unit; default_options%unit_warning = we_unit
   check = .true.

   write(*,"(a)")
   write(*,"(a)") "====================="
   write(*,"(a)") "Testing special cases"
   write(*,"(a)") "====================="

   do test = 1,2
      if (test == 1) then
         write(*,"(a)",advance="no") &
            " * Testing n = 0 (CSC)..................."
         a%n = 0
         a%ne = 0
         deallocate(a%ptr,a%row,a%val,stat=st)
         allocate(a%ptr(a%n+1),a%row(a%ne),a%val(a%ne))

         if (allocated(order)) deallocate(order)
         allocate(order(a%n))
         call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
            order=order)
         if(info%flag.ne.0) then
            write(*, "(a,i4)") &
               "Unexpected error during analyse. flag = ", info%flag
            errors = errors + 1
            call ssids_free(akeep, cuda_error)
            exit
         endif
      else
         write(*,"(a)",advance="no") &
            " * Testing n = 0 (coord)................."
         a%n = 0
         a%ne = 0
         deallocate(a%row,a%val)
         allocate(a%col(a%ne),a%row(a%ne),a%val(a%ne))

         if (allocated(order)) deallocate(order)
         allocate(order(a%n))
         call ssids_analyse_coord(a%n, a%ne, a%row, a%col, akeep, options, &
            info, order=order)
         if(info%flag.ne.0) then
            write(*, "(a,i4)") &
               "Unexpected error during analyse_coord. flag = ", info%flag
            errors = errors + 1
            call ssids_free(akeep, cuda_error)
            exit
         endif
      endif

      deallocate(scale,stat=st)
      allocate(scale(a%n))
      options%scaling = 1 ! MC64

      posdef = .true.

      call ssids_factor(posdef,a%val,akeep,fkeep,options,info,scale=scale)
      if(info%flag.ne.0) then
         write(*, "(a,i4)") &
            "Unexpected error during factor. flag = ", info%flag
         errors = errors + 1
         call ssids_free(akeep, cuda_error)
         call ssids_free(fkeep, cuda_error)
         exit
      endif

      if (allocated(x1)) deallocate(x1)
      allocate(x1(a%n))
      call ssids_solve(x1,akeep,fkeep,options,info)
      if(info%flag.ne.0) then
         write(*, "(a,i4)") &
            "Unexpected error during solve. flag = ", info%flag
         errors = errors + 1
         call ssids_free(akeep, cuda_error)
         call ssids_free(fkeep, cuda_error)
         exit
      endif

      call print_result(info%flag,0)
      call ssids_free(akeep, cuda_error)
      call ssids_free(fkeep, cuda_error)
   enddo
   deallocate(order, a%ptr, a%row, a%val)
   allocate(a%ptr(big_test_n+1), a%row(4*big_test_n), a%val(4*big_test_n))
   allocate(order(10))

   ! A matrix with entry in column 1 only, explicit zeroes on diagonal
   write(*,"(a)",advance="no") &
      " * Testing zero pivot code .............."
   a%n = 10
   a%ptr(1) = 1
   a%ptr(2) = a%n + 1
   do i = 1, a%n
      a%row(i) = i
      a%val(i) = i
      order(i) = i
   end do
   do i = 2, a%n
      a%ptr(i+1) = a%ptr(i) + 1
      a%row(a%ptr(i)) = i
      a%val(a%ptr(i)) = 0.0_wp
   end do
   options%ordering = 0
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call print_result(info%flag,SSIDS_SUCCESS)
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(.false., a, akeep, options, rhs, x, &
      res, SSIDS_WARNING_FACT_SINGULAR, fs=.true.)
   call ssids_free(akeep, cuda_error)

   ! Trigger block factor-solve code with zeroes
   ! (   0.0         )
   ! ( 1e-21 1.0     )
   ! (       2.0 3.0 )
   write(*,"(a)",advance="no") &
      " * Testing zero pivot code (block)......."
   a%n = 3
   a%ptr(1:4) = (/ 1, 3, 5, 6 /)
   a%row(1:5) = (/ 1, 2, 2, 3, 3 /)
   a%val(1:5) = (/ 0d0, 1d-2*options%small, 1d0, 2d0, 3d0 /)
   order(1:3) = (/ 1, 2, 3 /)
   options%ordering = 0
   options%nemin = 16
   options%scaling = 0
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call print_result(info%flag,SSIDS_SUCCESS)
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(.false., a, akeep, options, rhs, x, &
      res, SSIDS_WARNING_FACT_SINGULAR, fs=.true.)

   ! Trigger single column factor-solve code with zeroes
   ! (   0.0         )
   ! ( 1e-21 1.0     )
   ! (       2.0 3.0 )
   write(*,"(a)",advance="no") &
      " * Testing zero pivot code (column)......"
   a%n = 3
   a%ptr(1:4) = (/ 1, 3, 5, 6 /)
   a%row(1:5) = (/ 1, 2, 2, 3, 3 /)
   a%val(1:5) = (/ 0d0, 1d-2*options%small, 1d0, 2d0, 3d0 /)
   order(1:3) = (/ 1, 2, 3 /)
   options%ordering = 0
   options%nemin = 1
   options%scaling = 0
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
      order=order)
   call print_result(info%flag,SSIDS_SUCCESS)
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(.false., a, akeep, options, rhs, x, &
      res, SSIDS_WARNING_FACT_SINGULAR, fs=.true.)

   options = default_options
   write(*,"(a)",advance="no") &
      " * Testing n>1e5, ne<3.0*n, order=1......"
   a%n = big_test_n
   call gen_random_indef(a, 2_long*big_test_n, state)
   options%ordering = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info)
   call print_result(info%flag,SSIDS_SUCCESS)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") &
      " * Testing n>1e5, ne>3.0*n, order=1......"
   a%n = big_test_n
   call gen_random_indef(a, 4_long*big_test_n, state)
   options%ordering = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info)
   call print_result(info%flag,SSIDS_SUCCESS)
   call ssids_free(akeep, cuda_error)

   write(*,"(a)",advance="no") &
      " * Testing n>1e5, ne>3.0*n, order=2......"
   a%n = big_test_n
   call gen_random_indef(a, 4_long*big_test_n, state)
   options%ordering = 2
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, val=a%val)
   call print_result(info%flag,SSIDS_SUCCESS)
   call ssids_free(akeep, cuda_error)

   ! (     1.0 )
   ! ( 1.0     )
   write(*,"(a)",advance="no") &
      " * Testing n<1e5,oxo,m1<1.8*m2,order=1..."
   a%n = 2
   a%ptr(1:3) = (/ 1, 2, 2 /)
   a%row(1:1) = (/ 2 /)
   a%val(1:1) = (/ 1.0 /)
   options%ordering = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info)
   call print_result(info%flag,SSIDS_WARNING_MISSING_DIAGONAL)
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(.false., a, akeep, options, rhs, x, &
      res, SSIDS_WARNING_MISSING_DIAGONAL, fs=.true.)
   call ssids_free(akeep, cuda_error)

   ! (  x   x  1.0 )
   ! (  x   x  2.0 )
   ! ( 1.0 2.0  x  )
   write(*,"(a)",advance="no") &
      " * Testing n<1e5,oxo,m1>1.8*m2,order=1..."
   a%n = 3
   a%ptr(1:4) = (/ 1, 2, 3, 3 /)
   a%row(1:2) = (/ 3, 3 /)
   a%val(1:2) = (/ 1.0, 2.0 /)
   options%ordering = 1
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info)
   call print_result(info%flag,SSIDS_WARNING_MISSING_DIAGONAL)
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(.false., a, akeep, options, rhs, x, &
      res, SSIDS_WARNING_FACT_SINGULAR, fs=.true.)
   call ssids_free(akeep, cuda_error)

   ! (  x   x  1.0 )
   ! (  x   x  2.0 )
   ! ( 1.0 2.0  x  )
   write(*,"(a)",advance="no") &
      " * Testing n<1e5,oxo,m1>1.8*m2,order=2..."
   a%n = 3
   a%ptr(1:4) = (/ 1, 2, 3, 3 /)
   a%row(1:2) = (/ 3, 3 /)
   a%val(1:2) = (/ 1.0, 2.0 /)
   options%ordering = 2
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, val=a%val)
   call print_result(info%flag,SSIDS_WARNING_ANAL_SINGULAR)
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(.false., a, akeep, options, rhs, x, &
      res, SSIDS_WARNING_FACT_SINGULAR, fs=.true.)
   call ssids_free(akeep, cuda_error)

   ! Test posdef with node at least 384 and multiple nodes [for coverage]
   write(*,"(a)",advance="no") &
      " * Testing n=500, posdef, BBD............"
   options = default_options
   call gen_bordered_block_diag(.true., (/ 15, 455, 10 /), 20, a%n, a%ptr, &
      a%row, a%val, state)
   call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, val=a%val)
   call print_result(info%flag,SSIDS_SUCCESS)
   call gen_rhs(a, rhs, x1, x, res, 1)
   call chk_answer(.true., a, akeep, options, rhs, x, res, SSIDS_SUCCESS)
   call ssids_free(akeep, cuda_error)

end subroutine test_special

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Generates a bordered block diagonal form
! blocks have size and number given in dimn
! border is width of border - however final var is only included in final block
! to discourage merging with largest block during analysis
subroutine gen_bordered_block_diag(posdef, dimn, border, n, ptr, row, val, state)
   logical, intent(in) :: posdef
   integer, dimension(:), intent(in) :: dimn
   integer, intent(in) :: border
   integer, intent(out) :: n
   integer, dimension(:), allocatable :: ptr
   integer, dimension(:), allocatable :: row
   real(wp), dimension(:), allocatable :: val
   type(random_state), intent(inout) :: state

   integer :: i, j, k, blk, blk_sa
   integer :: nnz
   integer :: st

   ! Clear any previous allocs
   deallocate(ptr, stat=st)
   deallocate(row, stat=st)
   deallocate(val, stat=st)

   ! allocate arrays
   n = sum(dimn(:)) + border
   nnz = 0
   do blk = 1, size(dimn)
      j = dimn(blk)
      nnz = nnz + j*(j+1)/2 + j*border
   end do
   nnz = nnz + border*(border+1)/2
   allocate(ptr(n+1), row(nnz), val(nnz))

   ! Generate val = unif(-1,1)
   do i = 1, nnz
      val(i) = random_real(state)
   end do

   ! Generate ptr and row; make posdef if required
   j = 1
   blk_sa = 1
   do blk = 1, size(dimn)
      do i = blk_sa, blk_sa+dimn(blk)-1
         ptr(i) = j
         if(posdef) val(j) = abs(val(j)) + n ! make diagonally dominant
         do k = i, blk_sa+dimn(blk)-1
            row(j) = k
            j = j + 1
         end do
         do k = n-border+1, n-1
            row(j) = k
            j = j + 1
         end do
      end do
      blk_sa = blk_sa + dimn(blk)
   end do
   do i = n-border+1, n
      ptr(i) = j
      if(posdef) val(j) = abs(val(j)) + n ! make diagonally dominant
      do k = i, n
         row(j) = k
         j = j + 1
      end do
   end do
   ptr(n+1) = j
end subroutine gen_bordered_block_diag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_mat_lower(a,extra)
   ! simple pos def test matrix (lower triangular part only)
   type(matrix_type), intent(inout) :: a
   integer, optional, intent(in) :: extra

   integer :: myextra,st

   myextra = 0
   if(present(extra)) myextra = extra

   !
   ! Create the simple sparse matrix (lower triangular part only):
   !
   ! 10.0  2.0       3.0
   !  2.0 10.0
   !           10.0  4.0
   !  3.0       4.0 10.0
   !

   a%n = 4
   a%ne = 7
   deallocate(a%ptr, a%row, a%col, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 4, 5, 7, 8 /)
   allocate(a%col(a%ne+myextra))
   allocate(a%row(a%ne+myextra))
   allocate(a%val(a%ne+myextra))

   a%col(1:7) = (/ 1, 1, 1,     2,    3, 3,    4 /)
   a%row(1:7) = (/ 1, 2, 4,     2,    3, 4,    4 /)
   a%val(1:7) = (/   10.0, 2.0, 3.0, &
                     10.0, &
                     10.0, 4.0, &
                     10.0 /)

end subroutine simple_mat_lower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_mat_zero_diag(a)
   type(matrix_type), intent(inout) :: a

   integer :: st

   !
   ! Create the simple sparse matrix:
   !
   !  0.0  1.0  2.0
   !  1.0  0.0  1.0
   !  2.0  1.0  0.0

   a%n = 3
   a%ne = 3
   deallocate(a%ptr, a%row, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 3, 4, 4 /)
   allocate(a%row((a%ptr(a%n+1)-1+a%n)))
   allocate(a%val((a%ptr(a%n+1)-1+a%n)))

   a%row(1:3) = (/ 2,3,  3 /)
   a%val(1:3) = (/   1.0, 2.0, &
                     1.0 /)

end subroutine simple_mat_zero_diag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_sing_mat(a)
   type(matrix_type), intent(inout) :: a
   integer :: st

   !
   ! Create the simple singular sparse matrix:
   !
   !  0.0  2.0
   !  2.0  0.0  1.0
   !       1.0  0.0
   !
   ! we will not enter diagonal entries explicitly

   a%n = 3
   a%ne = 2
   deallocate(a%ptr, a%row, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 2, 3, 3 /)
   allocate(a%row(9))
   allocate(a%val(9))
   a%row(1:2) = (/ 2, 3 /)
   a%val(1:2) = (/   2.0, 1.0 /)

end subroutine simple_sing_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_sing_mat2(a)
   type(matrix_type), intent(inout) :: a
   integer :: st

   !
   ! Create the simple singular sparse matrix:
   !
   ! -1.0  2.0  0.0
   !  2.0  1.0  0.0
   !  0.0  0.0  0.0
   !
   ! by entering null column.

   a%n = 3
   a%ne = 3
   deallocate(a%ptr, a%row, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 3, 4, 4 /)
   allocate(a%row(9))
   allocate(a%val(9))
   a%row(1:3) = (/ 1, 2, 2 /)
   a%val(1:3) = (/  -1.0, 2.0, 1.0 /)

end subroutine simple_sing_mat2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_mat(a,extra)
   ! simple pos def test matrix (lower triangular part)
   type(matrix_type), intent(inout) :: a
   integer, optional, intent(in) :: extra

   integer :: myextra,st

   myextra = 0
   if(present(extra)) myextra = extra

   !
   ! Create the simple sparse matrix (lower and upper triangles):
   !
   ! 10.0  2.0       3.0
   !  2.0 10.0
   !           10.0  4.0
   !  3.0       4.0 10.0
   !

   a%n = 4
   a%ne = 7
   deallocate(a%ptr, a%row, a%col, a%val, stat=st)
   allocate(a%ptr(a%n+1))
   a%ptr = (/ 1, 4, 5, 7, 8 /)
   allocate(a%row(a%ne+myextra))
   allocate(a%col(a%ne+myextra))
   allocate(a%val(a%ne+myextra))

   a%row(1:7) = (/ 1,2,4,     2,    3, 4,  4 /)
   a%col(1:7) = (/ 1,1,1,     2,    3, 3,  4 /)
   a%val(1:7) = (/  10.0,  2.0, 3.0, &
                     10.0,      &
                     10.0,  4.0,      &
                      10.0 /)

end subroutine simple_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chk_answer(posdef, a, akeep, options, rhs, x, res, &
      expected_flag, fs)
   logical, intent(in) :: posdef
   type(matrix_type), intent(inout) :: a
   type(ssids_akeep), intent(inout) :: akeep
   type(ssids_options), intent(in) :: options
   real(wp), dimension(:,:), intent(inout) :: rhs
   real(wp), dimension(:,:), intent(inout) :: x
   real(wp), dimension(:,:), intent(inout) :: res
   integer, intent(in) :: expected_flag
   logical, optional, intent(in) :: fs

   type(ssids_fkeep) :: fkeep
   type(ssids_inform) :: info
   integer :: nrhs, cuda_error
   type(ssids_options) :: myoptions

   myoptions = options
   myoptions%unit_warning = -1 ! disable printing warnings

   write(*,"(a)",advance="no") " *    checking answer...................."

   nrhs = 1
   !if(present(fs) .and. .false.) then
   !   call ssids_factor_solve(posdef,a%val,nrhs,x,a%n,akeep,fkeep,myoptions, &
   !      info)
   !   if(info%flag .ne. expected_flag) then
   !      write(*, "(a,2i4)") "fail on factor_solve",info%flag,expected_flag
   !      errors = errors + 1
   !      go to 99
   !   endif
   !else
      call ssids_factor(posdef,a%val,akeep,fkeep,myoptions,info)
      if(info%flag .ne. expected_flag) then
         write(*, "(a,2i4)") "fail on factor",info%flag,expected_flag
         errors = errors + 1
         go to 99
      endif

      call ssids_solve(nrhs,x,a%n,akeep,fkeep,myoptions,info)
      if(info%flag .ne. expected_flag) then
         write(*, "(a,2i4)") "fail on solve", info%flag,expected_flag
         errors = errors + 1
         go to 99
      endif
   !endif

   ! Check residual
   call compute_resid(nrhs,a,x,a%n,rhs,a%n,res,a%n)

   if(maxval(abs(res(1:a%n,1:nrhs))) < err_tol) then
      write(*, "(a)") "ok"
   else
      write(*, "(a,es12.4)") "fail residual = ", &
         maxval(abs(res(1:a%n,1:nrhs)))
      errors = errors + 1
   endif

   ! remember: must call finalise
   99 call ssids_free(fkeep, cuda_error)

end subroutine chk_answer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generate rhs and copy into x
subroutine gen_rhs(a, rhs, x1, x, res, nrhs, state)
   type(matrix_type), intent(inout) :: a
   real(wp), dimension(:,:), allocatable, intent(inout) :: rhs
   real(wp), dimension(:,:), allocatable, intent(inout) :: x
   real(wp), dimension(:), allocatable, intent(inout) :: x1
   real(wp), dimension(:,:), allocatable, intent(inout) :: res
   integer, intent(in) :: nrhs
   type (random_state), optional :: state

   integer :: i, j, k, n
   real(wp) :: atemp, random

   n = a%n
   if (allocated(rhs)) deallocate(rhs)
   allocate(rhs(n, nrhs))

   if (allocated(x)) deallocate(x)
   allocate(x(n, nrhs))

   if (allocated(res)) deallocate(res)
   allocate(res(n, nrhs))

   if (allocated(x1)) deallocate(x1)
   allocate(x1(n))

   rhs = zero
   if (.not.present(state)) then
      ! Generate rhs assuming x = 1
      do k = 1, n
         do j = a%ptr(k), a%ptr(k+1)-1
            i = a%row(j)
            if(i < k .or. i > n) cycle
            rhs(i, 1:nrhs) = rhs(i, 1:nrhs) + a%val(j)
            if(i.eq.k) cycle
            rhs(k, 1:nrhs) = rhs(k, 1:nrhs) + a%val(j)
        end do
      end do
   else
      do i = 1,n
         random = random_real(state,.true.)
         x(i,1:nrhs) = random
      end do
      do k = 1, n
         do j = a%ptr(k), a%ptr(k+1)-1
            i = a%row(j)
            if(i < k .or. i > n) cycle
            atemp = a%val(j)
            rhs(i, 1:nrhs) = rhs(i, 1:nrhs) + atemp*x(k,1:nrhs)
            if(i.eq.k) cycle
            rhs(k, 1:nrhs) = rhs(k, 1:nrhs) + atemp*x(i,1:nrhs)
        end do
      end do
   endif

   x(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
   x1(1:a%n) = rhs(1:a%n,1)

end subroutine gen_rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_resid(nrhs,a,x,lx,rhs,lrhs,res,lres)
   integer, intent(in) :: nrhs, lrhs, lx, lres
   type(matrix_type), intent(in) :: a
   real(wp), intent(in) :: rhs(lrhs,nrhs)
   real(wp), intent(in) :: x(lx,nrhs)
   real(wp), intent(out) :: res(lres,nrhs)

   real(wp), dimension(:), allocatable :: work

   integer :: i, j, k
   real(wp) :: anorm, atemp, bnorm(1:nrhs), xnorm(1:nrhs)

   allocate(work(a%n))

   anorm = 0
   bnorm = 0
   xnorm = 0
   work = 0

   ! Check residual
   res(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
   do k = 1, a%n
      do j = a%ptr(k), a%ptr(k+1)-1
         i = a%row(j)
         if (i < 1 .or. i > a%n) cycle
         atemp = a%val(j)
         res(i, 1:nrhs) = res(i, 1:nrhs) - atemp*x(k,1:nrhs)
         work(i) = work(i) + abs(atemp)
         if(i.eq.k) cycle
         res(k, 1:nrhs) = res(k, 1:nrhs) - atemp*x(i,1:nrhs)
         work(k) = work(k) + abs(atemp)
      end do
   end do

   do k = 1, a%n
      anorm = max(anorm,work(k))
      do i = 1,nrhs
         bnorm(i) = max(bnorm(i),abs(rhs(k,i)))
         xnorm(i) = max(xnorm(i),abs(x(k,i)))
      end do
   end do

   do k = 1,a%n
      do i = 1,nrhs
         res(k,i) = res(k,i)/(anorm*xnorm(i) + bnorm(i))
      end do
   end do

end subroutine compute_resid

real(wp) function inf_norm(a)
   type(matrix_type), intent(in) :: a

   real(wp), dimension(:), allocatable :: work

   integer :: i, j, k
   real(wp) :: atemp

   allocate(work(a%n))

   inf_norm = 0
   work(:) = 0

   ! Check residual
   do k = 1, a%n
      do j = a%ptr(k), a%ptr(k+1)-1
         i = a%row(j)
         if (i < 1 .or. i > a%n) cycle
         atemp = a%val(j)
         work(i) = work(i) + abs(atemp)
         if(i.eq.k) cycle
         work(k) = work(k) + abs(atemp)
      end do
   end do

   do k = 1, a%n
      inf_norm = max(inf_norm,work(k))
   end do

end function inf_norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_random
   type(ssids_akeep) :: akeep
   type(ssids_fkeep) :: fkeep
   type(ssids_options) :: options
   type(ssids_inform) :: info

   !logical, parameter :: debug = .true.
   logical, parameter :: debug = .false.

   integer :: maxn = 1000
   integer :: maxnemin = 48
   integer :: maxnz =  1000000
   integer, parameter :: maxnrhs = 10
   integer, parameter :: nprob = 100
   type(random_state) :: state

   type(matrix_type) :: a
   real(wp), allocatable, dimension(:, :) :: rhs,x
   real(wp), allocatable, dimension(:) :: rhs1d,x1
   real(wp), allocatable, dimension(:, :) :: res

   logical :: posdef
   integer :: prblm, i, j, k, n1, nrhs, mt
   integer(long) :: ne, nza
   integer, dimension(:), allocatable :: order, piv_order
   integer, dimension(:), allocatable :: xindex, bindex
   logical, dimension(:), allocatable :: lflag
   real(wp), dimension(:,:), allocatable :: d
   real(wp), dimension(:), allocatable :: d1
   integer :: cuda_error
   logical :: check, coord
   real(wp) :: num_flops
   integer :: max_threads
   type(numa_region), dimension(:), allocatable :: fake_topology

   max_threads = 1
!$ max_threads = omp_get_max_threads()
   if(max_threads > 1) then
      allocate(fake_topology(2))
      do i = 1, 2
         fake_topology(i)%nproc = 2
         allocate(fake_topology(i)%gpus(0))
      end do
   else
      allocate(fake_topology(1))
      do i = 1, 1
         fake_topology(i)%nproc = 1
         allocate(fake_topology(i)%gpus(0))
      end do
   end if

   options%unit_error = we_unit
   options%unit_diagnostics = dl_unit

   write(*, "(a)")
   write(*, "(a)") "======================="
   write(*, "(a)") "Testing random matrices"
   write(*, "(a)") "======================="

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz), a%col(2*maxnz))
   allocate(order(maxn))
   allocate(piv_order(maxn))
   allocate(rhs(maxn,maxnrhs), res(maxn,maxnrhs), x(maxn,maxnrhs))
   allocate(x1(maxn), rhs1d(maxn))
   allocate(d(2,maxn), d1(maxn))
   allocate(bindex(maxn), xindex(maxn), lflag(maxn))

   options%multiplier = 1.0 ! Ensure we give reallocation a work out

   if(debug) options%print_level = 10000
   options%min_gpu_work = 0 ! alway allow some gpu working

   do prblm = 1, nprob
      !if(errors>0) stop
      !call random_set_seed(state, 1963535693)
      !print *, random_get_seed(state)

      ! decide whether pos. def. problem
      posdef = .false.
      mt = SPRAL_MATRIX_REAL_SYM_INDEF
      n1 = random_integer(state,maxn)
      if ((n1/2)*2 == n1) then
         posdef = .true.
         mt = SPRAL_MATRIX_REAL_SYM_PSDEF
      endif

      ! Generate parameters
      a%n = random_integer(state, maxn)
      if (prblm < 21) a%n = prblm ! check very small problems
      i = a%n**2/2 - a%n
      i = max(0,i)
      nza = random_integer(state, i)
      nza = nza + a%n

      options%nemin = random_integer(state,  maxnemin)
      options%nemin = 1 ! FIXME: remove

      if(nza.gt.maxnz .or. a%n.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5)") "n = ", a%n, " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      if(posdef) then
         call gen_random_posdef(a, nza, state)
      else
         call gen_random_indef(a, nza, state)
      endif
      if(debug) call print_matrix(6, -1, mt, a%n, a%n, a%ptr, a%row, a%val)

      options%ordering = random_integer(state, 2) - 1 ! user or metis

      if(random_logical(state)) then
         options%nstream = random_integer(state, 4)
      else
         options%nstream = 1
      endif

      if(options%ordering.eq.0) then
         ! Generate a pivot order
         call simple_metis_order(a, order)
      endif

      check = .false.
      n1 = random_integer(state, maxn)
      if ((n1/2)*2 == n1) check = .true.

      ! Peform analyse
      n1 = random_integer(state, maxn)

      if ( posdef ) then
         write(*, "(a, i3, a, i5, a, i7, a, i2, a)",advance="no") &
            " + no. ", prblm,  " n = ", a%n, " nza = ", nza, "..."
      else
         write(*, "(a, i3, a, i5, a, i7, a, i2, a)",advance="no") &
            " - no. ", prblm,  " n = ", a%n, " nza = ", nza, "..."
      end if

      if ((n1/3)*3 == n1) then
         ne = a%ptr(a%n+1) - 1
         if (allocated(a%col)) deallocate(a%col)
         allocate (a%col(ne))
         do i = 1,a%n
           do j = a%ptr(i),a%ptr(i+1)-1
              a%col(j) = i
           end do
         end do

         if(random_logical(state)) then
            ! Use fake_topology
            call ssids_analyse_coord(a%n, ne, a%row, a%col, &
                 akeep, options, info, order, topology=fake_topology)
         else
            ! Use real topology
            call ssids_analyse_coord(a%n, ne, a%row, a%col, &
                 akeep, options, info, order)
         endif
         coord = .true.
      else
         call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
            order)
         coord = .false.
      endif
      num_flops = info%num_flops
      if(info%flag .ne. SSIDS_SUCCESS) then
         write(*, "(a,i3)") "fail on analyse", info%flag
         call ssids_free(akeep, cuda_error)
         errors = errors + 1
         cycle
      endif

      nrhs = random_integer(state,  maxnrhs)

      ! Generate rhs assuming x(k) = k/maxn. remember we have only
      ! half matrix held.
      rhs(1:a%n, 1:nrhs) = zero
      do k = 1, a%n
         do j = a%ptr(k), a%ptr(k+1)-1
            i = a%row(j)
            rhs(i, 1:nrhs) = rhs(i, 1:nrhs) + a%val(j)*real(k)/real(maxn)
            if(i.eq.k) cycle
            rhs(k, 1:nrhs) = rhs(k, 1:nrhs) + a%val(j)*real(i)/real(maxn)
         end do
      end do
      rhs1d(1:a%n) = rhs(1:a%n, 1)

      !
      ! Perform straightforward factor then solve
      !

      ! set right-hand side
      x(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
      x1(1:a%n) = rhs1d(1:a%n)
      if (coord) then
         call ssids_factor(posdef,a%val, akeep, fkeep, options, info)
      else
         call ssids_factor(posdef,a%val, akeep, fkeep, options, info,  &
               ptr=a%ptr, row=a%row)
      endif

      if(info%flag .lt. SSIDS_SUCCESS) then
         write(*, "(a,i3)") "fail on factor", info%flag
         call ssids_free(akeep, fkeep, cuda_error)
         errors = errors + 1
         cycle
      endif
      write(*,'(a,f6.1,1x)',advance="no") ' num_flops:',num_flops*1e-6

      ! Perform solve
      if(posdef) then
         select case(mod(prblm,2))
         case(0)
            ! Vanilla solve
            call ssids_solve(x1, akeep, fkeep, options, info)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job absent", &
                  info%flag
               call ssids_free(akeep, fkeep, cuda_error)
               errors = errors + 1
               cycle
            endif
            call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job absent", &
                  info%flag
               call ssids_free(akeep, fkeep, cuda_error)
               errors = errors + 1
               cycle
            endif
         case(1)
            ! Fwd, then bwd solves
            call ssids_solve(x1, akeep, fkeep, options, info, job=1)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job = 1", &
                  info%flag
               call ssids_free(akeep, fkeep, cuda_error)
               errors = errors + 1
               cycle
            endif
            call ssids_solve(x1, akeep, fkeep, options, info, job=3)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job = 3", &
                  info%flag
               call ssids_free(akeep, fkeep, cuda_error)
               errors = errors + 1
               cycle
            endif
            call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=1)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job = 1", &
                  info%flag
               errors = errors + 1
               cycle
            endif
            call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=3)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job = 3", &
                  info%flag
               call ssids_free(akeep, fkeep, cuda_error)
               errors = errors + 1
               cycle
            endif
         end select
      else
         ! FIXME: Restore separate cases
         !select case(mod(prblm, 3))
         !case(0)
            ! Vanilla solve
            call ssids_solve(x1, akeep, fkeep, options, info)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job absent", &
                  info%flag
               call ssids_free(akeep, fkeep, cuda_error)
               errors = errors + 1
               cycle
            endif
            ! FIXME: restore multirhs
            !call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info)
            !if(info%flag .lt. SSIDS_SUCCESS) then
            !   write(*, "(a,i4)") " fail on 1d solve with job absent", &
            !      info%flag
            !   call ssids_free(akeep, fkeep, cuda_error)
            !   errors = errors + 1
            !   cycle
            !endif
         !case(1)
         !   ! Fwd, then (diag+bwd) solves
         !   call ssids_solve(x1, akeep, fkeep, options, info, job=1)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 1", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(x1, akeep, fkeep, options, info, job=4)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 4", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=1)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 1", &
         !         info%flag
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=4)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 4", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !case(2)
         !   ! Fwd, then diag, then bwd solves
         !   call ssids_solve(x1, akeep, fkeep, options, info, job=1)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 1", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(x1, akeep, fkeep, options, info, job=2)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 2", &
         !         info%flag
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(x1, akeep, fkeep, options, info, job=3)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 3", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=1)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 1", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=2)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 2", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !   call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=3)
         !   if(info%flag .lt. SSIDS_SUCCESS) then
         !      write(*, "(a,i4)") " fail on 1d solve with job = 3", &
         !         info%flag
         !      call ssids_free(akeep, fkeep, cuda_error)
         !      errors = errors + 1
         !      cycle
         !   endif
         !end select
      endif

      if(prblm.eq.nprob) then
         options%print_level = 2
         options%unit_diagnostics = dl_unit
      endif

      !print *, "x = "
      !write (6,'(6es12.4)') x(1:a%n,1:1)

      ! Check residuals
      call compute_resid(1,a,x1,maxn,rhs1d,maxn,res,maxn)
      if(maxval(abs(res(1:a%n,1))) < err_tol) then
         !write(*, "(a)", advance="no") "ok..." ! FIXME: restore
         write(*, "(a)") "ok..."
      else
         write(*, "(a,es12.4)") " f+s fail residual 1d = ", &
            maxval(abs(res(1:a%n,1)))
       !  write (6,'(6es12.4)') x1(1:a%n)
         errors = errors + 1
         cycle
      endif
      ! FIXME: restore multirhs
      !!call compute_resid(nrhs,a,x,maxn,rhs,maxn,res,maxn)
      !if(maxval(abs(res(1:a%n,1:nrhs))) < err_tol) then
      !   write(*, "(a)") "ok..."
      !else
      !   write(*, "(a)") " f+s fail residual 2d = "
      !   do i = 1, nrhs
      !      write(*, "(es12.4)") maxval(abs(res(1:a%n,i)))
      !   end do
      !   write(*, "()")
      ! !  write (6,'(6es12.4)') x1(1:a%n)
      !   errors = errors + 1
      !   call ssids_free(akeep, fkeep, cuda_error)
      !   cycle
      !endif

      if (.not.posdef) then
         call ssids_enquire_indef(akeep,fkeep,options,info,piv_order=piv_order, &
            d=d)
         if (info%flag < 0) then
           write (*,'(a)') ' Unexpected error from ssids_enquire_indef'
           call ssids_free(akeep, fkeep, cuda_error)
           errors = errors + 1
           cycle
         endif

         call ssids_alter(d,akeep,fkeep,options,info)
         if (info%flag < 0) then
            write (*,'(a)') ' Unexpected error from ssids_alter'
            call ssids_free(akeep, fkeep, cuda_error)
            errors = errors + 1
            cycle
         endif
      else
         call ssids_enquire_posdef(akeep,fkeep,options,info,d1)
         if (info%flag < 0 .and. info%flag.ne.SSIDS_ERROR_UNIMPLEMENTED) then
           ! NB: unknown error return as not yet implemented with gpu factors
           write (*,'(a)') ' Unexpected error from ssids_enquire_posdef'
           call ssids_free(akeep, fkeep, cuda_error)
           errors = errors + 1
           cycle
         endif
      endif

      !
      ! Cleanup ready for next iteration
      !
      if(random_logical(state)) &
         call ssids_free(akeep, fkeep, cuda_error)

   end do

   ! Ensure cleaned up at end
   call ssids_free(akeep, fkeep, cuda_error)

end subroutine test_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_big
   type(ssids_akeep) :: akeep
   type(ssids_fkeep) :: fkeep
   type(ssids_options) :: options
   type(ssids_inform) :: info

   type(random_state) :: state
   type(matrix_type) :: a
   real(wp), allocatable, dimension(:, :) :: rhs,x
   real(wp), allocatable, dimension(:, :) :: res

   logical :: posdef
   integer :: i, j, k, nrhs, cuda_error
   real(wp) :: num_flops

   write(*, "(a)")
   write(*, "(a)") "=================="
   write(*, "(a)") "Testing big matrix"
   write(*, "(a)") "=================="

   a%n = 2000
   a%ne = 5*a%n
   nrhs = 3

   allocate(a%ptr(a%n+1))
   allocate(a%row(2*a%ne), a%val(2*a%ne), a%col(2*a%ne))
   allocate(rhs(a%n,nrhs), res(a%n,nrhs), x(a%n,nrhs))

   posdef = .false.

   write(*, "(a, i9, a, i11, a, i2, a)",advance="no") &
      " * n = ", a%n, " nza = ", a%ne, "..."

   if(posdef) then
      call gen_random_posdef(a, a%ne, state)
   else
      call gen_random_indef(a, a%ne, state)
   endif

   ! Peform analyse
   call ssids_analyse(.false., a%n, a%ptr, a%row, akeep, options, info)
   if(info%flag .ne. SSIDS_SUCCESS) then
      write(*, "(a,i3)") "fail on analyse", info%flag
      call ssids_free(akeep, cuda_error)
      errors = errors + 1
      return
   endif
   num_flops = info%num_flops

   ! Generate rhs assuming x(k) = k/maxn. remember we have only
   ! half matrix held.
   rhs(1:a%n, 1:nrhs) = zero
   do k = 1, a%n
      do j = a%ptr(k), a%ptr(k+1)-1
         i = a%row(j)
         rhs(i, 1:nrhs) = rhs(i, 1:nrhs) + a%val(j)*real(k)/real(a%n)
         if(i.eq.k) cycle
         rhs(k, 1:nrhs) = rhs(k, 1:nrhs) + a%val(j)*real(i)/real(a%n)
      end do
   end do

   !
   ! Perform straightforward factor then solve
   !

   ! set right-hand side
   x(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
   call ssids_factor(posdef,a%val, akeep, fkeep, options, info,  &
         ptr=a%ptr, row=a%row)
   if(info%flag .lt. SSIDS_SUCCESS) then
      write(*, "(a,i3)") "fail on factor", info%flag
      call ssids_free(akeep, fkeep, cuda_error)
      errors = errors + 1
      return
   endif
   write(*,'(a,f6.1,1x)',advance="no") ' num_flops:',num_flops*1e-6

   ! Perform solve
   call ssids_solve(nrhs, x, a%n, akeep, fkeep, options, info)
   if(info%flag .lt. SSIDS_SUCCESS) then
      write(*, "(a,i4)") " fail on 1d solve with job absent", &
         info%flag
      call ssids_free(akeep, fkeep, cuda_error)
      errors = errors + 1
      return
   endif

   ! Check residuals
   call compute_resid(nrhs,a,x,a%n,rhs,a%n,res,a%n)
   if(maxval(abs(res(1:a%n,1:nrhs))) < err_tol) then
      write(*, "(a)", advance="no") "ok..."
   else
      write(*, "(a)") " f+s fail residual 2d = "
      do i = 1, nrhs
         write(*, "(es12.4)", advance="no") maxval(abs(res(1:a%n,i)))
      end do
      write(*, "()")
      errors = errors + 1
      call ssids_free(akeep, fkeep, cuda_error)
      return
   endif

   !
   ! Cleanup ready for next iteration
   !
   call ssids_free(akeep, fkeep, cuda_error)

end subroutine test_big

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_random_scale
   type(ssids_akeep) :: akeep
   type(ssids_fkeep) :: fkeep
   type(ssids_options) :: options
   type(ssids_inform) :: info

   integer :: maxn = 500
   integer :: maxnemin = 48
   integer :: maxnz =  1000000
   integer, parameter :: maxnrhs = 10
   integer, parameter :: nprob = 100
   type(random_state) :: state

   type(matrix_type) :: a
   real(wp), allocatable, dimension(:, :) :: rhs,x
   real(wp), allocatable, dimension(:) :: rhs1d,x1
   real(wp), allocatable, dimension(:, :) :: res
   real(wp), allocatable, dimension(:) :: scale
   integer, dimension(:), allocatable :: xindex, bindex
   logical, dimension(:), allocatable :: lflag

   logical :: posdef
   integer :: prblm, i, j, k, n1, nrhs
   integer(long) :: ne, nza
   integer, dimension(:), allocatable :: order
   logical :: check
   real(wp) :: num_flops

   type(hungarian_options) :: hoptions
   type(hungarian_inform) :: hinform

   integer :: cuda_error

   write(*, "(a)")
   write(*, "(a)") "================================"
   write(*, "(a)") "Testing random matrices (scaled)"
   write(*, "(a)") "================================"

   allocate(a%ptr(maxn+1))
   allocate(a%row(2*maxnz), a%val(2*maxnz), a%col(2*maxnz))
   allocate(order(maxn))
   allocate(rhs(maxn,maxnrhs), res(maxn,maxnrhs), x(maxn,maxnrhs))
   allocate(x1(maxn),rhs1d(maxn),scale(2*maxn))
   allocate(bindex(maxn), xindex(maxn), lflag(maxn))

   options%action = .false.
   options%nstream = 2

   do prblm = 1, nprob
      !call random_set_seed(state, 1221086619)
      !print *, random_get_seed(state)

      ! Generate parameters
      posdef = random_logical(state)
      a%n = random_integer(state,  maxn)
      if (prblm < 21) a%n = prblm ! check very small problems
      i = a%n**2/2 - a%n
      i = max(0,i)
      nza = random_integer(state,  i)
      nza = nza + a%n

      if(nza.gt.maxnz .or. a%n.gt.maxn) then
         write(*, "(a)") "bad random matrix."
         write(*, "(a,i5,a,i5)") "n = ", a%n, " > maxn = ", maxn
         write(*, "(a,i8,a,i8)") "or nza = ", nza, " > maxnz = ", maxnz
         cycle
      endif

      if(posdef) then
         call gen_random_posdef(a, nza, state)
      else
         call gen_random_indef(a, nza, state)
      endif

      ! Generate a pivot order
      call simple_metis_order(a, order)

      options%nemin = random_integer(state,  maxnemin)

      ! Pick a scaling that isn't 3 (ie doesn't require matching-based ordering)
      n1 = random_integer(state,  4)
      select case(n1-1)
      case(0:2)
         options%scaling = n1-1
      case(3)
         options%scaling = 4
      end select

      options%ordering = 1 ! METIS
      n1 = random_integer(state, maxn)
      if ((n1/2)*2 == n1) options%ordering = 0

      if(random_logical(state)) then
         options%nstream = random_integer(state, 4)
      else
         options%nstream = 1
      endif

      n1 = random_integer(state, maxn)
      if ((n1/3)*3 == n1) options%ordering = 2
      if (options%ordering.eq.2) options%scaling = 3

      check = .false.
      n1 = random_integer(state, maxn)
      if ((n1/2)*2 == n1) check = .true.

      write(*, "(a, i3, a, i5, a, i7, a, i2, a)",advance="no") &
         " * no. ", prblm,  " n = ", a%n, " nza = ", nza, " scal = ", &
         options%scaling, "..."

      ! Peform analyse
      if ((n1/3)*3 == n1) then
         ne = a%ptr(a%n+1) - 1
         if (allocated(a%col)) deallocate(a%col)
         allocate (a%col(ne))
         do i = 1,a%n
           do j = a%ptr(i),a%ptr(i+1)-1
              a%col(j) = i
           end do
         end do

         if (options%ordering.le.1) then
            call ssids_analyse_coord(a%n, ne, a%row, a%col, &
               akeep, options, info, order)
         else
            call ssids_analyse_coord(a%n, ne, a%row, a%col, &
               akeep, options, info, order, val=a%val)
         endif
      else
         if (options%ordering.le.1) then
            call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
               order)
         else
            call ssids_analyse(check, a%n, a%ptr, a%row, akeep, options, info, &
               order, val=a%val)
         endif
      endif
      num_flops = info%num_flops
      if(info%flag .lt. SSIDS_SUCCESS) then
         write(*, "(a,i3)") "fail on analyse", info%flag
         errors = errors + 1
         cycle
      endif

      nrhs = random_integer(state,  maxnrhs)

      if(prblm.eq.nprob) then
         options%print_level = 2
         options%unit_diagnostics = dl_unit
         nrhs = 3
      endif

      ! Generate rhs assuming x(k) = k/maxn. remember we have only
      ! half matrix held.
      rhs(1:a%n, 1:nrhs) = zero
      do k = 1, a%n
         do j = a%ptr(k), a%ptr(k+1)-1
            i = a%row(j)
            rhs(i, 1:nrhs) = rhs(i, 1:nrhs) + a%val(j)*real(k)/real(maxn)
            if(i.eq.k) cycle
            rhs(k, 1:nrhs) = rhs(k, 1:nrhs) + a%val(j)*real(i)/real(maxn)
         end do
      end do
      rhs1d(1:a%n) = rhs(1:a%n, 1)

      !
      ! Perform straight forward factor then solve
      !

      ! set right-hand side
      x(1:a%n,1:nrhs) = rhs(1:a%n,1:nrhs)
      x1(1:a%n) = rhs1d(1:a%n)

      if (options%scaling.eq.0) then
         call hungarian_scale_sym(a%n, a%ptr, a%row, a%val, scale, &
            hoptions, hinform)
         call ssids_factor(posdef, a%val, akeep, fkeep, options, info,  &
            ptr=a%ptr, row=a%row, scale=scale)
      else
         call ssids_factor(posdef, a%val, akeep, fkeep, options, info,  &
            ptr=a%ptr, row=a%row, scale=scale)
      endif

      if(info%flag .lt. SSIDS_SUCCESS) then
         if(info%flag .eq. SSIDS_ERROR_SINGULAR) then
            write(*, "(a)") "singular"
            cycle
         endif
         write(*, "(a,i3)") "fail on factor", info%flag
         errors = errors + 1
         cycle
      endif
      write(*,'(a,f6.1,1x)',advance="no") ' num_flops:',num_flops*1e-6

      ! Perform solve
      select case(mod(prblm,3))
      case(0)
         ! Vanilla solve
         call ssids_solve(x1, akeep, fkeep, options, info)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job absent", &
               info%flag
            errors = errors + 1
            cycle
         endif
         call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job absent", &
               info%flag
            errors = errors + 1
            cycle
         endif
      case(1)
         ! Fwd, then (diag+bwd) solves
         call ssids_solve(x1, akeep, fkeep, options, info, job=1)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 1", &
               info%flag
            errors = errors + 1
            cycle
         endif
         if(posdef) then
            call ssids_solve(x1, akeep, fkeep, options, info, job=3)
         else
            call ssids_solve(x1, akeep, fkeep, options, info, job=4)
         endif
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 4", &
               info%flag
            errors = errors + 1
            cycle
         endif
         call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=1)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 1", &
               info%flag
            errors = errors + 1
            cycle
         endif
         if(posdef) then
            call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=3)
         else
            call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=4)
         endif
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 4", &
               info%flag
            errors = errors + 1
            cycle
         endif
      case(2)
         ! Fwd, then diag, then bwd solves
         call ssids_solve(x1, akeep, fkeep, options, info, job=1)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 1", &
               info%flag
            errors = errors + 1
            cycle
         endif
         if(.not.posdef) then
            call ssids_solve(x1, akeep, fkeep, options, info, job=2)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job = 2", &
                  info%flag
               errors = errors + 1
               cycle
            endif
         endif
         call ssids_solve(x1, akeep, fkeep, options, info, job=3)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 3", &
               info%flag
            errors = errors + 1
            cycle
         endif
         call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=1)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 1", &
               info%flag
            errors = errors + 1
            cycle
         endif
         if(.not.posdef) then
            call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=2)
            if(info%flag .lt. SSIDS_SUCCESS) then
               write(*, "(a,i4)") " fail on 1d solve with job = 2", &
                  info%flag
               errors = errors + 1
               cycle
            endif
         endif
         call ssids_solve(nrhs, x, maxn, akeep, fkeep, options, info, job=3)
         if(info%flag .lt. SSIDS_SUCCESS) then
            write(*, "(a,i4)") " fail on 1d solve with job = 3", &
               info%flag
            errors = errors + 1
            cycle
         endif
      end select

      ! write (6,'(6es12.4)') x(1:a%n,1:1)

      ! Check residuals
      call compute_resid(1,a,x1,maxn,rhs1d,maxn,res,maxn)
      if(maxval(abs(res(1:a%n,1))) < err_tol_scale) then
         write(*, "(a)", advance="no") "ok..."
      else
         write(*, "(a,es12.4)") " f+s fail residual 1d = ", &
            maxval(abs(res(1:a%n,1)))
         !write (6,'(6es12.4)') x1(1:a%n)
         errors = errors + 1
      endif
      call compute_resid(nrhs,a,x,maxn,rhs,maxn,res,maxn)
      if(maxval(abs(res(1:a%n,1:nrhs))) < err_tol_scale) then
         write(*, "(a)") "ok..."
      else
         write(*, "(a)") " f+s fail residual 2d = "
         do i = 1, nrhs
            write(*, "(es12.4)") maxval(abs(res(1:a%n,i)))
         end do
         write(*, "()")
       !  write (6,'(6es12.4)') x1(1:a%n)
         errors = errors + 1
      endif

      !
      ! Cleanup for next iteration
      !
      call ssids_free(akeep, fkeep, cuda_error)

   end do

   ! Free any memory in the event of errors/singular systems on final itr
   call ssids_free(akeep, fkeep, cuda_error)

end subroutine test_random_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generate random symmetric indefinite matrix (lower part)
subroutine gen_random_indef_fred(a, nza, state)
   type(matrix_type), intent(inout) :: a
   integer(long), intent(in) :: nza
   type(random_state), intent(inout) :: state

   integer :: i, k, l, flag

   ! Generate a FIXME: move to 64-bit
   call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_INDEF, a%n, a%n, &
      int(nza), a%ptr, a%row, flag, val=a%val, sort=.true.)
   if(flag.ne.0) print *, "Bad flag from random_matrix_generate()"

   if (a%n.gt.3) then
      ! Put some zeros on diagonal, observing first entry in column
      ! is always the diagonal after sorting

      l = random_integer(state,  a%n/16)
      do k = 1, a%n, max(1,l)
         i = a%ptr(k)
         a%val(i) = zero
      end do
   end if


end subroutine gen_random_indef_fred

!********************************************************

! sparse matrix-vector multiplication y=a*x when only
! lower triangle held
subroutine matvec_lower(a,x,y)
   type(matrix_type), intent(in) :: a ! Holds the matrix
   real(wp), intent(in) :: x(*) ! input vector
   real(wp), intent(out) :: y(*) ! output vector y=a*x

   integer :: n,i,j,k,jstrt,jstop
   real(wp) :: sum

   n = a%n
   y(1:n) = zero

   jstrt = 1
   do i = 1,n
      jstop = a%ptr(i+1)
      sum = zero
      do j = jstrt,jstop-1
         k = a%row(j)
         sum = sum + a%val(j)*x(k)
         if (k.eq.i) cycle
         y(k) = y(k) + a%val(j)*x(i)
      end do
      y(i) = y(i) + sum
      jstrt = jstop
   end do

end subroutine matvec_lower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simple_metis_order(a, order)
   type(matrix_type), intent(in) :: a
   integer, dimension(:), allocatable :: order

   integer :: i, flag, stat
   integer, dimension(:), allocatable :: invp

   allocate(invp(a%n))

   ! Perform MeTiS
   ! FIXME: which way around should we have order and invp?
   call metis_order(a%n, a%ptr, a%row, order, invp, flag, stat)
   if(flag.ne.0) then
      ! Failed for some reason
      do i = 1, a%n
         order(i) = i
      end do
      return
   endif

end subroutine simple_metis_order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_random_indef(a, nza, state, zr)
   type(matrix_type), intent(inout) :: a
   integer(long), intent(in) :: nza
   type(random_state), intent(inout) :: state
   integer, optional, intent(in) :: zr ! if present, all entries in
     ! row zr are zero

   integer :: i, k, l, flag

   ! Generate a FIXME: move to 64-bit
   call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_INDEF, a%n, a%n, &
      int(nza), a%ptr, a%row, flag, val=a%val, nonsingular=.true., sort=.true.)
   if(flag.ne.0) print *, "Bad flag from random_matrix_generate()"

   if (present(zr)) then
      ! Scan along row
      do i = a%ptr(1), a%ptr(zr)-1
         if(a%row(i).eq.zr) a%val(i) = zero
      end do
      ! Scan along column
      do i = a%ptr(zr),a%ptr(zr+1)-1
         a%val(i) = zero
      end do
   elseif(a%n.gt.3) then
      ! Put some zeros on diagonal, observing first entry in column
      ! is always the diagonal after sorting
      ! but don't have all zeros in the col.
      l = random_integer(state,  a%n/2)
      do k = 1, a%n, max(1,l)
         if (a%ptr(k+1) > a%ptr(k) + 1) then
           i = a%ptr(k)
           a%val(i) = zero
         endif
      end do
      ! also make sure we have some large off diagonals
      do k = 1, a%n
         i = a%ptr(k+1) - 1
         a%val(i) = a%val(i)*1000
      end do
   endif

end subroutine gen_random_indef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_random_posdef(a, nza, state)
   type(matrix_type), intent(inout) :: a
   integer(long), intent(in) :: nza
   type(random_state), intent(inout) :: state

   integer :: i, j, k, flag
   real(wp) :: tempv

   ! Generate matrix FIXME: move to 64-bit
   call random_matrix_generate(state, SPRAL_MATRIX_REAL_SYM_PSDEF, a%n, a%n, &
      int(nza), a%ptr, a%row, flag, val=a%val, nonsingular=.true., sort=.true.)
   if(flag.ne.0) print *, "Bad flag from random_matrix_generate()"

   ! Make a diagonally dominant, observing first entry in column
   ! is always the diagonal after sorting
   do k = 1, a%n
      tempv = zero
      do j = a%ptr(k)+1, a%ptr(k+1)-1
         tempv = tempv + abs(a%val(j))
         i = a%ptr(a%row(j))
         a%val(i) = a%val(i) + abs(a%val(j))
      end do
      i = a%ptr(k)
      a%val(i) = one + a%val(i) + tempv
   end do
end subroutine gen_random_posdef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assert_flag(context, actual, expected)
   character(len=*), intent(in) :: context
   integer, intent(in) :: actual
   integer, intent(in) :: expected

   if(actual.eq.expected) return ! All is good

   ! Otherwise report/record error
   write(*, "(3a,i4,a,i4,a)") &
      "Unexpected error during ", context, ". flag = ", actual, &
      " (expected ", expected, ")"
   errors = errors + 1
   return
end subroutine assert_flag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_result(actual, expected, continued)
   integer :: actual
   integer :: expected
   logical, optional :: continued

   logical :: mycontinued

   mycontinued = .false.
   if(present(continued)) mycontinued = continued

   if(actual.eq.expected) then
      if(mycontinued) then
         write(*,"(a)", advance="no") "ok..."
      else
         write(*,"(a)") "ok"
      endif
      return
   endif

   write(*,"(a)") "fail"
   write(*,"(2(a,i4))") "returned ", actual, ", expected ", expected
   errors = errors + 1
end subroutine print_result


end program
