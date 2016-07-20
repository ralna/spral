!
! This module defines ssids_fkeep type and associated procedures (CPU version)
!
module spral_ssids_fkeep
   use spral_ssids_akeep, only : ssids_akeep_base
   use spral_ssids_datatypes
   use spral_ssids_inform, only : ssids_inform_base
   use spral_ssids_cpu_subtree, only : cpu_numeric_subtree, cpu_symbolic_subtree
   use spral_ssids_subtree, only : numeric_subtree_base
   use, intrinsic :: iso_c_binding
   implicit none

   private
   public :: ssids_fkeep_base

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type numeric_subtree_ptr
      class(numeric_subtree_base), pointer :: ptr
   end type numeric_subtree_ptr

   !
   ! Data type for data generated in factorise phase
   !
   type ssids_fkeep_base
      integer :: flag ! copy of error flag.
      real(wp), dimension(:), allocatable :: scaling ! Stores scaling for
         ! each entry (in original matrix order)
      logical :: pos_def ! set to true if user indicates matrix pos. definite

      ! Factored subtrees
      type(numeric_subtree_ptr), dimension(:), allocatable :: subtree

      ! Info components to copy on solve
      integer :: matrix_rank
      integer :: maxfront
      integer :: num_delay
      integer(long) :: num_factor
      integer(long) :: num_flops
      integer :: num_neg
      integer :: num_two
      integer :: not_first_pass
      integer :: not_second_pass

   contains
      procedure, pass(fkeep) :: inner_factor => inner_factor_cpu ! Do actual factorization
      procedure, pass(fkeep) :: inner_solve => inner_solve_cpu ! Do actual solve
      procedure, pass(fkeep) :: enquire_posdef => enquire_posdef_cpu
      procedure, pass(fkeep) :: enquire_indef => enquire_indef_cpu
      procedure, pass(fkeep) :: alter => alter_cpu ! Alter D values
      procedure, pass(fkeep) :: free => free_fkeep ! Frees memory
   end type ssids_fkeep_base

contains

subroutine inner_factor_cpu(fkeep, akeep, val, options, inform)
   class(ssids_akeep_base), intent(in) :: akeep
   class(ssids_fkeep_base), target, intent(inout) :: fkeep
   real(wp), dimension(*), target, intent(in) :: val
   class(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform

   logical(C_BOOL) :: fposdef
   class(numeric_subtree_base), pointer :: subtree

   fposdef = fkeep%pos_def

   ! Allocate space for subtrees
   allocate(fkeep%subtree(1), stat=inform%stat)
   if(inform%stat.ne.0) then
      inform%flag = SSIDS_ERROR_ALLOCATION
      return
   endif

   ! Call subtree factor routines
   if(allocated(fkeep%scaling)) then
      fkeep%subtree(1)%ptr => akeep%subtree(1)%ptr%factor( &
         fposdef, val, options, inform, scaling=fkeep%scaling &
         )
   else
      fkeep%subtree(1)%ptr => akeep%subtree(1)%ptr%factor( &
         fposdef, val, options, inform &
         )
   endif
end subroutine inner_factor_cpu

subroutine inner_solve_cpu(local_job, nrhs, x, ldx, akeep, fkeep, options, inform)
   class(ssids_akeep_base), intent(in) :: akeep
   class(ssids_fkeep_base), intent(inout) :: fkeep
   integer, intent(inout) :: local_job
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,nrhs), target, intent(inout) :: x
   type(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform

   integer :: i, r
   integer :: n
   real(wp), dimension(:,:), allocatable :: x2

   n = akeep%n

   allocate(x2(n, nrhs), stat=inform%stat)
   if(inform%stat.ne.0) goto 100

   ! Permute/scale
   if (allocated(fkeep%scaling) .and. (local_job == SSIDS_SOLVE_JOB_ALL .or. &
            local_job == SSIDS_SOLVE_JOB_FWD)) then
      ! Copy and scale
      do r = 1, nrhs
         do i = 1, n
            x2(i,r) = x(akeep%invp(i),r) * fkeep%scaling(i)
         end do
      end do
   else
      ! Just copy
      do r = 1, nrhs
         x2(1:n, r) = x(akeep%invp(1:n), r)
      end do
   end if


   ! Perform relevant solves
   if (local_job.eq.SSIDS_SOLVE_JOB_FWD .or. &
         local_job.eq.SSIDS_SOLVE_JOB_ALL) &
      call fkeep%subtree(1)%ptr%solve_fwd(nrhs, x2, n, inform)
   if (inform%stat .ne. 0) goto 100

   if (local_job.eq.SSIDS_SOLVE_JOB_DIAG) &
      call fkeep%subtree(1)%ptr%solve_diag(nrhs, x2, n, inform)
   if (inform%stat .ne. 0) goto 100

   if (local_job.eq.SSIDS_SOLVE_JOB_BWD) &
      call fkeep%subtree(1)%ptr%solve_bwd(nrhs, x2, n, inform)
   if (inform%stat .ne. 0) goto 100

   if (local_job.eq.SSIDS_SOLVE_JOB_DIAG_BWD .or. &
         local_job.eq.SSIDS_SOLVE_JOB_ALL) &
      call fkeep%subtree(1)%ptr%solve_diag_bwd(nrhs, x2, n, inform)
   if (inform%stat .ne. 0) goto 100

   ! Unscale/unpermute
   if (allocated(fkeep%scaling) .and. ( &
            local_job == SSIDS_SOLVE_JOB_ALL .or. &
            local_job == SSIDS_SOLVE_JOB_BWD .or. &
            local_job == SSIDS_SOLVE_JOB_DIAG_BWD)) then
      ! Copy and scale
      do r = 1, nrhs
         do i = 1, n
            x(akeep%invp(i),r) = x2(i,r) * fkeep%scaling(i)
         end do
      end do
   else
      ! Just copy
      do r = 1, nrhs
         x(akeep%invp(1:n), r) = x2(1:n, r)
      end do
   end if

   return

   100 continue
   inform%flag = SSIDS_ERROR_ALLOCATION
   return

end subroutine inner_solve_cpu

!****************************************************************************

subroutine enquire_posdef_cpu(akeep, fkeep, inform, d)
   class(ssids_akeep_base), intent(in) :: akeep
   class(ssids_fkeep_base), target, intent(in) :: fkeep
   class(ssids_inform_base), intent(inout) :: inform
   real(wp), dimension(*), intent(out) :: d

   integer :: n

   n = akeep%n
   ! ensure d is not returned undefined
   d(1:n) = 0.0 ! ensure do not returned with this undefined

   associate(subtree => fkeep%subtree(1)%ptr)
      select type(subtree)
      type is (cpu_numeric_subtree)
         call subtree%enquire_posdef(d)
      end select
   end associate
   
end subroutine enquire_posdef_cpu

!****************************************************************************

subroutine enquire_indef_cpu(akeep, fkeep, inform, piv_order, d)
   class(ssids_akeep_base), intent(in) :: akeep
   class(ssids_fkeep_base), target, intent(in) :: fkeep
   class(ssids_inform_base), intent(inout) :: inform
   integer, dimension(*), optional, intent(out) :: piv_order
      ! If i is used to index a variable, its position in the pivot sequence
      ! will be placed in piv_order(i), with its sign negative if it is
      ! part of a 2 x 2 pivot; otherwise, piv_order(i) will be set to zero.
   real(wp), dimension(2,*), optional, intent(out) :: d ! The diagonal
      ! entries of D^{-1} will be placed in d(1,:i) and the off-diagonal
      ! entries will be placed in d(2,:). The entries are held in pivot order.

   integer :: n

   n = akeep%n
   if(present(d)) then
      ! ensure d is not returned undefined
      d(1:2,1:n) = 0.0
   end if

   associate(subtree => fkeep%subtree(1)%ptr)
      select type(subtree)
      type is (cpu_numeric_subtree)
         call subtree%enquire_indef(akeep%invp, &
            piv_order=piv_order, d=d)
      end select
   end associate

end subroutine enquire_indef_cpu

! Alter D values
subroutine alter_cpu(d, akeep, fkeep, options, inform)
   real(wp), dimension(2,*), intent(in) :: d  ! The required diagonal entries
     ! of D^{-1} must be placed in d(1,i) (i = 1,...n)
     ! and the off-diagonal entries must be placed in d(2,i) (i = 1,...n-1).
   type(ssids_akeep_base), intent(in) :: akeep
   class(ssids_fkeep_base), target, intent(inout) :: fkeep
   type(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform

   associate(subtree => fkeep%subtree(1)%ptr)
      select type(subtree)
      type is (cpu_numeric_subtree)
         call subtree%alter(d)
      end select
   end associate
end subroutine alter_cpu

!****************************************************************************

subroutine free_fkeep(fkeep, flag)
   class(ssids_fkeep_base), intent(inout) :: fkeep
   integer, intent(out) :: flag ! not actually used for cpu version, set to 0

   integer :: i
   integer :: st

   flag = 0 ! Not used for basic SSIDS, just zet to zero

   if(allocated(fkeep%subtree)) then
      do i = 1, size(fkeep%subtree)
         if(associated(fkeep%subtree(i)%ptr)) deallocate(fkeep%subtree(i)%ptr)
      end do
      deallocate(fkeep%subtree)
   endif

   deallocate(fkeep%scaling, stat=st)

end subroutine free_fkeep

end module spral_ssids_fkeep
