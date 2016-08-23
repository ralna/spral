!> \file
!> \copyright 2016 The Science and Technology Facilities Council (STFC)
!> \licence   BSD licence, see LICENCE file for details
!> \author    Jonathan Hogg
!
!> \brief Define ssids_fkeep type and associated procedures (CPU version)
module spral_ssids_fkeep
   use, intrinsic :: iso_c_binding
   use :: omp_lib
   use spral_ssids_akeep, only : ssids_akeep
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_datatypes
   use spral_ssids_inform, only : ssids_inform
   use spral_ssids_subtree, only : numeric_subtree_base
   use spral_ssids_cpu_subtree, only : cpu_numeric_subtree
   use spral_ssids_cpu_profile, only : cpu_profile_begin, cpu_profile_end
   implicit none

   private
   public :: ssids_fkeep

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type numeric_subtree_ptr
      class(numeric_subtree_base), pointer :: ptr
   end type numeric_subtree_ptr

   !
   ! Data type for data generated in factorise phase
   !
   type ssids_fkeep
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
   end type ssids_fkeep

contains

subroutine inner_factor_cpu(fkeep, akeep, val, options, inform)
   type(ssids_akeep), intent(in) :: akeep
   class(ssids_fkeep), target, intent(inout) :: fkeep
   real(wp), dimension(*), target, intent(in) :: val
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(inout) :: inform

   integer :: i, numa_region, exec_loc
   integer :: total_threads
   logical :: all_region
   type(contrib_type), dimension(:), allocatable :: child_contrib

   ! Begin profile trace (noop if not enabled)
   call cpu_profile_begin()

   ! Allocate space for subtrees
   allocate(fkeep%subtree(akeep%nparts), stat=inform%stat)
   if(inform%stat.ne.0) goto 100

   ! Call subtree factor routines
   allocate(child_contrib(akeep%nparts), stat=inform%stat)
   if(inform%stat.ne.0) goto 100
   ! Split into numa regions; parallelism within a region is responsibility
   ! of subtrees.
   ! FIXME: do we not want to have within-node parallelism at a higher level?
   all_region = .false.
!$omp parallel proc_bind(spread) num_threads(size(akeep%topology)) &
!$omp    default(none) private(i, exec_loc, numa_region) &
!$omp    shared(akeep, fkeep, val, options, inform, child_contrib, all_region)
   numa_region = omp_get_thread_num()
   call omp_set_num_threads(akeep%topology(numa_region+1)%nproc)
   do i = 1, akeep%nparts
      exec_loc = akeep%subtree(i)%exec_loc
      if(numa_region.eq.0 .and. exec_loc.eq.-1) all_region = .true.
      if(mod(exec_loc,size(akeep%topology)).ne.numa_region) cycle
      if(allocated(fkeep%scaling)) then
         fkeep%subtree(i)%ptr => akeep%subtree(i)%ptr%factor( &
            fkeep%pos_def, val, &
            child_contrib(akeep%contrib_ptr(i):akeep%contrib_ptr(i+1)-1), &
            options, inform, scaling=fkeep%scaling &
            )
      else
         fkeep%subtree(i)%ptr => akeep%subtree(i)%ptr%factor( &
            fkeep%pos_def, val, &
            child_contrib(akeep%contrib_ptr(i):akeep%contrib_ptr(i+1)-1), &
            options, inform &
            )
      endif
      if(akeep%contrib_idx(i).gt.akeep%nparts) cycle ! part is a root
      child_contrib(akeep%contrib_idx(i)) = fkeep%subtree(i)%ptr%get_contrib()
!$omp flush
      child_contrib(akeep%contrib_idx(i))%ready = .true.
   end do
!$omp end parallel

   if(all_region) then
      ! At least some all region subtrees exist
      total_threads = 0
      do i = 1, size(akeep%topology)
         total_threads = total_threads + akeep%topology(i)%nproc
      end do
      call omp_set_num_threads(total_threads)
      do i = 1, akeep%nparts
         exec_loc = akeep%subtree(i)%exec_loc
         if(exec_loc.ne.-1) cycle
         if(allocated(fkeep%scaling)) then
            fkeep%subtree(i)%ptr => akeep%subtree(i)%ptr%factor( &
               fkeep%pos_def, val, &
               child_contrib(akeep%contrib_ptr(i):akeep%contrib_ptr(i+1)-1), &
               options, inform, scaling=fkeep%scaling &
               )
         else
            fkeep%subtree(i)%ptr => akeep%subtree(i)%ptr%factor( &
               fkeep%pos_def, val, &
               child_contrib(akeep%contrib_ptr(i):akeep%contrib_ptr(i+1)-1), &
               options, inform &
               )
         endif
         if(akeep%contrib_idx(i).gt.akeep%nparts) cycle ! part is a root
         child_contrib(akeep%contrib_idx(i)) = &
            fkeep%subtree(i)%ptr%get_contrib()
!$omp    flush
         child_contrib(akeep%contrib_idx(i))%ready = .true.
      end do
   end if

   ! End profile trace (noop if not enabled)
   call cpu_profile_end()

   return

   100 continue
   inform%flag = SSIDS_ERROR_ALLOCATION
   return
end subroutine inner_factor_cpu

subroutine inner_solve_cpu(local_job, nrhs, x, ldx, akeep, fkeep, options, inform)
   type(ssids_akeep), intent(in) :: akeep
   class(ssids_fkeep), intent(inout) :: fkeep
   integer, intent(inout) :: local_job
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,nrhs), target, intent(inout) :: x
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(inout) :: inform

   integer :: i, r, part
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
         local_job.eq.SSIDS_SOLVE_JOB_ALL) then
      do part = 1, akeep%nparts
         call fkeep%subtree(part)%ptr%solve_fwd(nrhs, x2, n, inform)
         if (inform%stat .ne. 0) goto 100
      end do
   endif

   if (local_job.eq.SSIDS_SOLVE_JOB_DIAG) then
      do part = 1, akeep%nparts
         call fkeep%subtree(part)%ptr%solve_diag(nrhs, x2, n, inform)
         if (inform%stat .ne. 0) goto 100
      end do
   endif

   if (local_job.eq.SSIDS_SOLVE_JOB_BWD) then
      do part = akeep%nparts, 1, -1
         call fkeep%subtree(part)%ptr%solve_bwd(nrhs, x2, n, inform)
         if (inform%stat .ne. 0) goto 100
      end do
   endif

   if (local_job.eq.SSIDS_SOLVE_JOB_DIAG_BWD .or. &
         local_job.eq.SSIDS_SOLVE_JOB_ALL) then
      do part = akeep%nparts, 1, -1
         call fkeep%subtree(part)%ptr%solve_diag_bwd(nrhs, x2, n, inform)
         if (inform%stat .ne. 0) goto 100
      end do
   endif

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
   type(ssids_akeep), intent(in) :: akeep
   class(ssids_fkeep), target, intent(in) :: fkeep
   type(ssids_inform), intent(inout) :: inform
   real(wp), dimension(*), intent(out) :: d

   integer :: n
   integer :: part, sa, en

   n = akeep%n
   ! ensure d is not returned undefined
   d(1:n) = 0.0 ! ensure do not returned with this undefined

   do part = 1, akeep%nparts
      sa = akeep%part(part)
      en = akeep%part(part+1)-1
      associate(subtree => fkeep%subtree(part)%ptr)
         select type(subtree)
         type is (cpu_numeric_subtree)
            call subtree%enquire_posdef(d(sa:en))
         end select
      end associate
   end do
   
end subroutine enquire_posdef_cpu

!****************************************************************************

subroutine enquire_indef_cpu(akeep, fkeep, inform, piv_order, d)
   type(ssids_akeep), intent(in) :: akeep
   class(ssids_fkeep), target, intent(in) :: fkeep
   type(ssids_inform), intent(inout) :: inform
   integer, dimension(akeep%n), optional, intent(out) :: piv_order
      ! If i is used to index a variable, its position in the pivot sequence
      ! will be placed in piv_order(i), with its sign negative if it is
      ! part of a 2 x 2 pivot; otherwise, piv_order(i) will be set to zero.
   real(wp), dimension(2,akeep%n), optional, intent(out) :: d ! The diagonal
      ! entries of D^{-1} will be placed in d(1,:i) and the off-diagonal
      ! entries will be placed in d(2,:). The entries are held in pivot order.

   integer :: part, sa
   integer :: i, n
   integer, dimension(:), allocatable :: po

   n = akeep%n
   if(present(d)) then
      ! ensure d is not returned undefined
      d(1:2,1:n) = 0.0
   end if

   ! We need to apply the invp externally to piv_order
   if(present(piv_order)) allocate(po(akeep%n)) ! FIXME: stat

   ! FIXME: should probably return nelim from each part, due to delays passing
   ! between them
   do part = 1, akeep%nparts
      sa = akeep%part(part)
      associate(subtree => fkeep%subtree(1)%ptr)
         select type(subtree)
         type is (cpu_numeric_subtree)
            if(present(d)) then
               if(present(piv_order)) then
                  call subtree%enquire_indef(piv_order=po(sa:), d=d(1:2,sa:))
               else
                  call subtree%enquire_indef(d=d(1:2,sa:))
               endif
            else
               if(present(piv_order)) then
                  call subtree%enquire_indef(piv_order=po(sa:))
               else
                  ! No-op
                  ! FIXME: should we report an error here? (or done higher up?)
               endif
            endif
         end select
      end associate
   end do

   ! Apply invp to piv_order
   if(present(piv_order)) then
      do i = 1, akeep%n
         piv_order( akeep%invp(i) ) = po(i)
      end do
   endif

end subroutine enquire_indef_cpu

! Alter D values
subroutine alter_cpu(d, akeep, fkeep, options, inform)
   real(wp), dimension(2,*), intent(in) :: d  ! The required diagonal entries
     ! of D^{-1} must be placed in d(1,i) (i = 1,...n)
     ! and the off-diagonal entries must be placed in d(2,i) (i = 1,...n-1).
   type(ssids_akeep), intent(in) :: akeep
   class(ssids_fkeep), target, intent(inout) :: fkeep
   type(ssids_options), intent(in) :: options
   type(ssids_inform), intent(inout) :: inform

   integer :: part

   do part = 1, akeep%nparts
      associate(subtree => fkeep%subtree(1)%ptr)
         select type(subtree)
         type is (cpu_numeric_subtree)
            call subtree%alter(d(1:2,akeep%part(part):akeep%part(part+1)-1))
         end select
      end associate
   end do
end subroutine alter_cpu

!****************************************************************************

subroutine free_fkeep(fkeep, flag)
   class(ssids_fkeep), intent(inout) :: fkeep
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
