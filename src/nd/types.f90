module spral_nd_types
   implicit none

   ! NB: Everything in here is public

   integer, parameter :: wp = kind(1.0D0)
   integer, parameter :: long = selected_int_kind(18)

   ! Error flags
   integer, parameter :: ND_ERR_MEMORY_ALLOC = -1,   &
                         ND_ERR_MEMORY_DEALLOC = -2, &
                         ND_ERR_N = -3                 ! n<1

   ! Partition flags
   integer, parameter :: ND_PART1_FLAG = 0, &
                         ND_PART2_FLAG = 2, &
                         ND_SEP_FLAG = 1

   ! *****************************************************************

   type :: nd_options
      integer :: print_level = 0 ! amount of informational output required
      integer :: unit_diagnostics = 6 ! stream number for diagnostic
         ! messages
      integer :: unit_error = 6 ! stream number for error messages

      integer :: amd_call = 5 ! call AMD if number of supervariables is
         ! less than or equal to amd_call  *UPDATE*
      integer :: amd_switch1 = 50 ! switch to (halo)AMD if matrix size is
         ! less than or equal to nd_switch
      integer :: amd_switch2 = 20 ! maximum number of ND levels allowed
      integer :: cost_function = 1 ! selects the cost function used for
         ! evaluating a partition
         ! <=1 : |P1||P2|/|S| + penalty for imbalance
         ! >=2 : |S|(1 +  0.5||P1|-|P2||)
      integer :: partition_method = 2 ! Are we allowed to use a multilevel
         ! strategy
         ! <= 0 : do not use multilevel
         ! == 1 : use multilevel
         ! >= 2 : automatic choice based on size of levelsets
      integer :: matching = 0 ! Which coarsening method to use
         ! > 0 : heavy-edge
         ! <= 0 : common neighbours
      integer :: coarse_partition_method = 1 ! Which partition method to use
         ! at coarsest level
         ! <=1 : Ashcraft method (half-level set)
         ! >=2 : Level-set method
      integer :: refinement = 7 ! Which sort of refinement to use
         ! <2 : trim + fm to increase weights of A1 and A2
         ! =2 : trim + fm to increase weight of smaller partition and
         !      decrease weight of larger partition
         ! =3 : Automatically choose between options 1 and 2
         ! =4 : maxflow + fm to increase weights of A1 and A2
         ! =5 : maxflow + fm to increase weight of smaller partition and
         !      decrease weight of larger partition
         ! =6 : Automatically choose between options 4 and 5
         ! >6 : Automatically choose between options 1 and 5
      integer :: refinement_band = 4 ! band width for FM refinement. Values
         ! less than 1 mean that no FM refinement is done
      logical :: remove_dense_rows = .true. ! test the input for dense rows
         ! and place them at the end of the ordering
      integer :: stop_coarsening2 = 20 ! Max number of levels in the
         ! multilevel grid
      integer :: stop_coarsening1 = 100 ! Stop coarsening once matrix has
         ! order at most stop_coarsening1
      integer :: ml_call = 12000 ! Stop coarsening once matrix has
         ! order at most stop_coarsening1

      ! minimum and maximum grid reduction factors that must be achieved
      ! during coarsening. If cgrid%size is greater than
      ! max_reduction*grid%size or cgrid%size is less than
      ! min_reduction*grid%size then carry on coarsening
      real(wp) :: min_reduction = 0.5 ! size of next multigrid
         ! matrix must be greater than min_reduction*(size of current
         ! multigrid matrix)
      real(wp) :: max_reduction = 0.9 ! size of next multigrid
         ! matrix must be less than max_reduction*(size of current multigrid
         ! matrix)
      real(wp) :: balance = 2.0 ! Try to make sure that
         ! max(P1,P2)/min(P1/P2) <= balance

      integer :: max_improve_cycles = 2 ! Having computed a minimal partition,
         ! expand and refine it at most max_improve_cycles times to improve
         ! the quality of the partition.
      logical :: find_supervariables = .true. ! If .true., after dense rows
         ! have been (optionally) removed, check for supervariables and
         ! compress matrix if supervariables found.
   end type nd_options

   ! *****************************************************************

   type :: nd_inform
      integer :: flag = 0 ! error/warning flag
      integer :: dense = 0 ! holds number of dense rows
      integer :: stat = 0 ! holds Fortran stat parameter
      integer :: nsuper = 0 ! holds number of supervariables + number of zero
         ! rows after dense rows removed
      integer :: nzsuper = 0 ! holds number of nonzeros in compressed graph
      integer :: num_components = 1 ! Number of independent components after
         ! dense rows (optionally) removed
      integer :: n_max_component = -1 ! holds number rows in largest indep
         ! component after dense rows removed and supervariables (optionally)
         ! compressed
      integer :: nz_max_component = -1 ! holds number nonzeros in largest indep
         ! component after dense rows removed and supervariables (optionally)
         ! compressed
      integer :: maxdeg_max_component = -1 ! holds number nonzeros in largest
         ! indep component after dense rows removed and supervariables
         ! (optionally) compressed
      real(wp) :: band = -1 ! holds L, where L is the size
         ! of the largest level set at the top level of nested dissection. If
         ! the matrix is reducible, then it holds the value for the largest
         ! of the irreducible components.
         ! Not returned if options%partition_method==1.
      real(wp) :: depth = -1 ! holds number of levels in level set structure
         ! at the top level of nested dissection. If the matrix is reducible,
         ! then it holds the value for the largest of the irreducible
         ! components.
         ! Not returned if options%partition_method==1.
   end type nd_inform

   ! *****************************************************************
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! General utility functions needed by all nd modules
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Prints out errors and warnings according to value of flag
!
subroutine nd_print_message(flag,unit,context)
   integer, intent(in) :: flag ! Error flag to print message for
   integer, intent(in) :: unit ! Fortran unit to print message to
   character(len=*), intent(in) :: context ! context to print with message

   if (unit.lt.0) return ! No output

   if (flag<0) write (unit,fmt='('' ERROR: '')')
   write (unit,advance='no',fmt='('' '', a,'':'')') &
      context(1:len_trim(context))

   select case (flag)
   case (0)
      write (unit,'(A)') ' successful completion'
   case (ND_ERR_MEMORY_ALLOC)
      write (unit,'(A)') ' memory allocation failure'
   case (ND_ERR_MEMORY_DEALLOC)
      write (unit,'(A)') ' memory deallocation failure'
   case (ND_ERR_N)
      write (unit,'(A)') ' n<1'
   end select
end subroutine nd_print_message

!
! Returns ptr(idx) if idx.le.n, or ne+1 otherwise
!
integer function nd_get_ptr(idx, n, ne, ptr)
   integer, intent(in) :: idx, n, ne
   integer, dimension(n), intent(in) :: ptr

   if(idx.le.n) then
      nd_get_ptr = ptr(idx)
   else
      nd_get_ptr = ne+1
   endif
end function nd_get_ptr

end module spral_nd_types
