module spral_nd_types
   implicit none

   ! NB: Everything in here is public

   integer, parameter :: wp = kind(1.0D0)
   integer, parameter :: long = selected_int_kind(18)

   ! Error flags
   integer, parameter :: ND_ERR_MEMORY_ALLOC = -1,   &
                         ND_ERR_MEMORY_DEALLOC = -2, &
                         ND_ERR_N = -3,              & ! n<1
                         ND_ERR_INTERNAL = -99         ! Shouldn't happen

   ! Partition flags
   integer, parameter :: ND_PART1_FLAG = 0,     &
                         ND_PART2_FLAG = 2,     &
                         ND_SEP_FLAG = 1,       &
                         ND_SPECIAL_FLAG = -1      ! Used by some algorithms

   ! Option values
   integer, parameter :: ND_PARTITION_HALF_LEVEL_SET = 1, &
                         ND_PARTITION_LEVEL_SET      = 2
   integer, parameter :: ND_REFINE_TRIM_FM_BOTH    = 1, &
                         ND_REFINE_TRIM_FM_SMALLER = 2, &
                         ND_REFINE_TRIM_FM_AUTO    = 3, &
                         ND_REFINE_MAXFLOW_BOTH    = 4, &
                         ND_REFINE_MAXFLOW_SMALLER = 5, &
                         ND_REFINE_MAXFLOW_AUTO    = 6, &
                         ND_REFINE_AUTO            = 7
   integer, parameter :: ND_MATCH_COMMON_NEIGHBOURS = 1, &
                         ND_MATCH_SHEM              = 2

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
         !  =2 : |S|(1 +  0.5||P1|-|P2||)
         ! >=3 : As 1, but with more convex imbalance penalty
      integer :: partition_method = 2 ! Are we allowed to use a multilevel
         ! strategy
         ! <= 0 : do not use multilevel
         ! == 1 : use multilevel
         ! >= 2 : automatic choice based on size of levelsets
      integer :: matching = ND_MATCH_SHEM ! Which coarsening method to use
         ! <= 1 : common neighbours matching (CNM)
         ! >= 2 : sorted heavy-edge matching (SHEM)
      integer :: coarse_partition_method = 1 ! Which partition method to use
         ! at coarsest level
         ! <=1 : Half-level set method (Ashcraft)
         !  =2 : Level-set method
         ! >=3 : Region-growing edge seperator vertex cover
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
      integer :: stop_coarsening1 = 100 ! Stop coarsening once matrix has
         ! order at most stop_coarsening1
      integer :: stop_coarsening2 = 20 ! Max number of levels in the
         ! multilevel grid
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
      real(wp) :: balance = 4.0 ! Try to make sure that
         ! max(P1,P2)/min(P1/P2) <= balance

      integer :: max_improve_cycles = 2 ! Having computed a minimal partition,
         ! expand and refine it at most max_improve_cycles times to improve
         ! the quality of the partition.
      logical :: find_supervariables = .true. ! If .true., after dense rows
         ! have been (optionally) removed, check for supervariables and
         ! compress matrix if supervariables found.

      ! FIXME: decide whether to keep?
      integer :: reord = 1 ! which reordering to use in preprocessing phase
         ! 1 = Jonathan's [literal copy]
         ! 2 = Sue's [slightly more ordered, but strange]
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

   type nd_matrix
      integer :: m ! number rows
      integer :: n ! number columns
      integer :: ne ! number entries in matrix
      integer, allocatable, dimension(:) :: ptr ! pointer into col array
      integer, allocatable, dimension(:) :: col ! column indices
      integer, allocatable, dimension(:) :: val ! values
   end type nd_matrix

   ! *****************************************************************

   type nd_multigrid
      integer :: size ! number of vertices in graph at this level
      type (nd_matrix) :: graph ! this level of matrix
      integer, allocatable, dimension(:) :: match ! matching
      integer, allocatable, dimension(:) :: where ! where each row of
         ! this level of matrix will go (ie ordering for this level)
      integer, allocatable, dimension(:) :: row_wgt ! number of
         ! vertices this vertex of the coarse graph matrix represents
      integer :: level = 0 ! the level
      integer :: part_div(2) ! number of vertices in each part
   end type nd_multigrid

end module spral_nd_types
