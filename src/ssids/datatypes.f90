module spral_ssids_datatypes
!$ use omp_lib
   use, intrinsic :: iso_c_binding
   use spral_cuda, only : cudaGetErrorString
   use spral_scaling, only : auction_options, auction_inform
   implicit none

   private
   public :: smalloc_type, stack_mem_type, ssids_akeep, &
      stack_type, thread_stats, real_ptr_type, &
      ssids_options, ssids_inform, node_type
   public :: ssids_print_flag

   integer, parameter, public :: wp = C_DOUBLE
   integer, parameter, public :: long = selected_int_kind(18)

   real(wp), parameter, public :: one = 1.0_wp
   real(wp), parameter, public :: zero = 0.0_wp

   integer, parameter, public :: nemin_default = 8 ! node amalgamation parameter

   ! Success flag
   integer, parameter, public :: SSIDS_SUCCESS               = 0

   ! Error flags
   integer, parameter, public :: SSIDS_ERROR_CALL_SEQUENCE     = -1
   integer, parameter, public :: SSIDS_ERROR_A_N_OOR           = -2
   integer, parameter, public :: SSIDS_ERROR_A_PTR             = -3
   integer, parameter, public :: SSIDS_ERROR_A_ALL_OOR         = -4
   integer, parameter, public :: SSIDS_ERROR_SINGULAR          = -5
   integer, parameter, public :: SSIDS_ERROR_NOT_POS_DEF       = -6
   integer, parameter, public :: SSIDS_ERROR_PTR_ROW           = -7
   integer, parameter, public :: SSIDS_ERROR_ORDER             = -8
   integer, parameter, public :: SSIDS_ERROR_VAL               = -9
   integer, parameter, public :: SSIDS_ERROR_X_SIZE            = -10
   integer, parameter, public :: SSIDS_ERROR_JOB_OOR           = -11
   integer, parameter, public :: SSIDS_ERROR_PRESOLVE_INCOMPAT = -12
   integer, parameter, public :: SSIDS_ERROR_NOT_LLT           = -13
   integer, parameter, public :: SSIDS_ERROR_NOT_LDLT          = -14
   integer, parameter, public :: SSIDS_ERROR_NO_SAVED_SCALING  = -15
   integer, parameter, public :: SSIDS_ERROR_ALLOCATION        = -50
   integer, parameter, public :: SSIDS_ERROR_CUDA_UNKNOWN      = -51
   integer, parameter, public :: SSIDS_ERROR_CUBLAS_UNKNOWN    = -52
   integer, parameter, public :: SSIDS_ERROR_UNKNOWN           = -99

   ! warning flags
   integer, parameter, public :: SSIDS_WARNING_IDX_OOR          = 1
   integer, parameter, public :: SSIDS_WARNING_DUP_IDX          = 2
   integer, parameter, public :: SSIDS_WARNING_DUP_AND_OOR      = 3
   integer, parameter, public :: SSIDS_WARNING_MISSING_DIAGONAL = 4
   integer, parameter, public :: SSIDS_WARNING_MISS_DIAG_OORDUP = 5
   integer, parameter, public :: SSIDS_WARNING_ANAL_SINGULAR    = 6
   integer, parameter, public :: SSIDS_WARNING_FACT_SINGULAR    = 7
   integer, parameter, public :: SSIDS_WARNING_MATCH_ORD_NO_SCALE=8

   ! solve job values
   integer, parameter, public :: SSIDS_SOLVE_JOB_ALL     = 0 !PLD(PL)^TX = B
   integer, parameter, public :: SSIDS_SOLVE_JOB_FWD     = 1 !PLX = B
   integer, parameter, public :: SSIDS_SOLVE_JOB_DIAG    = 2 !DX = B (indef)
   integer, parameter, public :: SSIDS_SOLVE_JOB_BWD     = 3 !(PL)^TX = B
   integer, parameter, public :: SSIDS_SOLVE_JOB_DIAG_BWD= 4 !D(PL)^TX=B (indef)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Note: below smalloc etc. types can't be in spral_ssids_alloc module as
   ! they are used as components of later datatypes.

   ! Type for custom allocator
   ! Used to aggregate many small allocations by doing a single big allocation
   ! and chopping it up.
   ! Note: Only supports freeall operation, not individual frees.
   type smalloc_type
      real(wp), dimension(:), allocatable :: rmem ! real memory
      integer(long) :: rmem_size ! needed as size(rmem,kind=long) is f2003
      integer(long) :: rhead = 0 ! last location containing useful information
         ! in rmem
      integer, dimension(:), allocatable :: imem ! integer memory
      integer(long) :: imem_size ! needed as size(imem,kind=long) is f2003
      integer(long) :: ihead = 0 ! last location containing useful information
         ! in imem
      type(smalloc_type), pointer :: next_alloc => null()
      type(smalloc_type), pointer :: top_real => null() ! Last page where real
         ! allocation was successful
      type(smalloc_type), pointer :: top_int => null() ! Last page where integer
         ! allocation was successful
!$    integer(omp_lock_kind) :: lock
   end type smalloc_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Stack memory allocation type
   type stack_mem_type
      real(wp), dimension(:), allocatable :: mem ! real memory
      integer(long) :: mem_size ! needed as size(mem,kind=long) is f2003
      integer(long) :: head = 0 ! last location containing useful information
      type(stack_mem_type), pointer :: below => null() ! next stack frame down
   end type stack_mem_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Data type for storing each node of the factors
   type node_type
      integer :: nelim
      integer :: ndelay
      integer(long) :: rdptr ! entry into (rebuilt) rlist_direct
      integer :: ncpdb ! #contrib. to parent's diag. block
      type(C_PTR) :: gpu_lcol
      real(wp), dimension(:), pointer :: lcol ! values in factors
         ! (will also include unneeded data for any columns delayed from this
         ! node)
      integer, dimension(:), pointer :: perm ! permutation of columns at this
         ! node: perm(i) is column index in expected global elimination order
         ! that is actually eliminated at local elimination index i
         ! Assuming no delays or permutation this will be
         ! sptr(node):sptr(node+1)-1
      ! Following components are used to index directly into contiguous arrays
      ! lcol and perm without taking performance hit for passing pointers
      type(smalloc_type), pointer :: rsmptr, ismptr
      integer(long) :: rsmsa, ismsa
   end type node_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for information generated in analyse phase
   !
   type ssids_akeep
      logical :: check ! copy of check as input to analyse phase
      integer :: flag ! copy of error flag.
      integer :: maxmn ! maximum value of blkm or blkn
      integer :: n ! Dimension of matrix
      integer :: ne ! Set to number of entries input by user.
      integer(long) :: nfactor 
      integer :: nnodes = -1 ! Number of nodes in assembly tree
      integer :: num_two ! This is set to 0 as we ignore any negative signs
         ! that indicate 2x2 pivot (analyse does not exploit them)

      ! child_list(child_ptr(node):child_ptr(node+1)-1) is list of children
      ! of node. Used to ensure we always sum contributions from children
      ! in the same order
      integer, dimension(:), allocatable :: child_ptr 
      integer, dimension(:), allocatable :: child_list
      integer(C_INT), dimension(:), allocatable :: invp ! inverse of pivot order
         ! that is passed to factorize phase
      integer, dimension(:), allocatable :: level ! level(i) of the assembly
         ! tree at which node i sits. root is level 1.
      integer, dimension(:,:), allocatable :: nlist ! map from A to factors
         ! For nodes i, the entries nlist(1:2, nptr(i):nptr(i+1)-1) define
         ! a relationship:
         ! nodes(node)%lcol( nlist(2,j) ) = val( nlist(1,j) )
     integer, dimension(:), allocatable :: nptr ! Entries into nlist for
         ! nodes of the assembly tree. Has length nnodes+1
      integer, dimension(:), allocatable :: rlist ! rlist(rptr(i):rptr(i+1)-1)
         ! contains the row indices for node i of the assembly tree. 
         ! At each node, the list
         ! is in elimination order. Allocated within mc78_analyse.
      integer(long), dimension(:), allocatable :: rptr ! Pointers into rlist
         ! for nodes of assembly tree. Has length nnodes+1. 
         ! Allocated within mc78_analyse.
      integer, dimension(:), allocatable :: sparent ! sparent(i) is parent
         ! of node i in assembly tree. sparent(i)=nnodes+1 if i is a root.
         ! The parent is always numbered higher than each of its children.
         ! Allocated within mc78_analyse.
      integer, dimension(:), allocatable :: sptr ! (super)node pointers.
         ! Supernode i consists of sptr(i) through sptr(i+1)-1.
         ! Allocated within mc78_analyse.
      integer(long), dimension(:), allocatable :: subtree_work ! For each node,
         ! the number of flops involved in the subtree rooted at that node.

      integer, dimension(:), allocatable :: rlist_direct ! List of rows
         ! for direct addressing of children (assuming no delays!)

      ! Following components are for cleaned up matrix data.
      ! LOWER triangle only. We have to retain these for factorize phase
      ! as used if the user wants to do scaling.
      ! These components are NOT used if check is set to .false.
      ! on call to ssids_analyse.
      integer, allocatable :: ptr(:) ! column pointers
      integer, allocatable :: row(:)! row indices
      integer :: lmap ! used by hsl_mc69
      integer, allocatable :: map(:) ! used by hsl_mc69

      ! Following components are cached members of inform
      integer :: matrix_dup
      integer :: matrix_outrange
      integer :: matrix_missing_diag
      integer :: maxdepth
      integer(long) :: num_flops ! not copied to inform in factor, but used to
         ! determine if parallelism should be used
      integer :: num_sup

      ! Scaling from matching-based ordering
      real(wp), dimension(:), allocatable :: scaling

      ! GPU pointers (copies of the CPU arrays of same name)
      type(C_PTR) :: gpu_nlist = C_NULL_PTR
      type(C_PTR) :: gpu_rlist = C_NULL_PTR
      type(C_PTR) :: gpu_rlist_direct = C_NULL_PTR
   end type ssids_akeep

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for temporary stack data that is only needed transiently during
   ! factorise phase
   ! Each instance represents a "page" of memory
   !
   type stack_type
      real(wp), dimension(:), pointer :: val => null() ! generated element
      ! Following components allow us to pass contiguous array val without
      ! taking performance hit for passing pointers
      type(stack_mem_type), pointer :: stptr => null()
      integer(long) :: stsa
   end type stack_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for per-thread stats. This is amalgamated after end of parallel
   ! section to get info parameters of same name.
   !
   type thread_stats
      integer :: flag = SSIDS_SUCCESS
      integer :: st = 0
      integer :: cuda_error = 0
      integer :: cublas_error = 0
      integer :: maxfront = 0 ! Maximum front size
      integer(long) :: num_factor = 0_long ! Number of entries in factors
      integer(long) :: num_flops = 0_long ! Number of floating point operations
      integer :: num_delay = 0 ! Number of delayed variables
      integer :: num_neg = 0 ! Number of negative pivots
      integer :: num_two = 0 ! Number of 2x2 pivots
      integer :: num_zero = 0 ! Number of zero pivots
   end type thread_stats

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! This type is used to pass buf around for each thread such that it can
   ! be reallocated independantly
   !
   type real_ptr_type
      real(wp), pointer :: chkptr => null()
      real(wp), dimension(:), allocatable :: val
   end type real_ptr_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for control parameters
   !
   type ssids_options
      !
      ! Printing options
      !
      integer :: print_level = 0 ! Controls diagnostic printing.
         ! Possible values are:
         !  < 0: no printing.
         !  0: error and warning messages only.
         !  1: as 0 plus basic diagnostic printing.
         !  > 1: as 1 plus some more detailed diagnostic messages.
      integer :: unit_diagnostics = 6 ! unit number for diagnostic printing.
         ! Printing is suppressed if unit_diagnostics  <  0.
      integer :: unit_error = 6 ! unit number for error messages.
         ! Printing is suppressed if unit_error  <  0.
      integer :: unit_warning = 6 ! unit number for warning messages.
         ! Printing is suppressed if unit_warning  <  0.

      !
      ! Options used ssids_analyse() and ssids_analyse_coord()
      !
      integer :: ordering = 1 ! controls choice of ordering
         ! 0 Order must be supplied by user
         ! 1 METIS ordering with default settings is used.
         ! 2 Matching with METIS on compressed matrix.
      integer :: nemin = nemin_default ! Min. number of eliminations at a tree
         ! node for amalgamation not to be considered.

      !
      ! Options used by ssids_factor() [both indef+posdef]
      !
      integer :: scaling = 0 ! controls use of scaling. 
         !  <=0: user supplied (or no) scaling
         !    1: Matching-based scaling by Hungarian Algorithm (MC64-like)
         !    2: Matching-based scaling by Auction Algorithm
         !    3: Scaling generated during analyse phase for matching-based order
         !  >=4: Norm equilibriation algorithm (MC77-like)

      !
      ! Options used by ssids_factor() with posdef=.false.
      !
      logical :: action = .true. ! used in indefinite case only.
         ! If true and the matrix is found to be
         ! singular, computation continues with a warning.
         ! Otherwise, terminates with error SSIDS_ERROR_SINGULAR.
      real(wp) :: small = 1e-20_wp ! Minimum pivot size (absolute value of a
         ! pivot must be of size at least small to be accepted).
      real(wp) :: u = 0.01

      !
      ! Options used by ssids_factor() and ssids_solve()
      !
      logical :: use_gpu_factor = .true. ! Use GPU for factor phase if true
         ! or CPU if false
      logical :: use_gpu_solve = .true. ! Use GPU for solve phase if true
         ! or CPU if false
      integer :: presolve = 0 ! If set to a non-zero level, triggers L-factor
         ! optimization for the sake of subsequent multiple solves.
         ! Future releases may offer different levels of optimization.

      !
      ! Undocumented
      !
      integer :: nstream = 1 ! Number of streams to use
      real(wp) :: multiplier = 1.1 ! size to multiply expected memory size by
         ! when doing initial memory allocation to allow for delays.
      type(auction_options) :: auction ! Auction algorithm parameters
      real :: min_loadbalance = 0.8 ! Minimum load balance required when
         ! finding level set used for multiple streams
   end type ssids_options

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for information returned by code
   !
   type ssids_inform_base
      integer :: flag ! Takes one of the enumerated flag values:
         ! SSIDS_SUCCESS
         ! SSIDS_ERROR_XXX
         ! SSIDS_WARNING_XXX
      integer :: matrix_dup = 0 ! Number of duplicated entries.
      integer :: matrix_missing_diag = 0 ! Number of missing diag. entries
      integer :: matrix_outrange = 0 ! Number of out-of-range entries.
      integer :: matrix_rank = 0 ! Rank of matrix (anal=structral, fact=actual)
      integer :: maxdepth ! Maximum depth of tree
      integer :: maxfront ! Maximum front size
      integer :: num_delay = 0 ! Number of delayed variables
      integer(long) :: num_factor = 0_long ! Number of entries in factors
      integer(long) :: num_flops = 0_long ! Number of floating point operations
      integer :: num_neg = 0 ! Number of negative pivots
      integer :: num_sup = 0 ! Number of supernodes
      integer :: num_two = 0 ! Number of 2x2 pivots used by factorization
      integer :: stat = 0 ! stat parameter
      type(auction_inform) :: auction
   contains
      procedure, pass(this) :: flagToCharacter
   end type ssids_inform_base

   type, extends(ssids_inform_base) :: ssids_inform
      integer :: cuda_error = 0 ! cuda error value
      integer :: cublas_error = 0 ! cublas error value
   contains
      procedure, pass(this) :: flagToCharacter => flagToCharacter_gpu
   end type ssids_inform

   integer, parameter, public :: BLOCK_SIZE = 8
   integer, parameter, public :: MNF_BLOCKS = 11
   integer, parameter, public :: HOGG_ASSEMBLE_TX = 128
   integer, parameter, public :: HOGG_ASSEMBLE_TY = 8

contains

!*************************************************
!
! Returns a string representation
! Member function inform%flagToCharacter
!
function flagToCharacter(this) result(msg)
   class(ssids_inform_base), intent(in) :: this
   character(len=200) :: msg ! return value

   select case(this%flag)
   !
   ! Success
   !
   case(SSIDS_SUCCESS)
      msg = 'Success'
   !
   ! Errors
   !
   case(SSIDS_ERROR_CALL_SEQUENCE)
      msg = 'Error in sequence of calls.'
   case(SSIDS_ERROR_A_N_OOR)
      msg = 'n or ne is out of range (or has changed)'
   case(SSIDS_ERROR_A_PTR)
      msg = 'Error in ptr'
   case(SSIDS_ERROR_A_ALL_OOR)
      msg = 'All entries in a column out-of-range (ssids_analyse) &
            &or all entries out-of-range (ssids_analyse_coord)'
   case(SSIDS_ERROR_SINGULAR)
      msg = 'Matrix found to be singular'
   case(SSIDS_ERROR_NOT_POS_DEF)
      msg = 'Matrix is not positive-definite'
   case(SSIDS_ERROR_PTR_ROW)
      msg = 'ptr and row should be present'
   case(SSIDS_ERROR_ORDER)
      msg = 'Either control%ordering out of range or error in user-supplied  &
            &elimination order'
   case(SSIDS_ERROR_X_SIZE)
      msg = 'Error in size of x or nrhs'
   case(SSIDS_ERROR_JOB_OOR)
      msg = 'job out of range'
   case(SSIDS_ERROR_NOT_LLT)
      msg = 'Not a LL^T factorization of a positive-definite matrix'
   case(SSIDS_ERROR_NOT_LDLT)
      msg = 'Not a LDL^T factorization of an indefinite matrix'
   case(SSIDS_ERROR_ALLOCATION)
      write (msg,'(a,i6)') 'Allocation error. stat parameter = ', this%stat
   case(SSIDS_ERROR_VAL)
      msg = 'Optional argument val not present when expected'
   case(SSIDS_ERROR_NO_SAVED_SCALING)
      msg = 'Requested use of scaling from matching-based &
            &ordering but matching-based ordering not used'
   case(SSIDS_ERROR_PRESOLVE_INCOMPAT)
      msg = 'Invalid combination of options%presolve, options%use_gpu_solve &
         &and requested operation - see documentation for legal combinations'

   !
   ! Warnings
   !
   case(SSIDS_WARNING_IDX_OOR)
      msg = 'out-of-range indices detected'
   case(SSIDS_WARNING_DUP_IDX)
      msg = 'duplicate entries detected'
   case(SSIDS_WARNING_DUP_AND_OOR)
      msg = 'out-of-range indices detected and duplicate entries detected'
   case(SSIDS_WARNING_MISSING_DIAGONAL)
      msg = 'one or more diagonal entries is missing'
   case(SSIDS_WARNING_MISS_DIAG_OORDUP)
      msg = 'one or more diagonal entries is missing and out-of-range and/or &
            &duplicate entries detected'
   case(SSIDS_WARNING_ANAL_SINGULAR)
      msg = 'Matrix found to be structually singular'
   case(SSIDS_WARNING_FACT_SINGULAR)
      msg = 'Matrix found to be singular'
   case(SSIDS_WARNING_MATCH_ORD_NO_SCALE)
      msg = 'Matching-based ordering used but associated scaling ignored'
   case default
      msg = 'SSIDS Internal Error'
   end select

end function flagToCharacter

!*************************************************
!
! Returns a string representation (gpu extension)
! Member function inform%flagToCharacter
!
function flagToCharacter_gpu(this) result(msg)
   class(ssids_inform), intent(in) :: this
   character(len=200) :: msg ! return value

   select case(this%flag)
   !
   ! Errors
   !
   case(SSIDS_ERROR_CUDA_UNKNOWN)
      write(msg,'(2a)') ' Unhandled CUDA error: ', &
         cudaGetErrorString(this%cuda_error)
   case(SSIDS_ERROR_CUBLAS_UNKNOWN)
      msg = 'Unhandled CUBLAS error:'
      ! FIXME?
   !
   ! Fall back to base case
   !
   case default
      msg = this%ssids_inform_base%flagToCharacter()
   end select

end function flagToCharacter_gpu

!*************************************************
!
! routine to print errors and warnings
!
subroutine ssids_print_flag(inform,nout,context)
   type(ssids_inform), intent(in) :: inform
   integer, intent(in) :: nout
   character (len=*), optional, intent(in) :: context

   character(len=200) :: msg

   if (nout < 0) return
   if (inform%flag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(context),&
         '. Error flag = ', inform%flag
   else
      write (nout,'(/3a,i3)') ' Warning from ',trim(context),&
         '. Warning flag = ', inform%flag
   end if
   msg = inform%flagToCharacter()
   write(nout, '(a)') msg

end subroutine ssids_print_flag

end module spral_ssids_datatypes
