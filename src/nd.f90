    MODULE spral_nd

      USE spral_amd

      IMPLICIT NONE
      PRIVATE

      ! ---------------------------------------------------
      ! Precision
      ! ---------------------------------------------------
      INTEGER, PARAMETER :: wp = KIND(1.0D0)
      INTEGER, PARAMETER :: long = SELECTED_INT_KIND(18)

      ! ---------------------------------------------------
      ! Error flags
      ! ---------------------------------------------------
      INTEGER, PARAMETER :: ND_ERR_MEMORY_ALLOC = -1, & ! memory alloc
      ! error
        ND_ERR_MEMORY_DEALLOC = -2, & ! memory dealloc error
        ND_ERR_N = -3 ! n<1

      ! ---------------------------------------------------
      ! Partition flags
      ! ---------------------------------------------------
      INTEGER, PARAMETER :: ND_PART1_FLAG = 0, & ! node in partition 1
        ND_PART2_FLAG = 2, & ! node in partition 2
        ND_SEP_FLAG = 1 ! node in separator

      ! ---------------------------------------------------
      ! Derived type definitions
      ! ---------------------------------------------------

      ! *****************************************************************
      TYPE, PUBLIC :: nd_options
        INTEGER :: print_level = 0 ! amount of informational output required
        INTEGER :: unit_diagnostics = 6 ! stream number for diagnostic
        ! messages
        INTEGER :: unit_error = 6 ! stream number for error messages

        INTEGER :: amd_call = 5 ! call AMD if number of supervariables is
        ! less than or equal to amd_call  *UPDATE*
        INTEGER :: amd_switch1 = 50 ! switch to (halo)AMD if matrix size is
        ! less
        ! than or equal to nd_switch
        INTEGER :: amd_switch2 = 20 ! maximum number of ND levels allowed
        INTEGER :: cost_function = 1 ! selects the cost function used for 
        ! evaluating a partition
        ! <=1 : |P1||P2|/|S| + penalty for imbalance
        ! >=2 : |S|(1 +  0.5||P1|-|P2||)
        INTEGER :: partition_method = 2 ! Are we allowed to use a multilevel
        ! strategy
        ! <= 0 : do not use multilevel
        ! == 1 : use multilevel
        ! >= 2 : automatic choice based on size of levelsets
        INTEGER :: matching = 0 ! Which coarsening method to use
        ! > 0 : heavy-edge
        ! <= 0 : common neighbours
        INTEGER :: coarse_partition_method = 1 ! Which partition method to use
        ! at coarsest level
        ! <=1 : Ashcraft method (half-level set)
        ! >=2 : Level-set method
        INTEGER :: refinement = 7 ! Which sort of refinement to use
        ! <2 : trim + fm to increase weights of A1 and A2
        ! =2 : trim + fm to increase weight of smaller partition and
        ! decrease weight of larger partition
        ! =3 : Automatically choose between options 1 and 2
        ! =4 : maxflow + fm to increase weights of A1 and A2
        ! =5 : maxflow + fm to increase weight of smaller partition and
        ! decrease weight of larger partition
        ! =6 : Automatically choose between options 4 and 5
        ! >6 : Automatically choose between options 1 and 5
        INTEGER :: refinement_band = 4 ! band width for FM refinement. Values
        ! less than 1 mean that no FM refinement is done
        LOGICAL :: remove_dense_rows = .TRUE. ! test the input for dense rows
        ! and place them at the end of the ordering
        INTEGER :: stop_coarsening2 = 20 ! Max number of levels in the
        ! multilevel grid
        INTEGER :: stop_coarsening1 = 100 ! Stop coarsening once matrix has
        ! order at most stop_coarsening1
        INTEGER :: ml_call = 12000 ! Stop coarsening once matrix has
        ! order at most stop_coarsening1

        ! minimum and maximum grid reduction factors that must be achieved
        ! during coarsening. If cgrid%size is greater than
        ! max_reduction*grid%size or cgrid%size is less than
        ! min_reduction*grid%size then carry on coarsening
        REAL (wp) :: min_reduction = 0.5 ! size of next
        ! multigrid
        ! matrix must be greater than min_reduction*(size of current
        ! multigrid matrix)
        REAL (wp) :: max_reduction = 0.9 ! size of next
        ! multigrid
        ! matrix must be less than max_reduction*(size of current multigrid
        ! matrix)
        REAL (wp) :: balance = 2.0 ! Try to make sure that
        ! max(P1,P2)/min(P1/P2) <= balance

        INTEGER :: max_improve_cycles = 2 ! Having computed a minimal
        ! partition,
        ! expand and refine it at most max_improve_cycles times to improve
        ! the quality of the partition.
        LOGICAL :: find_supervariables = .TRUE. ! If .TRUE., after dense rows
        ! have been (optionally) removed, check for supervariables and
        ! compress matrix if supervariables found.

        ! REAL (wp) :: ml_bandwidth = 1.0 ! Let B be matrix A after
        ! dense rows removed and compressed. If ml>1 and bandwidth of B after
        ! RCM
        ! ordering is larger than ml_bandwidth*order(B), continue as if ml=1;
        ! otherwise continue as if ml=0; If B is separable, perform test on
        ! each
        ! component separately use accordingly with each component. Note: RCM
        ! ordering is not computed.

      END TYPE nd_options

      ! *****************************************************************
      TYPE, PUBLIC :: nd_inform
        INTEGER :: flag = 0 ! error/warning flag
        INTEGER :: dense = 0 ! holds number of dense rows
        INTEGER :: stat = 0 ! holds Fortran stat parameter
        INTEGER :: nsuper = 0 ! holds number of supervariables + number
        ! of zero
        ! rows after dense rows removed
        INTEGER :: nzsuper = 0 ! holds number of nonzeros in compressed
        ! graph
        INTEGER :: num_components = 1 ! Number of independent components after
        ! dense rows (optionally) removed
        INTEGER :: n_max_component = -1 ! holds number rows in largest indep
        ! component after dense rows removed and supervariables (optionally)
        ! compressed
        INTEGER :: nz_max_component = -1 ! holds number nonzeros in largest
        ! indep
        ! component after dense rows removed and supervariables (optionally)
        ! compressed
        INTEGER :: maxdeg_max_component = -1 ! holds number nonzeros in
        ! largest indep
        ! component after dense rows removed and supervariables (optionally)
        ! compressed
        REAL (wp) :: band = -1 ! holds L, where L is the size
        ! of the largest level set at the top level of nested dissection. If
        ! the matrix is reducible, then it holds the value for the largest
        ! of the irreducible components.
        ! Not returned if control%partition_method==1.
        REAL (wp) :: depth = -1 ! holds number of levels in
        ! level set
        ! structure at the top level of nested dissection. If
        ! the matrix is reducible, then it holds the value for the largest
        ! of the irreducible components.
        ! Not returned if control%partition_method==1.
      END TYPE nd_inform

      ! *****************************************************************

      TYPE nd_multigrid
        INTEGER :: size ! size of this level (number of rows)
        TYPE (nd_matrix), POINTER :: graph => NULL() ! this level of matrix
        INTEGER, POINTER, DIMENSION (:) :: where => NULL() ! where each row of
        ! this
        ! level of matrix will go (ie ordering for this level)
        INTEGER, POINTER, DIMENSION (:) :: row_wgt => NULL() ! number of
        ! vertices
        ! this vertex of the coarse graph matrix represents
        INTEGER :: level = 0 ! the level
        INTEGER :: part_div(2) ! number of vertices in each part
        TYPE (nd_multigrid), POINTER :: coarse => NULL() ! pointer to the
        ! coarse grid
        TYPE (nd_multigrid), POINTER :: fine => NULL() ! pointer to the
        ! fine grid
        TYPE (nd_matrix), POINTER :: p => NULL() ! the prolongation
        ! operator
        TYPE (nd_matrix), POINTER :: r => NULL() ! the restriction operator
      END TYPE nd_multigrid
      ! *****************************************************************

      TYPE nd_matrix
        INTEGER :: m ! number rows
        INTEGER :: n ! number columns
        INTEGER :: ne ! number entries in matrix
        INTEGER, ALLOCATABLE, DIMENSION (:) :: ptr ! pointer into col array
        INTEGER, ALLOCATABLE, DIMENSION (:) :: col ! column indices
        INTEGER, ALLOCATABLE, DIMENSION (:) :: val ! values

      END TYPE nd_matrix

      ! This is a version of maxflow with no bandgraph and assumption that
      ! there
      ! are no duplicated edges.

      TYPE network
        INTEGER :: nnode ! number of nodes in the network
        INTEGER :: narc ! number of arcs in the network
        INTEGER :: source ! source node = 1
        INTEGER :: sink ! sink node = nnode
        INTEGER, ALLOCATABLE :: inheads(:) ! for each node u
        ! first incoming arc for u
        INTEGER, ALLOCATABLE :: outheads(:) ! for each node u
        ! first outgoing arc for u
        INTEGER, ALLOCATABLE :: firsts(:) ! for each arc e
        ! first node for arc e
        INTEGER, ALLOCATABLE :: seconds(:) ! for each arc e
        ! second node for arc e
        INTEGER, ALLOCATABLE :: capacities(:) ! for each arc e
        ! capacity of arc e
        INTEGER, ALLOCATABLE :: flows(:) ! for each arc e
        ! flow through arc e
        INTEGER, ALLOCATABLE :: nextin(:) ! for each arc e
        ! next incoming arc
        INTEGER, ALLOCATABLE :: nextout(:) ! for each arc e
        ! next outgoing arc
      END TYPE network


      ! ---------------------------------------------------
      ! Interfaces
      ! ---------------------------------------------------
      INTERFACE nd_order
        MODULE PROCEDURE nd_nested_order
      END INTERFACE

      PUBLIC nd_order, nd_refine_fm

      ! ---------------------------------------------------
      ! The main code
      ! ---------------------------------------------------

    CONTAINS

      ! ---------------------------------------------------
      ! nd_print_message
      ! ---------------------------------------------------
      SUBROUTINE nd_print_message(flag,unit,context)
        ! Prints out errors and warnings according to value of flag

        ! flag: is an integer scaler of intent(in). It is the information flag
        ! whose corresponding error message is printed
        INTEGER, INTENT (IN) :: flag

        ! unit: is an optional integer scaler of intent(in). It is the unit
        ! number the
        ! error/warning message should be printed on
        INTEGER, INTENT (IN) :: unit

        ! context: is an optional assumed size character array of intent(in).
        ! It describes the context under which the error occured
        CHARACTER (len=*), INTENT (IN) :: context
        INTEGER :: length, p_unit

        p_unit = unit

        IF (p_unit<=0) RETURN

        IF (flag<0) THEN
          WRITE (p_unit,advance='yes',fmt='('' ERROR: '')')
        END IF

        length = LEN_TRIM(context)
        WRITE (p_unit,advance='no',fmt='('' '', a,'':'')') context(1:length)

        SELECT CASE (flag)
        CASE (0)
          WRITE (p_unit,'(A)') ' successful completion'

        CASE (ND_ERR_MEMORY_ALLOC)
          WRITE (p_unit,'(A)') ' memory allocation failure'

        CASE (ND_ERR_MEMORY_DEALLOC)
          WRITE (p_unit,'(A)') ' memory deallocation failure'

        CASE (ND_ERR_N)
          WRITE (p_unit,'(A)') ' n<1'
        END SELECT
        RETURN

      END SUBROUTINE nd_print_message


      SUBROUTINE nd_nested_order(mtx,n,ptr,row,perm,control,info,seps)
        INTEGER, INTENT (IN) :: mtx ! 0 if lower triangular part matrix input,
        ! 1
        ! if both upper and lower parts input
        INTEGER, INTENT (IN) :: n ! number of rows in the matrix
        INTEGER, INTENT (IN) :: ptr(n+1) ! column pointers
        INTEGER, INTENT (IN) :: row(ptr(n+1)-1) ! row indices (lower triangle)
        INTEGER, INTENT (OUT) :: perm(n) ! permutation: row i becomes row
        ! perm(i)
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info
        INTEGER, INTENT (OUT), OPTIONAL :: seps(n)
        ! seps(i) is -1 if vertex is not in a separator; otherwise it
        ! is equal to l, where l is the nested dissection level at
        ! which it became part of the separator

        IF (mtx<1) THEN
          IF (PRESENT(seps)) THEN
            CALL nd_nested_lower(n,ptr,row,perm,control,info,seps)
          ELSE
            CALL nd_nested_lower(n,ptr,row,perm,control,info)
          END IF
        ELSE
          IF (PRESENT(seps)) THEN
            CALL nd_nested_full(n,ptr,row,perm,control,info,seps)
          ELSE
            CALL nd_nested_full(n,ptr,row,perm,control,info)
          END IF
        END IF

      END SUBROUTINE nd_nested_order

      ! ---------------------------------------------------
      ! nd_nested_lower
      ! ---------------------------------------------------
      SUBROUTINE nd_nested_lower(n,ptr,row,perm,control,info,seps)
        INTEGER, INTENT (IN) :: n ! number of rows in the matrix
        INTEGER, INTENT (IN) :: ptr(n+1) ! column pointers
        INTEGER, INTENT (IN) :: row(ptr(n+1)-1) ! row indices (lower triangle)
        INTEGER, INTENT (OUT) :: perm(n) ! permutation: row i becomes row
        ! perm(i)
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info
        INTEGER, INTENT (OUT), OPTIONAL :: seps(n)
        ! seps(i) is -1 if vertex is not in a separator; otherwise it
        ! is equal to l, where l is the nested dissection level at
        ! which it became part of the separator

        ! ---------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: a_n ! dimension of the expanded matrix
        INTEGER :: a_ne ! number off-diagonal entries stored in
        ! expanded matrix
        INTEGER :: i, j, k
        INTEGER :: unit_error ! unit on which to print errors
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER, ALLOCATABLE, DIMENSION (:) :: a_ptr ! a_ptr(i) will contain
        ! the
        ! position in a_row that column i ends for the expanded matrix
        INTEGER, ALLOCATABLE, DIMENSION (:) :: a_row ! contains for row
        ! indices
        ! for the pattern of the expanded matrix. The row indices of
        ! column j are stored before those of column j+1, j=1,..,n-1.
        LOGICAL :: printe, printi, printd

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        unit_error = control%unit_error
        printe = (control%print_level>=0 .AND. unit_error>=0)
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
        ! ---------------------------------------------------

        ! Expand the matrix to hold all of its pattern
        ! Throughout the code we will use a modified compressed sparse
        ! column/row
        ! format a_ptr(i) will contain the position in a_row that column i
        ! begins. This
        ! is to avoid lots of allocations during the nested part of the
        ! algorithm. a_ne
        ! contains the number of entries stored for the expanded matrix.

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'nd_nested_lower:'
        END IF



        ! Set the dimension of the expanded matrix
        a_n = n
        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a,i10)') 'n = ', n
        END IF
        IF (n<1) THEN
          info%flag = ND_ERR_N
          IF (printe) CALL nd_print_message(info%flag,unit_error, &
            'nd_nested_lower')
          RETURN
        END IF


        ! Allocate space to store pointers for expanded matrix
        ALLOCATE (a_ptr(a_n),STAT=info%stat)
        IF (info%stat/=0) GO TO 10

        ! Fill a_col and a_ptr removing any diagonal entries
        a_ptr(:) = 0

        ! Set a_ptr(j) to hold no. nonzeros in column j
        DO j = 1, n
          DO k = ptr(j), ptr(j+1) - 1
            i = row(k)
            IF (j/=i) THEN
              a_ptr(i) = a_ptr(i) + 1
              a_ptr(j) = a_ptr(j) + 1
            END IF
          END DO
        END DO

        ! Set a_ptr(j) to point to where row indices will end in a_row
        DO j = 2, n
          a_ptr(j) = a_ptr(j-1) + a_ptr(j)
        END DO
        a_ne = a_ptr(a_n)
        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a,i10)') &
            'entries in expanded matrix with diags removed = ', a_ne
        END IF

        ! Allocate space to store row indices of expanded matrix
        ALLOCATE (a_row(a_ne),STAT=info%stat)
        IF (info%stat/=0) GO TO 10

        ! Initialise all of a_row to 0
        a_row(:) = 0

        ! Fill a_row and a_ptr
        DO j = 1, n
          DO k = ptr(j), ptr(j+1) - 1
            i = row(k)
            IF (j/=i) THEN
              a_row(a_ptr(i)) = j
              a_row(a_ptr(j)) = i
              a_ptr(i) = a_ptr(i) - 1
              a_ptr(j) = a_ptr(j) - 1
            END IF
          END DO
        END DO

        ! Reset a_ptr to point to where column starts
        DO j = 1, a_n
          a_ptr(j) = a_ptr(j) + 1
        END DO

        IF (printd) THEN
          ! Print out a_ptr and a_row
          WRITE (unit_diagnostics,'(a8)') 'a_ptr = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n)
          WRITE (unit_diagnostics,'(a8)') 'a_row = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne)
        ELSE IF (printi) THEN
          ! Print out first few entries of a_ptr and a_row
          WRITE (unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,MIN(5,a_n))
          WRITE (unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,MIN(5,a_ne))
        END IF
        IF (PRESENT(seps)) THEN
          CALL nd_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info,seps)
        ELSE
          CALL nd_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info)
        END IF

        IF (printd) THEN
          ! Print out perm
          WRITE (unit_diagnostics,'(a8)') 'perm = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,a_n)
        ELSE IF (printi) THEN
          ! Print out first few entries of perm
          WRITE (unit_diagnostics,'(a21)') 'perm(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,MIN(5,a_n))
        END IF

        ! Deallocate arrays
        DEALLOCATE (a_ptr,STAT=info%stat)
        IF (info%stat/=0) GO TO 20
        DEALLOCATE (a_row,STAT=info%stat)
        IF (info%stat/=0) GO TO 20

        info%flag = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info%flag,unit_diagnostics,'nd_nested_lower')
        END IF
        RETURN

10      info%flag = ND_ERR_MEMORY_ALLOC
        IF (printe) CALL nd_print_message(info%flag,unit_error, &
          'nd_nested_lower')
        RETURN

20      info%flag = ND_ERR_MEMORY_DEALLOC
        IF (printe) CALL nd_print_message(info%flag,unit_error, &
          'nd_nested_lower')
        RETURN

      END SUBROUTINE nd_nested_lower


      ! ---------------------------------------------------
      ! nd_nested
      ! ---------------------------------------------------
      SUBROUTINE nd_nested_full(n,ptr,row,perm,control,info,seps)
        INTEGER, INTENT (IN) :: n ! number of rows in the matrix
        INTEGER, INTENT (IN) :: ptr(n+1) ! column pointers
        INTEGER, INTENT (IN) :: row(ptr(n+1)-1) ! row indices (lower triangle)
        INTEGER, INTENT (OUT) :: perm(n) ! permutation: row i becomes row
        ! perm(i)
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info
        INTEGER, INTENT (OUT), OPTIONAL :: seps(n)
        ! seps(i) is -1 if vertex is not in a separator; otherwise it
        ! is equal to l, where l is the nested dissection level at
        ! which it became part of the separator

        ! ---------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: a_n ! dimension of the expanded matrix
        INTEGER :: a_ne ! number off-diagonal entries stored in
        ! expanded matrix
        INTEGER :: i, j, k, l, ndiags
        INTEGER :: unit_error ! unit on which to print errors
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER, ALLOCATABLE, DIMENSION (:) :: a_ptr ! a_ptr(i) will contain
        ! the
        ! position in a_row that column i ends for the expanded matrix
        INTEGER, ALLOCATABLE, DIMENSION (:) :: a_row ! contains for row
        ! indices
        ! for the pattern of the expanded matrix. The row indices of
        ! column j are stored before those of column j+1, j=1,..,n-1.
        LOGICAL :: printe, printi, printd

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        unit_error = control%unit_error
        printe = (control%print_level>=0 .AND. unit_error>=0)
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
        ! ---------------------------------------------------

        ! Expand the matrix to hold all of its pattern
        ! Throughout the code we will use a modified compressed sparse
        ! column/row
        ! format a_ptr(i) will contain the position in a_row that column i
        ! begins. This
        ! is to avoid lots of allocations during the nested part of the
        ! algorithm. a_ne
        ! contains the number of entries stored for the expanded matrix.

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'nd_nested_full:'
        END IF

        ! Set the dimension of the expanded matrix
        a_n = n
        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a,i10)') 'n = ', n
        END IF
        IF (n<1) THEN
          info%flag = ND_ERR_N
          IF (printe) CALL nd_print_message(info%flag,unit_error, &
            'nd_nested_full')
          RETURN
        END IF


        ! Work out how many diagonal entries need removing
        ndiags = 0
        DO j = 1, n
          DO l = ptr(j), ptr(j+1) - 1
            i = row(l)
            IF (i==j) ndiags = ndiags + 1
          END DO
        END DO
        a_ne = ptr(n+1) - 1 - ndiags

        ! Allocate space to store pointers and rows for expanded matrix
        ALLOCATE (a_ptr(a_n),STAT=info%stat)
        IF (info%stat/=0) GO TO 10
        ALLOCATE (a_row(a_ne),STAT=info%stat)
        IF (info%stat/=0) GO TO 10

        IF (ndiags==0) THEN
          ! No diagonal entries so do direct copy
          a_ptr(1:a_n) = ptr(1:n)
          a_row(1:a_ne) = row(1:a_ne)

        ELSE
          ! Diagonal entries present
          k = 1
          DO i = 1, n
            a_ptr(i) = k
            DO l = ptr(i), ptr(i+1) - 1
              j = row(l)
              IF (i/=j) THEN
                a_row(k) = j
                k = k + 1
              END IF
            END DO
          END DO
        END IF


        IF (printd) THEN
          ! Print out a_ptr and a_row
          WRITE (unit_diagnostics,'(a8)') 'a_ptr = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n)
          WRITE (unit_diagnostics,'(a8)') 'a_row = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne)
        ELSE IF (printi) THEN
          ! Print out first few entries of a_ptr and a_row
          WRITE (unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,MIN(5,a_n))
          WRITE (unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,MIN(5,a_ne))
        END IF
        IF (PRESENT(seps)) THEN
          CALL nd_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info,seps)
        ELSE
          CALL nd_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info)
        END IF

        IF (printd) THEN
          ! Print out perm
          WRITE (unit_diagnostics,'(a8)') 'perm = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,a_n)
        ELSE IF (printi) THEN
          ! Print out first few entries of perm
          WRITE (unit_diagnostics,'(a21)') 'perm(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,MIN(5,a_n))
        END IF

        ! Deallocate arrays
        DEALLOCATE (a_ptr,STAT=info%stat)
        IF (info%stat/=0) GO TO 20
        DEALLOCATE (a_row,STAT=info%stat)
        IF (info%stat/=0) GO TO 20

        info%flag = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info%flag,unit_diagnostics, &
            'nd_nested_full')
        END IF
        RETURN

10      info%flag = ND_ERR_MEMORY_ALLOC
        IF (printe) CALL nd_print_message(info%flag,unit_error, &
          'nd_nested_full')
        RETURN

20      info%flag = ND_ERR_MEMORY_DEALLOC
        IF (printe) CALL nd_print_message(info%flag,unit_error, &
          'nd_nested_full')
        RETURN

      END SUBROUTINE nd_nested_full

      ! ---------------------------------------------------
      ! nd_nested
      ! ---------------------------------------------------
      SUBROUTINE nd_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info, &
          seps)
        INTEGER, INTENT (IN) :: a_n ! number of rows in the matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in the matrix
        INTEGER, INTENT (INOUT) :: a_ptr(a_n) ! column pointers
        INTEGER, INTENT (INOUT) :: a_row(a_ne) ! row indices (lower and upper)
        INTEGER, INTENT (OUT) :: perm(a_n) ! permutation: row i becomes row
        ! perm(i)
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info
        INTEGER, INTENT (OUT), OPTIONAL :: seps(a_n)
        ! seps(i) is -1 if vertex is not in a separator; otherwise it
        ! is equal to l, where l is the nested dissection level at
        ! which it became part of the separator

        ! ---------------------------------------------------
        ! LOCAL VARIABLES
        INTEGER :: i, j, k, l, ll, a_ne_new, a_n_new, lwork, lirn
        INTEGER :: a_n_curr, a_ne_curr, num_zero_row
        INTEGER :: amd_order_irn, amd_order_ip, amd_order_sep, amd_order_perm, amd_order_work, &
          amd_order_iperm
        INTEGER :: work_iperm, work_seps
        INTEGER :: unit_error ! unit on which to print errors
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: st, nsvar, svinfo
        INTEGER :: sv_ptr, sv_perm, sv_invp, sv_svar, sv_ptr2, sv_row2, &
          sumweight
        INTEGER, ALLOCATABLE, DIMENSION (:) :: a_weight ! a_weight(i) will
        ! contain the weight of variable (column) i ends for the
        ! expanded matrix
        INTEGER, ALLOCATABLE, DIMENSION (:) :: iperm ! row iperm(i) will
        ! become row i when matrix reordered
        INTEGER, ALLOCATABLE, DIMENSION (:) :: work ! space for doing work
        INTEGER, ALLOCATABLE, DIMENSION (:) :: svwork ! supervariable work
        ! space
        LOGICAL :: printe, printi, printd
        LOGICAL :: use_multilevel
        TYPE (nd_multigrid) :: grid

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        unit_error = control%unit_error
        printe = (control%print_level>=0 .AND. unit_error>=0)
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
        ! ---------------------------------------------------

        ! Expand the matrix to hold all of its pattern
        ! Throughout the code we will use a modified compressed sparse
        ! column/row
        ! format a_ptr(i) will contain the position in a_row that column i
        ! begins. This
        ! is to avoid lots of allocations during the nested part of the
        ! algorithm. a_ne
        ! contains the number of entries stored for the expanded matrix.

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'nd_nested_both:'
        END IF
        work_seps = 0
        work_iperm = 0

        ! Allocate iperm to have length n
        ALLOCATE (iperm(a_n),STAT=info%stat)
        IF (info%stat/=0) GO TO 10

        ! Initialise iperm to identity permutation
        iperm(:) = (/ (i,i=1,a_n) /)
        IF (PRESENT(seps)) THEN
          seps(:) = -1
        END IF

        ! Check for dense rows
        IF (control%remove_dense_rows) THEN
          ! Initially allocate work to have length 4*a_n
          ALLOCATE (work(4*a_n),STAT=info%stat)
          IF (info%stat/=0) GO TO 10

          CALL nd_dense_rows(a_n,a_ne,a_ptr,a_row,i,j,iperm,work(1:4*a_n), &
            control,info)
          a_n_new = i
          a_ne_new = j

          IF (printd) THEN
            WRITE (unit_diagnostics,'(a,i10)') ' No. dense rows removed: ', &
              a_n - a_n_new
          END IF

          DEALLOCATE (work,STAT=info%stat)
          IF (info%stat/=0) GO TO 20

        ELSE
          a_n_new = a_n
          a_ne_new = a_ne

        END IF

        ! Check whether matrix is diagonal
        IF (a_ne_new==0) THEN

          IF (printi .OR. printd) THEN
            WRITE (unit_diagnostics,'(a)') ' '
            WRITE (unit_diagnostics,'(a)') 'Matrix is diagonal'
          END IF

          info%nsuper = a_n_new
          info%nzsuper = 0
          ! Create perm from iperm
          DO i = 1, a_n
            j = iperm(i)
            perm(j) = i
          END DO

          RETURN
        END IF

        IF (control%find_supervariables) THEN
          ! Check for supervariables
          ALLOCATE (svwork(5*a_n_new+a_ne_new+2),STAT=info%stat)
          IF (info%stat/=0) GO TO 10
          sv_ptr = 0
          sv_perm = sv_ptr + a_n_new + 1
          sv_invp = sv_perm + a_n_new
          sv_svar = sv_invp + a_n_new
          sv_ptr2 = sv_svar + a_n_new
          sv_row2 = sv_ptr2 + a_n_new + 1


          svwork(sv_ptr+1:sv_ptr+a_n_new) = a_ptr(1:a_n_new)
          svwork(sv_ptr+a_n_new+1) = a_ne_new + 1
          svwork(sv_perm+1:sv_perm+a_n_new) = (/ (i,i=1,a_n_new) /)
          svwork(sv_invp+1:sv_invp+a_n_new) = (/ (i,i=1,a_n_new) /)
          i = a_n_new
          CALL nd_supervars(i,svwork(sv_ptr+1:sv_ptr+a_n_new+1), &
            a_row(1:a_ne_new),svwork(sv_perm+1:sv_perm+a_n_new), &
            svwork(sv_invp+1:sv_invp+a_n_new),nsvar, &
            svwork(sv_svar+1:sv_svar+a_n_new),svinfo,st)
          IF (svinfo==ND_ERR_MEMORY_ALLOC) GO TO 10
          IF (svinfo==ND_ERR_MEMORY_DEALLOC) GO TO 20

          num_zero_row = a_n_new - i
          IF (printd) THEN
            WRITE (unit_diagnostics,'(a,i10)') 'Number supervariables: ', &
              nsvar + num_zero_row
          END IF
          IF (nsvar+num_zero_row==a_n_new) THEN
            DEALLOCATE (svwork,STAT=info%stat)
            IF (info%stat/=0) GO TO 20
            a_n_curr = a_n_new
            a_ne_curr = a_ne_new
            ALLOCATE (a_weight(a_n_curr),STAT=info%stat)
            IF (info%stat/=0) GO TO 10

            ! Initialise a_weight
            a_weight(:) = 1

          ELSE
            CALL nd_compress_by_svar(a_n_new,a_ne_new, &
              svwork(sv_ptr+1:sv_ptr+a_n_new+1),a_row(1:a_ne_new), &
              svwork(sv_invp+1:sv_invp+a_n_new),nsvar, &
              svwork(sv_svar+1:sv_svar+a_n_new),svwork(sv_ptr2+1:sv_ptr2+ &
              a_n_new+1),svwork(sv_row2+1:sv_row2+a_ne_new),svinfo,st)
            IF (svinfo==ND_ERR_MEMORY_ALLOC) GO TO 10
            IF (svinfo==ND_ERR_MEMORY_DEALLOC) GO TO 20

            a_n_curr = nsvar
            ! Fill a_ptr removing any diagonal entries
            a_ptr(:) = 0

            ! Set a_ptr(j) to hold no. nonzeros in column j
            DO j = 1, a_n_curr
              DO k = svwork(sv_ptr2+j), svwork(sv_ptr2+j+1) - 1
                i = svwork(sv_row2+k)
                IF (j<i) THEN
                  a_ptr(i) = a_ptr(i) + 1
                  a_ptr(j) = a_ptr(j) + 1
                END IF
              END DO
            END DO

            ! Set a_ptr(j) to point to where row indices will end in a_row
            DO j = 2, a_n_curr
              a_ptr(j) = a_ptr(j-1) + a_ptr(j)
            END DO
            a_ne_curr = a_ptr(a_n_curr)
            ! Initialise all of a_row to 0
            a_row(1:a_ne_curr) = 0

            ! Fill a_row and a_ptr
            DO j = 1, a_n_curr
              DO k = svwork(sv_ptr2+j), svwork(sv_ptr2+j+1) - 1
                i = svwork(sv_row2+k)
                IF (j<i) THEN
                  a_row(a_ptr(i)) = j
                  a_row(a_ptr(j)) = i
                  a_ptr(i) = a_ptr(i) - 1
                  a_ptr(j) = a_ptr(j) - 1
                END IF
              END DO
            END DO

            ! Reset a_ptr to point to where column starts
            DO j = 1, a_n_curr
              a_ptr(j) = a_ptr(j) + 1
            END DO

            ALLOCATE (a_weight(a_n_curr+num_zero_row),STAT=info%stat)
            IF (info%stat/=0) GO TO 10

            ! Initialise a_weight
            a_weight(1:a_n_curr) = svwork(sv_svar+1:sv_svar+a_n_curr)
            a_weight(a_n_curr+1:a_n_curr+num_zero_row) = 1

            ! Add zero rows/cols to matrix
            a_ptr(a_n_curr+1:a_n_curr+num_zero_row) = a_ne_curr + 1
            a_n_curr = a_n_curr + num_zero_row

            ! set svwork(sv_svar+1:sv_svar+a_n_new) such that
            ! svwork(sv_svar+i)
            ! points to the end of the list of variables in sv_invp for
            ! supervariable i
            DO i = 2, nsvar
              svwork(sv_svar+i) = svwork(sv_svar+i) + svwork(sv_svar+i-1)
            END DO
            j = svwork(sv_svar+nsvar)
            DO i = 1, num_zero_row
              svwork(sv_svar+nsvar+i) = j + 1
              j = j + 1
            END DO

          END IF

        ELSE
          a_n_curr = a_n_new
          a_ne_curr = a_ne_new
          nsvar = a_n_new
          num_zero_row = 0
          ALLOCATE (a_weight(a_n_curr+num_zero_row),STAT=info%stat)
          IF (info%stat/=0) GO TO 10

          ! Initialise a_weight
          a_weight(:) = 1
        END IF

        ! Carryout nested dissection on matrix once dense rows removed
        info%nzsuper = a_ne_curr
        info%nsuper = a_n_curr

        IF (control%amd_switch2<=0 .OR. a_n_curr<=MAX(2,MAX(control%amd_call, &
            control%amd_switch1))) THEN
          ! Apply AMD to matrix
          ! Allocate work to have length 5*a_n_curr+a_ne_curr
          IF (printd) THEN
            WRITE (unit_diagnostics,'(a)') ' Form AMD ordering'
          END IF
          ALLOCATE (work(11*a_n_curr+a_ne_curr+a_n_new),STAT=info%stat)
          IF (info%stat/=0) GO TO 10
          lirn = a_ne_curr + a_n_curr
          amd_order_irn = 0
          amd_order_ip = amd_order_irn + lirn
          amd_order_sep = amd_order_ip + a_n_curr
          amd_order_perm = amd_order_sep + a_n_curr
          amd_order_work = amd_order_perm + a_n_curr
          amd_order_iperm = amd_order_work + 7*a_n_curr

          work(amd_order_irn+1:amd_order_irn+a_ne_curr) = a_row(1:a_ne_curr)
          work(amd_order_irn+a_ne_curr+1:amd_order_irn+lirn) = 0
          work(amd_order_ip+1:amd_order_ip+a_n_curr) = a_ptr(1:a_n_curr)
          work(amd_order_sep+1:amd_order_sep+a_n_curr) = ND_SEP_FLAG - 1
          CALL amd_order(a_n_curr,a_ne_curr,lirn,work(amd_order_irn+1:amd_order_irn+lirn), &
            work(amd_order_ip+1:amd_order_ip+a_n_curr),work(amd_order_sep+1:amd_order_sep+a_n_curr &
            ),work(amd_order_perm+1:amd_order_perm+a_n_curr),work(amd_order_work+1:amd_order_work+ &
            7*a_n_curr))

          ! Extract perm from amd_order to apply to iperm
          IF (nsvar+num_zero_row==a_n_new) THEN
            DO i = 1, a_n_curr
              j = work(amd_order_perm+i)
              work(amd_order_work+i) = iperm(j)
            END DO
            iperm(1:a_n_curr) = work(amd_order_work+1:amd_order_work+a_n_curr)
          ELSE
            DO i = 1, a_n_curr
              j = work(amd_order_perm+i)
              work(amd_order_work+i) = j
            END DO
            ! Expand to matrix before supervariables detected
            k = 1
            DO i = 1, a_n_curr
              j = work(amd_order_work+i)
              IF (j==1) THEN
                ll = 1
              ELSE
                ll = svwork(sv_svar+j-1) + 1
              END IF
              DO l = ll, svwork(sv_svar+j)
                work(amd_order_iperm+k) = iperm(svwork(sv_invp+l))
                k = k + 1
              END DO
            END DO

            iperm(1:a_n_new) = work(amd_order_iperm+1:amd_order_iperm+a_n_new)
          END IF

        ELSE
          ! Apply ND to matrix

          IF (printd) THEN
            WRITE (unit_diagnostics,'(a)') ' Form ND ordering'
          END IF

          IF (nsvar+num_zero_row==a_n_new) THEN
            ! Allocate work to have length 14*a_n_curr+a_ne_curr
            ALLOCATE (work(a_n_new+14*a_n_curr+a_ne_curr),STAT=info%stat)
          ELSE
            ! Allocate work to have length a_n_new+14*a_n_curr+a_ne_curr
            ALLOCATE (work(3*a_n_new+14*a_n_curr+a_ne_curr),STAT=info%stat)
            work_iperm = 14*a_n_curr + a_ne_curr + a_n_new
            work(work_iperm+1:work_iperm+a_n_curr) = (/ (i,i=1,a_n_curr) /)
            IF (PRESENT(seps)) THEN
              work_seps = work_iperm + a_n_curr
              work(work_seps+1:work_seps+a_n_curr) = -1
            END IF
          END IF
          IF (info%stat/=0) GO TO 10

          use_multilevel = .TRUE.
          IF (nsvar+num_zero_row==a_n_new) THEN
            sumweight = SUM(a_weight(1:a_n_curr))
            lwork = 12*a_n_curr + sumweight + a_ne_curr
            IF (PRESENT(seps)) THEN
              CALL nd_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr), &
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight, &
                iperm(1:a_n_curr),work(1:lwork),work(lwork+1:lwork+a_n_curr), &
                work(lwork+a_n_curr+1:lwork+2*a_n_curr),0,control,info, &
                .FALSE.,use_multilevel,grid,seps(1:a_n_curr))
            ELSE
              CALL nd_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr), &
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight, &
                iperm(1:a_n_curr),work(1:lwork),work(lwork+1:lwork+a_n_curr), &
                work(lwork+a_n_curr+1:lwork+2*a_n_curr),0,control,info, &
                .FALSE.,use_multilevel,grid)

            END IF
          ELSE
            sumweight = SUM(a_weight(1:a_n_curr))
            lwork = 12*a_n_curr + sumweight + a_ne_curr
            IF (PRESENT(seps)) THEN
              CALL nd_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr), &
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight, &
                work(work_iperm+1:work_iperm+a_n_curr),work(1:lwork), &
                work(lwork+1:lwork+a_n_curr),work(lwork+a_n_curr+1:lwork+2* &
                a_n_curr),0,control,info,.FALSE.,use_multilevel,grid, &
                work(work_seps+1:work_seps+a_n_curr))
            ELSE
              CALL nd_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr), &
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight, &
                work(work_iperm+1:work_iperm+a_n_curr),work(1:lwork), &
                work(lwork+1:lwork+a_n_curr),work(lwork+a_n_curr+1:lwork+2* &
                a_n_curr),0,control,info,.FALSE.,use_multilevel,grid)

            END IF
          END IF

          IF (grid%level==1) CALL mg_grid_destroy(grid,info%flag)

          IF (nsvar+num_zero_row<a_n_new) THEN
            IF (PRESENT(seps)) THEN
              ! Expand and reorder seps
              DO i = 1, a_n_curr
                j = work(work_iperm+i)
                IF (j==1) THEN
                  ll = 1
                ELSE
                  ll = svwork(sv_svar+j-1) + 1
                END IF
                DO l = ll, svwork(sv_svar+j)
                  seps(svwork(sv_invp+l)) = work(work_seps+i)
                END DO
              END DO
            END IF

            ! Expand iperm to matrix before supervariables detected
            k = a_n_new
            DO i = a_n_curr, 1, -1
              j = work(work_iperm+i)
              IF (j==1) THEN
                ll = 1
              ELSE
                ll = svwork(sv_svar+j-1) + 1
              END IF
              DO l = ll, svwork(sv_svar+j)
                work(work_iperm+k) = iperm(svwork(sv_invp+l))
                k = k - 1
              END DO
            END DO

            iperm(1:a_n_new) = work(work_iperm+1:work_iperm+a_n_new)

          ELSE
            IF (PRESENT(seps)) THEN
              ! reorder seps!
              DO i = 1, a_n_curr
                j = iperm(i)
                work(j) = seps(i)
              END DO
              seps(1:a_n_curr) = work(1:a_n_curr)
            END IF
          END IF
        END IF
        ! Create perm from iperm
        DO i = 1, a_n
          j = iperm(i)
          perm(j) = i
        END DO

        ! Deallocate arrays
        IF (nsvar+num_zero_row<a_n_new) THEN
          DEALLOCATE (svwork,STAT=info%stat)
          IF (info%stat/=0) GO TO 20
        END IF
        DEALLOCATE (a_weight,STAT=info%stat)
        IF (info%stat/=0) GO TO 20
        DEALLOCATE (iperm,STAT=info%stat)
        IF (info%stat/=0) GO TO 20
        DEALLOCATE (work,STAT=info%stat)
        IF (info%stat/=0) GO TO 20

        info%flag = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info%flag,unit_diagnostics, &
            'nd_nested_both')
        END IF
        RETURN

10      info%flag = ND_ERR_MEMORY_ALLOC
        IF (printe) CALL nd_print_message(info%flag,unit_error, &
          'nd_nested_both')
        RETURN

20      info%flag = ND_ERR_MEMORY_DEALLOC
        IF (printe) CALL nd_print_message(info%flag,unit_error, &
          'nd_nested_both')
        RETURN

      END SUBROUTINE nd_nested_both

      ! ---------------------------------------------------
      ! nd_dense_rows
      ! ---------------------------------------------------
      ! Identifies and removes dense rows
      SUBROUTINE nd_dense_rows(a_n_in,a_ne_in,a_ptr,a_row,a_n_out,a_ne_out, &
          iperm,work,control,info)
        INTEGER, INTENT (IN) :: a_n_in ! dimension of subproblem before dense
        ! rows removed
        INTEGER, INTENT (IN) :: a_ne_in ! no. nonzeros of subproblem before
        ! dense rows removed
        INTEGER, INTENT (INOUT) :: a_ptr(a_n_in) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start. This is then
        ! used to hold positions for submatrices after dense row removed
        INTEGER, INTENT (INOUT) :: a_row(a_ne_in) ! On input a_row contains
        ! row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.This is then used to hold row indices for
        ! submatrices after partitioning
        INTEGER, INTENT (INOUT) :: iperm(a_n_in) ! On input, iperm(i) contains
        ! the row in the original matrix (when nd_nested was called) that
        ! row i in this sub problem maps to. On output, this is updated to
        ! reflect the computed permutation.
        INTEGER, INTENT (OUT) :: a_n_out ! dimension of subproblem after dense
        ! rows removed
        INTEGER, INTENT (OUT) :: a_ne_out ! no. nonzeros of subproblem after
        ! dense rows removed
        INTEGER, INTENT (OUT) :: work(4*a_n_in) ! Used during the algorithm to
        ! reduce need for allocations. The output is garbage.
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info

        ! ---------------------------------------------
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: deg, prev, next, dense ! pointers into work array
        INTEGER :: ndense ! number of dense rows found
        INTEGER :: max_deg ! maximum degree
        INTEGER :: degree, i, j, k, l, inext, ilast, l1, l2, m, m1
        LOGICAL :: printi, printd

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
        ! ---------------------------------------------------
        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Find and remove dense rows'
        END IF

        ! Set pointers into work array
        deg = 0
        prev = deg + a_n_in
        next = prev + a_n_in
        dense = next + a_n_in

        ! By the end of this loop work(dense+i) will be
        ! 0 if row is not dense
        ! <0 otherwise. The larger the number, the earlier the row was
        ! was determined to be dense.
        ndense = 0
        max_deg = 0
        work(deg+1:deg+a_n_in) = 0

        ! Calculate degree of each row before anything removed
        DO i = 1, a_n_in
          k = a_ptr(i)
          IF (i<a_n_in) THEN
            degree = a_ptr(i+1) - k
          ELSE
            degree = a_ne_in - a_ptr(a_n_in) + 1
          END IF
          work(dense+i) = degree
          IF (degree/=0) THEN
            max_deg = MAX(max_deg,degree)
            CALL dense_add_to_list(a_n_in,work(next+1:next+a_n_in),&
                work(prev+1:prev+a_n_in),work(deg+1:deg+a_n_in),i,degree)
          END IF
        END DO
        degree = max_deg
        a_n_out = a_n_in
        a_ne_out = a_ne_in

        DO WHILE (REAL(degree)-REAL(a_ne_out)/REAL(a_n_out)>=40*(REAL(a_n_out- &
            1)/REAL(a_n_out))*LOG(REAL(a_n_out)) .AND. degree>0)
          ! DO WHILE (real(degree) - real(a_ne_out)/real(a_n_out)>= &
          ! 300*(real(a_n_out-1)/real(a_n_out)) &
          ! .AND. degree>0)
          i = work(deg+degree)
          ndense = ndense + 1
          work(dense+i) = -ndense
          CALL dense_remove_from_list(a_n_in,work(next+1:next+a_n_in),&
                work(prev+1:prev+a_n_in),work(deg+1:deg+a_n_in),i,degree)
          ! update degrees of adjacent vertices
          IF (i<a_n_in) THEN
            l = a_ptr(i+1) - 1
          ELSE
            l = a_ne_in
          END IF
          DO k = a_ptr(i), l
            j = a_row(k)
            IF (work(dense+j)>0) THEN
              CALL dense_remove_from_list(a_n_in,work(next+1:next+a_n_in),&
                work(prev+1:prev+a_n_in),work(deg+1:deg+a_n_in),j,work(dense+j))
              work(dense+j) = work(dense+j) - 1
              IF (work(dense+j)>0) THEN
                CALL dense_add_to_list(a_n_in,work(next+1:next+a_n_in),&
                work(prev+1:prev+a_n_in),work(deg+1:deg+a_n_in),j,work(dense+j))
              END IF
            END IF
          END DO
          a_n_out = a_n_out - 1
          a_ne_out = a_ne_out - 2*degree
          IF (work(deg+degree)==0) THEN
            ! Find next largest degree
            degree = degree - 1
            DO
              IF (degree==0) EXIT
              IF (work(deg+degree)>0) EXIT
              degree = degree - 1
            END DO
          END IF
        END DO

        ! By the end of this loop work(dense+i) will be
        ! >=0 if row is not dense
        ! <0 otherwise. The larger the number, the earlier the row was
        ! was determined to be dense.
        ! !!!!

        IF (ndense>0) THEN
          IF (printi .OR. printd) THEN
            WRITE (unit_diagnostics,'(a)') ' '
            WRITE (unit_diagnostics,'(i10,a)') ndense, ' dense rows detected'
          END IF
          info%dense = ndense

          a_n_out = 0
          l = a_n_in + 1
          DO i = 1, a_n_in
            k = work(dense+i)
            IF (k>=0) THEN
              a_n_out = a_n_out + 1
              work(dense+i) = a_n_out
              work(next+a_n_out) = i
            ELSE
              work(next+l+k) = i
            END IF
          END DO

          k = 1
          j = 1

          DO i = 1, a_n_in
            l1 = a_ptr(i)
            IF (i<a_n_in) THEN
              l2 = a_ptr(i+1) - 1
            ELSE
              l2 = a_ne_in
            END IF
            IF (work(dense+i)>=0) THEN
              a_ptr(j) = k
              DO l = l1, l2
                m = a_row(l)
                m1 = work(dense+m)
                IF (m1>=0) THEN
                  a_row(k) = m1
                  k = k + 1
                END IF
              END DO
              j = j + 1
            END IF
          END DO
          a_ptr(j) = k
          IF (printd) THEN
            ! Print out a_ptr and a_row
            WRITE (unit_diagnostics,'(a11)') 'a_n_out = '
            WRITE (unit_diagnostics,'(i15)') a_n_out
            WRITE (unit_diagnostics,'(a11)') 'a_ne_out = '
            WRITE (unit_diagnostics,'(i15)') a_ne_out
            WRITE (unit_diagnostics,'(a8)') 'a_ptr = '
            WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n_out)
            WRITE (unit_diagnostics,'(a8)') 'a_row = '
            WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne_out)
          ELSE IF (printi) THEN
            ! Print out first few entries of a_ptr and a_row
            WRITE (unit_diagnostics,'(a11)') 'a_n_out = '
            WRITE (unit_diagnostics,'(i15)') a_n_out
            WRITE (unit_diagnostics,'(a11)') 'a_ne_out = '
            WRITE (unit_diagnostics,'(i15)') a_ne_out
            WRITE (unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n_out)) = '
            WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,MIN(5,a_n_out))
            WRITE (unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne_out)) = '
            WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,MIN(5,a_ne_out))
          END IF
        ELSE

          a_n_out = a_n_in
          a_ne_out = a_ne_in
          work(next+1:next+a_n_in) = (/ (i,i=1,a_n_in) /)
        END IF

        DO i = 1, a_n_in
          j = work(next+i)
          work(next+i) = iperm(j)
        END DO

        DO i = 1, a_n_in
          iperm(i) = work(next+i)
        END DO

        info%flag = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info%flag,unit_diagnostics, &
            'nd_dense_rows')
        END IF

      END SUBROUTINE nd_dense_rows

        SUBROUTINE dense_remove_from_list(n,next,prev,deg,irm,ig)
          INTEGER, INTENT(IN) :: n ! order matrix
          INTEGER, INTENT(INOUT) :: next(n),prev(n),deg(n)
          INTEGER, INTENT(IN) :: irm, ig
          INTEGER :: inext, ilast

          inext = next(irm)
          ilast = prev(irm)
          IF (ilast==0) THEN
            deg(ig) = inext
            IF (inext/=0) prev(inext) = 0
          ELSE
            next(ilast) = inext
            IF (inext/=0) prev(inext) = ilast
          END IF
        END SUBROUTINE dense_remove_from_list

        SUBROUTINE dense_add_to_list(n,next,prev,deg,irm,ig)
          INTEGER, INTENT(IN) :: n ! order matrix
          INTEGER, INTENT(INOUT) :: next(n),prev(n),deg(n)
          INTEGER, INTENT(IN) :: irm, ig
          INTEGER :: inext

          inext = deg(ig)
          deg(ig) = irm
          next(irm) = inext
          IF (inext/=0) prev(inext) = irm
          prev(irm) = 0
        END SUBROUTINE dense_add_to_list

      ! ---------------------------------------------------
      ! nd_nested_internal
      ! ---------------------------------------------------
      ! Does the recursive nested dissection

      RECURSIVE SUBROUTINE nd_nested_internal(a_n,a_ne,a_ptr,a_row, &
          a_weight,sumweight,iperm,work,work_comp_n,work_comp_nz,level, &
          control,info,use_amd,use_multilevel,grid,seps)

        INTEGER, INTENT (IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT (INOUT) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start. This is then
        ! used to hold positions for submatrices after partitioning
        INTEGER, INTENT (INOUT) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.This is then used to hold row indices for
        ! submatrices after partitioning
        INTEGER, INTENT (INOUT) :: a_weight(a_n) ! On input a_weight(i)
        ! contains
        ! weight of column i. This is then
        ! used to hold the weights for submatrices after partitioning
        INTEGER, INTENT (IN) :: sumweight ! sum entries in a_weight
        ! (unchanged)
        INTEGER, INTENT (INOUT) :: iperm(a_n) ! On input, iperm(i) contains
        ! the
        ! row in the original matrix (when nd_nested was called) that
        ! row i in this sub problem maps to. On output, this is updated to
        ! reflect the computed permutation.
        INTEGER, INTENT (OUT) :: work_comp_n(a_n)
        INTEGER, INTENT (OUT) :: work_comp_nz(a_n)
        INTEGER, INTENT (OUT) :: work(12*a_n+sumweight+a_ne) ! Used during the
        ! algorithm to reduce need for allocations. The output is garbage.
        INTEGER, INTENT (IN) :: level ! which level of nested dissection is
        ! this
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info
        LOGICAL, INTENT (IN) :: use_amd
        LOGICAL, INTENT (INOUT) :: use_multilevel
        INTEGER, INTENT (INOUT), OPTIONAL :: seps(a_n)
        ! seps(i) is -1 if vertex i of permuted submatrix is not in a
        ! separator; otherwise it is equal to l, where l is the nested
        ! dissection level at which it became part of the separator
        TYPE (nd_multigrid), INTENT (INOUT) :: grid

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: i, j, k, l, m, s
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: num_components ! Number of independent components found
        INTEGER :: compwork, lwork
        INTEGER :: sumweight_sub, maxdeg_max_component
        INTEGER :: offset_ptr, offset_row
        INTEGER :: a_n1, a_n2, a_ne1, a_ne2
        INTEGER :: amd_order_irn, amd_order_ip, amd_order_sep, amd_order_perm, amd_order_work, lirn
        LOGICAL :: printi, printd, use_amdi

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
        use_amdi = .FALSE.

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a,i6)') 'Nested dissection level ', level
        END IF

        ! Check whether matrix is diagonal and act accordingly
        IF (a_ne==0) THEN
          IF (printi .OR. printd) THEN
            WRITE (unit_diagnostics,'(a)') ' '
            WRITE (unit_diagnostics,'(a)') 'Submatrix is diagonal'
          END IF
          RETURN
        END IF

        IF (level==0) THEN
          maxdeg_max_component = a_ne + 1 - a_ptr(a_n)
          DO i = 1, a_n - 1
            IF (maxdeg_max_component<a_ptr(i+1)-a_ptr(i)) &
              maxdeg_max_component = a_ptr(i+1) - a_ptr(i)
          END DO

        END IF


        ! Check whether max number of levels has been reached or if matrix
        ! size is below
        ! nd_switch
        IF (level>=control%amd_switch2 .OR. a_n<=MAX(2,control%amd_switch1) &
          .OR. use_amd) &
          GO TO 10
        lwork = 12*a_n + sumweight + a_ne
        CALL nd_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,level, &
          a_n1,a_n2,a_ne1,a_ne2,iperm,work(1:lwork),control,info, &
          use_multilevel,grid)


        IF (a_n1==a_n) THEN
          GO TO 10
        END IF


        IF (a_n1/=0 .AND. a_n2/=0 .AND. a_n1+a_n2==a_n) THEN
          ! matrix is reducible
          IF (printi .OR. printd) THEN
            WRITE (unit_diagnostics,'(a)') ' '
            WRITE (unit_diagnostics,'(a)') 'Matrix is reducible'
          END IF
          compwork = 0
          ! work array needs to be total length 5*a_n+a_ne
          CALL nd_find_indep_comps(a_n,a_ne,a_ptr,a_row,a_weight,iperm, &
            num_components,work_comp_n(1:a_n),work_comp_nz(1:a_n), &
            work(compwork+1:compwork+3*a_n+a_ne),control,info)
          IF (num_components==1) THEN
            k = control%amd_switch2 ! Should never be reached - indep. comps.
            ! only found if it has been detected that they exist
          ELSE
            k = level
            IF (k==0) info%num_components = num_components
          END IF

          ! Apply the ND to each component - do not test for indep comps
          offset_ptr = a_n + 1
          offset_row = a_ne + 1
          i = num_components
          ! work from last component so that space at end of work_comp_n and
          ! work_comp_nz can be reused without worrying about overwriting
          ! important data
          j = a_n
          DO WHILE (i>=1)
            IF (level==0) use_multilevel = .TRUE.
            l = work_comp_n(i)
            IF (level==0 .AND. l .LE. control%amd_call) THEN
               use_amdi = .TRUE.
            ELSE
               use_amdi = .FALSE.
            END IF
            m = work_comp_nz(i)
            s = SUM(a_weight(offset_ptr-l:offset_ptr-1))
            IF (m>0) THEN
              ! Matrix not diagonal
              IF (PRESENT(seps)) THEN
                CALL nd_nested_internal(l,m,a_ptr(offset_ptr-l:offset_ptr-1 &
                  ),a_row(offset_row-m:offset_row-1), &
                  a_weight(offset_ptr-l:offset_ptr-1),s, &
                  iperm(offset_ptr-l:offset_ptr-1),work(compwork+1:compwork+12 &
                  *l+s+m),work_comp_n(j-l+1:j),work_comp_nz(j-l+1:j),k, &
                  control,info,use_amdi,use_multilevel,grid,seps(offset_ptr-l: &
                  offset_ptr-1))
              ELSE
                CALL nd_nested_internal(l,m,a_ptr(offset_ptr-l:offset_ptr-1 &
                  ),a_row(offset_row-m:offset_row-1), &
                  a_weight(offset_ptr-l:offset_ptr-1),s, &
                  iperm(offset_ptr-l:offset_ptr-1),work(compwork+1:compwork+12 &
                  *l+s+m),work_comp_n(j-l+1:j),work_comp_nz(j-l+1:j),k, &
                  control,info,use_amdi,use_multilevel,grid)

              END IF
            END IF
            offset_ptr = offset_ptr - l
            offset_row = offset_row - m
            j = j - l
            i = i - 1
          END DO
          RETURN
        ELSE
          IF (level==0 .AND. a_n>info%n_max_component) THEN
            info%n_max_component = a_n
            info%nz_max_component = a_ne
            info%maxdeg_max_component = maxdeg_max_component
          END IF
        END IF

        IF (PRESENT(seps)) THEN
          seps(a_n1+a_n2+1:a_n) = level
        END IF
        IF (a_n1>MAX(2,control%amd_switch1)) THEN
          sumweight_sub = SUM(a_weight(1:a_n1))
          IF (PRESENT(seps)) THEN
            CALL nd_nested_internal(a_n1,a_ne1,a_ptr(1:a_n1), &
              a_row(1:a_ne1),a_weight(1:a_n1),sumweight_sub,iperm(1:a_n1), &
              work(1:12*a_n1+sumweight_sub+a_ne1),work_comp_n(1:a_n1), &
              work_comp_nz(1:a_n1),level+1,control,info,use_amdi,&
              use_multilevel,grid, &
              seps(1:a_n1))
          ELSE
            CALL nd_nested_internal(a_n1,a_ne1,a_ptr(1:a_n1), &
              a_row(1:a_ne1),a_weight(1:a_n1),sumweight_sub,iperm(1:a_n1), &
              work(1:12*a_n1+sumweight_sub+a_ne1),work_comp_n(1:a_n1), &
              work_comp_nz(1:a_n1),level+1,control,info,use_amdi,&
              use_multilevel,grid)
          END IF

        END IF

        IF (a_n2>MAX(2,control%amd_switch1)) THEN
          IF (a_n1>MAX(2,control%amd_switch1)) THEN
            sumweight_sub = SUM(a_weight(a_n1+1:a_n1+a_n2))
            IF (PRESENT(seps)) THEN
              CALL nd_nested_internal(a_n2,a_ne2,a_ptr(a_n1+1:a_n1+a_n2), &
                a_row(a_ne1+1:a_ne1+a_ne2),a_weight(a_n1+1:a_n1+a_n2), &
                sumweight_sub,iperm(a_n1+1:a_n1+a_n2), &
                work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2), &
                work_comp_nz(1:a_n2),level+1,control,info,use_amdi,&
                use_multilevel,grid, &
                seps(a_n1+1:a_n1+a_n2))
            ELSE
              CALL nd_nested_internal(a_n2,a_ne2,a_ptr(a_n1+1:a_n1+a_n2), &
                a_row(a_ne1+1:a_ne1+a_ne2),a_weight(a_n1+1:a_n1+a_n2), &
                sumweight_sub,iperm(a_n1+1:a_n1+a_n2), &
                work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2), &
                work_comp_nz(1:a_n2),level+1,control,info,use_amdi,&
                use_multilevel,grid)
            END IF
          ELSE
            sumweight_sub = SUM(a_weight(a_n1+1:a_n1+a_n2))
            IF (PRESENT(seps)) THEN
              CALL nd_nested_internal(a_n2,a_ne2,a_ptr(1:a_n2), &
                a_row(1:a_ne2),a_weight(a_n1+1:a_n1+a_n2),sumweight_sub, &
                iperm(a_n1+1:a_n1+a_n2),work(1:12*a_n2+sumweight_sub+a_ne2), &
                work_comp_n(1:a_n2),work_comp_nz(1:a_n2),level+1,control,info, &
                use_amdi,use_multilevel,grid,seps(a_n1+1:a_n1+a_n2))
            ELSE
              CALL nd_nested_internal(a_n2,a_ne2,a_ptr(1:a_n2), &
                a_row(1:a_ne2),a_weight(a_n1+1:a_n1+a_n2),sumweight_sub, &
                iperm(a_n1+1:a_n1+a_n2),work(1:12*a_n2+sumweight_sub+a_ne2), &
                work_comp_n(1:a_n2),work_comp_nz(1:a_n2),level+1,control,info, &
                use_amdi,use_multilevel,grid)

            END IF

          END IF

        END IF
        GO TO 20
        RETURN

        ! No partition found or max number of levels have been reached
10      CONTINUE
        ! Apply AMD to matrix (note: a_n always greater than 1)
        lirn = a_ne + a_n
        amd_order_irn = 0
        amd_order_ip = amd_order_irn + lirn
        amd_order_sep = amd_order_ip + a_n
        amd_order_perm = amd_order_sep + a_n
        amd_order_work = amd_order_perm + a_n

        work(amd_order_irn+1:amd_order_irn+a_ne) = a_row(1:a_ne)
        work(amd_order_irn+a_ne+1:amd_order_irn+lirn) = 0
        work(amd_order_ip+1:amd_order_ip+a_n) = a_ptr(1:a_n)
        work(amd_order_sep+1:amd_order_sep+a_n) = ND_SEP_FLAG - 1
        CALL amd_order(a_n,a_ne,lirn,work(amd_order_irn+1:amd_order_irn+lirn), &
          work(amd_order_ip+1:amd_order_ip+a_n),work(amd_order_sep+1:amd_order_sep+a_n), &
          work(amd_order_perm+1:amd_order_perm+a_n),work(amd_order_work+1:amd_order_work+7*a_n))

        ! Extract perm from amd_order to apply to iperm
        DO i = 1, a_n
          j = work(amd_order_perm+i)
          work(amd_order_work+i) = iperm(j)
        END DO
        iperm(1:a_n) = work(amd_order_work+1:amd_order_work+a_n)

20      info%flag = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info%flag,unit_diagnostics, &
            'nd_nested_internal')
        END IF
        RETURN

      END SUBROUTINE nd_nested_internal

      ! ---------------------------------------------------
      ! nd_find_indep_comps
      ! ---------------------------------------------------
      ! Finds and forms independent components in a matrix
      SUBROUTINE nd_find_indep_comps(a_n,a_ne,a_ptr,a_row,a_weight,iperm, &
          comp_num,compsizes,compnzs,work,control,info)
        INTEGER, INTENT (IN) :: a_n ! size of matrix
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros in matrix
        INTEGER, INTENT (INOUT) :: a_ptr(a_n) ! On entry, column ptrs for
        ! input
        ! matrix.
        ! On exit, a_ptr(1:compsizes(1)) contains column ptrs for compontent
        ! 1;
        ! a_ptr(compsizes(1)+1:compsizes(1)+compsizes(2)) contains column ptrs
        ! for compontent 2; etc.
        INTEGER, INTENT (INOUT) :: a_row(a_ne) ! On entry, row indices for
        ! input
        ! matrix.
        ! On exit, a_row(1:compnzs(1)) contains row indices for compontent 1;
        ! a_ptr(compnzs(1)+1:compnzs(1)+compnzs(2)) contains row indices
        ! for compontent 2; etc.
        INTEGER, INTENT (INOUT) :: a_weight(a_n) ! On entry, a_weight(i)
        ! contains
        ! weight of column i for input matrix.
        ! On exit, a_weight(1:compsizes(1)) contains column weights for
        ! compontent 1;
        ! a_weight(compsizes(1)+1:compsizes(1)+compsizes(2)) contains column
        ! weights or compontent 2; etc.
        INTEGER, INTENT (INOUT) :: iperm(a_n) ! On input, iperm(i) contains
        ! the
        ! row in the original matrix (when nd_nested was called) that
        ! row i in this sub problem maps to. On output, this is updated to
        ! reflect the computed permutation.
        INTEGER, INTENT (OUT) :: comp_num ! number independent components
        ! found
        INTEGER, INTENT (OUT) :: compsizes(a_n) ! compsizes(i) will contain
        ! the
        ! size of compontent i
        INTEGER, INTENT (OUT) :: compnzs(a_n) ! compnzs(i) will contain the
        ! number of nonzeros in compontent i
        INTEGER, INTENT (OUT) :: work(3*a_n+a_ne) ! used as work arrays during
        ! computation
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: mask, ptr_temp, row_temp, front ! pointers into work array
        INTEGER :: root
        INTEGER :: front_sta, front_sto ! start and end of front list
        INTEGER :: num_assigned ! number of cols that have been assigned to
        ! a
        ! component
        INTEGER :: i, j, l, u, v, p
        INTEGER :: ptr_temp_pos, row_temp_pos, front_pos
        LOGICAL :: printi, printd

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Start independent component search'
        END IF

        mask = 0
        front = mask + a_n
        ptr_temp = front + a_n
        row_temp = ptr_temp + a_n ! total length of work array 3*a_n+a_ne

        ! Initialise arrays
        work(1:3*a_n) = 0
        compsizes(1:a_n) = 0
        compnzs(1:a_n) = 0

        ! Check from all possible roots
        num_assigned = 0
        comp_num = 0
        DO root = 1, a_n
          IF (work(mask+root)==0) THEN
            comp_num = comp_num + 1
            front_sta = num_assigned + 1
            front_sto = front_sta
            work(front+front_sta) = root
            compsizes(comp_num) = 1
            compnzs(comp_num) = 0
            work(mask+root) = compsizes(comp_num)
            num_assigned = num_assigned + 1

            DO WHILE (front_sto-front_sta>=0)
              DO i = front_sta, front_sto
                ! pick vertex from front
                v = work(front+i)
                ! update compnzs
                IF (v<a_n) THEN
                  l = a_ptr(v+1)
                ELSE
                  l = a_ne + 1
                END IF
                compnzs(comp_num) = compnzs(comp_num) + l - a_ptr(v)
                DO j = a_ptr(v), l - 1
                  ! pick a neighbour
                  u = a_row(j)
                  IF (work(mask+u)/=0) CYCLE
                  ! found unmasked vertex
                  compsizes(comp_num) = compsizes(comp_num) + 1
                  num_assigned = num_assigned + 1
                  ! mask this vertex
                  work(mask+u) = compsizes(comp_num)
                  ! add vertex to component
                  work(front+num_assigned) = u
                END DO
              END DO
              front_sta = front_sto + 1
              front_sto = num_assigned
            END DO
          END IF
        END DO

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(i5,a)') comp_num, ' components found'
        END IF

        IF (comp_num>1) THEN
          ! Reorder matrix into block diagonal form
          ! Indexing will be local to each block
          ptr_temp_pos = 1
          row_temp_pos = 1
          front_pos = 1
          DO l = 1, comp_num ! for each component
            work(ptr_temp+ptr_temp_pos) = 1
            ptr_temp_pos = ptr_temp_pos + 1
            DO u = 1, compsizes(l)
              v = work(front+front_pos) ! for each column in the component
              front_pos = front_pos + 1
              IF (v==a_n) THEN
                p = a_ne
              ELSE
                p = a_ptr(v+1) - 1
              END IF
              DO i = a_ptr(v), p ! for each nonzero in the column
                work(row_temp+row_temp_pos) = work(mask+a_row(i))
                row_temp_pos = row_temp_pos + 1
              END DO
              IF (u<compsizes(l)) THEN
                work(ptr_temp+ptr_temp_pos) = work(ptr_temp+ptr_temp_pos-1) + &
                  p + 1 - a_ptr(v)
                ptr_temp_pos = ptr_temp_pos + 1
              END IF
            END DO
          END DO

          a_ptr(1:a_n) = work(ptr_temp+1:ptr_temp+a_n)
          a_row(1:a_ne) = work(row_temp+1:row_temp+a_ne)

          ! Reorder iperm and a_weight
          DO front_pos = 1, a_n
            j = work(front+front_pos)
            work(front+front_pos) = iperm(j)
            work(mask+front_pos) = a_weight(j)
          END DO
          iperm(1:a_n) = work(front+1:front+a_n)
          a_weight(1:a_n) = work(mask+1:mask+a_n)
        END IF

        info%flag = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info%flag,unit_diagnostics, &
            'nd_find_indep_comps')
        END IF
        RETURN

      END SUBROUTINE nd_find_indep_comps

      ! ---------------------------------------------------
      ! nd_partition matrix
      ! ---------------------------------------------------
      ! Partition the matrix and if one (or more) of the generated submatrices
      ! is
      ! small enough, apply halo amd

      SUBROUTINE nd_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
          level,a_n1,a_n2,a_ne1,a_ne2,iperm,work,control,info,use_multilevel, &
          grid)

        INTEGER, INTENT (IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT (INOUT) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start. This is then
        ! used to hold positions for submatrices after partitioning
        INTEGER, INTENT (INOUT) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.This is then used to hold row indices for
        ! submatrices after partitioning
        INTEGER, INTENT (INOUT) :: a_weight(a_n) ! On input a_weight(i)
        ! contains
        ! the weight of column i. This is then used to hold the weights for
        ! the submatrices after partitioning.
        INTEGER, INTENT (IN) :: sumweight ! Sum entries in a_weight.
        ! Unchanged.
        INTEGER, INTENT (IN) :: level ! Current nested dissection level
        INTEGER, INTENT (OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT (OUT) :: a_ne1, a_ne2 ! no. nonzeros in two
        ! submatrices
        INTEGER, INTENT (INOUT) :: iperm(a_n) ! On input, iperm(i) contains
        ! the
        ! row in the original matrix (when nd_nested was called) that
        ! row i in this sub problem maps to. On output, this is updated to
        ! reflect the computed permutation.
        LOGICAL, INTENT (INOUT) :: use_multilevel
        INTEGER, INTENT (OUT) :: work(12*a_n+sumweight+a_ne)
        TYPE (nd_options), INTENT (IN) :: control
        TYPE (nd_inform), INTENT (INOUT) :: info
        TYPE (nd_multigrid), INTENT (INOUT) :: grid

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        LOGICAL :: printi, printd
        INTEGER :: partition_ptr ! pointer into work array
        INTEGER :: work_ptr ! pointer into work array
        INTEGER :: part_ptr ! pointer into work array
        INTEGER :: a_ptr_sub_ptr ! pointer into work array
        INTEGER :: a_row_sub_ptr ! pointer into work array
        INTEGER :: partition_method
        INTEGER :: i, j, k
        INTEGER :: a_weight_1, a_weight_2, a_weight_sep
        INTEGER :: ref_method, ref_control
        INTEGER :: a_n1_new, a_n2_new, a_weight_1_new, a_weight_2_new, &
          a_weight_sep_new
        REAL (wp) :: ratio, tau_best, tau, band, depth
        LOGICAL :: imbal


        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Start finding a partition'
        END IF

        ! If matrix is full, then don't partition
        IF (a_ne==a_n*(a_n-1)) THEN
          a_n1 = a_n
          a_n2 = 0
          a_ne1 = a_ne
          a_ne2 = 0
          GO TO 10
        END IF

        ! Find the partition
        IF (control%coarse_partition_method<=1) THEN
          partition_method = 1
        ELSE
          partition_method = 2
        END IF


        partition_ptr = 0 ! length a_n
        work_ptr = partition_ptr + a_n ! max length 9*a_n+sumweight+a_ne/2
        ! DO p = 1,a_n
        ! DO q = p+1,a_n
        ! IF (control%partition_method .GE. 2) use_multilevel = .TRUE.
        IF (partition_method==1) THEN
          ! Ashcraft method
          CALL nd_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,level, &
            a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
            work(partition_ptr+1:partition_ptr+a_n), &
            work(work_ptr+1:work_ptr+9*a_n+sumweight),control,info%flag,band, &
            depth,use_multilevel,grid)
        ELSE
          CALL nd_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
            level,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
            work(partition_ptr+1:partition_ptr+a_n), &
            work(work_ptr+1:work_ptr+9*a_n+sumweight),control,info%flag, &
            band,depth,use_multilevel,grid)

        END IF

        IF (a_n1+a_n2==a_n) THEN
          RETURN
        END IF

        IF (level==0) THEN
          IF (a_n>info%n_max_component .OR. (a_n==info%n_max_component .AND. &
              band>info%band)) THEN
            info%band = band
            info%depth = depth
          END IF

        END IF

        IF (a_n1/=0 .AND. a_n2/=0 .AND. a_n>=3) THEN
          IF ( .NOT. use_multilevel) THEN

            IF (control%refinement>6) THEN
              ref_control = 3
            ELSE
              IF (control%refinement<1) THEN
                ref_control = 1
              ELSE
                ref_control = control%refinement
              END IF
            END IF

            SELECT CASE (ref_control)
            CASE (1)
              ref_method = 1

            CASE (2)
              ref_method = 2

            CASE (3)
              IF (REAL(MAX(a_weight_1,a_weight_2))/REAL(MIN(a_weight_1, &
                  a_weight_2)+a_weight_sep)>=MAX(REAL(1.0, &
                  wp),control%balance)) THEN
                ref_method = 2
              ELSE
                ref_method = 1
              END IF

            CASE (4)
              ref_method = 0

            CASE (5)
              ref_method = 2

            CASE (6)
              IF (REAL(MAX(a_weight_1,a_weight_2))/REAL(MIN(a_weight_1, &
                  a_weight_2)+a_weight_sep)>=MAX(REAL(1.0, &
                  wp),control%balance)) THEN
                ref_method = 2
              ELSE
                ref_method = 0
              END IF

            END SELECT
            IF (printd) THEN
              WRITE (unit_diagnostics,'(a)') 'Partition before refinement'
              WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', a_n1, &
                ',  a_n2=', a_n2, ',  a_n_sep=', a_n - a_n1 - a_n2
            END IF

            SELECT CASE (ref_method)

            CASE (0)
              CALL nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1, &
                a_n2,a_weight_1,a_weight_2,a_weight_sep, &
                work(partition_ptr+1:partition_ptr+a_n), &
                work(work_ptr+1:work_ptr+8),control)

            CASE (1)
              IF (MIN(a_weight_1,a_weight_2)+a_weight_sep< &
                  MAX(a_weight_1,a_weight_2)) THEN
                CALL nd_refine_block_trim(a_n,a_ne,a_ptr,a_row, &
                  a_weight,sumweight,a_n1,a_n2,a_weight_1,a_weight_2, &
                  a_weight_sep,work(partition_ptr+1:partition_ptr+a_n), &
                  work(work_ptr+1:work_ptr+5*a_n),control)
              ELSE
                CALL nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
                  sumweight,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
                  work(partition_ptr+1:partition_ptr+a_n), &
                  work(work_ptr+1:work_ptr+3*a_n),control)
              END IF


            CASE (2)
              CALL nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
                a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
                work(partition_ptr+1:partition_ptr+a_n), &
                work(work_ptr+1:work_ptr+3*a_n),control)


            END SELECT


            IF (printd) THEN
              WRITE (unit_diagnostics,'(a)') 'Partition after refinement'
              WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', a_n1, &
                ',  a_n2=', a_n2, ',  a_n_sep=', a_n - a_n1 - a_n2
            END IF

            IF (control%max_improve_cycles>0) THEN
              ratio = MAX(REAL(1.0,wp),control%balance)
              IF (ratio>REAL(sumweight-2)) THEN
                imbal = .FALSE.
              ELSE
                imbal = .TRUE.
              END IF
              CALL cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight, &
                ratio,imbal,control%cost_function,tau_best)
              a_n1_new = a_n1
              a_n2_new = a_n2
              a_weight_1_new = a_weight_1
              a_weight_2_new = a_weight_2
              a_weight_sep_new = a_weight_sep
            END IF

            part_ptr = work_ptr + 5*a_n
            work(part_ptr+1:part_ptr+a_n) = work(partition_ptr+1:partition_ptr &
              +a_n)

            k = control%max_improve_cycles
            DO i = 1, k

              CALL expand_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1_new, &
                a_n2_new,a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
                work(part_ptr+1:part_ptr+a_n),work(work_ptr+1:work_ptr+5*a_n))


              IF (printd) THEN
                WRITE (unit_diagnostics,'(a)') &
                  'Partition sizes after expansion'
                WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', &
                  a_n1_new, ',  a_n2=', a_n2_new, ',  a_n_sep=', &
                  a_n - a_n1_new - a_n2_new
              END IF

              ! call
              ! check_partition1(a_n,a_ne,a_ptr,a_row,a_n1_new,a_n2_new,work(p
              ! art_ptr+1:part_ptr+a_n))

              SELECT CASE (ref_control)

              CASE (3)
                IF (REAL(MAX(a_weight_1_new,a_weight_2_new))/REAL(MIN( &
                    a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>MAX(REAL( &
                    1.0,wp),control%balance)) THEN
                  ref_method = 2
                ELSE
                  ref_method = 1
                END IF

              CASE (6)
                IF (REAL(MAX(a_weight_1_new,a_weight_2_new))/REAL(MIN( &
                    a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>MAX(REAL( &
                    1.0,wp),control%balance)) THEN
                  ref_method = 2
                ELSE
                  ref_method = 0
                END IF
              END SELECT


              SELECT CASE (ref_method)

              CASE (0)
                CALL nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight, &
                  a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
                  a_weight_sep_new,work(part_ptr+1:part_ptr+a_n), &
                  work(work_ptr+1:work_ptr+8),control)

              CASE (1)
                IF (MIN(a_weight_1,a_weight_2)+a_weight_sep< &
                    MAX(a_weight_1,a_weight_2)) THEN
                  CALL nd_refine_block_trim(a_n,a_ne,a_ptr,a_row, &
                    a_weight,sumweight,a_n1_new,a_n2_new,a_weight_1_new, &
                    a_weight_2_new,a_weight_sep_new, &
                    work(part_ptr+1:part_ptr+a_n),work(work_ptr+1:work_ptr+5* &
                    a_n),control)
                ELSE
                  CALL nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
                    sumweight,a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
                    a_weight_sep_new,work(part_ptr+1:part_ptr+a_n), &
                    work(work_ptr+1:work_ptr+3*a_n),control)
                END IF


              CASE (2)
                CALL nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
                  a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
                  a_weight_sep_new,work(part_ptr+1:part_ptr+a_n), &
                  work(work_ptr+1:work_ptr+3*a_n),control)

              END SELECT


              IF (printd) THEN
                WRITE (unit_diagnostics,'(a)') &
                  'Partition sizes after refinement'
                WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', &
                  a_n1_new, ',  a_n2=', a_n2_new, ',  a_n_sep=', &
                  a_n - a_n1_new - a_n2_new
              END IF


              CALL cost_function(a_weight_1_new,a_weight_2_new, &
                a_weight_sep_new,sumweight,ratio,imbal,control%cost_function,tau)
              IF (tau<tau_best) THEN
                tau_best = tau
                work(partition_ptr+1:partition_ptr+a_n) &
                  = work(part_ptr+1:part_ptr+a_n)
                a_n1 = a_n1_new
                a_n2 = a_n2_new
                a_weight_1 = a_weight_1_new
                a_weight_2 = a_weight_2_new
                a_weight_sep = a_weight_sep_new
              ELSE
                EXIT
              END IF
            END DO

            CALL nd_refine_fm(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
              a_n2,a_weight_1,a_weight_2,a_weight_sep, &
              work(partition_ptr+1:partition_ptr+a_n), &
              work(work_ptr+1:work_ptr+8*a_n+sumweight),control)


          END IF
        ELSE
          GO TO 10
        END IF

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') 'Partition found:'
          WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1=', a_n1, &
            ', a_n2=', a_n2, ', a_nsep=', a_n - a_n1 - a_n2
        END IF

        IF ((a_n1<=MAX(2,control%amd_switch1) .AND. a_n2<=MAX(2, &
            control%amd_switch1))) THEN
          ! apply halo amd to submatrics
          CALL amd_order_both(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
            work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight, &
            work(work_ptr+1:work_ptr+12*a_n+a_ne))
          a_ne1 = 0
          a_ne2 = 0

        ELSE IF (a_n1<=MAX(2,control%amd_switch1)) THEN
          ! apply halo amd to [A1, B1'; B1, I] using two levels
          ! return A2 and apply ND to it
          CALL amd_order_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
            work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight,a_ne2, &
            work(work_ptr+1:work_ptr+12*a_n+a_ne))
          a_ne1 = 0


        ELSE IF (a_n2<=MAX(2,control%amd_switch1)) THEN
          ! apply halo amd to [A2, B2'; B2, I] using two levels
          ! return A1 and apply ND to it
          CALL amd_order_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
            work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight,a_ne1, &
            work(work_ptr+1:work_ptr+12*a_n+a_ne))
          a_ne2 = 0

        ELSE
          ! return A1 and A2 and apply ND to them
          a_ptr_sub_ptr = work_ptr + a_n ! length a_n
          a_row_sub_ptr = a_ptr_sub_ptr + a_n ! length a_ne

          CALL extract_both_matrices(a_n,a_ne,a_ptr,a_row,a_n1,a_n2, &
            work(partition_ptr+1:partition_ptr+a_n1+a_n2),a_ne1,a_ne2, &
            work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2),a_ne, &
            work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne), &
            work(work_ptr+1:work_ptr+a_n))

          ! Copy extracted matrices
          a_ptr(1:a_n1+a_n2) = work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2)
          a_row(1:a_ne1+a_ne2) = work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne1+ &
            a_ne2)

          ! Update iperm and a_weight
          DO i = 1, a_n
            j = work(partition_ptr+i)
            work(a_ptr_sub_ptr+i) = iperm(j)
            work(work_ptr+i) = a_weight(j)
          END DO
          DO i = 1, a_n
            iperm(i) = work(a_ptr_sub_ptr+i)
            a_weight(i) = work(work_ptr+i)
          END DO

        END IF
        GO TO 20

10      IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') 'No partition found'
        END IF

20      info%flag = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info%flag,unit_diagnostics, &
            'nd_partition')
        END IF
        RETURN

      END SUBROUTINE nd_partition

      ! ---------------------------------------------------
      ! nd_ashcraft
      ! ---------------------------------------------------
      ! Partition the matrix using the Ashcraft method
      RECURSIVE SUBROUTINE nd_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight, &
          sumweight,level,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
          partition,work,control,info,band,depth,use_multilevel,grid)

        INTEGER, INTENT (IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i.
        INTEGER, INTENT (IN) :: sumweight ! sum of entries in a_weight
        INTEGER, INTENT (IN) :: level ! current level of nested dissection
        INTEGER, INTENT (OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT (OUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weight
        ! ed
        ! size of partitions and separator
        INTEGER, INTENT (OUT) :: partition(a_n) ! First a_n1 entries will
        ! contain
        ! list of (local) indices in partition 1; next a_n2 entries will
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end
        INTEGER, INTENT (OUT) :: work(9*a_n+sumweight) ! used as work array
        TYPE (nd_options), INTENT (IN) :: control
        INTEGER, INTENT (INOUT) :: info
        REAL (wp), INTENT (OUT) :: band ! If level = 0, then on
        ! output band = 100*L/a_n, where L is the size of the
        ! largest levelset
        REAL (wp), INTENT (OUT) :: depth ! If level = 0, then
        ! on
        ! output depth = num_levels_nend
        LOGICAL, INTENT (INOUT) :: use_multilevel ! are we allowed to use a
        ! multilevel
        ! partitioning strategy
        TYPE (nd_multigrid), INTENT (INOUT) :: grid

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: nstrt, nend
        INTEGER :: i, j, dptr, p1sz, p2sz, sepsz, lwork, k
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: mask_p, level_p, level_ptr_p, level2_p, level2_ptr_p, &
          work_p
        INTEGER :: num_levels_nend ! no. levels in structure rooted at nend
        INTEGER :: num_levels_nstrt ! no. levels in structure rooted at nstrt
        INTEGER :: num_entries ! no. entries in level set structure
        INTEGER :: best_sep_start
        INTEGER :: distance
        INTEGER :: distance_ptr
        INTEGER :: lwidth, mindeg, degree, max_search
        INTEGER :: ww
        INTEGER :: stop_coarsening2 ! max no. multigrid levels
        REAL (wp) :: bestval
        REAL (wp) :: val
        REAL (wp) :: ratio
        LOGICAL :: printi, printd
        LOGICAL :: imbal, use_multilevel_copy

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Use two-sided level set method'
        END IF
        ratio = MAX(REAL(1.0,wp),control%balance)
        IF (ratio>REAL(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF
        p2sz = 0
        sepsz = 0
        use_multilevel_copy = use_multilevel


        band = -1
        depth = -1

        IF (control%partition_method==1 .AND. use_multilevel) GO TO 10

        IF (control%partition_method>1 .AND. level>0 .AND. use_multilevel) &
          GO TO 10

        ! Find pseudoperipheral nodes nstart and nend, and the level structure
        ! rooted at nend
        mask_p = 0 ! size a_n
        level_ptr_p = mask_p + a_n ! size a_n
        level_p = level_ptr_p + a_n ! size a_n
        level2_ptr_p = level_p + a_n ! size a_n
        level2_p = level2_ptr_p + a_n ! size a_n
        work_p = level2_p + a_n ! size 2*a_n

        nend = -1

        ! Choose nstrt
        ! node with minimum degree
        mindeg = sumweight + 1
        DO i = 1, a_n
          IF (i<a_n) THEN
            k = a_ptr(i+1) - 1
          ELSE
            k = a_ne
          END IF
          degree = 0
          DO j = a_ptr(i), k
            degree = degree + a_weight(a_row(j))
          END DO
          IF (degree<mindeg) THEN
            mindeg = degree
            nstrt = i
          END IF
        END DO

        max_search = 5

        CALL nd_find_pseudo(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
          work(level_ptr_p+1:level_ptr_p+a_n),work(level_p+1:level_p+a_n), &
          nstrt,nend,max_search,work(work_p+1:work_p+2*a_n),num_levels_nend, &
          num_entries,lwidth)

        IF (num_entries<a_n) THEN
          ! matrix is separable
          a_n1 = num_entries
          a_n2 = a_n - a_n1
          a_weight_sep = 0
          a_weight_1 = 0
          work(work_p+1:work_p+a_n) = 0
          DO i = 1, num_entries
            j = work(level_p+i)
            partition(i) = j
            work(work_p+j) = 1
            a_weight_1 = a_weight_1 + a_weight(j)
          END DO
          a_weight_2 = sumweight - a_weight_1
          j = num_entries + 1
          DO i = 1, a_n
            IF (work(work_p+i)==0) THEN
              partition(j) = i
              j = j + 1
            END IF

          END DO
          IF (level==0) THEN
            band = -REAL(lwidth,wp)
          END IF

          RETURN
        END IF

        ! ********************************************************************
        ! ***
        IF (level==0) THEN
          band = 100.0*REAL(lwidth,wp)/REAL(a_n,wp)
          ! band = max(band,real(lwidth,wp))
          ! write(*,*) sqrt(real(lwidth))/real(num_levels_nend), lwidth, &
          ! real(sumweight)/real(num_levels_nend), &
          ! real(num_levels_nend)*(real(lwidth)/real(sumweight))
          depth = 100.0*REAL(num_levels_nend,wp)/ &
            REAL(a_n,wp)
        END IF
        IF (control%stop_coarsening2<=0 .OR. control%partition_method<1) THEN
          use_multilevel = .FALSE.
        END IF
        IF (control%partition_method>=2 .AND. use_multilevel) THEN
          ! IF (real(lwidth,wp) .LE. 2.0* &
          ! real(sumweight,wp)/real(num_levels_nend,wp) )
          ! THEN
          IF (100.0*REAL(lwidth,wp)/REAL(sumweight,wp)<=3.0 .OR. &
             a_n.LT. control%ml_call) &
              THEN
            use_multilevel = .FALSE.
          ELSE
            use_multilevel = .TRUE.
          END IF
        END IF

        IF (use_multilevel) GO TO 10


        ! Find level structure rooted at nstrt
        work(mask_p+1:mask_p+a_n) = 1
        CALL nd_level_struct(nstrt,a_n,a_ne,a_ptr,a_row, &
          work(mask_p+1:mask_p+a_n),work(level2_ptr_p+1:level2_ptr_p+a_n), &
          work(level2_p+1:level2_p+a_n),num_levels_nstrt,lwidth,num_entries)

        ! Calculate difference in distances from nstart and nend, and set up
        ! lists D_i of nodes with same distance
        distance_ptr = work_p
        distance = mask_p
        CALL nd_distance(a_n,num_levels_nend,work(level_ptr_p+1:level_ptr_p &
          +a_n),work(level_p+1:level_p+a_n),num_levels_nstrt, &
          work(level2_ptr_p+1:level2_ptr_p+a_n),work(level2_p+1:level2_p+a_n), &
          work(distance_ptr+1:distance_ptr+2*a_n-1), &
          work(distance+1:distance+a_n))

        ! Do not need the information in work(level_ptr_p+1:level2_ptr_p)
        ! Calculate total weight in each distance level
        dptr = level_ptr_p + a_n
        work(dptr+1-a_n:dptr+a_n-1) = 0
        DO i = 1 - num_levels_nstrt, num_levels_nend - 2
          DO j = work(distance_ptr+a_n+i), work(distance_ptr+a_n+i+1) - 1
            work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
          END DO
        END DO
        i = num_levels_nend - 1
        DO j = work(distance_ptr+a_n+i), a_n
          work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
        END DO

        ! Find first possible separator
        ww = 2
        p1sz = 0
        DO i = 1 - num_levels_nstrt, num_levels_nend - ww - 2
          p1sz = p1sz + work(dptr+i)
          sepsz = work(dptr+i+1)
          DO j = 1, ww - 1
            sepsz = sepsz + work(dptr+i+1+j)
          END DO
          p2sz = sumweight - p1sz - sepsz
          IF (p1sz>0 .AND. sepsz>0) THEN
            EXIT
          END IF
        END DO

        IF (i+ww>=num_levels_nend-1 .OR. p2sz==0) THEN
          ! Not possible to find separator
          ! This can only be invoked for a fully connected graph. The
          ! partition
          ! subroutine checks for this case so it should never be called
          a_n1 = a_n
          a_n2 = 0
          partition(1:a_n) = (/ (i,i=1,a_n) /)
          RETURN
        END IF

        ! Entries in level i will form first possible partition
        best_sep_start = i + 1
        a_n1 = work(distance_ptr+a_n+i+1) - 1
        a_n2 = a_n - work(distance_ptr+a_n+i+1+ww) + 1
        a_weight_1 = p1sz
        a_weight_2 = p2sz
        a_weight_sep = sepsz
        CALL cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
            control%cost_function,bestval)

        ! Search for best separator using tau

        DO j = i + 1, num_levels_nend - ww - 2
          p1sz = p1sz + work(dptr+j)
          sepsz = work(dptr+j+1)
          DO k = 1, ww - 1
            sepsz = sepsz + work(dptr+j+1+k)
          END DO
          p2sz = sumweight - sepsz - p1sz
          IF (p2sz==0) EXIT
          CALL cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
             control%cost_function,val)
          IF (val<bestval) THEN
            bestval = val
            best_sep_start = j + 1
            a_n1 = work(distance_ptr+a_n+j+1) - 1
            a_n2 = a_n - work(distance_ptr+a_n+j+1+ww) + 1
            a_weight_1 = p1sz
            a_weight_2 = p2sz
            a_weight_sep = sepsz
          END IF
        END DO

        IF (imbal .AND. use_multilevel_copy .AND. control%partition_method>=2) &
            THEN
          IF (REAL(MAX(a_weight_1,a_weight_2))/REAL(MIN(a_weight_1, &
              a_weight_2))>ratio) THEN
            use_multilevel = .TRUE.
            GO TO 10

          END IF
        END IF

        ! Rearrange partition
        ! Entries in partition 1
        j = 1
        DO i = 1, work(distance_ptr+a_n+best_sep_start) - 1
          partition(j) = work(distance+i)
          j = j + 1
        END DO

        ! Entries in partition 2
        DO i = work(distance_ptr+a_n+best_sep_start+ww), a_n
          partition(j) = work(distance+i)
          j = j + 1
        END DO

        ! Entries in separator
        DO i = work(distance_ptr+a_n+best_sep_start), &
            work(distance_ptr+a_n+best_sep_start+ww) - 1
          partition(j) = work(distance+i)
          j = j + 1
        END DO
        GO TO 20

10      CONTINUE

        IF (use_multilevel) THEN
          stop_coarsening2 = control%stop_coarsening2
          lwork = 9*a_n + sumweight
          CALL multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
            partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,control, &
            info,lwork,work(1:lwork),stop_coarsening2,grid)
          RETURN
        END IF


20      info = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info,unit_diagnostics,'nd_ashcraft')
        END IF
        RETURN

      END SUBROUTINE nd_ashcraft


      ! ---------------------------------------------------
      ! nd_level_set
      ! ---------------------------------------------------
      ! Partition the matrix using the level set method
      RECURSIVE SUBROUTINE nd_level_set(a_n,a_ne,a_ptr,a_row,a_weight, &
          sumweight,level,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
          partition,work,control,info,band,depth,use_multilevel,grid)

        INTEGER, INTENT (IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i.
        INTEGER, INTENT (IN) :: sumweight ! sum of entries in a_weight
        INTEGER, INTENT (IN) :: level ! current nested dissection level
        INTEGER, INTENT (OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT (OUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weight
        ! ed
        ! size of partitions and separator
        INTEGER, INTENT (OUT) :: partition(a_n) ! First a_n1 entries will
        ! contain
        ! list of (local) indices in partition 1; next a_n2 entries will
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end
        INTEGER, INTENT (OUT) :: work(9*a_n+sumweight)
        TYPE (nd_options), INTENT (IN) :: control
        INTEGER, INTENT (INOUT) :: info
        REAL (wp), INTENT (OUT) :: band ! If level = 0, then on
        ! output band = 100*L/a_n, where L is the size of the
        ! largest levelset
        REAL (wp), INTENT (OUT) :: depth ! If level = 0, then
        ! on
        ! output band = num_levels_nend
        LOGICAL, INTENT (INOUT) :: use_multilevel ! are we allowed to use a
        ! multilevel
        ! partitioning strategy
        TYPE (nd_multigrid), INTENT (INOUT) :: grid

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: nstrt, nend
        LOGICAL :: printi, printd
        INTEGER :: level_p, level_ptr_p, work_p
        INTEGER :: num_levels_nend ! no. levels in structure rooted at nend
        INTEGER :: num_entries ! no. entries in level structure rooted at
        ! nend
        INTEGER :: best_sep_start
        INTEGER :: i, j, k, p1sz, p2sz, sepsz, lwidth
        INTEGER :: stop_coarsening2, lwork
        INTEGER :: mindeg, degree, max_search
        REAL (wp) :: bestval
        REAL (wp) :: val
        REAL (wp) :: ratio
        LOGICAL :: imbal

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Use one-sided level set method'
        END IF

        ratio = MAX(REAL(1.0,wp),control%balance)
        IF (ratio>REAL(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF

        band = -1
        depth = -1

        IF (control%partition_method==1 .AND. use_multilevel) THEN
          use_multilevel = .TRUE.
          GO TO 10
        END IF

        IF (control%partition_method>1 .AND. level>0 .AND. use_multilevel) &
          GO TO 10

        ! Find pseudoperipheral nodes nstart and nend, and the level structure
        ! rooted at nend
        level_ptr_p = 0 ! size a_n
        level_p = level_ptr_p + a_n ! size a_n
        work_p = level_p + a_n ! size 2*a_n
        mindeg = sumweight + 1
        DO i = 1, a_n
          IF (i<a_n) THEN
            k = a_ptr(i+1) - 1
          ELSE
            k = a_ne
          END IF
          degree = 0
          DO j = a_ptr(i), k
            degree = degree + a_weight(a_row(j))
          END DO
          IF (degree<mindeg) THEN
            mindeg = degree
            nstrt = i
          END IF
        END DO
        max_search = 5

        CALL nd_find_pseudo(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
          work(level_ptr_p+1:level_ptr_p+a_n),work(level_p+1:level_p+a_n), &
          nstrt,nend,max_search,work(work_p+1:work_p+2*a_n),num_levels_nend, &
          num_entries,lwidth)

        IF (num_entries<a_n) THEN
          ! matrix is separable
          a_n1 = num_entries
          a_n2 = a_n - a_n1
          work(work_p+1:work_p+a_n) = 0
          a_weight_1 = 0
          DO i = 1, num_entries
            j = work(level_p+i)
            partition(i) = j
            work(work_p+j) = 1
            a_weight_1 = a_weight_1 + a_weight(j)
          END DO
          j = num_entries + 1
          a_weight_2 = 0
          DO i = 1, a_n
            IF (work(work_p+i)==0) THEN
              partition(j) = i
              a_weight_2 = a_weight_2 + a_weight(j)
              j = j + 1
            END IF
          END DO
          a_weight_sep = 0
          IF (level==0) THEN
            band = -REAL(lwidth,wp)
          END IF
          RETURN
        END IF

        IF (level==0) THEN
          band = 100.0*REAL(lwidth,wp)/REAL(sumweight,wp)
          depth = 100.0*REAL(num_levels_nend,wp)/ &
            REAL(sumweight,wp)
          ! band = max(band,real(lwidth,wp))
        END IF

        IF ((control%partition_method<=0) .OR. (use_multilevel .AND. control% &
            stop_coarsening2<=0)) THEN
          use_multilevel = .FALSE.
        END IF
        IF (control%partition_method>=2 .AND. use_multilevel) THEN
          IF (100.0*REAL(lwidth,wp)/REAL(sumweight,wp)<=3.0 .OR. &
             a_n.LT. control%ml_call) &
              THEN
            use_multilevel = .FALSE.
          ELSE
            use_multilevel = .TRUE.
          END IF
        END IF

10      CONTINUE

        IF (use_multilevel) THEN
          stop_coarsening2 = control%stop_coarsening2
          lwork = 9*a_n + sumweight
          
          CALL multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
            partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,control, &
            info,lwork,work(1:lwork),stop_coarsening2,grid)
          RETURN
        END IF

        IF (num_levels_nend<=2) THEN
          ! Not possible to find separator
          ! This can only be invoked for a full connected graph. The partition
          ! subroutine checks for this case so it should never be called
          a_n1 = a_n
          a_n2 = 0
          partition(1:a_n) = (/ (i,i=1,a_n) /)
          RETURN
        END IF

        ! Calculate total weight in each level set
        work(work_p+1:work_p+num_levels_nend) = 0
        DO i = 1, num_levels_nend - 1
          DO j = work(level_ptr_p+i), work(level_ptr_p+i+1) - 1
            work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
          END DO
        END DO
        i = num_levels_nend
        DO j = work(level_ptr_p+i), a_n
          work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
        END DO


        ! First possible separator contains all of the nodes in level 2
        p1sz = work(work_p+1)
        sepsz = work(work_p+2)
        p2sz = SUM(work(work_p+3:work_p+num_levels_nend))
        a_weight_1 = p1sz
        a_weight_2 = p2sz
        a_weight_sep = sepsz
        best_sep_start = 2
        a_n1 = work(level_ptr_p+2) - 1
        a_n2 = a_n - work(level_ptr_p+3) + 1
        CALL cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
         control%cost_function,bestval)

        ! Search for best separator using tau
        DO j = 2, num_levels_nend - 4
          p1sz = p1sz + work(work_p+j)
          sepsz = work(work_p+j+1)
          p2sz = p2sz - work(work_p+j+1)
          CALL cost_function(p1sz,p2sz,sepsz,sumweight,ratio,imbal,&
             control%cost_function,val)
          IF (val<bestval) THEN
            bestval = val
            best_sep_start = j + 1
            a_n1 = work(level_ptr_p+j+1) - 1
            a_n2 = a_n - work(level_ptr_p+j+2) + 1
            a_weight_1 = p1sz
            a_weight_2 = p2sz
            a_weight_sep = sepsz
          END IF
        END DO


        ! Rearrange partition
        ! Entries in partition 1
        j = 1
        DO i = 1, work(level_ptr_p+best_sep_start) - 1
          partition(j) = work(level_p+i)
          j = j + 1
        END DO

        ! Entries in partition 2
        DO i = work(level_ptr_p+best_sep_start+1), a_n
          partition(j) = work(level_p+i)
          j = j + 1
        END DO

        ! Entries in separator
        DO i = work(level_ptr_p+best_sep_start), work(level_ptr_p+ &
            best_sep_start+1) - 1
          partition(j) = work(level_p+i)
          j = j + 1
        END DO

        info = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info,unit_diagnostics,'nd_level_set')
        END IF
        RETURN

      END SUBROUTINE nd_level_set

      ! ---------------------------------------------------
      ! nd_find_pseudo
      ! ---------------------------------------------------
      ! Find pseudoperipheral pairs of nodes for irreducible graph
      SUBROUTINE nd_find_pseudo(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
          level_ptr,level,nstrt,nend,max_search,work,num_levels,num_entries, &
          lwidth)
        INTEGER, INTENT (IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! weight of vertex i
        INTEGER, INTENT (IN) :: sumweight ! sum of entries in a_weight
        INTEGER, INTENT (OUT) :: level_ptr(a_n) ! On output level_ptr(i)
        ! contains
        ! position in level that entries for level i start.
        INTEGER, INTENT (OUT) :: level(a_n) ! On output level contains lists
        ! of
        ! rows according to the level set that they are in
        INTEGER, INTENT (INOUT) :: nstrt ! Starting pseudoperipheral node
        INTEGER, INTENT (OUT) :: nend ! End pseudoperipheral node
        INTEGER, INTENT (IN) :: max_search
        INTEGER, INTENT (OUT) :: work(2*a_n)
        INTEGER, INTENT (OUT) :: num_levels
        INTEGER, INTENT (OUT) :: num_entries ! number of entries in level
        ! structure
        INTEGER, INTENT (OUT) :: lwidth
        ! Based on MC60HD

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: i, j, k, l, ll
        INTEGER :: mindeg, maxdep, main, lwidth1
        INTEGER :: mask, list
        INTEGER :: nstop ! ending pseudoperipheral node
        INTEGER :: node ! Index of graph node
        INTEGER :: lsize ! size of final levelset
        INTEGER :: nlsize ! no. nodes of differing degrees in final
        ! level set
        INTEGER :: minwid ! Minimum levelset width
        INTEGER :: nlvl ! number of levels in level structure
        INTEGER :: minwid1

        j = 0
        mask = 0
        list = mask + a_n
        ! Initialise work(mask+1:mask+a_n) and work(list+1:list+a_n)
        work(mask+1:mask+a_n) = 1
        work(list+1:list+a_n) = 0
        level_ptr(:) = 0


        ! First guess for starting node is input nstrt

        ! Generate level structure for node nstrt
        CALL nd_level_struct(nstrt,a_n,a_ne,a_ptr,a_row, &
          work(mask+1:mask+a_n),level_ptr,level,maxdep,lwidth,num_entries)
        IF (num_entries<a_n) THEN
          ! matrix is separable
          num_levels = maxdep
          RETURN

        END IF

        nstop = level(a_n)
        DO main = 1, MIN(a_n,10)
          ! Store nodes in the final level set and their (weighted) degrees
          lsize = 0
          j = level_ptr(maxdep)
          DO i = j, a_n
            node = level(i)
            lsize = lsize + 1
            work(list+lsize) = node
            level_ptr(node) = 0
            IF (node==a_n) THEN
              k = a_ne
            ELSE
              k = a_ptr(node+1) - 1
            END IF
            DO l = a_ptr(node), k
              level_ptr(node) = level_ptr(node) + a_weight(a_row(l))
            END DO
          END DO

          ! Choose at most max_search nodes
          DO nlsize = 1, max_search
            ! Look for candiate with least degree
            mindeg = sumweight + 1
            ! mindeg = -1
            DO i = nlsize, lsize
              IF (level_ptr(work(list+i))<mindeg) THEN
                j = i
                mindeg = level_ptr(work(list+i))
              END IF
            END DO
            ! Jump out of loop if no candidates left
            IF (mindeg==sumweight+1) GO TO 10
            ! IF (mindeg .EQ. -1) GO TO 55
            ! Swap chose candidate to next position
            node = work(list+j)
            work(list+j) = work(list+nlsize)
            work(list+nlsize) = node
            ! Rule out the neighbours of the chosen node
            IF (node==a_n) THEN
              k = a_ne
            ELSE
              k = a_ptr(node+1) - 1
            END IF
            DO i = a_ptr(node), k
              level_ptr(a_row(i)) = sumweight + 1
            END DO
          END DO
10        nlsize = nlsize - 1

          ! Loop over nodes in list
          minwid = HUGE(a_n)
          minwid1 = HUGE(a_n)

          DO i = 1, nlsize
            node = work(list+i)

            ! Form rooted level structures for node
            CALL nd_level_struct(node,a_n,a_ne,a_ptr,a_row, &
              work(mask+1:mask+a_n),level_ptr,level,nlvl,lwidth,num_entries)
            ! IF (lwidth .LE. minwid) THEN
            lwidth1 = 0
            DO k = 1, nlvl - 1
              DO l = level_ptr(k), level_ptr(k+1) - 1
                ll = level(l)
                lwidth1 = lwidth1 + k*a_weight(ll)
              END DO
            END DO
            DO l = level_ptr(nlvl), a_n
              ll = level(l)
              lwidth1 = lwidth1 + nlvl*a_weight(ll)
            END DO

            IF (nlvl>maxdep) THEN
              ! Level structure of greater depth. Begin a new iteration.
              nstrt = node
              maxdep = nlvl

              GO TO 20
            ELSE
              IF (lwidth1<minwid1) THEN

                nstop = node
                minwid = lwidth
                minwid1 = lwidth1
              END IF
            END IF
            ! END IF
          END DO
          GO TO 30
20        CONTINUE
        END DO
30      IF (nstop/=node) THEN
          CALL nd_level_struct(node,a_n,a_ne,a_ptr,a_row, &
            work(mask+1:mask+a_n),level_ptr,level,nlvl,lwidth,num_entries)
        END IF
        num_levels = maxdep
        nend = nstop

        RETURN
      END SUBROUTINE nd_find_pseudo


      ! ---------------------------------------------------
      ! nd_level_struct
      ! ---------------------------------------------------
      ! Given a root, calculate the level structure of a given graph from that
      ! root

      SUBROUTINE nd_level_struct(root,a_n,a_ne,a_ptr,a_row,mask,level_ptr, &
          level,num_levels,lwidth,num_entries)

        INTEGER, INTENT (IN) :: root ! Root node for level structure
        INTEGER, INTENT (IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (INOUT) :: mask(a_n) ! Always restored to input value
        ! at
        ! end of the call. mask(node) > 0 for all visible nodes
        INTEGER, INTENT (OUT) :: level_ptr(a_n) ! On output level_ptr(i)
        ! contains
        ! position in level that entries for level i start.
        INTEGER, INTENT (OUT) :: level(a_n) ! On output level contains lists
        ! of
        ! rows according to the level set that they are in
        INTEGER, INTENT (OUT) :: num_levels ! On output num_levels contains
        ! the
        ! number of levels
        INTEGER, INTENT (OUT) :: lwidth ! On output, contains the width of the
        ! structure
        INTEGER, INTENT (OUT) :: num_entries ! On output, contains number of
        ! entries in the tree structure containing root

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: lvlend ! End of previous level in level
        INTEGER :: lnbr ! Next position in level
        INTEGER :: lbegin ! Beginning of present level in level
        INTEGER :: lw ! Level width
        INTEGER :: nbr ! neighbour node
        INTEGER :: node ! current node
        INTEGER :: nvars
        INTEGER :: i, j, k

        ! Based on mc60ld
        ! Establish level 1
        level_ptr(:) = 0
        mask(root) = -mask(root)
        level(1) = root
        lvlend = 0
        nvars = 0
        lnbr = 1
        lwidth = 1
        DO num_levels = 1, a_n
          ! Generate next level by finding all unmasked neighbours of nodes in
          ! present level
          lbegin = lvlend + 1
          lvlend = lnbr
          level_ptr(num_levels) = lbegin
          lw = 0
          DO i = lbegin, lvlend
            node = level(i)
            IF (node==a_n) THEN
              k = a_ne
            ELSE
              k = a_ptr(node+1) - 1
            END IF
            DO j = a_ptr(node), k
              nbr = a_row(j)
              IF (mask(nbr)>0) THEN
                lnbr = lnbr + 1
                level(lnbr) = nbr
                mask(nbr) = -mask(nbr)
                ! lw = lw + a_weight(nbr)
                lw = lw + 1
              END IF
            END DO
          END DO
          lwidth = MAX(lwidth,lw)
          nvars = nvars + lw
          ! If no neighbours found, we are done
          IF (lnbr==lvlend) EXIT
          ! Abort construnction if level structure too wide
          ! IF (lwidth .GT. maxwid) EXIT
        END DO
        ! Reset mask
        DO i = 1, lnbr
          mask(level(i)) = ABS(mask(level(i)))
        END DO
        num_entries = lnbr

      END SUBROUTINE nd_level_struct


      ! ---------------------------------------------------
      ! nd_distance
      ! ---------------------------------------------------
      ! Given two level structures, calculate the difference in each nodes
      ! distance
      ! from the start and end node
      SUBROUTINE nd_distance(a_n,num_levels_nend,level_ptr_nend,level_nend, &
          num_levels_nstrt,level_ptr_nstrt,level_nstrt,distance_ptr,distance)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: num_levels_nend ! number of levels with root
        ! nend
        INTEGER, INTENT (IN) :: level_ptr_nend(a_n) ! level_ptr(i) contains
        ! position in level that entries for level i start (root = nend)
        INTEGER, INTENT (IN) :: level_nend(a_n) ! Contains lists of rows
        ! according to the level set that they are in (root = nend)
        INTEGER, INTENT (IN) :: num_levels_nstrt ! no. of levels with root
        ! nstrt
        INTEGER, INTENT (INOUT) :: level_ptr_nstrt(a_n) ! level_ptr(i)
        ! contains
        ! position in level that entries for level i start (root = nstrt)
        ! Reused during subroutine
        INTEGER, INTENT (IN) :: level_nstrt(a_n) ! Contains lists of rows
        ! according to the level set that they are in (root = nstrt)
        INTEGER, INTENT (OUT) :: distance_ptr(2*a_n-1) ! distance(i) contains
        ! position in distance where entries with distance i-a_n
        INTEGER, INTENT (OUT) :: distance(a_n) ! Contains lists of rows
        ! ordered
        ! according to their distance

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: j, k
        INTEGER :: dptr ! used as offset into distance_ptr
        INTEGER :: lev ! stores current level

        dptr = a_n
        distance_ptr(dptr+1-a_n:dptr+a_n-1) = 0
        distance(1:a_n) = 0

        ! set distance(i) to hold the level that i belongs to in the structure
        ! rooted at nend (= distance from end node + 1)
        DO lev = 1, num_levels_nend - 1
          DO j = level_ptr_nend(lev), level_ptr_nend(lev+1) - 1
            distance(level_nend(j)) = lev
          END DO
        END DO
        DO j = level_ptr_nend(num_levels_nend), a_n
          distance(level_nend(j)) = num_levels_nend
        END DO

        ! now consider level structure rooted at start node
        DO lev = 1, num_levels_nstrt - 1
          DO j = level_ptr_nstrt(lev), level_ptr_nstrt(lev+1) - 1
            distance(level_nstrt(j)) = distance(level_nstrt(j)) - lev
            k = distance(level_nstrt(j))
            distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
          END DO
        END DO
        DO j = level_ptr_nstrt(num_levels_nstrt), a_n
          distance(level_nstrt(j)) = distance(level_nstrt(j)) - lev
          k = distance(level_nstrt(j))
          distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
        END DO

        ! Copy distance into level_ptr_nstrt to save memory
        level_ptr_nstrt(1:a_n) = distance(1:a_n)

        ! Set distance_ptr(i) to point one place to the right of where entries
        ! with distance i will be
        DO j = 1 - a_n, a_n - 1
          IF (distance_ptr(dptr+j)>0) THEN
            EXIT
          ELSE
            distance_ptr(dptr+j) = 1
          END IF
        END DO
        distance_ptr(dptr+j) = distance_ptr(dptr+j) + 1
        DO k = j + 1, a_n - 1
          distance_ptr(dptr+k) = distance_ptr(dptr+k) + distance_ptr(dptr+k-1)
        END DO

        ! Set up lists of rows with same value of distance
        DO j = 1, a_n
          k = level_ptr_nstrt(j)
          distance_ptr(dptr+k) = distance_ptr(dptr+k) - 1
          distance(distance_ptr(dptr+k)) = j
        END DO


      END SUBROUTINE nd_distance

      ! ---------------------------------------------------
      ! nd_convert_partition_flags
      ! ---------------------------------------------------
      ! Given a partition array, convert the partition into a flag array
      SUBROUTINE nd_convert_partition_flags(a_n,a_n1,a_n2,partition,flag_1, &
          flag_2,flag_sep,flags)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (IN) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (IN) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end.
        INTEGER, INTENT (IN) :: flag_1 ! flag for rows in partition 1
        INTEGER, INTENT (IN) :: flag_2 ! flag for rows in partition 2
        INTEGER, INTENT (IN) :: flag_sep ! flag for rows in separator
        INTEGER, INTENT (OUT) :: flags(a_n) ! flags(i) contains flag for row i
        ! and indicates which partition it is in

        INTEGER :: j, k

        DO j = 1, a_n1
          k = partition(j)
          flags(k) = flag_1
        END DO
        DO j = a_n1 + 1, a_n1 + a_n2
          k = partition(j)
          flags(k) = flag_2
        END DO
        DO j = a_n1 + a_n2 + 1, a_n
          k = partition(j)
          flags(k) = flag_sep
        END DO

      END SUBROUTINE nd_convert_partition_flags


      ! ---------------------------------------------------
      ! nd_flags_partition
      ! ---------------------------------------------------
      ! Given a partition array, convert the partition into a flag array
      SUBROUTINE nd_convert_flags_partition(a_n,a_n1,a_n2,flags,flag_1, &
          flag_2,partition)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (IN) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (IN) :: flags(a_n) ! flags(i) contains flag for row i
        ! and indicates which partition it is in
        INTEGER, INTENT (IN) :: flag_1 ! flag for rows in partition 1
        INTEGER, INTENT (IN) :: flag_2 ! flag for rows in partition 2
        INTEGER, INTENT (OUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end.

        INTEGER :: i, j, k, l

        i = 1
        j = a_n1 + 1
        k = a_n1 + a_n2 + 1
        DO l = 1, a_n
          IF (flags(l)==flag_1) THEN
            partition(i) = l
            i = i + 1
          ELSE
            IF (flags(l)==flag_2) THEN
              partition(j) = l
              j = j + 1
            ELSE
              partition(k) = l
              k = k + 1
            END IF
          END IF
        END DO

      END SUBROUTINE nd_convert_flags_partition

      ! ---------------------------------------------------
      ! nd_move_partition
      ! ---------------------------------------------------
      ! Given a flag array, move the separator by forming an edge separator
      ! between the input separator and the larger of P1 and P2
      SUBROUTINE nd_move_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
          a_weight_1,a_weight_2,a_weight_sep,flag_1,flag_2,flag_sep,flags)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (IN) :: flag_1 ! flag for rows in partition 1
        INTEGER, INTENT (IN) :: flag_2 ! flag for rows in partition 2
        INTEGER, INTENT (IN) :: flag_sep ! flag for rows in separator
        INTEGER, INTENT (INOUT) :: flags(a_n) ! flags(i) contains flag for row
        ! i
        ! and indicates which partition it is in. This is updated

        INTEGER :: part_no, a_n1_temp, a_n2_temp
        INTEGER :: a_weight_1_temp, a_weight_2_temp, a_weight_sep_temp
        INTEGER :: j, k, l
        LOGICAL :: stay_sep

        ! Decide whether new (large) separator is formed from
        ! intersection of the separator with partition 1 or the intersection
        ! of
        ! the separator with partition 2
        IF (a_weight_1>a_weight_2) THEN
          ! Apply to partition 1
          part_no = flag_1
        ELSE
          ! Apply to partition 2
          part_no = flag_2
        END IF

        a_n1_temp = a_n1
        a_n2_temp = a_n2
        a_weight_1_temp = a_weight_1
        a_weight_2_temp = a_weight_2
        a_weight_sep_temp = a_weight_sep

        DO j = 1, a_n
          IF (flags(j)/=flag_sep) CYCLE
          ! j is in initial separator
          IF (j==a_n) THEN
            k = a_ne
          ELSE
            k = a_ptr(j+1) - 1
          END IF
          stay_sep = .FALSE.
          DO l = a_ptr(j), k
            IF (flags(a_row(l))==part_no) THEN
              stay_sep = .TRUE.
              IF (part_no==flag_1) THEN
                IF (a_n1_temp>1) THEN
                  a_n1_temp = a_n1_temp - 1
                  flags(a_row(l)) = -1
                  a_weight_sep_temp = a_weight_sep_temp + a_weight(a_row(l))
                  a_weight_1_temp = a_weight_1_temp - a_weight(a_row(l))
                END IF
              ELSE
                IF (a_n2_temp>1) THEN
                  a_n2_temp = a_n2_temp - 1
                  flags(a_row(l)) = -1
                  a_weight_sep_temp = a_weight_sep_temp + a_weight(a_row(l))
                  a_weight_2_temp = a_weight_2_temp - a_weight(a_row(l))
                END IF
              END IF
            ELSE IF (flags(a_row(l))==-1) THEN
              stay_sep = .TRUE.
            END IF
          END DO
          IF ( .NOT. stay_sep) THEN
            IF (part_no==flag_1) THEN
              flags(j) = flag_2
            ELSE
              flags(j) = flag_1
            END IF
            a_weight_sep_temp = a_weight_sep_temp - a_weight(j)
            IF (part_no==flag_1) THEN
              a_n2_temp = a_n2_temp + 1
              a_weight_2_temp = a_weight_2_temp + a_weight(j)
            ELSE
              a_n1_temp = a_n1_temp + 1
              a_weight_1_temp = a_weight_1_temp + a_weight(j)
            END IF
          END IF
        END DO

        DO j = 1, a_n
          IF (flags(j)==-1) THEN
            flags(j) = flag_sep
          END IF
        END DO

        a_n1 = a_n1_temp
        a_n2 = a_n2_temp
        a_weight_1 = a_weight_1_temp
        a_weight_2 = a_weight_2_temp
        a_weight_sep = a_weight_sep_temp

      END SUBROUTINE nd_move_partition


      ! ---------------------------------------------------
      ! nd_refine_edge
      ! ---------------------------------------------------
      ! Given a partition, refine the partition to improve the (weighted) value 
      ! of the cost function. An edge separator is formed between the input 
      ! separator and the larger partition, and this is then minimal using 
      ! trimming or max flow
      SUBROUTINE nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
          a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end. This is updated to the new
        ! partition
        INTEGER, INTENT (OUT) :: work(3*a_n) ! Work array
        TYPE (nd_options), INTENT (IN) :: control

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: fm_flags ! pointer into work array for start of
        ! flags from FM

        fm_flags = 0 ! length a_n

        ! Initialise work(fm_flags+1:fm_flags+a_n)
        CALL nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
          ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
          work(fm_flags+1:fm_flags+a_n))

        ! Create new separator by forming edge separator between input
        ! separator and largest of P1 and P2

        CALL nd_move_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
          a_weight_1,a_weight_2,a_weight_sep,ND_PART1_FLAG, &
          ND_PART2_FLAG,ND_SEP_FLAG,work(fm_flags+1:fm_flags+a_n))

        ! Update partition
        CALL nd_convert_flags_partition(a_n,a_n1,a_n2, &
          work(fm_flags+1:fm_flags+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
          partition(1:a_n))

        IF (control%refinement>3) THEN
          CALL nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
            a_weight_1,a_weight_2,a_weight_sep,partition,work(1:8),control)
        ELSE
          CALL nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
            a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work(1:3*a_n), &
            control)
        END IF



      END SUBROUTINE nd_refine_edge

      ! ---------------------------------------------------
      ! nd_refine_fm
      ! ---------------------------------------------------
      ! Given a partition, refine the partition using FM refinement. Wrapper
      ! for code
      SUBROUTINE nd_refine_fm(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
          a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end. This is updated to the new
        ! partition
        INTEGER, INTENT (OUT) :: work(8*a_n+sumweight) ! Work array
        TYPE (nd_options), INTENT (IN) :: control

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: fm_flags ! pointer into work array for start of
        ! flags from FM
        INTEGER :: fm_ipart ! pointer into work array for start of
        ! ipart from FM
        INTEGER :: fm_next ! pointer into work array for start of next
        ! from FM
        INTEGER :: fm_last ! pointer into work array for start of last
        ! from FM
        INTEGER :: fm_gain1 ! pointer into work array for start of
        ! gain1 from FM
        INTEGER :: fm_gain2 ! pointer into work array for start of
        ! gain2 from FM
        INTEGER :: fm_done ! pointer into work array for start of done
        ! from FM
        INTEGER :: fm_head ! pointer into work array for start of head
        ! from FM
        INTEGER :: fm_distance ! pointer into work array for start of head
        ! from FM
        INTEGER :: icut, mult ! Used within FM refinement
        INTEGER :: band
        REAL (wp) :: ratio
        LOGICAL :: imbal ! Should we check for imbalance?

        IF (control%refinement_band .LT. 1) RETURN

        ratio = MAX(REAL(1.0,wp),control%balance)
        IF (ratio>REAL(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF
        fm_flags = 0 ! length a_n

        ! Initialise work(fm_flags+1:fm_flags+a_n)
        CALL nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
          ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
          work(fm_flags+1:fm_flags+a_n))

        fm_ipart = fm_flags + a_n ! length a_n
        fm_next = fm_ipart + a_n ! length a_n
        fm_last = fm_next + a_n ! length a_n
        fm_gain1 = fm_last + a_n ! length a_n
        fm_gain2 = fm_gain1 + a_n ! length a_n
        fm_done = fm_gain2 + a_n ! length a_n
        fm_head = fm_done + a_n ! length icut+mult+1
        icut = MIN(sumweight-1,3*(sumweight/a_n))
        icut = MIN(icut,5*MAXVAL(a_weight))
        ! icut = sumweight/2
        ! mult = min(sumweight/20,10*sumweight/a_n) - 1
        mult = sumweight - icut - 1
        mult = MIN(mult,icut)
        ! mult = sumweight/2-1
        fm_distance = fm_head + icut + mult + 1

        band = MIN(control%refinement_band,a_n)

        CALL nd_fm_refinement(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,icut, &
          mult,control%cost_function,a_n1,a_n2,a_weight_1,a_weight_2,&
          a_weight_sep,band,ratio, &
          work(fm_flags+1:fm_flags+a_n),work(fm_ipart+1:fm_ipart+a_n), &
          work(fm_next+1:fm_next+a_n),work(fm_last+1:fm_last+a_n), &
          work(fm_gain1+1:fm_gain1+a_n),work(fm_gain2+1:fm_gain2+a_n), &
          work(fm_done+1:fm_done+a_n),work(fm_head+1:fm_head+icut+mult+1), &
          work(fm_distance+1:fm_distance+a_n))

        ! Update partition
        CALL nd_convert_flags_partition(a_n,a_n1,a_n2, &
          work(fm_flags+1:fm_flags+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
          partition(1:a_n))

      END SUBROUTINE nd_refine_fm




      ! The subroutine nd_fm_refinement uses a version of the Fiduccia-
      ! Mattheyses refinement algorithm on a tripartite partitioing of the
      ! nodes of a
      ! graph where a node in the first partition is not connected to any node
      ! in the second partition and any path between nodes in partition 1 and
      ! partition 2 must go through a node in the cutset (partition 3).
      ! The intention of the algorithm is to reduce f(P1,P2,P3), where
      ! f(P1,P2,P3) =   |P3|/(|P1||P2|) if min(|P1|,|P2|)/max(|P1|,|P2|) >=
      ! ratio
      ! =  sumweight - 2 + max(|P1|,|P2|)/min(|P1|,|P2|), otherwise
      ! This is a banded version so only nodes a distance of at most band from
      ! the
      ! input separator can be moved into the new separator

      SUBROUTINE nd_fm_refinement(n,a_ne,ptr,col,weight,sumweight,icut, &
          mult,costf,a_n1,a_n2,wnv1,wnv2,wns,band,ratio,flags,ipart,next,last,gain1, &
          gain2,done,head,dist)

        ! Matrix is held in matrix using compressed column scheme
        INTEGER, INTENT (IN) :: n ! size of matrix
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros in matrix
        INTEGER, INTENT (IN) :: ptr(n) ! row pointers
        INTEGER, INTENT (IN) :: col(a_ne) ! column indices
        ! TYPE (nd_matrix), INTENT (INOUT) :: matrix
        ! The array weight is used to hold a weight on the vertices indicating
        ! how many vertices from the finer graphs have been combined into the
        ! current coarse graph vertex.
        INTEGER, INTENT (IN) :: weight(n)
        INTEGER, INTENT (IN) :: sumweight
        INTEGER, INTENT (IN) :: icut ! Used to limit search
        INTEGER, INTENT (IN) :: mult ! Used to bound search
        INTEGER, INTENT (IN) :: costf ! Determines which cost function used
        INTEGER, INTENT (INOUT) :: a_n1 ! No. vertices partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! No. vertices partition 2
        INTEGER, INTENT (INOUT) :: wnv1 ! Weighted sum of vertices partition 1
        INTEGER, INTENT (INOUT) :: wnv2 ! Weighted sum of vertices partition 2
        INTEGER, INTENT (INOUT) :: wns ! Weighted sum of vertices separator
        INTEGER, INTENT (IN) :: band ! width of band around initial separator
        ! that the separator can lie in
        REAL (wp), INTENT (IN) :: ratio ! ratio to determine
        ! whether
        ! partition is balanced

        ! flags holds a list of nodes stating which partition node i is in.
        ! The whole point of this routine is to return a revised partition
        ! with better properties.  Normally less nodes in the cutset while
        ! maintaining
        ! a balance between the number of nodes in the two components.
        ! flags(i) == ND_PART1_FLAG : i is in partition 1
        ! flags(i) == ND_PART2_FLAG : i is in partition 2
        ! flags(i) == ND_SEP_FLAG   : i is in separator/cutset
        INTEGER, INTENT (INOUT) :: flags(n)
        ! info holds parameters giving information about the performance of
        ! the
        ! subroutine
        INTEGER, INTENT (OUT) :: ipart(n), next(n), last(n)
        INTEGER, INTENT (OUT) :: gain1(n), gain2(n), done(n), head(-mult:icut)
        INTEGER, INTENT (OUT) :: dist(n)

        ! Number nodes in each partition
        INTEGER :: nv1, ns, nv2, inv1, inv2, ins
        ! Weighted nodes in each partition
        INTEGER :: winv1, winv2, wins
        INTEGER :: i, j, ii, jj, eye, k, l
        INTEGER :: inn, outer
        INTEGER :: move, ming, gain, old_gain, inext, ilast, idummy
        INTEGER :: first, tail
        REAL (wp) :: eval, evalc, evalo, eval1, eval2
        LOGICAL :: imbal


        IF (ratio>REAL(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF

        ! Set-up distance to hold (min) distance of node from separator. Only
        ! find
        ! distance for nodes within band distance from separator
        first = 0
        tail = 0
        DO i = 1, n
          IF (flags(i)==ND_SEP_FLAG) THEN
            dist(i) = 0
            IF (first==0) THEN
              first = i
              tail = i
            ELSE
              next(tail) = i
              tail = i
            END IF
          ELSE
            dist(i) = -2
          END IF
        END DO

        DO WHILE (first/=0)
          j = dist(first)
          IF (j==band-1) THEN
            IF (first==n) THEN
              l = a_ne
            ELSE
              l = ptr(first+1) - 1
            END IF
            DO i = ptr(first), l
              IF (dist(col(i))==-2) THEN
                k = col(i)
                dist(k) = j + 1
              END IF
            END DO

          ELSE
            IF (first==n) THEN
              l = a_ne
            ELSE
              l = ptr(first+1) - 1
            END IF
            DO i = ptr(first), l
              IF (dist(col(i))==-2) THEN
                k = col(i)
                dist(k) = j + 1
                next(tail) = k
                tail = k
              END IF
            END DO

          END IF
          IF (first==tail) THEN
            first = 0
          ELSE
            k = next(first)
            first = k
          END IF
        END DO
        next(1:n) = 0

        ! nv1,nv2,ns are the number of nodes in partitions 1, 2 and the cutset
        ! in
        ! the current partition
        ! inv1,inv2,ins,ipart are the equivalent quantities within the inner
        ! loop
        ! The same identifiers prefixed by w refer to weighted counts
        ! inner and outer are the two main loop indices

        ! Initialize nv1,nv2,ns
        nv1 = a_n1
        nv2 = a_n2
        ns = n - (a_n1+a_n2)
        ii = 1
        jj = ii + nv1

        ! Initialize ipart
        ipart(1:n) = flags(1:n)

        ! Initialize array done that flags that a node has been considered in
        ! an inner loop pass
        done = 0

        ! Compute evaluation function for current partitioning

        CALL cost_function(wnv1+1,wnv2+1,wns,sumweight,ratio,imbal,costf,evalc)

        ! icut is set to limit search in inner loop .. may later be a
        ! parameter
        ! we allow gains of up to max(weight)*5

        head(-mult:icut) = 0

        ! Set up doubly linked list linking nodes with same gain and headers
        ! to starts (up to cut off value icut)

        ! Compute gains for nodes in cutset
        ming = sumweight
        DO i = 1, n

          IF (flags(i)==ND_SEP_FLAG) THEN
            ! Node i is in cutset
            ! gain1(i) is change to cutset size if node i is moved to
            ! partition 1.
            ! gain2(i) is change to cutset size if node i is moved to
            ! partition 2.
            ! Run through all neighbours of node i to see what loss/gain is if
            ! node
            ! i is moved

            CALL compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,flags)
            gain = MAX(-mult,MIN(gain1(i),gain2(i)))

            IF (gain<ming) ming = gain
            IF (gain>icut) CYCLE
            ! New node is put at head of list
            CALL add_to_list(n,mult,icut,next,last,head,i,gain)
          END IF
        END DO

        ! Initilialization finished.  Now perform F-M algorithm in two loops.
        ! In each inner loop we choose the best obtained and if the evaluation
        ! function for this is better than previous best we perform another
        ! inner
        ! loop; otherwise we terminate.
        evalo = evalc
        DO outer = 1, n
          ! Set partition that will be altered in inner loop
          inv1 = nv1
          inv2 = nv2
          ins = ns
          winv1 = wnv1
          winv2 = wnv2
          wins = wns
          ipart(1:n) = flags(1:n)
INNER:    DO inn = 1, n

            ! Choose best eligible move
            DO idummy = 1, n

              DO gain = ming, icut
                IF (head(gain)/=0) EXIT
              END DO
              IF (gain>icut) EXIT INNER

              ! Now cycle through nodes of least gain
              ! Currently inefficient because of re-searching linked list
              inext = head(gain)
              k = 0
10            i = inext
              IF (i==0) CYCLE
              ! Use node if it has not been considered already
              IF (done(i)<outer) GO TO 20
              inext = next(i)
              ! !! Extra statements to trap infinite loop
              k = k + 1
              IF (k>ins) THEN
                ! WRITE (*,*) 'Bug in code because of infinite loop'
                ! !! You may wish to change this to a stop
                EXIT
              END IF
              GO TO 10
            END DO
            EXIT INNER
            ! Node i has been selected as the best eligible node
            ! Set flag so only considered once in this pass
20          done(i) = outer
            ! As i will not be chosen again in this pass, remove from list
            CALL remove_from_list(n,mult,icut,next,last,head,i,gain)
            ! Move the node to the appropriate partition and reset partition
            ! information
            ! We will try both weighted and unweighted

            IF (wnv1==0 .AND. wnv2>0) THEN

              ! Move node i to partition 1
              move = ND_PART1_FLAG
              inv1 = inv1 + 1
              winv1 = winv1 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            ELSE IF (wnv2==0 .AND. wnv1>0) THEN
              ! Move node i to partition 2
              move = ND_PART2_FLAG
              inv2 = inv2 + 1
              winv2 = winv2 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)

            ELSE
              CALL cost_function(winv1+weight(i),winv2+1-gain1(i)-weight(i), &
                wins+gain1(i)-1,sumweight,ratio,imbal,costf,eval1)

              CALL cost_function(winv1+1-gain2(i)-weight(i),winv2+weight(i), &
                wins+gain2(i)-1,sumweight,ratio,imbal,costf,eval2)
              IF ((eval1<eval2) .OR. ((eval1==eval2) .AND. (wnv1<wnv2))) THEN
                ! Move node i to partition 1
                move = ND_PART1_FLAG
                inv1 = inv1 + 1
                winv1 = winv1 + weight(i)
                ins = ins - 1
                wins = wins - weight(i)
              ELSE
                ! Move node i to partition 2
                move = ND_PART2_FLAG
                inv2 = inv2 + 1
                winv2 = winv2 + weight(i)
                ins = ins - 1
                wins = wins - weight(i)
              END IF
            END IF
            ! Set new partition for node i
            ipart(i) = move
            ! Run through neigbours of node i to update data
            IF (i==n) THEN
              l = a_ne
            ELSE
              l = ptr(i+1) - 1
            END IF
            DO jj = ptr(i), l
              j = col(jj)
              ! Check which partition node j is in and take appropriate action
              IF (ipart(j)==move) CYCLE
              ! If node j is in cutset, update its gain value
              IF (ipart(j)==ND_SEP_FLAG) THEN
                ! If it has already been chosen in this pass just skip it
                IF (done(j)==outer .OR. dist(j)==-2) CYCLE
                ! old_gain is present gain

                old_gain = MAX(-mult,MIN(gain1(j),gain2(j)))
                ! old_gain = min(gain1(j),gain2(j))

                IF (move==ND_PART1_FLAG) gain2(j) = gain2(j) + weight(i)
                IF (move==ND_PART2_FLAG) gain1(j) = gain1(j) + weight(i)
                gain = MAX(-mult,MIN(gain1(j),gain2(j)))
                ! gain = min(gain1(j),gain2(j))

                IF (old_gain==gain) CYCLE
                ! Remove from old list
                IF (old_gain<=icut) THEN
                  CALL remove_from_list(n,mult,icut,next,last,head,j,old_gain)
                END IF
                ! gain has changed so move to new linked list if less than
                ! icut
                IF (gain<=icut) THEN
                  ! Reset ming if necessary
                  IF (gain<ming) ming = gain
                  CALL add_to_list(n,mult,icut,next,last,head,j,gain)
                END IF
              END IF
              IF (ipart(j)==2-move) THEN
                ! We have a new node in the cutset.
                ipart(j) = ND_SEP_FLAG
                ! Compute gains for this new node in the cutset and place in
                ! linked list
                ! We intentionally did not do this earlier but we do now
                ! [maybe not since won't access this node again in this pass]
                ! We use done array to record this but probably not necessary
                ! as not put
                ! in head linked list so won't be accessed
                ! First check that it was not earlier moved from cutset
                IF (done(j)/=outer .AND. dist(j)/=-2) THEN
                  ! Compute gain
                  CALL compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,j,ipart)
                  gain = MAX(-mult,MIN(gain1(j),gain2(j)))
                  ! gain = min(gain1(j),gain2(j))
                  ! !! Just added this
                  IF (gain<ming) ming = gain
                  ! Add to  list
                  IF (gain<=icut) THEN
                    CALL add_to_list(n,mult,icut,next,last,head,j,gain)
                  END IF
                END IF
                ! Update partition and gain of any nodes in cutset connected
                ! to node j
                ins = ins + 1
                wins = wins + weight(j)
                IF (move==ND_PART1_FLAG) THEN
                  inv2 = inv2 - 1
                  winv2 = winv2 - weight(j)
                END IF
                IF (move==ND_PART2_FLAG) THEN
                  inv1 = inv1 - 1
                  winv1 = winv1 - weight(j)
                END IF
                ! Check neighbours of j since any in cut set will have gain
                ! changed
                IF (j==n) THEN
                  l = a_ne
                ELSE
                  l = ptr(j+1) - 1
                END IF
                DO ii = ptr(j), l
                  eye = col(ii)
                  IF (ipart(eye)/=ND_SEP_FLAG) CYCLE
                  IF (dist(eye)==-2) CYCLE
                  ! Neighbour is in cutset. Recompute gain and insert in
                  ! linked list.
                  IF (done(eye)==outer) CYCLE
                  ! old_gain is present gain
                  old_gain = MAX(-mult,MIN(gain1(eye),gain2(eye)))
                  ! old_gain = min(gain1(eye),gain2(eye))


                  IF (move==ND_PART1_FLAG) THEN
                    gain1(eye) = gain1(eye) - weight(j)
                  END IF
                  IF (move==ND_PART2_FLAG) THEN
                    gain2(eye) = gain2(eye) - weight(j)
                  END IF
                  ! gain is new gain
                  gain = MAX(-mult,MIN(gain1(eye),gain2(eye)))
                  ! gain = min(gain1(eye),gain2(eye))
                  IF (old_gain==gain) CYCLE
                  ! Remove from old list
                  IF (old_gain<=icut) THEN
                    CALL remove_from_list(n,mult,icut,next,last,head,eye,old_gain)
                  END IF
                  ! gain has changed so move to new linked list if less than
                  ! icut
                  IF (gain<=icut) THEN
                    ! Reset ming if necessary
                    IF (gain<ming) ming = gain
                    CALL add_to_list(n,mult,icut,next,last,head,eye,gain)
                  END IF
                END DO
              END IF
              ! end of neighbours loop
            END DO

            ! ii = 0
            ! do i = 1,n
            ! if (ipart(i) == 2) ii = ii + 1
            ! enddo
            ! if (ii .ne. inv2) write(6,*) 'problem in partition',ii,inv2

            ! Evaluate new partition
            CALL cost_function(winv1+1,winv2+1,wins,sumweight,ratio,imbal, &
              costf,eval)
            ! Compare this with best so far in inner loop and store partition
            ! information if it is the best
            IF (inv1*inv2>0 .AND. nv1*nv2==0) THEN
              ! Might have to store gains and who is in the cutset
              evalc = eval
              nv1 = inv1
              nv2 = inv2
              ns = ins
              wnv1 = winv1
              wnv2 = winv2
              wns = wins
              flags = ipart

            ELSE IF (eval<evalc .AND. (inv1*inv2>0)) THEN
              ! Might have to store gains and who is in the cutset
              evalc = eval
              nv1 = inv1
              nv2 = inv2
              ns = ins
              wnv1 = winv1
              wnv2 = winv2
              wns = wins
              flags(1:n) = ipart(1:n)
            END IF
            ! End inner loop
          END DO INNER
          ! Leave loop if inner loop has not found better partition
          IF (evalc>=(1.0-1.0/(LOG(REAL(sumweight))**2.3))*evalo) EXIT
          ! Otherwise we reset evalo and go back to inner loop
          evalo = evalc
          ! Recompute gains for this new partition
          ! Compute gains for nodes in cutset
          ! This is very inefficient but is in now to test functionality
          head(-mult:icut) = 0
          ming = icut + 1
          DO i = 1, n
            IF (flags(i)/=ND_SEP_FLAG) CYCLE
            IF (dist(i)==-2) CYCLE
            ! Node i is in cutset
            ! gain1(i) is change to cutset size if node i is moved to
            ! partition 1.
            ! gain2(i) is change to cutset size if node i is moved to
            ! partition 2.
            ! Run through all neighbours of node i to see what loss/gain is if
            ! node
            ! i is moved
            CALL compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,flags)
            ! Recalculate doubly linked list linking nodes with same gain and
            ! headers
            ! to starts (up to cut off value icut)
            ! Initialize array done that flags that a node has been considered
            ! in
            ! an inner loop pass
            gain = MAX(-mult,MIN(gain1(i),gain2(i)))
            ! gain = min(gain1(i),gain2(i))
            IF (gain>icut) CYCLE
            IF (gain<ming) ming = gain
            ! New node is put at head of list
            CALL add_to_list(n,mult,icut,next,last,head,i,gain)
          END DO
          ! End of outer loop
        END DO
        a_n1 = nv1
        a_n2 = nv2
        RETURN

      END SUBROUTINE nd_fm_refinement


        SUBROUTINE remove_from_list(n,mult,icut,next,last,head,irm,ig)
          INTEGER, INTENT(IN) :: n,mult,icut ! order matrix
          INTEGER, INTENT(INOUT) :: next(n),last(n),head(-mult:icut)
          INTEGER, INTENT(IN) :: irm, ig
          INTEGER :: inext, ilast

          inext = next(irm)
          ilast = last(irm)
          IF (ilast==0) THEN
            head(ig) = inext
            IF (inext/=0) last(inext) = 0
          ELSE
            next(ilast) = inext
            IF (inext/=0) last(inext) = ilast
          END IF
        END SUBROUTINE remove_from_list


        SUBROUTINE add_to_list(n,mult,icut,next,last,head,irm,ig)
          INTEGER, INTENT(IN) :: n,mult,icut ! order matrix
          INTEGER, INTENT(INOUT) :: next(n),last(n),head(-mult:icut)
          INTEGER, INTENT(IN) :: irm, ig
          INTEGER :: inext, ilast

          inext = head(ig)
          head(ig) = irm
          next(irm) = inext
          IF (inext/=0) last(inext) = irm
          last(irm) = 0
        END SUBROUTINE add_to_list


        SUBROUTINE compute_gain(n,a_ne,ptr,col,gain1,gain2,weight,i,partit)
          INTEGER, INTENT(IN) :: n,a_ne,ptr(n),col(a_ne),weight(n)
          INTEGER, INTENT(INOUT) :: gain1(n), gain2(n)  
          INTEGER, INTENT(IN) :: i, partit(:)
          INTEGER :: j, jj, l
          ! Initialize gain ... knowing node i will be removed from cutset
          ! The +1 is to give identical result to previous code when unit
          ! weights
          gain1(i) = -weight(i) + 1
          gain2(i) = -weight(i) + 1
          IF (i==n) THEN
            l = a_ne
          ELSE
            l = ptr(i+1) - 1
          END IF
          DO jj = ptr(i), l
            j = col(jj)
            ! Check which partition node j is in and adjust gain array
            ! appropriately
            IF (partit(j)==ND_PART1_FLAG) THEN
              gain2(i) = gain2(i) + weight(j)
            END IF
            IF (partit(j)==ND_PART2_FLAG) THEN
              gain1(i) = gain1(i) + weight(j)
            END IF
          END DO
        END SUBROUTINE compute_gain

      SUBROUTINE amd_order_both(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
          iperm,a_weight,work)
        INTEGER, INTENT (IN) :: a_n ! order of matrix being partitioned
        INTEGER, INTENT (IN) :: a_ne ! no. entries in matrix being partitioned
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! col ptrs for matrix being
        ! partitioned
        INTEGER, INTENT (IN) :: a_row(a_ne) ! row indices for matrix
        ! being partitioned.
        INTEGER, INTENT (IN) :: a_n1 ! no. rows in partition 1
        INTEGER, INTENT (IN) :: a_n2 ! no. rows in partition 2
        INTEGER, INTENT (IN) :: partition(a_n) ! the partitions
        INTEGER, INTENT (INOUT) :: iperm(a_n) ! maps current permutation to
        ! the
        ! column indices of the matrix whose ordering is being computed
        INTEGER, INTENT (INOUT) :: a_weight(a_n) ! weights of vertices
        INTEGER, INTENT (OUT) :: work(12*a_n+a_ne)


        ! Local variables
        INTEGER :: i, j
        INTEGER :: extract_work ! pointers into work array for mask arrays
        INTEGER :: amd_order_perm ! pointer into work array for perm array
        INTEGER :: amd_order_work ! pointer into work array for amd_order work
        ! array
        INTEGER :: a_ptr_sub ! pointer into work for col ptrs of
        ! submatrix
        INTEGER :: a_irn_sub ! pointer into work for irn array of
        ! submatrix
        INTEGER :: rows_sub ! pointer into work for rows_sub array
        INTEGER :: a_lirn_sub ! length of irn array of submatrix
        INTEGER :: a_n_1 ! order of submatrix 1
        INTEGER :: a_n_2 ! order of submatrix 2
        INTEGER :: a_n_sep ! number entries in separator
        INTEGER :: a_ne_sub ! number entries in submatrix
        INTEGER :: len_a_row_sub ! used when extracting submatrices


        ! Set orders of submatrices
        a_n_1 = a_n1
        a_n_2 = a_n2
        a_n_sep = 0


        ! Set pointers into work array
        amd_order_perm = 0 ! length a_n
        a_ptr_sub = amd_order_perm + a_n ! length a_n
        a_lirn_sub = a_ne + MAX(a_n_1,a_n_2) + 1
        ! max(a_n_1,a_n_2) + 1 .le. a_n
        a_irn_sub = a_ptr_sub + a_n ! length a_lirn_sub
        rows_sub = a_irn_sub + a_lirn_sub ! length a_n
        extract_work = rows_sub + a_n ! length a_n
        amd_order_work = rows_sub + a_n ! length 7*a_n


        ! Form submatrix 1
        ! IF (a_n1 .NE. 1) THEN
        work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
        len_a_row_sub = a_ne
        CALL extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,a_n_sep, &
          work(rows_sub+1:rows_sub+a_n_1),a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+ &
          a_n_1),len_a_row_sub,work(a_irn_sub+1:a_irn_sub+a_ne), &
          work(extract_work+1:extract_work+a_n))

        ! Apply amd_order
        work(rows_sub+1:rows_sub+a_n1) = ND_PART1_FLAG
        CALL amd_order(a_n_1,a_ne_sub,a_lirn_sub,work(a_irn_sub+1:a_irn_sub+ &
          a_lirn_sub),work(a_ptr_sub+1:a_ptr_sub+a_n_1), &
          work(rows_sub+1:rows_sub+a_n_1),work(amd_order_perm+1:amd_order_perm+a_n_1), &
          work(amd_order_work+1:amd_order_work+7*a_n_1))

        ! Overwrite first a_n1 entries of amd_order_perm with first a_n1 entries
        ! that will form new iperm. Similarly, overwrite first a_n1 entries of
        ! rows_sub with first a_n1 entries that will form new a_weight
        ! no longer need info in a_ptr
        DO i = 1, a_n1
          j = work(amd_order_perm+i)
          work(a_ptr_sub+i) = iperm(partition(j))
          work(rows_sub+i) = a_weight(partition(j))
        END DO
        work(amd_order_perm+1:amd_order_perm+a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)

        ! Build second submatrix
        ! IF (a_n2 .NE. 1) THEN
        work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2) = partition(a_n1+1:a_n1+a_n2 &
          )
        len_a_row_sub = a_ne
        CALL extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,a_n_sep, &
          work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2),a_ne_sub, &
          work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub, &
          work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
          a_n))
        ! Apply amd_order
        work(rows_sub+a_n1+1:rows_sub+a_n1+a_n2) = ND_PART1_FLAG
        CALL amd_order(a_n_2,a_ne_sub,a_lirn_sub,work(a_irn_sub+1:a_irn_sub+ &
          a_lirn_sub),work(a_ptr_sub+1:a_ptr_sub+a_n_2), &
          work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2), &
          work(amd_order_perm+a_n1+1:amd_order_perm+a_n),work(amd_order_work+1:amd_order_work+7* &
          a_n_2))
        ! ELSE
        ! work(amd_order_perm+a_n1+1) = 1
        ! END IF
        DO i = 1, a_n_2
          j = work(amd_order_perm+a_n1+i)
          work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
          work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))
        END DO
        DO i = a_n_2 + 1, a_n_2 + (a_n-a_n1-a_n2)
          j = i
          work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
          work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))

        END DO
        iperm(a_n1+1:a_n) = work(a_ptr_sub+1+a_n1:a_ptr_sub+a_n)
        iperm(1:a_n1) = work(amd_order_perm+1:amd_order_perm+a_n1)
        a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

      END SUBROUTINE amd_order_both


      SUBROUTINE extract_matrix(a_n,a_ne,a_ptr,a_row,a_n_part,a_n_sep, &
          rows_sub,a_ne_sub,a_ptr_sub,len_a_row_sub,a_row_sub,work)
        INTEGER, INTENT (IN) :: a_n ! order of matrix being partitioned
        INTEGER, INTENT (IN) :: a_ne ! no. entries in matrix being partitioned
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! col ptrs for matrix being
        ! partitioned
        INTEGER, INTENT (IN) :: a_row(a_ne) ! row indices for matrix
        ! being partitioned.
        INTEGER, INTENT (IN) :: a_n_part ! no. rows in partition
        INTEGER, INTENT (IN) :: a_n_sep ! no. rows in partition
        INTEGER, INTENT (IN) :: rows_sub(a_n_part+a_n_sep) ! rows/cols of
        ! matrix
        ! to be extracted. Intersecting rows/cols of separator will be
        ! replaced
        ! by matrix of all zeros
        INTEGER, INTENT (OUT) :: a_ne_sub ! no. entries stored in extracted
        ! matrix
        INTEGER, INTENT (OUT) :: a_ptr_sub(a_n_part+a_n_sep) ! col ptrs for
        ! extracted matrix
        INTEGER, INTENT (IN) :: len_a_row_sub ! length of a_row_sub
        INTEGER, INTENT (OUT) :: a_row_sub(len_a_row_sub) ! row indices for
        ! extracted matrix
        INTEGER, INTENT (OUT) :: work(a_n)

        ! Local variables
        INTEGER :: i, j, k, l, m, p
        INTEGER :: a_n_sub ! Order of extracted matrix
        INTEGER :: mask ! pointer into work array for mask arrays

        ! Set pointers into work array
        mask = 0 ! length a_n

        ! Set mask
        a_n_sub = a_n_part + a_n_sep
        work(mask+1:mask+a_n) = 0
        DO i = 1, a_n_part
          j = rows_sub(i)
          work(mask+j) = i
        END DO
        DO i = a_n_part + 1, a_n_sub
          j = rows_sub(i)
          work(mask+j) = -i
        END DO
        a_row_sub(:) = 0

        ! Count number of entries in  submatrix and set-up column ptrs
        a_ptr_sub(1:a_n_sub) = 0
        DO j = 1, a_n_part
          a_ptr_sub(j) = 0
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)/=0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) + 1

            END IF
          END DO
        END DO
        DO j = a_n_part + 1, a_n_sub
          a_ptr_sub(j) = 0
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)>0) a_ptr_sub(j) = a_ptr_sub(j) + 1
          END DO
        END DO
        a_ptr_sub(1) = a_ptr_sub(1) + 1
        DO j = 2, a_n_sub
          a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
        END DO
        a_ne_sub = a_ptr_sub(a_n_sub) - 1

        ! Form a_row_sub
        DO j = 1, a_n_part
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)/=0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) - 1
              p = a_ptr_sub(j)
              a_row_sub(p) = ABS(work(mask+m))
            END IF
          END DO
        END DO
        DO j = a_n_part + 1, a_n_sub
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)>0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) - 1
              p = a_ptr_sub(j)
              a_row_sub(p) = work(mask+m)
            END IF
          END DO
        END DO
      END SUBROUTINE extract_matrix

      SUBROUTINE extract_both_matrices(a_n,a_ne,a_ptr,a_row,a_n_part1, &
          a_n_part2,rows_sub,a_ne_sub1,a_ne_sub2,a_ptr_sub,len_a_row_sub, &
          a_row_sub,work)
        INTEGER, INTENT (IN) :: a_n ! order of matrix being partitioned
        INTEGER, INTENT (IN) :: a_ne ! no. entries in matrix being partitioned
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! col ptrs for matrix being
        ! partitioned
        INTEGER, INTENT (IN) :: a_row(a_ne) ! row indices for matrix
        ! being partitioned.
        INTEGER, INTENT (IN) :: a_n_part1 ! no. rows in partition 1
        INTEGER, INTENT (IN) :: a_n_part2 ! no. rows in partition 2
        INTEGER, INTENT (IN) :: rows_sub(a_n_part1+a_n_part2) ! rows/cols of
        ! matrices to be extracted. First a_n_part1 entries contain the
        ! rows/cols forming the first matrix to be extracted.
        INTEGER, INTENT (OUT) :: a_ne_sub1 ! no. entries in extracted matrix 1
        INTEGER, INTENT (OUT) :: a_ne_sub2 ! no. entries in extracted matrix 2
        INTEGER, INTENT (OUT) :: a_ptr_sub(a_n_part1+a_n_part2) ! col ptrs for
        ! extracted matrices. First a_n_part1 are for first matrix, etc
        INTEGER, INTENT (IN) :: len_a_row_sub ! length of a_row_sub
        INTEGER, INTENT (OUT) :: a_row_sub(len_a_row_sub) ! row indices for
        ! extracted matrices. First a_ne_part1 entries are for first matrix;
        ! the immediately following a_ne_part2 entries are for second matrix
        INTEGER, INTENT (OUT) :: work(a_n)

        ! Local variables
        INTEGER :: i, j, k, l, m, p
        INTEGER :: mask ! pointer into work array for mask arrays


        ! Set pointers into work array
        mask = 0 ! length a_n

        ! Set mask
        work(mask+1:mask+a_n) = 0
        DO i = 1, a_n_part1
          j = rows_sub(i)
          work(mask+j) = i
        END DO
        DO i = a_n_part1 + 1, a_n_part1 + a_n_part2
          j = rows_sub(i)
          work(mask+j) = -i + a_n_part1
        END DO
        a_row_sub(:) = 0

        ! Count number of entries in each submatrix and set-up column ptrs
        a_ptr_sub(1:a_n_part1+a_n_part2) = 0
        DO j = 1, a_n_part1
          a_ptr_sub(j) = 0
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)>0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) + 1

            END IF
          END DO
        END DO

        DO j = a_n_part1 + 1, a_n_part1 + a_n_part2
          a_ptr_sub(j) = 0
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)<0) a_ptr_sub(j) = a_ptr_sub(j) + 1
          END DO
        END DO

        a_ptr_sub(1) = a_ptr_sub(1) + 1
        DO j = 2, a_n_part1
          a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
        END DO
        a_ne_sub1 = a_ptr_sub(a_n_part1) - 1

        a_ptr_sub(a_n_part1+1) = a_ptr_sub(a_n_part1+1) + 1
        DO j = a_n_part1 + 2, a_n_part1 + a_n_part2
          a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
        END DO
        a_ne_sub2 = a_ptr_sub(a_n_part1+a_n_part2) - 1

        ! Form a_row_sub
        DO j = 1, a_n_part1
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)>0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) - 1
              p = a_ptr_sub(j)
              a_row_sub(p) = ABS(work(mask+m))
            END IF
          END DO
        END DO

        DO j = a_n_part1 + 1, a_n_part1 + a_n_part2
          i = rows_sub(j)
          IF (i==a_n) THEN
            l = a_ne
          ELSE
            l = a_ptr(i+1) - 1
          END IF
          DO k = a_ptr(i), l
            m = a_row(k)
            IF (work(mask+m)<0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) - 1
              p = a_ptr_sub(j)
              a_row_sub(p+a_ne_sub1) = -work(mask+m)
            END IF
          END DO
        END DO
      END SUBROUTINE extract_both_matrices

      SUBROUTINE amd_order_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,iperm, &
          a_weight,a_ne_sub,work)
        ! Apply amd_order to the smaller partition and extract the other matrix
        ! into
        ! a_ptr and a_row
        INTEGER, INTENT (IN) :: a_n ! order of matrix being partitioned
        INTEGER, INTENT (IN) :: a_ne ! no. entries in matrix being partitioned
        INTEGER, INTENT (INOUT) :: a_ptr(a_n) ! col ptrs for matrix being
        ! partitioned and, on return, for the extracted submatrix
        INTEGER, INTENT (INOUT) :: a_row(a_ne) ! row indices for matrix
        ! being partitioned and, on return, for the extracted submatrix
        INTEGER, INTENT (IN) :: a_n1 ! no. rows in partition 1
        INTEGER, INTENT (IN) :: a_n2 ! no. rows in partition 2
        INTEGER, INTENT (IN) :: partition(a_n) ! the partitions
        INTEGER, INTENT (INOUT) :: iperm(a_n) ! maps current permuation to the
        ! column indices of the matrix whose ordering is being computed
        INTEGER, INTENT (INOUT) :: a_weight(a_n) ! weights of vertices
        INTEGER, INTENT (OUT) :: a_ne_sub ! number entries in returned
        ! submatrix
        INTEGER, INTENT (OUT) :: work(12*a_n+a_ne)

        ! Local variables
        INTEGER :: i, j
        INTEGER :: extract_work ! pointers into work array for mask arrays
        INTEGER :: amd_order_perm ! pointer into work array for perm array
        INTEGER :: amd_order_work ! pointer into work array for amd_order work
        ! array
        INTEGER :: a_ptr_sub ! pointer into work for col ptrs of
        ! submatrix
        INTEGER :: a_irn_sub ! pointer into work for irn array of
        ! submatrix
        INTEGER :: rows_sub ! pointer into work for rows_sub array
        INTEGER :: a_lirn_sub ! length of irn array of submatrix
        INTEGER :: a_n_1 ! order of submatrix 1
        INTEGER :: a_n_2 ! order of submatrix 2
        INTEGER :: a_n_sep ! number entries in separator
        INTEGER :: len_a_row_sub ! used when extracting submatrices

        ! Set orders of submatrices
        IF (a_n1<a_n2) THEN
          ! Applying amd_order to first partition
          a_n_1 = a_n - a_n2
          a_n_2 = a_n2
          a_n_sep = a_n - a_n1 - a_n2
        ELSE
          ! Applying amd_order to second partition
          a_n_2 = a_n - a_n1
          a_n_1 = a_n1
          a_n_sep = a_n - a_n1 - a_n2

        END IF


        ! Set pointers into work array
        amd_order_perm = 0 ! length a_n
        a_ptr_sub = amd_order_perm + a_n ! length a_n
        IF (a_n1<a_n2) THEN
          a_lirn_sub = a_ne + a_n_1 + 1
        ELSE
          a_lirn_sub = a_ne + a_n_2 + 1
        END IF
        ! max(a_n_1,a_n_2) + 1 .le. a_n
        a_irn_sub = a_ptr_sub + a_n ! length a_lirn_sub
        rows_sub = a_irn_sub + a_lirn_sub ! length a_n
        extract_work = rows_sub + a_n ! length a_n
        amd_order_work = rows_sub + a_n ! length 7*a_n

        ! Extract matrix that amd_order is applied to
        IF (a_n1<a_n2) THEN

          ! Start by updating iperm and a_weight
          DO i = 1, a_n
            work(amd_order_perm+i) = iperm(partition(i))
            work(rows_sub+i) = a_weight(partition(i))
          END DO
          iperm(1:a_n) = work(amd_order_perm+1:amd_order_perm+a_n)
          a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

          ! Form submatrix 1
          work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
          work(rows_sub+a_n1+1:rows_sub+a_n_1) = partition(a_n1+a_n2+1:a_n)
          len_a_row_sub = a_ne
          CALL extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,a_n_sep, &
            work(rows_sub+1:rows_sub+a_n_1),a_ne_sub, &
            work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub, &
            work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
            a_n))


          ! Apply amd_order
          work(rows_sub+1:rows_sub+a_n1) = ND_PART1_FLAG
          work(rows_sub+a_n1+1:rows_sub+a_n_1) = ND_SEP_FLAG
          CALL amd_order(a_n_1,a_ne_sub,a_lirn_sub,work(a_irn_sub+1:a_irn_sub+ &
            a_lirn_sub),work(a_ptr_sub+1:a_ptr_sub+a_n_1), &
            work(rows_sub+1:rows_sub+a_n_1),work(amd_order_perm+1:amd_order_perm+a_n_1), &
            work(amd_order_work+1:amd_order_work+7*a_n_1))

          ! Overwrite first a_n1 entries of amd_order_perm with first a_n1 entries
          ! that will form new iperm. Similarly, overwrite first a_n1 entries
          ! of
          ! rows_sub with first a_n1 entries that will form new a_weight
          ! no longer need info in a_ptr
          ! no longer need info in a_ptr
          DO i = 1, a_n1
            j = work(amd_order_perm+i)
            work(a_ptr_sub+i) = iperm(j)
            work(rows_sub+i) = a_weight(j)
          END DO
          ! work(amd_order_perm+1:amd_order_perm+a_n1) =
          ! work(a_ptr_sub+1:a_ptr_sub+a_n1)


          iperm(1:a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)
          a_weight(1:a_n1) = work(rows_sub+1:rows_sub+a_n1)

          ! Build second submatrix
          work(rows_sub+1:rows_sub+a_n_2) = partition(a_n1+1:a_n1+a_n2)
          len_a_row_sub = a_ne
          CALL extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,0, &
            work(rows_sub+1:rows_sub+a_n_2),a_ne_sub, &
            work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub, &
            work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
            a_n))

          ! Copy matrix to a_ptr and a_row
          a_ptr(1:a_n_2) = work(a_ptr_sub+1:a_ptr_sub+a_n_2)
          a_row(1:a_ne_sub) = work(a_irn_sub+1:a_irn_sub+a_ne_sub)

        ELSE


          ! Form submatrix 2
          work(rows_sub+1:rows_sub+a_n_2) = partition(a_n1+1:a_n)
          len_a_row_sub = a_ne
          CALL extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,a_n_sep, &
            work(rows_sub+1:rows_sub+a_n_2),a_ne_sub, &
            work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub, &
            work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
            a_n))

          ! Apply amd_order
          work(rows_sub+1:rows_sub+a_n2) = ND_PART2_FLAG
          work(rows_sub+a_n2+1:rows_sub+a_n_2) = ND_SEP_FLAG
          CALL amd_order(a_n_2,a_ne_sub,a_lirn_sub,work(a_irn_sub+1:a_irn_sub+ &
            a_lirn_sub),work(a_ptr_sub+1:a_ptr_sub+a_n_2), &
            work(rows_sub+1:rows_sub+a_n_2),work(amd_order_perm+a_n1+1:amd_order_perm+ &
            a_n),work(amd_order_work+1:amd_order_work+7*a_n_2))

          DO i = 1, a_n_2
            j = work(amd_order_perm+a_n1+i)
            work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
            work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))
          END DO

          DO i = 1, a_n_1
            work(a_ptr_sub+i) = iperm(partition(i))
            work(rows_sub+i) = a_weight(partition(i))
          END DO
          iperm(1:a_n) = work(a_ptr_sub+1:a_ptr_sub+a_n)
          a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

          ! Form submatrix 1
          work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
          len_a_row_sub = a_ne
          CALL extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,0, &
            work(rows_sub+1:rows_sub+a_n_1),a_ne_sub, &
            work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub, &
            work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+ &
            a_n))

          ! Copy matrix to a_ptr and a_row
          a_ptr(1:a_n_1) = work(a_ptr_sub+1:a_ptr_sub+a_n_1)
          a_row(1:a_ne_sub) = work(a_irn_sub+1:a_irn_sub+a_ne_sub)

        END IF

      END SUBROUTINE amd_order_one

      SUBROUTINE multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
          partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,control, &
          info1,lwork,work,stop_coarsening2,grid)

        INTEGER, INTENT (IN) :: a_n ! order of matrix being partitioned
        INTEGER, INTENT (IN) :: a_ne ! no. entries in matrix being partitioned
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! col ptrs for matrix being
        ! partitioned and, on return, for the extracted submatrix
        INTEGER, INTENT (IN) :: a_row(a_ne) ! row indices for matrix
        INTEGER, INTENT (IN) :: a_weight(a_n) ! weights associated with rows
        ! of matrix (useful if matrix has already been compressed)
        INTEGER, INTENT (IN) :: sumweight ! sum of entries in a_weight
        INTEGER, INTENT (OUT) :: partition(a_n) ! computed partition
        INTEGER, INTENT (OUT) :: a_n1 ! number of entries in partition 1
        INTEGER, INTENT (OUT) :: a_n2 ! number of entries in partition 2
        INTEGER, INTENT (OUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weight
        ! ed
        ! size of partitions and separator
        TYPE (nd_options), INTENT (IN) :: control
        INTEGER, INTENT (IN) :: lwork ! length of work array: must be atleast
        ! 9a_n + sumweight
        INTEGER, INTENT (OUT) :: work(lwork) ! work array
        INTEGER, INTENT (INOUT) :: info1
        INTEGER, INTENT (IN) :: stop_coarsening2 ! no. levels in the
        ! multilevel grid (default
        ! 10)

        TYPE (nd_multigrid), INTENT (INOUT) :: grid ! the multilevel of
        ! graphs (matrices)

        INTEGER :: i, j, k, inv1, inv2, ins
        INTEGER :: mp
        INTEGER :: mglevel_cur ! current level
        INTEGER :: err, print_level ! printing
        LOGICAL :: lerr

        info1 = 0
        ! Set up printing
        IF (control%print_level<0) print_level = 0
        ! The default is control%print_level = 0
        IF (control%print_level==0) print_level = 1
        IF (control%print_level==1) print_level = 2
        IF (control%print_level>1) print_level = 3
        mp = control%unit_diagnostics
        IF (mp<0) print_level = 0
        ! Set error controls
        lerr = control%unit_error >= 0 .AND. print_level > 0
        err = control%unit_error


        IF (print_level>1) THEN
          WRITE (mp,'(a)') 'Start multilevel_partition:'
        END IF

        ! construct the multigrid at this level

        IF ( .NOT. ASSOCIATED(grid%graph)) ALLOCATE (grid%graph)

        CALL nd_matrix_construct(grid%graph,a_n,a_n,a_ne,info1)
        IF (info1<0) THEN
          IF (lerr) CALL nd_print_message(info1,err, &
            ' multilevel_partition')
          RETURN
        END IF

        grid%graph%ptr(1:a_n) = a_ptr(1:a_n)
        grid%graph%ptr(a_n+1) = a_ne + 1
        grid%graph%col(1:a_ne) = a_row(1:a_ne)

        DO i = 1, a_n - 1
          DO j = a_ptr(i), a_ptr(i+1) - 1
            k = a_row(j)
            grid%graph%val(j) = a_weight(i)*a_weight(k)
          END DO
        END DO
        DO j = a_ptr(a_n), a_ne
          k = a_row(j)
          grid%graph%val(j) = a_weight(a_n)*a_weight(k)
        END DO

        grid%size = a_n
        grid%level = 1
        ! NULLIFY (grid%p)

        CALL nd_assoc(grid%where,a_n,info1)
        IF (info1<0) THEN
          IF (lerr) CALL nd_print_message(info1,err, &
            ' multilevel_partition')
          RETURN
        END IF

        CALL nd_assoc(grid%row_wgt,a_n,info1)
        IF (info1<0) THEN
          IF (lerr) CALL nd_print_message(info1,err, &
            ' multilevel_partition')
          RETURN
        END IF

        ! Initialise row weights
        grid%row_wgt(1:a_n) = a_weight(1:a_n)

        ! initialise mglevel_cur to the maximum number of levels
        ! allowed for this bisection
        mglevel_cur = stop_coarsening2
        CALL multilevel(grid,control,sumweight,mglevel_cur,mp,print_level, &
          lwork,work,info1)

        IF (info1/=0) THEN
          IF (lerr) CALL nd_print_message(info1,err, &
            ' multilevel_partition')
          RETURN
        END IF

        inv1 = 1
        inv2 = grid%part_div(1) + 1
        ins = grid%part_div(1) + grid%part_div(2) + 1

        a_weight_1 = 0
        a_weight_2 = 0
        a_weight_sep = 0
        DO i = 1, a_n
          SELECT CASE (grid%where(i))
          CASE (ND_PART1_FLAG)
            partition(inv1) = i
            inv1 = inv1 + 1
            a_weight_1 = a_weight_1 + a_weight(i)
          CASE (ND_PART2_FLAG)
            partition(inv2) = i
            inv2 = inv2 + 1
            a_weight_2 = a_weight_2 + a_weight(i)
          CASE DEFAULT
            partition(ins) = i
            ins = ins + 1
            a_weight_sep = a_weight_sep + a_weight(i)
          END SELECT
        END DO

        a_n1 = grid%part_div(1)
        a_n2 = grid%part_div(2)

        IF (.FALSE.) THEN
          WRITE (*,'(a)') ' '
          WRITE (*,'(a)') 'Multilevel partition found'
          WRITE (*,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, ', a_n2=', a_n2, &
            ', a_n_sep=', a_n - a_n1 - a_n2
          WRITE (*,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', a_weight_1, &
            ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
            sumweight - a_weight_1 - a_weight_2
        END IF

        ! deallocate the finest level
        ! CALL multigrid_deallocate_first(a_n,a_n,grid,info1)
        IF (info1/=0) THEN
          IF (lerr) CALL nd_print_message(info1,err, &
            ' multilevel_partition')
          RETURN
        END IF

        ! DEALLOCATE (matrix%ptr,matrix%col,matrix%val,STAT=st)
        ! IF (st/=0) info1 = ND_ERR_MEMORY_DEALLOC
        IF (info1<0) THEN
          IF (lerr) CALL nd_print_message(info1,err, &
            ' multilevel_partition')
          RETURN
        END IF

        IF (print_level>2) THEN
          WRITE (mp,'(a)') 'multilevel_partition: successful completion'
        END IF

      END SUBROUTINE multilevel_partition

      ! ********************************************************

      ! main subroutine for computing multilevel structure.
      ! Offers heavy-edge collapsing and maximal independent vertex
      ! set for coarsening. We will need to test out to see
      ! which is better.

      RECURSIVE SUBROUTINE multilevel(grid,control,sumweight,mglevel_cur,mp, &
          print_level,lwork,work,info)

        REAL (wp), PARAMETER :: half = 0.5_wp
        REAL (wp), PARAMETER :: one = 1.0_wp

        ! Arguments
        TYPE (nd_multigrid), INTENT (INOUT), TARGET :: grid ! this level
        ! of matrix (grid)
        TYPE (nd_options), INTENT (IN) :: control
        INTEGER, INTENT (IN) :: sumweight ! sum of weights (unchanged between
        ! coarse and fine grid
        INTEGER, INTENT (INOUT) :: mglevel_cur ! current grid level
        INTEGER, INTENT (IN) :: mp, print_level ! diagnostic printing
        INTEGER, INTENT (IN) :: lwork ! length of work array
        ! (.GE.9*grid%graph%n +sumweight)
        INTEGER, INTENT (OUT) :: work(lwork) ! work array
        INTEGER, INTENT (INOUT) :: info ! Error flag

        ! Local variables
        TYPE (nd_multigrid), POINTER :: cgrid ! the coarse level grid
        INTEGER :: cnvtx ! number of vertices (rows) in the coarse
        ! matrix
        TYPE (nd_matrix), POINTER :: p ! the coarse grid prolongator
        TYPE (nd_matrix), POINTER :: r ! the coarse grid restrictor (= p')

        INTEGER, DIMENSION (:), POINTER :: fwhere ! partition on fine grid
        INTEGER, DIMENSION (:), POINTER :: cwhere ! partition on coarse grid
        TYPE (nd_matrix), POINTER :: cgraph ! the coarse graph
        TYPE (nd_matrix), POINTER :: graph ! the fine graph
        INTEGER, DIMENSION (:), POINTER :: row_wgt ! fine
        ! graph vertex weights
        INTEGER, DIMENSION (:), POINTER :: crow_wgt ! coarse
        ! graph vertex weights
        REAL (wp) :: grid_rdc_fac_min ! min grid reduction
        ! factor
        REAL (wp) :: grid_rdc_fac_max ! max grid reduction
        ! factor
        REAL (wp) :: one1
        INTEGER :: stop_coarsening1 ! controls when to stop coarsening
        INTEGER :: partition_ptr, part_ptr, work_ptr, a_ne, ref_control, &
          clwork
        INTEGER :: i, j, k, l, a_weight_1, a_weight_2, a_weight_sep, &
          ref_method, lwk
        INTEGER :: a_n1_new, a_n2_new, a_weight_1_new, a_weight_2_new, &
          a_weight_sep_new
        LOGICAL :: imbal
        REAL (wp) :: tau, ratio, tau_best
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!
        info = 0
        one1 = 1.0

        stop_coarsening1 = MAX(2,control%stop_coarsening1)
        IF (print_level.GE.2) CALL level_print(mp,'size of grid on level ', &
          grid%level,' is ',REAL(grid%size,wp))

        grid_rdc_fac_min = MAX(0.01_wp,control%min_reduction)
        ! max grid reduction factor must be at least half and at most one
        grid_rdc_fac_max = MAX(half,control%max_reduction)
        grid_rdc_fac_max = MIN(one,grid_rdc_fac_max)

        ! Test to see if this is either the last level or
        ! if the matrix size too small
        IF (grid%level>=mglevel_cur .OR. grid%size<=stop_coarsening1) THEN
          IF (print_level.GE.2) CALL level_print(mp,'end of level ',grid%level)

          ! coarsest level in multilevel so compute separator
          a_ne = grid%graph%ptr(grid%graph%n+1) - 1
          CALL nd_coarse_partition(grid%graph%n,a_ne,grid%graph%ptr, &
            grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
            grid%part_div(2),grid%where,lwork,work,control,info)
          RETURN
        END IF

        ! Coarsest level not yet reached so carry on coarsening
        IF (control%matching==1) THEN
          lwk = grid%size
          CALL coarsen_hec(grid,lwk,work(1:lwk),info)
        ELSE
          IF (control%matching>1) THEN
            lwk = 3*grid%size
            CALL coarsen_best(grid,lwk,work(1:lwk),info)
          ELSE
            lwk = 2*grid%size
            CALL coarsen_cn(grid,lwk,work(1:lwk))
          END IF
        END IF
        IF (info<0) RETURN

        cgrid => grid%coarse
        cnvtx = cgrid%size
        ! allocate coarse grid quantities
        CALL nd_assoc(cgrid%where,cnvtx,info)
        IF (info/=0) THEN
          RETURN
        END IF

        CALL nd_assoc(cgrid%row_wgt,cnvtx,info)
        IF (info/=0) THEN
          RETURN
        END IF

        ! see if the grid reduction is achieved, if not, set the allowed
        ! maximum level to current level and partition this level
        ! deallocate the coarse grid quantities that haves been allocated so
        ! far
        IF (REAL(cgrid%size)/REAL(grid%size)>grid_rdc_fac_max .OR. &
            REAL(cgrid%size)/REAL(grid%size)<grid_rdc_fac_min .OR. &
            cgrid%size<4) THEN

          IF (print_level.GE.2) THEN
            ! IF (.true.) THEN
            WRITE (mp,'(a,i10,a,f12.4,i4)') 'at level ', grid%level, &
              ' further coarsening gives reduction factor', &
              cgrid%size/REAL(grid%size)
            WRITE (mp,'(a,i10)') 'current size = ', grid%size
          END IF

          ! set current grid level and recurse
          mglevel_cur = grid%level

          CALL multilevel(grid,control,sumweight,mglevel_cur,mp,print_level, &
            lwork,work,info)
          IF (info<0) RETURN

          RETURN
        END IF

        ! restriction ================

        ! form the coarse grid graph and matrix
        ! cmatrix = P^T*matrix = R*matrix
        p => cgrid%p
        r => cgrid%r
        graph => grid%graph
        cgraph => cgrid%graph

        ! get the coarse matrix
        lwk = 3*grid%size
        CALL galerkin_graph(graph,p,r,cgraph,info,lwk,work(1:lwk))
        IF (info<0) RETURN

        ! check if matrix is full
        IF (REAL(cgrid%graph%ptr(cgrid%graph%n+1)-1)/REAL(cgrid%graph%n)>=REAL &
            (cgrid%graph%n-1)) THEN
          IF (print_level.GE.2) THEN
            WRITE (mp,'(a,i10,a)') 'at level ', grid%level, &
              ' further coarsening gives full matrix'
          END IF

          ! set current grid level and recurse
          mglevel_cur = grid%level - 1
          CALL multilevel(grid,control,sumweight,mglevel_cur,mp,print_level, &
            lwork,work,info)
          IF (info<0) RETURN

          RETURN
        END IF

        ! row weight cw = R*w
        row_wgt => grid%row_wgt(1:grid%size)
        crow_wgt => cgrid%row_wgt(1:cgrid%size)
        CALL nd_matrix_multiply_vec(r,row_wgt,crow_wgt)
        clwork = 9*cgrid%graph%n + sumweight
        CALL multilevel(cgrid,control,sumweight,mglevel_cur,mp,print_level, &
          clwork,work(1:clwork),info)


        ! check if partition is returned
        IF (cgrid%part_div(1)==0 .OR. cgrid%part_div(2)==0) THEN
          ! Unlikely to be called because 99.999% of cases caught in full 
          ! matrix check above. Follows same procedure as when full matrix found
          IF (print_level.GE.2) THEN
            WRITE (mp,'(a,i10,a)') 'at level ', grid%level, &
              ' no partition found'
          END IF

          ! set current grid level and recurse
          mglevel_cur = grid%level - 1
          CALL multilevel(grid,control,sumweight,mglevel_cur,mp,print_level, &
            lwork,work,info)
          IF (info<0) RETURN

          RETURN
        END IF

        ! prolongation ================

        ! injection of the order from coarse grid to the
        ! fine grid, since cwhere(i) is the index of the
        ! i-th vertex in the new ordering, the order
        ! of this vertex should be where(i)
        ! grid%where = P*order_on_coarse_grid
        ! here P is a special matrix with only one non-zero entry per row
        fwhere => grid%where(1:grid%size)
        cwhere => cgrid%where(1:cgrid%size)
        grid%part_div(1:2) = 0
        CALL nd_matrix_multiply_vec(p,cwhere,fwhere)

        DO i = 1, grid%size
          IF (fwhere(i)==ND_PART1_FLAG) THEN
            grid%part_div(1) = grid%part_div(1) + 1
          ELSE
            IF (fwhere(i)==ND_PART2_FLAG) THEN
              grid%part_div(2) = grid%part_div(2) + 1
            END IF
          END IF
        END DO
        a_weight_1 = 0
        a_weight_2 = 0
        a_weight_sep = 0

        ! Set partition
        partition_ptr = 0
        work_ptr = partition_ptr + grid%graph%n
        i = 1
        j = grid%part_div(1) + 1
        k = grid%part_div(1) + grid%part_div(2) + 1
        DO l = 1, grid%size
          SELECT CASE (grid%where(l))
          CASE (ND_PART1_FLAG)
            work(partition_ptr+i) = l
            a_weight_1 = a_weight_1 + grid%row_wgt(l)
            i = i + 1
          CASE (ND_PART2_FLAG)
            work(partition_ptr+j) = l
            a_weight_2 = a_weight_2 + grid%row_wgt(l)
            j = j + 1
          CASE (ND_SEP_FLAG)
            work(partition_ptr+k) = l
            a_weight_sep = a_weight_sep + grid%row_wgt(l)
            k = k + 1
          END SELECT
        END DO
        a_ne = grid%graph%ptr(grid%graph%n+1) - 1

        IF (a_weight_sep>0) THEN
          ! Do not refine if separable graph

          IF (control%refinement>6) THEN
            ref_control = 3
          ELSE
            IF (control%refinement<1) THEN
              ref_control = 1
            ELSE
              ref_control = control%refinement
            END IF
          END IF

          SELECT CASE (ref_control)
          CASE (1)
            ref_method = 1

          CASE (2)
            ref_method = 2

          CASE (3)
            IF (REAL(MAX(a_weight_1,a_weight_2))/REAL(MIN(a_weight_1, &
                a_weight_2)+a_weight_sep)>MAX(REAL(1.0, &
                wp),control%balance)) THEN
              ref_method = 2
            ELSE
              ref_method = 1
            END IF

          CASE (4)
            ref_method = 0

          CASE (5)
            ref_method = 2

          CASE (6)
            IF (REAL(MAX(a_weight_1,a_weight_2))/REAL(MIN(a_weight_1, &
                a_weight_2)+a_weight_sep)>MAX(REAL(1.0, &
                wp),control%balance)) THEN
              ref_method = 2
            ELSE
              ref_method = 0
            END IF
          END SELECT

          SELECT CASE (ref_method)
          CASE (0)
            CALL nd_refine_max_flow(grid%graph%n,a_ne,grid%graph%ptr, &
              grid%graph%col,grid%row_wgt,grid%part_div(1),grid%part_div(2), &
              a_weight_1,a_weight_2,a_weight_sep,work(partition_ptr+1: &
              partition_ptr+grid%graph%n),work(work_ptr+1:work_ptr+8),control)
          CASE (1)
            IF (MIN(a_weight_1,a_weight_2)+a_weight_sep< &
                MAX(a_weight_1,a_weight_2)) THEN
              CALL nd_refine_block_trim(grid%graph%n,a_ne, &
                grid%graph%ptr,grid%graph%col,grid%row_wgt,sumweight, &
                grid%part_div(1),grid%part_div(2),a_weight_1,a_weight_2, &
                a_weight_sep,work(partition_ptr+1:partition_ptr+grid%graph%n), &
                work(work_ptr+1:work_ptr+5*grid%graph%n),control)
            ELSE
              CALL nd_refine_trim(grid%graph%n,a_ne,grid%graph%ptr, &
                grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
                grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep, &
                work(partition_ptr+1:partition_ptr+grid%graph%n), &
                work(work_ptr+1:work_ptr+3*grid%graph%n),control)

            END IF
          CASE (2)
            CALL nd_refine_edge(grid%graph%n,a_ne,grid%graph%ptr, &
              grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
              grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep, &
              work(partition_ptr+1:partition_ptr+grid%graph%n), &
              work(work_ptr+1:work_ptr+3*grid%graph%n),control)
          END SELECT

          IF (control%max_improve_cycles>0) THEN
            ratio = MAX(REAL(1.0,wp),control%balance)
            IF (ratio>REAL(sumweight-2)) THEN
              imbal = .FALSE.
            ELSE
              imbal = .TRUE.
            END IF
            CALL cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight, &
              ratio,imbal,control%cost_function,tau_best)
            a_n1_new = grid%part_div(1)
            a_n2_new = grid%part_div(2)
            a_weight_1_new = a_weight_1
            a_weight_2_new = a_weight_2
            a_weight_sep_new = a_weight_sep
          END IF

          part_ptr = work_ptr + 5*grid%graph%n
          work(part_ptr+1:part_ptr+grid%graph%n) = work(partition_ptr+1: &
            partition_ptr+grid%graph%n)

          k = control%max_improve_cycles
          DO i = 1, k

            CALL expand_partition(grid%graph%n,a_ne,grid%graph%ptr, &
              grid%graph%col,grid%row_wgt,a_n1_new,a_n2_new,a_weight_1_new, &
              a_weight_2_new,a_weight_sep_new,work(part_ptr+1:part_ptr+grid% &
              graph%n),work(work_ptr+1:work_ptr+5*grid%graph%n))


            ! call
            ! check_partition1(a_n,a_ne,a_ptr,a_row,a_n1_new,a_n2_new,work(par
            ! t_ptr+1:part_ptr+a_n))

            SELECT CASE (ref_control)

            CASE (3)
              IF (REAL(MAX(a_weight_1_new,a_weight_2_new))/REAL(MIN( &
                  a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>MAX(REAL( &
                  1.0,wp),control%balance)) THEN
                ref_method = 2
              ELSE
                ref_method = 1
              END IF

            CASE (6)
              IF (REAL(MAX(a_weight_1_new,a_weight_2_new))/REAL(MIN( &
                  a_weight_1_new,a_weight_2_new)+a_weight_sep_new)>MAX(REAL( &
                  1.0,wp),control%balance)) THEN
                ref_method = 2
              ELSE
                ref_method = 0
              END IF
            END SELECT


            SELECT CASE (ref_method)

            CASE (0)
              CALL nd_refine_max_flow(grid%graph%n,a_ne,grid%graph%ptr, &
                grid%graph%col,grid%row_wgt,a_n1_new,a_n2_new,a_weight_1_new, &
                a_weight_2_new,a_weight_sep_new,work(part_ptr+1:part_ptr+grid% &
                graph%n),work(work_ptr+1:work_ptr+8),control)

            CASE (1)
              IF (MIN(a_weight_1,a_weight_2)+a_weight_sep< &
                  MAX(a_weight_1,a_weight_2)) THEN
                CALL nd_refine_block_trim(grid%graph%n,a_ne, &
                  grid%graph%ptr,grid%graph%col,grid%row_wgt,sumweight, &
                  a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new, &
                  a_weight_sep_new,work(part_ptr+1:part_ptr+grid%graph%n), &
                  work(work_ptr+1:work_ptr+5*grid%graph%n),control)
              ELSE
                CALL nd_refine_trim(grid%graph%n,a_ne,grid%graph%ptr, &
                  grid%graph%col,grid%row_wgt,sumweight,a_n1_new,a_n2_new, &
                  a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
                  work(part_ptr+1:part_ptr+grid%graph%n), &
                  work(work_ptr+1:work_ptr+3*grid%graph%n),control)
              END IF


            CASE (2)
              CALL nd_refine_edge(grid%graph%n,a_ne,grid%graph%ptr, &
                grid%graph%col,grid%row_wgt,sumweight,a_n1_new,a_n2_new, &
                a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
                work(part_ptr+1:part_ptr+grid%graph%n), &
                work(work_ptr+1:work_ptr+3*grid%graph%n),control)

            END SELECT


            CALL cost_function(a_weight_1_new,a_weight_2_new,a_weight_sep_new, &
              sumweight,ratio,imbal,control%cost_function,tau)
            IF (tau<tau_best) THEN
              tau_best = tau
              work(partition_ptr+1:partition_ptr+grid%graph%n) &
                = work(part_ptr+1:part_ptr+grid%graph%n)
              grid%part_div(1) = a_n1_new
              grid%part_div(2) = a_n2_new
              a_weight_1 = a_weight_1_new
              a_weight_2 = a_weight_2_new
              a_weight_sep = a_weight_sep_new
            ELSE
              EXIT
            END IF
          END DO




          ! IF (grid%level .LE.2) THEN
          CALL nd_refine_fm(grid%graph%n,a_ne,grid%graph%ptr, &
            grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1), &
            grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep, &
            work(partition_ptr+1:partition_ptr+grid%graph%n), &
            work(work_ptr+1:work_ptr+8*grid%graph%n+sumweight),control)
          ! END IF

        END IF

        DO i = 1, grid%part_div(1)
          j = work(partition_ptr+i)
          grid%where(j) = ND_PART1_FLAG
        END DO
        DO i = grid%part_div(1) + 1, grid%part_div(1) + grid%part_div(2)
          j = work(partition_ptr+i)
          grid%where(j) = ND_PART2_FLAG
        END DO
        DO i = grid%part_div(1) + grid%part_div(2) + 1, grid%graph%n
          j = work(partition_ptr+i)
          grid%where(j) = ND_SEP_FLAG
        END DO

        IF (info<0) RETURN

        IF (print_level==3) CALL level_print(mp,' after post smoothing ', &
          grid%level)

        ! deallocate the previous level
        ! CALL multigrid_deallocate(cgrid,info)

      END SUBROUTINE multilevel

      ! ***************************************************************
      ! ---------------------------------------------------
      ! nd_partition matrix
      ! ---------------------------------------------------
      ! Partition the matrix and if one (or more) of the generated submatrices
      ! is
      ! small enough, apply halo amd

      SUBROUTINE nd_coarse_partition(a_n,a_ne,a_ptr,a_row,a_weight, &
          sumweight,a_n1,a_n2,where1,lwork,work,control,info)

        INTEGER, INTENT (IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT (IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT (INOUT) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start. This is then
        ! used to hold positions for submatrices after partitioning
        INTEGER, INTENT (INOUT) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.This is then used to hold row indices for
        ! submatrices after partitioning
        INTEGER, INTENT (INOUT) :: a_weight(a_n) ! On input a_weight(i)
        ! contains
        ! the weight of column i. This is then used to hold the weights for
        ! the submatrices after partitioning.
        INTEGER, INTENT (IN) :: sumweight ! Sum entries in a_weight.
        ! Unchanged.
        INTEGER, INTENT (OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT (OUT) :: where1(a_n) ! Computed partition
        INTEGER, INTENT (IN) :: lwork ! .GE. 9*a_n+sumweight
        INTEGER, INTENT (OUT) :: work(lwork)
        TYPE (nd_options), INTENT (IN) :: control
        INTEGER, INTENT (INOUT) :: info
        ! REAL (wp), OPTIONAL, INTENT(OUT) :: real_work(a_n)

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        LOGICAL :: printi, printd, use_multilevel
        INTEGER :: partition_ptr ! pointer into work array
        INTEGER :: work_ptr ! pointer into work array
        INTEGER :: partition_method
        INTEGER :: st
        INTEGER :: a_weight_1, a_weight_2, a_weight_sep, ref_method, &
          ref_control
        INTEGER, ALLOCATABLE :: work1(:)
        REAL (wp) :: dummy, dummy1
        TYPE (nd_multigrid) :: gridtemp

        ! ---------------------------------------------
        ! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)

        IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Start finding a coarse partition'
          WRITE (unit_diagnostics,'(a,i10,a,i10)') 'a_n=', a_n, ', a_ne=', &
            a_ne
        END IF

        ! Find the partition
        IF (control%coarse_partition_method<=1) THEN
          partition_method = 1
        ELSE
          partition_method = 2
        END IF

        partition_ptr = 0 ! length a_n
        work_ptr = partition_ptr + a_n ! max length needed 9*a_n+a_ne

        ALLOCATE (work1(a_n),STAT=st)
        IF (st/=0) THEN
          info = ND_ERR_MEMORY_ALLOC
          RETURN
        END IF


        SELECT CASE (partition_method)
        CASE (1)
          ! Ashcraft method
          use_multilevel = .FALSE.
          CALL nd_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
            a_n2,a_weight_1,a_weight_2,a_weight_sep, &
            work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
            control,info,dummy,dummy1,use_multilevel,gridtemp)

          IF (printi .OR. printd) THEN
            WRITE (unit_diagnostics,'(a)') ' '
            WRITE (unit_diagnostics,'(a)') 'Initial half-level set partition'
            WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
              ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
            WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
              a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
              sumweight - a_weight_1 - a_weight_2
          END IF

        CASE (2)
          ! Level set method
          use_multilevel = .FALSE.
          CALL nd_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
            a_n2,a_weight_1,a_weight_2,a_weight_sep, &
            work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
            control,info,dummy,dummy1,use_multilevel,gridtemp)


          IF (printi .OR. printd) THEN
            WRITE (unit_diagnostics,'(a)') ' '
            WRITE (unit_diagnostics,'(a)') 'Initial level set partition'
            WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
              ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
            WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
              a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
              sumweight - a_weight_1 - a_weight_2
          END IF

        END SELECT



        IF (a_n1/=0 .AND. a_n2/=0 .AND. a_n>=3) THEN
          IF (a_n1+a_n2<a_n) THEN
            ! Refine the partition
            IF (control%refinement>6) THEN
              ref_control = 3
            ELSE
              IF (control%refinement<1) THEN
                ref_control = 1
              ELSE
                ref_control = control%refinement
              END IF
            END IF

            SELECT CASE (ref_control)
            CASE (1)
              ref_method = 1

            CASE (2)
              ref_method = 2

            CASE (3)
              IF (REAL(MAX(a_weight_1,a_weight_2))/REAL(MIN(a_weight_1, &
                  a_weight_2)+a_weight_sep)>MAX(REAL(1.0, &
                  wp),control%balance)) THEN
                ref_method = 2
              ELSE
                ref_method = 1
              END IF

            CASE (4)
              ref_method = 0

            CASE (5)
              ref_method = 2

            CASE (6)
              IF (REAL(MAX(a_weight_1,a_weight_2))/REAL(MIN(a_weight_1, &
                  a_weight_2)+a_weight_sep)>MAX(REAL(1.0, &
                  wp),control%balance)) THEN
                ref_method = 2
              ELSE
                ref_method = 0
              END IF
            END SELECT

            SELECT CASE (ref_method)

            CASE (0)
              CALL nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1, &
                a_n2,a_weight_1,a_weight_2,a_weight_sep, &
                work1(partition_ptr+1:partition_ptr+a_n),work(1:8),control)

            CASE (1)
              IF (MIN(a_weight_1,a_weight_2)+a_weight_sep< &
                  MAX(a_weight_1,a_weight_2)) THEN
                CALL nd_refine_block_trim(a_n,a_ne,a_ptr,a_row, &
                  a_weight,sumweight,a_n1,a_n2,a_weight_1,a_weight_2, &
                  a_weight_sep,work1(partition_ptr+1:partition_ptr+a_n), &
                  work(1:5*a_n),control)
              ELSE
                CALL nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
                  sumweight,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
                  work1(partition_ptr+1:partition_ptr+a_n),work(1:3*a_n), &
                  control)
              END IF

            CASE (2)
              CALL nd_refine_edge(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
                a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
                work1(partition_ptr+1:partition_ptr+a_n),work(1:3*a_n), &
                control)


            END SELECT

            IF (printi .OR. printd) THEN
              WRITE (unit_diagnostics,'(a)') ' '
              WRITE (unit_diagnostics,'(a)') 'Trimmed partition found'
              WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
                ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
              WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
                a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
                sumweight - a_weight_1 - a_weight_2
            END IF

            CALL nd_refine_fm(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1, &
              a_n2,a_weight_1,a_weight_2,a_weight_sep, &
              work1(partition_ptr+1:partition_ptr+a_n), &
              work(1:8*a_n+sumweight),control)

          END IF
        ELSE
          GO TO 10
        END IF


        CALL nd_convert_partition_flags(a_n,a_n1,a_n2, &
          work1(partition_ptr+1:partition_ptr+a_n),ND_PART1_FLAG, &
          ND_PART2_FLAG,ND_SEP_FLAG,where1(1:a_n))

        DEALLOCATE (work1,STAT=st)
        IF (st/=0) THEN
          info = ND_ERR_MEMORY_ALLOC
          RETURN
        END IF



        IF (printi .OR. printd .OR. .FALSE.) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Initial coarse partition found'
          WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
            ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
          WRITE (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
            a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
            sumweight - a_weight_1 - a_weight_2
        END IF
        GO TO 20

10      IF (printi .OR. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'No partition found'
        END IF

20      info = 0
        IF (printi .OR. printd) THEN
          CALL nd_print_message(info,unit_diagnostics, &
            'nd_coarse_partition')
        END IF
        RETURN

      END SUBROUTINE nd_coarse_partition


      ! *****************************************************************

      RECURSIVE SUBROUTINE mg_grid_destroy(grid,info)
        ! deallocate a grid structure
        TYPE (nd_multigrid) :: grid
        INTEGER :: info

        IF (ASSOCIATED(grid%coarse)) THEN

          CALL mg_grid_destroy(grid%coarse,info)

          IF (grid%level/=1) THEN

            CALL multigrid_deallocate(grid,info)

          ELSE

            CALL multigrid_deallocate_first(grid,info)

          END IF

        ELSE

          IF (grid%level/=1) THEN

            CALL multigrid_deallocate_last(grid,info)

          ELSE

            CALL multigrid_deallocate_first(grid,info)

          END IF

        END IF

      END SUBROUTINE mg_grid_destroy


      ! *****************************************************************
      SUBROUTINE multigrid_deallocate(grid,info)
        ! deallocate a grid (at given level between last and first)
        TYPE (nd_multigrid) :: grid
        INTEGER :: st, info

        CALL nd_matrix_destruct(grid%graph,info)
        IF (info/=0) THEN
          RETURN
        END IF

        CALL nd_matrix_destruct(grid%p,info)
        IF (info/=0) THEN
          RETURN
        END IF



        CALL nd_matrix_destruct(grid%r,info)
        IF (info/=0) THEN
          RETURN
        END IF

        DEALLOCATE (grid%graph,grid%p,grid%r,grid%where,grid%row_wgt,STAT=st)
        IF (st/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF
        NULLIFY (grid%graph,grid%coarse)

      END SUBROUTINE multigrid_deallocate

      ! *****************************************************************
      SUBROUTINE multigrid_deallocate_last(grid,info)

        ! deallocate a grid (at the last level). In this case the matrix
        ! grid%graph
        ! has not been formed yet
        TYPE (nd_multigrid) :: grid
        INTEGER, INTENT (INOUT) :: info


        INTEGER :: ierr

        CALL nd_matrix_destruct(grid%p,ierr)
        IF (ierr/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF

        CALL nd_matrix_destruct(grid%r,ierr)
        IF (ierr/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF
        DEALLOCATE (grid%graph,grid%p,grid%r,grid%where,grid%row_wgt, &
          STAT=ierr)
        IF (ierr/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF
        NULLIFY (grid%graph,grid%coarse)

      END SUBROUTINE multigrid_deallocate_last
      ! *****************************************************************
      SUBROUTINE multigrid_deallocate_first(grid,info)
        ! deallocate a grid (at the first level). In this case the matrix
        ! grid%p
        ! does not exist
        TYPE (nd_multigrid) :: grid
        INTEGER, INTENT (INOUT) :: info
        INTEGER :: ierr

        IF (ASSOCIATED(grid%graph)) THEN
          CALL nd_matrix_destruct(grid%graph,ierr)
          IF (ierr/=0) THEN
            info = ND_ERR_MEMORY_DEALLOC
            RETURN
          END IF
        END IF

        ! in subroutine front, grid%graph is not allocated but is pointed to
        ! the
        ! finest level graph. so no need to deallocate
        NULLIFY (grid%graph)
        DEALLOCATE (grid%where,grid%row_wgt,STAT=ierr)
        IF (ierr/=0) info = ND_ERR_MEMORY_DEALLOC

      END SUBROUTINE multigrid_deallocate_first

      ! ***************************************************************
      SUBROUTINE coarsen_hec(grid,lwork,work,info)
        ! coarsen the grid using heavy-edge collapsing and set up the
        ! coarse grid equation, the prolongator and restrictor

        TYPE (nd_multigrid), INTENT (INOUT), TARGET :: grid
        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)
        INTEGER, INTENT (INOUT) :: info


        IF ( .NOT. ASSOCIATED(grid%coarse)) ALLOCATE (grid%coarse)

        grid%coarse%fine => grid

        ! find the prolongator
        CALL prolng_heavy_edge(grid,lwork,work,info)


        grid%coarse%level = grid%level + 1

      END SUBROUTINE coarsen_hec


      ! ***************************************************************
      SUBROUTINE coarsen_cn(grid,lwork,work)
        ! coarsen the grid using common neighbours collapsing and set up the
        ! coarse grid equation, the prolongator and restrictor

        TYPE (nd_multigrid), INTENT (INOUT), TARGET :: grid

        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)

        IF ( .NOT. ASSOCIATED(grid%coarse)) ALLOCATE (grid%coarse)

        grid%coarse%fine => grid

        ! find the prolongator

        CALL prolng_common_neigh(grid,lwork,work)

        grid%coarse%level = grid%level + 1

      END SUBROUTINE coarsen_cn

      ! ***************************************************************
      SUBROUTINE coarsen_best(grid,lwork,work,info)
        ! coarsen the grid using common neighbours collapsing and set up the
        ! coarse grid equation, the prolongator and restrictor

        INTEGER, INTENT (INOUT) :: info
        TYPE (nd_multigrid), INTENT (INOUT), TARGET :: grid

        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)

        IF ( .NOT. ASSOCIATED(grid%coarse)) ALLOCATE (grid%coarse)

        grid%coarse%fine => grid

        ! find the prolongator

        CALL prolng_best(grid,lwork,work,info)

        grid%coarse%level = grid%level + 1

      END SUBROUTINE coarsen_best

      ! ***************************************************************
      SUBROUTINE nd_matrix_multiply_vec(matrix,x,y)
        ! subroutine nd_matrix_multiply_vec(matrix,x,y)

        ! y = matrix*x where x and y are integer vectors. Entries of
        ! matrix is assumed to be one. Dimension of y
        ! is checked and returned if it is smaller than the row dimension
        ! of x    !

        ! matrix: of the derived type nd_matrix, INTENT (IN),
        ! the sparse matrix in compressed sparse row format
        TYPE (nd_matrix), INTENT (IN) :: matrix

        ! x: integer array of intent (IN), a vector to be
        ! multiplied with the matrix
        INTEGER, INTENT (IN), DIMENSION (*) :: x

        ! y: integer array of intent (OUT), the result of
        ! matrix*x or matrix^T*x
        INTEGER, INTENT (OUT), DIMENSION (*) :: y

        ! local ==========
        INTEGER :: m, n, i, l1, l2

        m = matrix%m
        n = matrix%n

        DO i = 1, m
          l1 = matrix%ptr(i)
          l2 = matrix%ptr(i+1) - 1
          y(i) = SUM(x(matrix%col(l1:l2)))
        END DO
      END SUBROUTINE nd_matrix_multiply_vec


      ! ***************************************************************
      SUBROUTINE nd_matrix_destruct(matrix,info,stat)
        ! subroutine nd_matrix_destruct(matrix,info):

        ! destruct the matrix object by deallocating all
        ! space occupied by
        ! matrix. including matrix%ptr, matrix%col and matrix%val.

        ! matrix: is of the derived type nd_matrix,
        ! with INTENT (INOUT). It
        ! the sparse matrix object to be destroyed.
        TYPE (nd_matrix), INTENT (INOUT) :: matrix

        ! info: is an integer scaler of INTENT (OUT).
        ! = 0 if successful
        ! = nd_ERR_MEMORY_DEALLOC if memory deallocation failed
        INTEGER, INTENT (OUT) :: info

        ! stat: is an integer scaler of INTENT (OUT). If supplied,
        ! on exit it holds the error tag for memory allocation
        INTEGER, OPTIONAL, INTENT (OUT) :: stat

        ! ===================== local variables =============
        ! ierr: error tag for deallocation
        INTEGER :: ierr

        info = 0
        IF (PRESENT(stat)) stat = 0
        DEALLOCATE (matrix%col,matrix%ptr,STAT=ierr)
        IF (PRESENT(stat)) stat = ierr
        IF (ierr/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF

        DEALLOCATE (matrix%val,STAT=ierr)
        IF (PRESENT(stat)) stat = ierr
        IF (ierr/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF

      END SUBROUTINE nd_matrix_destruct


      ! ***************************************************************

      SUBROUTINE nd_matrix_construct(p,m,n,ne,info)
        ! Construct data structure for storing sparse matrix
        ! Arrays in nd_matrix will only be (re)allocated if they are not
        ! long
        ! enough. On exit,
        ! size(p%val) <-  max(ne, size(p%val)
        ! size(p%col) <-  max(ne, size(p%col)
        ! size(p%ptr) <-  max(m+1, size(p%ptr)
        TYPE (nd_matrix), INTENT (INOUT) :: p ! matrix being formed using
        ! CSR
        INTEGER, INTENT (IN) :: m ! number of rows
        INTEGER, INTENT (IN) :: n ! number of columns
        INTEGER, INTENT (IN) :: ne ! number entries
        INTEGER, INTENT (OUT) :: info

        info = 0

        p%m = m
        p%n = n
        p%ne = ne

        CALL nd_alloc(p%ptr,m+1,info)
        IF (info<0) THEN
          RETURN
        END IF

        CALL nd_alloc(p%col,ne,info)
        IF (info<0) THEN
          RETURN
        END IF

        CALL nd_alloc(p%val,ne,info)
        IF (info<0) THEN
          RETURN
        END IF

      END SUBROUTINE nd_matrix_construct

      ! ***************************************************************

      SUBROUTINE nd_alloc(v,n,info)
        INTEGER, INTENT (INOUT), ALLOCATABLE :: v(:)
        INTEGER, INTENT (IN) :: n
        INTEGER, INTENT (OUT) :: info

        INTEGER :: st

        info = 0

        IF (ALLOCATED(v)) THEN
          IF (SIZE(v)<n) THEN
            DEALLOCATE (v,STAT=st)
            IF (st<0) THEN
              info = ND_ERR_MEMORY_ALLOC
              RETURN
            END IF
          ELSE
            RETURN
          END IF
        END IF

        ALLOCATE (v(n),STAT=st)
        IF (st<0) THEN
          info = ND_ERR_MEMORY_DEALLOC
        END IF


      END SUBROUTINE nd_alloc


      ! ********************************************************

      SUBROUTINE nd_assoc(arr,sz,info1)
        ! If arr has size at least sz, do nothing. Otherwise, create array arr
        ! of size
        ! sz.
        INTEGER, POINTER, INTENT (INOUT) :: arr(:)
        INTEGER, INTENT (IN) :: sz
        INTEGER, INTENT (INOUT) :: info1

        INTEGER :: st

        info1 = 0

        IF (ASSOCIATED(arr)) THEN
          IF (SIZE(arr)<sz) THEN
            DEALLOCATE (arr,STAT=st)
            IF (st/=0) THEN
              info1 = ND_ERR_MEMORY_DEALLOC
              RETURN
            END IF
          END IF
        END IF

        IF ( .NOT. ASSOCIATED(arr)) THEN
          ALLOCATE (arr(sz),STAT=st)
          IF (st/=0) info1 = ND_ERR_MEMORY_ALLOC
        END IF

      END SUBROUTINE nd_assoc

      ! ********************************************
      SUBROUTINE prolng_heavy_edge(grid,lwork,work,info)

        ! calculate the prolongator for heavy-edge collapsing:
        ! match the vertices of the heaviest edges

        INTEGER, INTENT (INOUT) :: info
        ! input fine grid
        TYPE (nd_multigrid), INTENT (INOUT) :: grid
        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)

        ! coarse grid based on the fine grid
        TYPE (nd_multigrid), POINTER :: cgrid

        ! the fine grid row connectivity graph
        TYPE (nd_matrix), POINTER :: graph

        ! the coarse grid prolongator
        TYPE (nd_matrix), POINTER :: p

        ! the coarse grid restrictor
        TYPE (nd_matrix), POINTER :: r

        ! the number of fine and coarse grid vertices
        INTEGER :: nvtx, cnvtx

        ! working variables
        INTEGER :: v, u, j, i, k
        INTEGER :: nz

        ! whether a vertex is matched already
        INTEGER, PARAMETER :: unmatched = -1

        ! matching status of each vertex
        INTEGER :: ptr_match

        ! maximum weight and index of edges connected to the current vertex
        INTEGER :: maxwgt
        INTEGER :: maxind

        ! allocate the prolongation matrix pointers
        cgrid => grid%coarse
        graph => grid%graph

        ! allocate the graph and matrix pointer and the mincut pointer
        ! so that everything is defined
        IF ( .NOT. ASSOCIATED(cgrid%graph)) ALLOCATE (cgrid%graph)

        nvtx = graph%n

        ! prolongator start here ================================

        ! initialise the matching status and randomly permute the vertex order
        ptr_match = 0

        work(ptr_match+1:ptr_match+nvtx) = unmatched

        ! loop over each vertex and match along the heaviest edge
        cnvtx = 0
        nz = 0
        DO i = 1, nvtx
          v = i
          ! If already matched, next vertex please
          IF (work(ptr_match+v)/=unmatched) CYCLE
          maxwgt = -HUGE(0)
          ! in the case no match is found then match itself
          maxind = v
          ! Loop over entries in row v
          DO j = graph%ptr(v), graph%ptr(v+1) - 1
            ! u is col index of en entry in row v (so u is neighbor of v)
            u = graph%col(j)
            ! heavy edge matching
            ! if u is unmatched and value of the entry in col. u is greater
            ! than maxwgt, select u as the matching.
            IF (work(ptr_match+u)==unmatched .AND. maxwgt<ABS(graph%val(j))) &
                THEN
              maxwgt = ABS(graph%val(j))
              maxind = u
            END IF
          END DO
          ! NOTE: maxind .GE. v
          ! the neighbor with heaviest weight
          work(ptr_match+v) = maxind
          ! mark maxind as having been matched
          work(ptr_match+maxind) = v
          ! increase number of vertices in coarse graph by 1
          cnvtx = cnvtx + 1
          ! construct the prolongation matrix: find vertex v and maxind is
          ! linked
          ! with the coarse grid vertex cnvtx
          nz = nz + 1
          IF (maxind/=v) THEN
            nz = nz + 1
          END IF
        END DO

        ! storage allocation for col. indices and values of prolongation
        ! matrix P (nvtx * cnvtx)
        IF ( .NOT. ASSOCIATED(cgrid%p)) THEN
          ALLOCATE (cgrid%p)
          p => cgrid%p
          CALL nd_matrix_construct(p,nvtx,cnvtx,nz,info)
        ELSE
          p => cgrid%p
          CALL nd_matrix_construct(p,nvtx,cnvtx,nz,info)
        END IF


        ! storage allocation for col. indices and values of restiction
        ! matrix R (cnvtx * nvtx)
        IF ( .NOT. ASSOCIATED(cgrid%r)) THEN
          ALLOCATE (cgrid%r)
          r => cgrid%r
          CALL nd_matrix_construct(r,cnvtx,nvtx,nz,info)
        ELSE
          r => cgrid%r
          CALL nd_matrix_construct(r,cnvtx,nvtx,nz,info)
        END IF

        r%val(1:nz) = 1

        ! store restriction matrix
        r%ptr(cnvtx+1) = nz + 1

        j = 1
        k = 1
        DO i = 1, nvtx
          IF (work(ptr_match+i)==i) THEN
            r%ptr(k) = j
            r%col(j) = i
            j = j + 1
            k = k + 1
          ELSE
            IF (work(ptr_match+i)>i) THEN
              r%ptr(k) = j
              r%col(j) = i
              r%col(j+1) = work(ptr_match+i)
              j = j + 2
              k = k + 1
            END IF
          END IF
        END DO

        ! store prolongation matrix

        p%ptr(1) = 1
        DO i = 1, nvtx
          p%ptr(i+1) = p%ptr(i) + 1
        END DO

        p%val(1:nz) = 1

        j = 1
        DO i = 1, nvtx
          k = work(ptr_match+i)
          IF (k==i) THEN
            p%col(p%ptr(i)) = j
            j = j + 1
          ELSE
            IF (k>i) THEN
              p%col(p%ptr(i)) = j
              p%col(p%ptr(k)) = j
              j = j + 1
            END IF
          END IF
        END DO

        ! size of coarse grid
        cgrid%size = cnvtx

      END SUBROUTINE prolng_heavy_edge

      ! *******************************************************************

      SUBROUTINE prolng_common_neigh(grid,lwork,work)

        ! calculate the prolongator:
        ! match the vertices of with most neighbours in common

        ! input fine grid
        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)
        TYPE (nd_multigrid), INTENT (INOUT) :: grid

        ! coarse grid based on the fine grid
        TYPE (nd_multigrid), POINTER :: cgrid

        ! the fine grid row connectivity graph
        TYPE (nd_matrix), POINTER :: graph

        ! the coarse grid prolongator
        TYPE (nd_matrix), POINTER :: p

        ! the fine grid restrictor
        TYPE (nd_matrix), POINTER :: r

        ! the number of fine and coarse grid vertices
        INTEGER :: nvtx, cnvtx

        ! working variables
        INTEGER :: v, u, w, j, i, k

        ! whether a vertex is matched already
        INTEGER, PARAMETER :: unmatched = -1

        ! matching status of each vertex
        INTEGER :: ptr_match
        ! flag array to flag up neighbours of a  node
        INTEGER :: ptr_flag

        ! maximum no. neighbours and index of edges connected to the current
        ! vertex
        INTEGER :: max_neigh, maxind, num

        INTEGER :: info
        INTEGER :: nz

        ! allocate the prolongation matrix pointers
        cgrid => grid%coarse
        graph => grid%graph

        ! allocate the graph and matrix pointer and the mincut pointer
        ! so that everything is defined
        IF ( .NOT. ASSOCIATED(cgrid%graph)) ALLOCATE (cgrid%graph)

        nvtx = graph%n

        ! prolongator start here ================================

        ! initialise the matching status
        ptr_match = 0
        ptr_flag = ptr_match + nvtx

        work(ptr_match+1:ptr_match+nvtx) = unmatched

        work(ptr_flag+1:ptr_flag+nvtx) = 0

        ! loop over each vertex and match based on number of neighbours in
        ! common
        cnvtx = 0
        nz = 0
        DO i = 1, nvtx
          v = i
          ! If already matched, next vertex please
          IF (work(ptr_match+v)/=unmatched) CYCLE
          ! access the col. indices of row v

          ! in the case no match is found then match itself
          maxind = v
          ! Loop over entries in row v and set flag for each entry
          work(ptr_flag+v) = i
          DO j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
            ! u is col index of en entry in row v (so u is neighbor of v)
            u = grid%graph%col(j)
            work(ptr_flag+u) = i
          END DO
          ! For each unmatched neighbour of v, count the number of
          ! neighbours it has in common with v
          max_neigh = 0
          DO j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
            u = grid%graph%col(j)
            ! cycle is u is already matched
            IF (work(ptr_match+u)/=unmatched) CYCLE
            num = 0
            DO k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
              w = grid%graph%col(k)
              IF (work(ptr_flag+w)==i) num = num + 1
            END DO
            IF (num>max_neigh) THEN
              max_neigh = num
              maxind = u
            END IF
          END DO

          ! the neighbor with largest number of neighbours in common with v
          work(ptr_match+v) = maxind
          ! mark maxind as having been matched
          work(ptr_match+maxind) = v
          ! increase number of vertices in coarse graph by 1
          cnvtx = cnvtx + 1
          ! construct the prolongation matrix: find vertex v and maxind is
          ! linked
          ! with the coarse grid vertex cnvtx
          nz = nz + 1
          IF (maxind/=v) THEN
            nz = nz + 1
          END IF
        END DO

        ! storage allocation for col. indices and values of prolongation
        ! matrix P (order nvtx * cnvtx)
        IF ( .NOT. ASSOCIATED(cgrid%p)) THEN
          ALLOCATE (cgrid%p)
          p => cgrid%p
          CALL nd_matrix_construct(p,nvtx,cnvtx,nz,info)
        ELSE
          p => cgrid%p
          CALL nd_matrix_construct(p,nvtx,cnvtx,nz,info)
        END IF
        p%val(1:nz) = 0

        ! storage allocation for col. indices and values of restiction
        ! matrix R (cnvtx * nvtx)
        IF ( .NOT. ASSOCIATED(cgrid%r)) THEN
          ALLOCATE (cgrid%r)
          r => cgrid%r
          CALL nd_matrix_construct(r,cnvtx,nvtx,nz,info)
        ELSE
          r => cgrid%r
          CALL nd_matrix_construct(r,cnvtx,nvtx,nz,info)
        END IF

        r%val(1:nz) = 1

        ! store restriction matrix
        r%ptr(cnvtx+1) = nz + 1

        j = 1
        k = 1
        DO i = 1, nvtx
          IF (work(ptr_match+i)==i) THEN
            r%ptr(k) = j
            r%col(j) = i
            j = j + 1
            k = k + 1
          ELSE
            IF (work(ptr_match+i)>i) THEN
              r%ptr(k) = j
              r%col(j) = i
              r%col(j+1) = work(ptr_match+i)
              j = j + 2
              k = k + 1
            END IF
          END IF
        END DO


        ! store prolongation matrix

        p%ptr(1) = 1
        DO i = 1, nvtx
          p%ptr(i+1) = p%ptr(i) + 1
        END DO

        p%val(1:nz) = 1

        j = 1
        DO i = 1, nvtx
          k = work(ptr_match+i)
          IF (k==i) THEN
            p%col(p%ptr(i)) = j
            j = j + 1
          ELSE
            IF (k>i) THEN
              p%col(p%ptr(i)) = j
              p%col(p%ptr(k)) = j
              j = j + 1
            END IF
          END IF
        END DO


        ! size of coarse grid
        cgrid%size = cnvtx

      END SUBROUTINE prolng_common_neigh

      ! ********************************************
      SUBROUTINE prolng_best(grid,lwork,work,info)

        ! calculate the prolongator for heavy-edge collapsing:
        ! match the vertices of the heaviest edges

        INTEGER, INTENT (INOUT) :: info
        ! input fine grid
        TYPE (nd_multigrid), INTENT (INOUT) :: grid
        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT), TARGET :: work(lwork)

        ! coarse grid based on the fine grid
        TYPE (nd_multigrid), POINTER :: cgrid

        ! the fine grid row connectivity graph
        TYPE (nd_matrix), POINTER :: graph

        ! the coarse grid prolongator
        TYPE (nd_matrix), POINTER :: p

        ! the coarse grid restrictor
        TYPE (nd_matrix), POINTER :: r

        ! the number of fine and coarse grid vertices
        INTEGER :: nvtx, cnvtx, cnvtx1

        ! working variables
        INTEGER :: v, u, j, i, k
        INTEGER :: nz

        ! whether a vertex is matched already
        INTEGER, PARAMETER :: unmatched = -1

        ! matching status of each vertex
        INTEGER :: ptr_match, ptr_match1, ptr_flag, max_neigh, num, w

        ! maximum weight and index of edges connected to the current vertex
        INTEGER :: maxwgt
        INTEGER :: maxind

        INTEGER, POINTER, DIMENSION (:) :: matching

        ! allocate the prolongation matrix pointers
        cgrid => grid%coarse
        graph => grid%graph

        ! allocate the graph and matrix pointer and the mincut pointer
        ! so that everything is defined
        IF ( .NOT. ASSOCIATED(cgrid%graph)) ALLOCATE (cgrid%graph)

        nvtx = graph%n

        ptr_match = 0
        ptr_match1 = ptr_match + nvtx

        ! -----------------------------------------------------------------
        ! Find heavy-edge matching

        ! initialise the matching status and randomly permute the vertex order

        work(ptr_match+1:ptr_match+nvtx) = unmatched

        ! loop over each vertex and match along the heaviest edge
        cnvtx = 0
        DO i = 1, nvtx
          v = i
          ! If already matched, next vertex please
          IF (work(ptr_match+v)/=unmatched) CYCLE
          maxwgt = -HUGE(0)
          ! in the case no match is found then match itself
          maxind = v
          ! Loop over entries in row v
          DO j = graph%ptr(v), graph%ptr(v+1) - 1
            ! u is col index of en entry in row v (so u is neighbor of v)
            u = graph%col(j)
            ! heavy edge matching
            ! if u is unmatched and value of the entry in col. u is greater
            ! than maxwgt, select u as the matching.
            IF (work(ptr_match+u)==unmatched .AND. maxwgt<ABS(graph%val(j))) &
                THEN
              maxwgt = ABS(graph%val(j))
              maxind = u
            END IF
          END DO
          ! NOTE: maxind .GE. v
          ! the neighbor with heaviest weight
          work(ptr_match+v) = maxind
          ! mark maxind as having been matched
          work(ptr_match+maxind) = v
          ! increase number of vertices in coarse graph by 1
          cnvtx = cnvtx + 1
          ! construct the prolongation matrix: find vertex v and maxind is
          ! linked
          ! with the coarse grid vertex cnvtx
        END DO
        nz = nvtx

        ! -----------------------------------------------------------------
        ! Find common neighbours matching

        ! initialise the matching status
        ptr_match1 = ptr_match + nvtx
        ptr_flag = ptr_match1 + nvtx

        work(ptr_match1+1:ptr_match1+nvtx) = unmatched

        work(ptr_flag+1:ptr_flag+nvtx) = 0

        ! loop over each vertex and match based on number of neighbours in
        ! common
        cnvtx1 = 0
        nz = 0
        DO i = 1, nvtx
          v = i
          ! If already matched, next vertex please
          IF (work(ptr_match1+v)/=unmatched) CYCLE
          ! access the col. indices of row v

          ! in the case no match is found then match itself
          maxind = v
          ! Loop over entries in row v and set flag for each entry
          work(ptr_flag+v) = i
          DO j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
            ! u is col index of en entry in row v (so u is neighbor of v)
            u = grid%graph%col(j)
            work(ptr_flag+u) = i
          END DO
          ! For each unmatched neighbour of v, count the number of
          ! neighbours it has in common with v
          max_neigh = 0
          DO j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
            u = grid%graph%col(j)
            ! cycle is u is already matched
            IF (work(ptr_match1+u)/=unmatched) CYCLE
            num = 0
            DO k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
              w = grid%graph%col(k)
              IF (work(ptr_flag+w)==i) num = num + 1
            END DO
            IF (num>max_neigh) THEN
              max_neigh = num
              maxind = u
            END IF
          END DO

          ! the neighbor with largest number of neighbours in common with v
          work(ptr_match1+v) = maxind
          ! mark maxind as having been matched
          work(ptr_match1+maxind) = v
          ! increase number of vertices in coarse graph by 1
          cnvtx1 = cnvtx1 + 1
          ! construct the prolongation matrix: find vertex v and maxind is
          ! linked
          ! with the coarse grid vertex cnvtx1
          nz = nz + 1
          IF (maxind/=v) THEN
            nz = nz + 1
          END IF
        END DO

        ! --------------------------------------------------------------------
        ! -
        IF (cnvtx<=cnvtx1) THEN
          ! use heavy-edge matching
          matching => work(ptr_match+1:ptr_match+nvtx)
        ELSE
          ! use common neighbours matching
          matching => work(ptr_match1+1:ptr_match1+nvtx)
          cnvtx = cnvtx1
        END IF


        ! storage allocation for col. indices and values of prolongation
        ! matrix P (nvtx * cnvtx)
        IF ( .NOT. ASSOCIATED(cgrid%p)) THEN
          ALLOCATE (cgrid%p)
          p => cgrid%p
          CALL nd_matrix_construct(p,nvtx,cnvtx,nz,info)
        ELSE
          p => cgrid%p
          CALL nd_matrix_construct(p,nvtx,cnvtx,nz,info)
        END IF


        ! storage allocation for col. indices and values of restiction
        ! matrix R (cnvtx * nvtx)
        IF ( .NOT. ASSOCIATED(cgrid%r)) THEN
          ALLOCATE (cgrid%r)
          r => cgrid%r
          CALL nd_matrix_construct(r,cnvtx,nvtx,nz,info)
        ELSE
          r => cgrid%r
          CALL nd_matrix_construct(r,cnvtx,nvtx,nz,info)
        END IF

        r%val(1:nz) = 1

        ! store restriction matrix
        r%ptr(cnvtx+1) = nz + 1

        j = 1
        k = 1
        DO i = 1, nvtx
          IF (matching(i)==i) THEN
            r%ptr(k) = j
            r%col(j) = i
            j = j + 1
            k = k + 1
          ELSE
            IF (matching(i)>i) THEN
              r%ptr(k) = j
              r%col(j) = i
              r%col(j+1) = matching(i)
              j = j + 2
              k = k + 1
            END IF
          END IF
        END DO

        ! store prolongation matrix

        p%ptr(1) = 1
        DO i = 1, nvtx
          p%ptr(i+1) = p%ptr(i) + 1
        END DO

        p%val(1:nz) = 1

        j = 1
        DO i = 1, nvtx
          k = matching(i)
          IF (k==i) THEN
            p%col(p%ptr(i)) = j
            j = j + 1
          ELSE
            IF (k>i) THEN
              p%col(p%ptr(i)) = j
              p%col(p%ptr(k)) = j
              j = j + 1
            END IF
          END IF
        END DO

        ! size of coarse grid
        cgrid%size = cnvtx
        NULLIFY (matching)

      END SUBROUTINE prolng_best


      ! *******************************************************************
      SUBROUTINE level_print(mp,title1,level,title2,res)

        CHARACTER (len=*), INTENT (IN) :: title1
        INTEGER, INTENT (IN) :: mp, level
        REAL (wp), OPTIONAL, INTENT (IN) :: res
        CHARACTER (len=*), OPTIONAL, INTENT (IN) :: title2
        INTEGER :: char_len1, char_len2

        char_len1 = LEN_TRIM(title1)

        IF (PRESENT(res) .AND. PRESENT(title2)) THEN
          char_len2 = LEN_TRIM(title2)
          WRITE (mp,'(a,i4,a,g14.3)') title1, level, title2, res
        ELSE
          WRITE (mp,'(a,i4)') title1, level
        END IF

      END SUBROUTINE level_print



      ! *************************************************
      SUBROUTINE galerkin_graph(matrix,p,r,cmatrix,info,lwork,work)

        ! Given matrix on fine grid and a prolongation operator p,
        ! find the coarse matrix R*A*P

        ! matrix: fine grid matrix
        TYPE (nd_matrix), INTENT (IN) :: matrix
        ! p: prolongation operator
        TYPE (nd_matrix), INTENT (IN) :: p
        ! r: restriction operator
        TYPE (nd_matrix), INTENT (IN) :: r
        ! cmatrix: coarse grid matrix
        TYPE (nd_matrix), INTENT (INOUT) :: cmatrix
        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)

        ! nvtx,cnvtx: size of fine and coarse grid
        INTEGER :: nvtx, cnvtx
        INTEGER :: nz

        INTEGER, INTENT (INOUT) :: info

        ! CALL mc65_matrix_transpose(p,r,info65)
        ! IF (info65<0) THEN
        ! info = info65
        ! RETURN
        ! END IF
        nvtx = matrix%n
        cnvtx = p%n

        ! get the size of the coarse matrix first
        CALL galerkin_graph_rap_size(nvtx,cnvtx,nz,p%ptr(nvtx+1)-1,p%col, &
          p%ptr,matrix%ptr(nvtx+1)-1,matrix%col,matrix%ptr,r%ptr(cnvtx+1)-1, &
          r%col,r%ptr,lwork,work(1:lwork))
        IF (info<0) RETURN

        CALL nd_matrix_construct(cmatrix,cnvtx,cnvtx,nz,info)
        IF (info<0) THEN
          RETURN
        END IF

        CALL galerkin_graph_rap(nvtx,cnvtx,p%ptr(nvtx+1)-1,p%val,p%col,p%ptr, &
          matrix%ptr(nvtx+1)-1,matrix%val,matrix%col,matrix%ptr, &
          r%ptr(cnvtx+1)-1,r%val,r%col,r%ptr,nz,cmatrix%val,cmatrix%col, &
          cmatrix%ptr,lwork,work(1:lwork))
        IF (info<0) RETURN

      END SUBROUTINE galerkin_graph

      ! *************************************************

      SUBROUTINE galerkin_graph_rap_size(nvtx,cnvtx,nz,nzp,pcol,pptr,nzaa, &
          acol,aptr,nzr,rcol,rptr,lwork,work)
        ! get the number of nonzeros in R*A*P
        ! nvtx: size of aa matrix
        ! cnvtx: size of ca matrix
        INTEGER, INTENT (IN) :: nvtx, cnvtx
        ! nz: number of nonzeros in R*A*P
        INTEGER, INTENT (OUT) :: nz

        ! P: matrix
        INTEGER, INTENT (IN) :: nzp
        INTEGER, INTENT (IN), DIMENSION (nzp) :: pcol
        INTEGER, INTENT (IN), DIMENSION (nvtx+1) :: pptr
        ! aa: matrix
        INTEGER, INTENT (IN) :: nzaa
        INTEGER, INTENT (IN), DIMENSION (nzaa) :: acol
        INTEGER, INTENT (IN), DIMENSION (nvtx+1) :: aptr
        ! R: matrix
        INTEGER, INTENT (IN) :: nzr
        INTEGER, INTENT (IN), DIMENSION (nzr) :: rcol
        INTEGER, INTENT (IN), DIMENSION (cnvtx+1) :: rptr

        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)

        ! mask: masking array to see if an entry has been seen before
        INTEGER :: ptr_mask
        ! i,j,k: loop index
        INTEGER :: i, j, k
        ! nz: number of nonzeros so far in ca
        INTEGER :: nz1
        ! various neighbors
        INTEGER :: neigh, neighneigh

        ! col: column index of a row of r*matrix
        INTEGER :: ptr_col

        ptr_mask = 0
        ptr_col = ptr_mask + nvtx
        work(ptr_mask+1:ptr_mask+nvtx) = 0
        nz = 0
        ! loop over coarse grid points
        DO i = 1, cnvtx
          ! first form row i of (r*matrix)
          nz1 = 0
          ! for each vertex D that restricts to C (including itself).
          DO j = rptr(i), rptr(i+1) - 1
            neigh = rcol(j)
            ! find D's neighbor
            DO k = aptr(neigh), aptr(neigh+1) - 1
              neighneigh = acol(k)
              IF (work(ptr_mask+neighneigh)/=i) THEN
                nz1 = nz1 + 1
                work(ptr_col+nz1) = neighneigh
                work(ptr_mask+neighneigh) = i
              END IF
            END DO
          END DO
          ! form row i of (r*matrix)*p
          DO j = 1, nz1
            neigh = work(ptr_col+j)
            DO k = pptr(neigh), pptr(neigh+1) - 1
              neighneigh = pcol(k)
              IF (work(ptr_mask+neighneigh)/=-i .AND. neighneigh/=i) THEN
                nz = nz + 1
                work(ptr_mask+neighneigh) = -i
              END IF
            END DO
          END DO
        END DO

      END SUBROUTINE galerkin_graph_rap_size
      ! ******************************************************
      SUBROUTINE galerkin_graph_rap(nvtx,cnvtx,nzp,pa,pcol,pptr,nzaa,aa,acol, &
          aptr,nzr,ra,rcol,rptr,nzca,ca,ccol,cptr,lwork,work)
        ! multiply R*A*P to get CA
        ! nvtx: size of aa matrix
        ! cnvtx: size of ca matrix
        INTEGER, INTENT (IN) :: nvtx, cnvtx
        ! p: matrix
        INTEGER, INTENT (IN) :: nzp
        INTEGER, INTENT (IN), DIMENSION (nzp) :: pa
        INTEGER, INTENT (IN), DIMENSION (nzp) :: pcol
        INTEGER, INTENT (IN), DIMENSION (nvtx+1) :: pptr
        ! aa: matrix
        INTEGER, INTENT (IN) :: nzaa
        INTEGER, INTENT (IN), DIMENSION (:) :: aa
        INTEGER, INTENT (IN), DIMENSION (nzaa) :: acol
        INTEGER, INTENT (IN), DIMENSION (:) :: aptr
        ! r: matrix
        INTEGER, INTENT (IN) :: nzr
        INTEGER, INTENT (IN), DIMENSION (nzr) :: ra
        INTEGER, INTENT (IN), DIMENSION (nzr) :: rcol
        INTEGER, INTENT (IN), DIMENSION (cnvtx+1) :: rptr
        ! ca: matrix
        INTEGER, INTENT (IN) :: nzca
        INTEGER, INTENT (INOUT), DIMENSION (nzca) :: ca
        INTEGER, INTENT (INOUT), DIMENSION (nzca) :: ccol
        INTEGER, INTENT (INOUT), DIMENSION (cnvtx+1) :: cptr

        INTEGER, INTENT (IN) :: lwork
        INTEGER, INTENT (OUT) :: work(lwork)


        ! mask: masking array to see if an entry has been seen before
        INTEGER :: ptr_mask
        ! i,j,k,l: loop index
        INTEGER :: i, j, k
        ! nz: number of nonzeros so far in ca
        INTEGER :: nz, nzz, nz1
        ! various neighbors
        INTEGER :: neigh, neighneigh
        ! r_ij: (i,j) element of r
        INTEGER :: r_ij
        ! col: column index of a row of r*matrix
        ! a: values of a row of r*matrix
        INTEGER :: ptr_col, ptr_a

        ptr_mask = 0
        ptr_col = ptr_mask + nvtx
        ptr_a = ptr_col + nvtx
        ! now get the entries of the coarse matrix
        cptr(1) = 1
        work(ptr_mask+1:ptr_mask+nvtx) = 0
        nz = 0
        ! loop over every coarse grid point
        DO i = 1, cnvtx
          ! first form row i of (r*matrix)
          nz1 = 0
          ! foreach each vertex D that restricts to C (including itself).
          DO j = rptr(i), rptr(i+1) - 1
            neigh = rcol(j)
            r_ij = ra(j)
            ! find D's neighbor
            DO k = aptr(neigh), aptr(neigh+1) - 1
              neighneigh = acol(k)
              nzz = work(ptr_mask+neighneigh)
              IF (nzz==0) THEN
                nz1 = nz1 + 1
                work(ptr_col+nz1) = neighneigh
                work(ptr_a+nz1) = r_ij*aa(k)
                work(ptr_mask+neighneigh) = nz1
              ELSE
                work(ptr_a+nzz) = work(ptr_a+nzz) + r_ij*aa(k)
              END IF
            END DO
          END DO
          DO j = 1, nz1
            work(ptr_mask+work(ptr_col+j)) = 0
          END DO

          ! form row i of (r*matrix)*p
          DO j = 1, nz1
            neigh = work(ptr_col+j)
            r_ij = work(ptr_a+j)
            DO k = pptr(neigh), pptr(neigh+1) - 1
              neighneigh = pcol(k)
              IF (neighneigh==i) CYCLE
              nzz = work(ptr_mask+neighneigh)
              IF (nzz==0) THEN
                nz = nz + 1
                work(ptr_mask+neighneigh) = nz
                ca(nz) = r_ij*pa(k)
                ccol(nz) = neighneigh
              ELSE
                ca(nzz) = ca(nzz) + r_ij*pa(k)
              END IF
            END DO
          END DO

          DO j = cptr(i), nz
            work(ptr_mask+ccol(j)) = 0
          END DO
          cptr(i+1) = nz + 1
        END DO


      END SUBROUTINE galerkin_graph_rap

      ! ---------------------------------------------------
      ! nd_refine_trim
      ! ---------------------------------------------------
      ! Given a partition, trim the partition to make it minimal
      SUBROUTINE nd_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight, &
          a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end. This is updated to the new
        ! partition
        INTEGER, INTENT (OUT) :: work(3*a_n) ! Work array
        TYPE (nd_options), INTENT (IN) :: control

        ! ---------------------------------------------
        ! Local variables
        ! ---------------------------------------------
        INTEGER :: work_part, work_next, work_prev, a_n1_orig, a_n2_orig
        INTEGER :: head1, head2, tail1, tail2
        INTEGER, PARAMETER :: sep1 = -1
        INTEGER, PARAMETER :: sep2 = -2
        INTEGER, PARAMETER :: sep3 = -3
        INTEGER :: i, j, k, l, m, p, q, w1, w2
        LOGICAL :: next1, next2, imbal
        REAL (wp) :: t1, t2
        REAL (wp) :: ratio

        ratio = MAX(REAL(1.0,wp),control%balance)
        IF (ratio>REAL(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF

        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what
        ! part of the partition the nodes are in
        work_part = 0
        CALL nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
          ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
          work(work_part+1:work_part+a_n))

        a_n1_orig = a_n1
        a_n2_orig = a_n2


        ! Create two lists
        work_next = work_part + a_n
        work_prev = work_next + a_n
        head1 = 0
        head2 = 0
        tail1 = 0
        tail2 = 0
        DO i = a_n1 + a_n2 + 1, a_n
          next1 = .FALSE.
          next2 = .FALSE.
          j = partition(i)
          IF (j<a_n) THEN
            k = a_ptr(j+1) - 1
          ELSE
            k = a_ne
          END IF
          DO l = a_ptr(j), k
            m = a_row(l)
            IF (work(work_part+m)==ND_PART1_FLAG) THEN
              next1 = .TRUE.
            ELSE IF (work(work_part+m)==ND_PART2_FLAG) THEN
              next2 = .TRUE.
            END IF
          END DO
          IF ((next1 .AND. .NOT. next2) .OR. ( .NOT. next1 .AND. .NOT. next2 &
              .AND. a_n1==0)) THEN
            ! Add to list 1
            IF (head1==0) THEN
              head1 = j
              work(work_next+j) = 0
              work(work_prev+j) = 0
              tail1 = j
            ELSE
              work(work_next+tail1) = j
              work(work_prev+j) = tail1
              work(work_next+j) = 0
              tail1 = j
            END IF
            work(work_part+j) = sep1
          ELSE IF ((next2 .AND. .NOT. next1) .OR. ( .NOT. next1 .AND. .NOT. &
              next2 .AND. a_n2==0)) THEN
            ! Add to list 2
            IF (head2==0) THEN
              head2 = j
              work(work_next+j) = 0
              work(work_prev+j) = 0
              tail2 = j
            ELSE
              work(work_next+tail2) = j
              work(work_prev+j) = tail2
              work(work_next+j) = 0
              tail2 = j
            END IF
            work(work_part+j) = sep2
          ELSE IF (next1 .AND. next2) THEN
            work(work_part+j) = sep3
          ELSE
            CONTINUE
          END IF
        END DO

        DO WHILE (head1>0 .OR. head2>0)
          IF (head1>0 .AND. head2>0) THEN
            w1 = a_weight(head1)
            w2 = a_weight(head2)
            CALL cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1, &
              sumweight,ratio,imbal,control%cost_function,t1)
            CALL cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2, &
              sumweight,ratio,imbal,control%cost_function,t2)

            IF (t1<t2) THEN
              GO TO 10
            ELSE
              GO TO 20
            END IF

          ELSE IF (head1>0) THEN
            GO TO 10
          ELSE
            GO TO 20
          END IF


          ! move entry from separator to partition1
10        i = head1
          work(work_part+i) = ND_PART1_FLAG
          head1 = work(work_next+i)
          work(work_next+i) = 0
          a_n1 = a_n1 + 1
          a_weight_1 = a_weight_1 + a_weight(i)
          a_weight_sep = a_weight_sep - a_weight(i)
          ! update list
          IF (i<a_n) THEN
            k = a_ptr(i+1) - 1
          ELSE
            k = a_ne
          END IF
          DO l = a_ptr(i), k
            j = a_row(l)
            m = work(work_part+j)
            SELECT CASE (m)
            CASE (ND_SEP_FLAG)
              ! Add to list 1
              work(work_next+tail1) = j
              work(work_prev+j) = tail1
              work(work_next+j) = 0
              tail1 = j
              work(work_part+j) = sep1
              IF (head1==0) THEN
                head1 = j
              END IF

            CASE (sep2)
              ! Remove from list 2
              p = work(work_prev+j)
              q = work(work_next+j)

              IF (j/=head2 .AND. j/=tail2) THEN
                work(work_prev+q) = p
                work(work_next+p) = q
                work(work_prev+j) = 0
                work(work_next+j) = 0
              ELSE IF (j/=head2 .AND. j==tail2) THEN
                work(work_next+p) = 0
                work(work_prev+j) = 0
                tail2 = p
              ELSE IF (j/=tail2 .AND. j==head2) THEN
                work(work_prev+q) = p
                work(work_next+j) = 0
                head2 = q
              ELSE
                head2 = 0
                tail2 = 0

              END IF
              work(work_part+j) = sep3

            END SELECT
          END DO
          GO TO 30

          ! move entry from separator to partition 2
20        i = head2
          work(work_part+i) = ND_PART2_FLAG
          head2 = work(work_next+i)
          work(work_next+i) = 0
          a_n2 = a_n2 + 1
          a_weight_2 = a_weight_2 + a_weight(i)
          a_weight_sep = a_weight_sep - a_weight(i)
          ! update list
          IF (i<a_n) THEN
            k = a_ptr(i+1) - 1
          ELSE
            k = a_ne
          END IF
          DO l = a_ptr(i), k
            j = a_row(l)
            m = work(work_part+j)
            SELECT CASE (m)
            CASE (ND_SEP_FLAG)
              ! Add to list 2
              work(work_next+tail2) = j
              work(work_prev+j) = tail2
              work(work_next+j) = 0
              tail2 = j
              work(work_part+j) = sep2
              IF (head2==0) THEN
                head2 = j
              END IF

            CASE (sep1)
              ! Remove from list 1
              p = work(work_prev+j)
              q = work(work_next+j)

              IF (j/=head1 .AND. j/=tail1) THEN
                work(work_prev+q) = p
                work(work_next+p) = q
                work(work_prev+j) = 0
                work(work_next+j) = 0
              ELSE IF (j/=head1 .AND. j==tail1) THEN
                work(work_next+p) = 0
                work(work_prev+j) = 0
                tail1 = p
              ELSE IF (j/=tail1 .AND. j==head1) THEN
                work(work_prev+q) = p
                work(work_next+j) = 0
                head1 = q
              ELSE
                head1 = 0
                tail1 = 0
              END IF
              work(work_part+j) = sep3
            END SELECT
          END DO

30        CONTINUE

        END DO

        ! Check for any entries in separator that are still inside boundary
        ! and
        ! move into a partition
        work(work_next+a_n1_orig+a_n2_orig+1:work_next+a_n) = 0
        DO i = a_n1_orig + a_n2_orig + 1, a_n
          j = partition(i)
          IF (work(work_part+j)==ND_SEP_FLAG) THEN
            ! j is not on the boundary
            IF (a_weight_1<a_weight_2) THEN
              ! Move j into partition 1
              work(work_part+j) = ND_PART1_FLAG
              a_n1 = a_n1 + 1
              a_weight_1 = a_weight_1 + a_weight(j)
              a_weight_sep = a_weight_sep - a_weight(j)

              head1 = j
              tail1 = j
              DO WHILE (head1>0)
                q = head1
                IF (q<a_n) THEN
                  k = a_ptr(q+1) - 1
                ELSE
                  k = a_ne
                END IF
                DO l = a_ptr(q), k
                  p = a_row(l)
                  IF (work(work_part+p)==ND_SEP_FLAG) THEN
                    work(work_part+p) = ND_PART1_FLAG
                    a_n1 = a_n1 + 1
                    a_weight_1 = a_weight_1 + a_weight(p)
                    a_weight_sep = a_weight_sep - a_weight(p)
                    work(work_next+tail1) = p
                    tail1 = p
                  END IF
                END DO
                IF (head1==tail1) THEN
                  head1 = 0
                  tail1 = 0
                ELSE
                  head1 = work(work_next+q)
                  work(work_next+q) = 0
                END IF

              END DO

            ELSE
              ! Move j into partition 2
              work(work_part+j) = ND_PART2_FLAG
              a_n2 = a_n2 + 1
              a_weight_2 = a_weight_2 + a_weight(j)
              a_weight_sep = a_weight_sep - a_weight(j)
              head2 = j
              tail2 = j

              DO WHILE (head2>0)
                q = head2
                IF (q<a_n) THEN
                  k = a_ptr(q+1) - 1
                ELSE
                  k = a_ne
                END IF
                DO l = a_ptr(q), k
                  p = a_row(l)
                  IF (work(work_part+p)==ND_SEP_FLAG) THEN
                    work(work_part+p) = ND_PART2_FLAG
                    a_n2 = a_n2 + 1
                    a_weight_2 = a_weight_2 + a_weight(p)
                    a_weight_sep = a_weight_sep - a_weight(p)
                    work(work_next+tail2) = p
                    tail2 = p
                  END IF
                END DO
                IF (head2==tail2) THEN
                  head2 = 0
                  tail2 = 0
                ELSE
                  head2 = work(work_next+q)
                  work(work_next+q) = 0
                END IF
              END DO
            END IF
          END IF
        END DO

        a_weight_sep = sumweight - a_weight_1 - a_weight_2

        ! Reset partition matrix
        CALL nd_convert_flags_partition(a_n,a_n1,a_n2, &
          work(work_part+1:work_part+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
          partition(1:a_n))

        ! call
        ! check_partition1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,work(work_part+1:wor
        ! k_part+a_n),a_weight_1,a_weight_2,a_weight)

      END SUBROUTINE nd_refine_trim


      ! ---------------------------------------------------
      ! nd_refine_block_trim
      ! ---------------------------------------------------
      ! Given a partition, trim the partition using blocks to make it minimal
      SUBROUTINE nd_refine_block_trim(a_n,a_ne,a_ptr,a_row,a_weight, &
          sumweight,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,partition, &
          work,control)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end. This is updated to the new
        ! partition
        INTEGER, INTENT (OUT) :: work(5*a_n) ! Work array
        TYPE (nd_options), INTENT (IN) :: control

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: work_part, work_next1, work_next2, work_level1, work_level2
        INTEGER :: a_n1_orig, a_n2_orig
        INTEGER :: head1, head2, tail1, tail2, maxlevel1, maxlevel2
        INTEGER :: currlevel1, currlevel2
        INTEGER :: i, j, k, l, m, w1, w2, l1, l2
        LOGICAL :: next1, next2, imbal
        REAL (wp) :: t1, t2
        REAL (wp) :: ratio

        ratio = MAX(REAL(1.0,wp),control%balance)
        IF (ratio>REAL(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF

        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what
        ! part of the partition the nodes are in
        work_part = 0
        CALL nd_convert_partition_flags(a_n,a_n1,a_n2,partition, &
          ND_PART1_FLAG,ND_PART2_FLAG,ND_SEP_FLAG, &
          work(work_part+1:work_part+a_n))
        a_n1_orig = a_n1
        a_n2_orig = a_n2
        a_weight_sep = sumweight - a_weight_1 - a_weight_2


        ! Set work(work_next1+1:work_next1+a_n) to hold list pointers
        ! Set work(work_next2+1:work_next2+a_n) to hold list pointers
        ! Set work(work_level1+1:work_level1+a_n) to hold distance from orig
        ! partition 1
        ! Set work(work_level2+1:work_level2+a_n) to hold distance from orig
        ! partition 2
        work_next1 = work_part + a_n
        work_next2 = work_next1 + a_n
        work_level1 = work_next2 + a_n
        work_level2 = work_level1 + a_n
        work(work_next1+1:work_next1+a_n) = 0
        work(work_next2+1:work_next2+a_n) = 0
        work(work_level1+1:work_level1+a_n) = 0
        work(work_level2+1:work_level2+a_n) = 0

        ! Create two lists
        head1 = 0
        head2 = 0
        DO i = a_n1 + a_n2 + 1, a_n
          next1 = .FALSE.
          next2 = .FALSE.
          j = partition(i)
          IF (j<a_n) THEN
            k = a_ptr(j+1) - 1
          ELSE
            k = a_ne
          END IF
          DO l = a_ptr(j), k
            m = a_row(l)
            IF (work(work_part+m)==ND_PART1_FLAG) THEN
              next1 = .TRUE.
            ELSE IF (work(work_part+m)==ND_PART2_FLAG) THEN
              next2 = .TRUE.
            END IF
          END DO
          IF (next1) THEN
            ! Add to list 1
            IF (head1==0) THEN
              head1 = j
            ELSE
              work(work_next1+tail1) = j
            END IF
            tail1 = j
            work(work_level1+j) = 1
          END IF
          IF (next2) THEN
            ! Add to list 2
            IF (head2==0) THEN
              head2 = j
            ELSE
              work(work_next2+tail2) = j
            END IF
            tail2 = j
            work(work_level2+j) = 1
          END IF
        END DO

        ! Breadth first search of separator from entries adjacent to partition
        ! 1
        l1 = head1
        DO WHILE (l1>0)
          IF (l1<a_n) THEN
            k = a_ptr(l1+1) - 1
          ELSE
            k = a_ne
          END IF
          DO l = a_ptr(l1), k
            m = a_row(l)
            IF (work(work_part+m)==ND_SEP_FLAG .AND. &
                work(work_level1+m)==0) THEN
              ! Add to list (note list is non-empty)
              work(work_next1+tail1) = m
              tail1 = m
              work(work_level1+m) = work(work_level1+l1) + 1
            END IF
          END DO
          l1 = work(work_next1+l1)
        END DO
        maxlevel1 = work(work_level1+tail1)

        ! Breadth first search of separator from entries adjacent to partition
        ! 2
        l1 = head2
        DO WHILE (l1>0)
          IF (l1<a_n) THEN
            k = a_ptr(l1+1) - 1
          ELSE
            k = a_ne
          END IF
          DO l = a_ptr(l1), k
            m = a_row(l)
            IF (work(work_part+m)==ND_SEP_FLAG .AND. &
                work(work_level2+m)==0) THEN
              ! Add to list (note list is non-empty)
              work(work_next2+tail2) = m
              tail2 = m
              work(work_level2+m) = work(work_level2+l1) + 1
            END IF
          END DO
          l1 = work(work_next2+l1)
        END DO
        maxlevel2 = work(work_level2+tail2)

        ! Check for any entries in separator only reachable from one partition
        DO i = a_n1 + a_n2 + 1, a_n
          j = partition(i)
          IF (work(work_level2+j)==0) THEN
            work(work_level2+j) = maxlevel2 + 1
          ELSE
            IF (work(work_level1+j)==0) THEN
              work(work_level1+j) = maxlevel1 + 1
            END IF
          END IF
        END DO

        ! Trim the separator
        currlevel1 = 1
        currlevel2 = 1
        l1 = head1
        l2 = head2
        DO WHILE (currlevel1<=maxlevel1 .OR. currlevel2<=maxlevel2)
          IF (currlevel1>maxlevel1) THEN
            t1 = HUGE(1.0_wp)
          ELSE
            w1 = 0
            j = l1
            DO WHILE (work(work_level1+j)==currlevel1)
              IF (work(work_level2+j)>currlevel2) THEN
                w1 = w1 + a_weight(j)
              END IF
              j = work(work_next1+j)
              IF (j==0) EXIT
            END DO
            IF (w1==0) THEN
              currlevel1 = currlevel1 + 1
              l1 = j
              CYCLE
            ELSE
              CALL cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1, &
                sumweight,ratio,imbal,control%cost_function,t1)
            END IF
          END IF

          IF (currlevel2>maxlevel2) THEN
            t2 = HUGE(1.0_wp)
          ELSE
            w2 = 0
            j = l2
            DO WHILE (work(work_level2+j)==currlevel2)
              IF (work(work_level1+j)>currlevel1) THEN
                w2 = w2 + a_weight(j)
              END IF
              j = work(work_next2+j)
              IF (j==0) EXIT
            END DO
            IF (w2==0) THEN
              currlevel2 = currlevel2 + 1
              l2 = j
              CYCLE
            ELSE
              CALL cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2, &
                sumweight,ratio,imbal,control%cost_function,t2)
            END IF
          END IF

          ! Add entries to relevant partition and update a_n1, a_n2 etc
          IF (t1<t2) THEN
            j = l1
            DO WHILE (work(work_level1+j)==currlevel1)
              IF (work(work_level2+j)>currlevel2) THEN
                work(work_part+j) = ND_PART1_FLAG
                a_n1 = a_n1 + 1
              END IF
              j = work(work_next1+j)
              IF (j==0) EXIT
            END DO
            a_weight_1 = a_weight_1 + w1
            a_weight_sep = a_weight_sep - w1
            l1 = j
            IF (j==0) THEN
              currlevel1 = maxlevel1 + 1
            ELSE
              currlevel1 = (work(work_level1+l1))
            END IF

          ELSE
            j = l2
            DO WHILE (work(work_level2+j)==currlevel2)
              IF (work(work_level1+j)>currlevel1) THEN
                work(work_part+j) = ND_PART2_FLAG
                a_n2 = a_n2 + 1
              END IF
              j = work(work_next2+j)
              IF (j==0) EXIT
            END DO
            a_weight_2 = a_weight_2 + w2
            a_weight_sep = a_weight_sep - w2
            l2 = j
            IF (j==0) THEN
              currlevel2 = maxlevel2 + 1
            ELSE
              currlevel2 = (work(work_level2+l2))
            END IF
          END IF
        END DO

        ! Reset partition matrix
        CALL nd_convert_flags_partition(a_n,a_n1,a_n2, &
          work(work_part+1:work_part+a_n),ND_PART1_FLAG,ND_PART2_FLAG, &
          partition(1:a_n))
        a_weight_sep = sumweight - a_weight_1 - a_weight_2
        ! call check_partition1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition)

      END SUBROUTINE nd_refine_block_trim


      ! ---------------------------------------------------
      ! nd_refine_max_flow
      ! ---------------------------------------------------
      ! Given a partition, trim the partition using blocks to make it minimal
      SUBROUTINE nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1, &
          a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end. This is updated to the new
        ! partition
        INTEGER, INTENT (OUT) :: work(8) ! Work array
        TYPE (nd_options), INTENT (IN) :: control

        ! ---------------------------------------------
        ! Local variables
        INTEGER :: msglvl
        REAL (wp) :: cost, ratio

        msglvl = 0
        IF (control%print_level==1 .AND. control%unit_diagnostics>=0) &
          msglvl = 1
        IF (control%print_level>=2 .AND. control%unit_diagnostics>=0) &
          msglvl = 3


        IF (a_n-a_n1-a_n2>1) THEN
          ratio = MAX(REAL(1.0,wp),control%balance)
          CALL nd_maxflow(a_n,a_ne,a_ptr,a_row,a_weight,control%cost_function,&
            a_n1,a_n2, &
            a_weight_1,a_weight_2,a_weight_sep,partition,ratio,msglvl, &
            work(1:8),cost)
        END IF


      END SUBROUTINE nd_refine_max_flow

      SUBROUTINE expand_partition(a_n,a_ne,a_ptr,a_row,a_weight,a_n1,a_n2, &
          a_weight_1,a_weight_2,a_weight_sep,partition,work)
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT (IN) :: a_ptr(a_n) ! On input a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(a_ne) ! On input a_row contains row
        ! indices of the non-zero rows. Diagonal entries have been removed
        ! and the matrix expanded.
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end. This is updated to the new
        ! partition
        INTEGER, INTENT (OUT) :: work(a_n) ! Work array

        ! Local variables
        INTEGER :: i, j, k, l, m, w
        INTEGER :: work_part, a_weight_sep_orig
        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what
        ! part of the partition the nodes are in
        work_part = 0
        DO i = 1, a_n1
          j = partition(i)
          work(work_part+j) = ND_PART1_FLAG
        END DO
        DO i = a_n1 + 1, a_n1 + a_n2
          j = partition(i)
          work(work_part+j) = ND_PART2_FLAG
        END DO
        DO i = a_n1 + a_n2 + 1, a_n
          j = partition(i)
          work(work_part+j) = ND_SEP_FLAG
        END DO

        ! IF (a_weight_1 .LT. a_weight_2) THEN
        ! side = ND_PART2_FLAG
        ! ELSE IF (a_weight_1 .GT. a_weight_2) THEN
        ! side = ND_PART1_FLAG
        ! ELSE
        ! side = ND_SEP_FLAG
        ! END IF
        a_weight_sep_orig = a_weight_sep

        DO i = a_n1 + a_n2 + 1, a_n
          j = partition(i)
          ! search neighbours of j and add to separator
          IF (j==a_n) THEN
            k = a_ne
          ELSE
            k = a_ptr(j+1) - 1
          END IF
          DO l = a_ptr(j), k
            m = a_row(l)
            IF (work(work_part+m)==ND_PART1_FLAG .AND. a_n1>1) THEN
              ! IF (side .EQ. ND_PART1_FLAG .OR. side .EQ. ND_SEP_FLAG)
              ! THEN
              work(work_part+m) = ND_SEP_FLAG
              a_n1 = a_n1 - 1
              w = a_weight(m)
              a_weight_1 = a_weight_1 - w
              a_weight_sep = a_weight_sep + w
              ! END IF
            ELSE IF (work(work_part+m)==ND_PART2_FLAG .AND. a_n2>1) THEN
              ! IF (side .EQ. ND_PART2_FLAG .OR. side .EQ. ND_SEP_FLAG)
              ! THEN
              work(work_part+m) = ND_SEP_FLAG
              a_n2 = a_n2 - 1
              w = a_weight(m)
              a_weight_2 = a_weight_2 - w
              a_weight_sep = a_weight_sep + w
              ! END IF
            END IF
          END DO
        END DO
        j = 1
        k = j + a_n1
        l = k + a_n2
        DO i = 1, a_n
          m = work(work_part+i)
          SELECT CASE (m)
          CASE (ND_PART1_FLAG)
            partition(j) = i
            j = j + 1
          CASE (ND_PART2_FLAG)
            partition(k) = i
            k = k + 1
          CASE DEFAULT
            partition(l) = i
            l = l + 1
          END SELECT
        END DO

      END SUBROUTINE expand_partition


      SUBROUTINE cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight, &
          ratio,imbal,costf,tau)

        INTEGER, INTENT (IN) :: a_weight_1, a_weight_2, a_weight_sep ! Weighte
        ! d
        ! size of partitions and separator
        INTEGER, INTENT (IN) :: sumweight
        REAL (wp), INTENT (IN) :: ratio
        LOGICAL, INTENT (IN) :: imbal ! Use penalty function?
        INTEGER, INTENT (IN) :: costf
        REAL (wp), INTENT (OUT) :: tau
        REAL (wp) :: beta, a_wgt1, a_wgt2

        beta = 0.5
        a_wgt1 = MAX(1,a_weight_1)
        a_wgt2 = MAX(1,a_weight_2)

        IF (costf.LE.1) THEN
          tau = ((REAL(a_weight_sep)**1.0)/REAL(a_wgt1))/REAL(a_wgt2)
          IF (imbal .AND. REAL(MAX(a_wgt1,a_wgt2))/REAL(MIN(a_wgt1, &
              a_wgt2))>=ratio) THEN
            tau = REAL(sumweight-2) + tau
          END IF
        ELSE
          IF (imbal .AND. REAL(MAX(a_wgt1,a_wgt2))/REAL(MIN(a_wgt1, &
              a_wgt2))>=ratio) THEN
            tau = REAL(sumweight)*(1.0+beta) + REAL(a_weight_sep)*( &
              1.0_wp+beta*REAL(ABS(a_wgt1-a_wgt2))/REAL(sumweight))
          ELSE
            tau = REAL(a_weight_sep)*(1.0_wp+beta*REAL(ABS(a_wgt1- &
              a_wgt2))/REAL(sumweight))
          END IF

        END IF

      END SUBROUTINE cost_function

      ! ---------------------------------------------------
      ! nd_maxflow
      ! ---------------------------------------------------
      ! Given a partition, get better partition using maxflow algorithm
      SUBROUTINE nd_maxflow(a_n,a_ne,a_ptr,a_row,a_weight,costf,a_n1,a_n2, &
          a_weight_1,a_weight_2,a_weight_sep,partition,alpha,msglvl,stats, &
          cost)
        IMPLICIT NONE

        ! Input matrix: a_n, a_ne, a_ptr, a_row
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix (lower and
        ! upper triangle)
        INTEGER, INTENT (IN) :: a_ptr(:) ! On input, a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(:) ! On input, a_row contains row
        ! indices of the nonzero entries. Diagonal entries have been
        ! removed and the matrix expanded.
        ! At the moment weights are not used at all
        INTEGER, INTENT (IN) :: a_weight(a_n) ! On input, a_weight(i) contains
        ! the weight of column i
        INTEGER, INTENT (IN) :: costf ! Determines which cost function is used
        ! Data on partition a_n1, a_n2, partition ... will be updated
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1 (ie B)
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2 (ie W)
        INTEGER, INTENT (INOUT) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
        ! hted
        ! size of partitions and separator
        INTEGER, INTENT (INOUT) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end. This is updated to the new
        ! partition.

        ! Parameters alpha (for balance) for cost function
        REAL (wp), INTENT (IN) :: alpha
        INTEGER, INTENT (IN) :: msglvl
        ! output --
        ! stats[1] -- weight of vertices in S
        ! stats[2] -- weight of vertices in B
        ! stats[3] -- weight of vertices in W
        ! stats[4] -- weight of edges in A_{S,S}
        ! stats[5] -- weight of edges in A_{S,B}
        ! stats[6] -- weight of edges in A_{S,W}
        ! stats[7] -- weight of edges in A_{B,B}
        ! stats[8] -- weight of edges in A_{B,W}
        ! cost     -- cost of new partition
        INTEGER, INTENT (OUT) :: stats(8)
        REAL (wp), INTENT (OUT) :: cost

        TYPE (network) :: netw

        ! Work arrays
        ! map,mapL,mapR of length a_n
        ! dmapL,dmapR,vwts of length a_ns
        ! mark1,mark2,pred,list of length number of nodes in network that
        ! is bounded by 2*a_ns+2

        ! Eventually these arrays will be allocated higher up the chain and
        ! will be passed as parameters
        INTEGER, ALLOCATABLE :: map(:), mapl(:), mapr(:)
        INTEGER, ALLOCATABLE :: dmapl(:), dmapr(:)
        INTEGER, ALLOCATABLE :: vwts(:)
        INTEGER, ALLOCATABLE :: sedge(:,:)
        INTEGER, ALLOCATABLE :: mark1(:), mark2(:), pred(:), list(:)
        REAL (wp), ALLOCATABLE :: imb(:)

        ! Local variables
        INTEGER :: a_ns, i, istart_s, j1, k, lp, wtw, wtb, statsr(9), &
          statsl(9)
        INTEGER nedge, matsiz
        REAL (wp) :: costr, costl

        lp = 6

        ! write(9,*) 'Entering maxflow'
        ! write(0,*) 'Entering maxflow'
        ALLOCATE (map(a_n),mapl(a_n),mapr(a_n))


        ! Number vertices in separator
        a_ns = a_n - a_n1 - a_n2

        ALLOCATE (dmapl(2*a_ns),dmapr(2*a_ns),vwts(a_ns))

        ! Allocate network work arrays.  Length is upper bound
        ! ALLOCATE (mark1(2*a_ns+2),mark2(2*a_ns+2),pred(2*a_ns+2), &
        ! list(2*a_ns+2))
        matsiz = MAX(a_n,2*a_ns+2)
        ALLOCATE (mark1(matsiz),mark2(matsiz),pred(2*a_ns+2),list(matsiz))
        ALLOCATE (imb(2*a_ns))

        ! Set up map array to define in what partition each vertex lies
        ! At same time set weights for partition (can check with Sue's input)

        CALL nd_convert_partition_flags(a_n,a_n1,a_n2,partition,1,2,0, &
          map(1:a_n))

        wtb = a_weight_1
        wtw = a_weight_2

        DO i = a_n1 + a_n2 + 1, a_n
          k = partition(i)
          vwts(i-a_n1-a_n2) = a_weight(k)
        END DO


        ! Count edges to get upper bound on size of sedge array
        nedge = 0
        DO k = 1, a_ns
          i = partition(a_n1+a_n2+k)
          j1 = a_ptr(i)
          IF (i==a_n) THEN
            nedge = nedge + a_ne + 1 - a_ptr(i)
          ELSE
            nedge = nedge + a_ptr(i+1) - a_ptr(i)
          END IF
        END DO
        ALLOCATE (sedge(nedge,2))


        ! Generate network graph.  The structure for our maxflow algorithm

        ! Work arrays dmapL and dmapR used to hold isadjsource and isadjsink
        ! Work array  mapL used to hold sep_map
        ! Source is associated with partition B (size a_n1)
        ! Sink is associated with partition W   (size a_n2)
        ! write(9,*) 'Calling mk_network'
        ! write(0,*) 'Calling mk_network'
        CALL mk_network(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,map,a_ns, &
          msglvl,netw,vwts,wtb,wtw,sedge,mapl,dmapl,dmapr,pred,list,mark1, &
          mark2,imb)
        ! write(9,*) 'Leaving mk_network'
        ! write(0,*) 'Leaving mk_network'


        ! solve a max flow problem to find the two new maps dmapL and dmapR

        CALL solvemaxflow(netw,dmapl,dmapr,mark1,mark2,pred,list)


        IF (msglvl>2) THEN
          WRITE (lp,*) 'dmapL ...'
          WRITE (lp,'(10I4)') dmapl
          WRITE (lp,*) 'dmapR ...'
          WRITE (lp,'(10I4)') dmapr
        END IF

        mapl = map
        mapr = map
        istart_s = a_n1 + a_n2
        DO i = 1, a_ns
          mapl(partition(istart_s+i)) = dmapl(i)
          mapr(partition(istart_s+i)) = dmapr(i)
        END DO

        IF (msglvl>2) THEN
          WRITE (lp,*) 'mapL ...'
          WRITE (lp,'(10I4)') mapl
          WRITE (lp,*) 'mapR ...'
          WRITE (lp,'(10I4)') mapr
        END IF

        ! Use evaluation function to choose best partition from among these
        ! two
        ! Use Sue's weighted code
        CALL evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,mapl,alpha,costf,statsl,costl)
        CALL evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,mapr,alpha,costf,statsr,costr)


        ! Find the better of the two partitions


        IF (statsl(9)==1 .AND. statsr(9)==1) THEN
          IF (msglvl>0) WRITE (lp,'(A)') 'both maps are acceptable'
          IF (costl<=costr) THEN
            map = mapl
            stats = statsl(1:8)
            cost = costl
            IF (msglvl>0) WRITE (lp,'(A)') 'left map accepted'
          ELSE
            map = mapr
            stats = statsr(1:8)
            cost = costr
            IF (msglvl>0) WRITE (lp,'(A)') 'right map accepted'
          END IF
        ELSE IF (statsl(9)==1) THEN
          map = mapl
          stats = statsl(1:8)
          cost = costl
          IF (msglvl>0) WRITE (lp,'(A)') &
            'right map NOT acceptable, left map accepted'
        ELSE IF (statsr(9)==1) THEN
          map = mapr
          stats = statsr(1:8)
          cost = costr
          IF (msglvl>0) WRITE (lp,'(A)') &
            'left map NOT acceptable, right map accepted'
        ELSE
          IF (msglvl>0) WRITE (lp,'(A)') 'NEITHER map acceptable'
          IF (costl<=costr) THEN
            map = mapl
            stats = statsl(1:8)
            cost = costl
            IF (msglvl>0) WRITE (lp,'(A)') 'left map accepted'
          ELSE
            map = mapr
            stats = statsr(1:8)
            cost = costr
            IF (msglvl>0) WRITE (lp,'(A)') 'right map accepted'
          END IF
        END IF
        a_weight_1 = stats(2)
        a_weight_2 = stats(3)
        a_weight_sep = stats(1)

        ! Now update partition
        ! First count number of vertices in each part
        a_n1 = 0
        a_n2 = 0
        a_ns = 0
        DO i = 1, a_n
          IF (map(i)==1) a_n1 = a_n1 + 1
          IF (map(i)==2) a_n2 = a_n2 + 1
          IF (map(i)==0) a_ns = a_ns + 1
        END DO


        CALL nd_convert_flags_partition(a_n,a_n1,a_n2,map(1:a_n),1,2, &
          partition(1:a_n))

        DEALLOCATE (map)
        DEALLOCATE (dmapl,dmapr)
        DEALLOCATE (mapl,mapr)
        DEALLOCATE (mark1,mark2,pred,list,vwts)
        DEALLOCATE (sedge)


      END SUBROUTINE nd_maxflow

      SUBROUTINE solvemaxflow(netw,maps1,maps2,mark1,mark2,pred,list)
        ! Find two partitions of a wide separator graph by solving a max flow
        ! problem




        ! output --

        ! mapS1[n_S] -- first map from wide separator to {0,1,2} = {S,B,W}
        ! mapS2[n_S] -- second map from wide separator to {0,1,2} = {S,B,W}


        ! input/output --

        ! network graph

        ! work arrays

        ! mark1, mark2, pred, list of length nnode


        IMPLICIT NONE
        TYPE (network), INTENT (INOUT) :: netw
        INTEGER, INTENT (OUT) :: maps1(:), maps2(:)

        INTEGER :: mark1(:), mark2(:), pred(:), list(:)

        ! Local variables
        INTEGER ii, lp, narc, nnode, u

        lp = 6


        nnode = netw%nnode
        narc = netw%narc


        ! Find maxflow through the network using the Ford-Fulkerson algorithm

        CALL findmaxflow(netw,pred,list,mark1,mark2)



        ! Find the two mincuts

        CALL findmincut(netw,mark1,mark2,list)


        ! Use mark1 and mark2 to generate maps1 and maps2


        maps1 = -1
        maps2 = -1

        DO ii = 2, nnode - 1, 2
          u = ii/2
          IF (mark1(ii)==1) THEN
            IF (mark1(ii+1)==1) THEN
              maps1(u) = 1
            ELSE
              maps1(u) = 0
            END IF
          END IF
          IF (mark1(ii)==2) THEN
            IF (mark1(ii+1)==2) THEN
              maps1(u) = 2
            ELSE
              maps1(u) = 0
            END IF
          END IF
        END DO

        DO ii = 2, nnode - 1, 2
          u = ii/2
          IF (mark2(ii)==1) THEN
            IF (mark2(ii+1)==1) THEN
              maps2(u) = 1
            ELSE
              maps2(u) = 0
            END IF
          END IF
          IF (mark2(ii)==2) THEN
            IF (mark2(ii+1)==2) THEN
              maps2(u) = 2
            ELSE
              maps2(u) = 0
            END IF
          END IF
        END DO

      END SUBROUTINE solvemaxflow


      SUBROUTINE mk_network(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,map,nvtx, &
          msglvl,netw,vwts,wtb,wtw,sedge,sep_map,isadjtosource,isadjtosink, &
          count,list,mark1,mark2,imb)
        ! Create and return a network structure

        IMPLICIT NONE

        ! Input matrix: a_n, a_ne, a_ptr, a_row
        INTEGER, INTENT (IN) :: a_n ! order of matrix
        INTEGER, INTENT (IN) :: a_ne ! number of entries in matrix (lower and
        ! upper triangle)
        INTEGER, INTENT (IN) :: a_ptr(:) ! On input, a_ptr(i) contains
        ! position in a_row that entries for column i start.
        INTEGER, INTENT (IN) :: a_row(:) ! On input, a_row contains row
        ! indices of the nonzero entries. Diagonal entries have been
        ! removed and the matrix expanded.

        INTEGER, INTENT (IN) :: partition(a_n) ! First a_n1 entries contain
        ! list of (local) indices in partition 1; next a_n2 entries
        ! contain list of (local) entries in partition 2; entries in
        ! separator are listed at the end.
        INTEGER, INTENT (IN) :: map(a_n) ! First a_n1 entries contain
        INTEGER, INTENT (IN) :: msglvl, nvtx
        INTEGER, INTENT (IN) :: vwts(:)
        INTEGER, INTENT (INOUT) :: wtb, wtw
        INTEGER, INTENT (INOUT) :: a_n1 ! Size of partition 1 (ie B)
        INTEGER, INTENT (INOUT) :: a_n2 ! Size of partition 2 (ie W)

        ! Note that we still do the allocations here.  Doing it further up the
        ! call tree would need to allocate much more storage.
        TYPE (network), INTENT (OUT) :: netw

        INTEGER i, iarc, ii, jj, lp, narc1, narc2, narc3, narc4, narc, nnode, &
          nedge, u, v, wts
        INTEGER j1, j2, j, k

        ! Work arrays sedge(nedge,2),sep_map(a_n),isAdjToSource(nvtx),
        ! isAdjToSink(nvtx)
        INTEGER :: sedge(:,:)
        INTEGER :: sep_map(:)
        INTEGER :: isadjtosource(:), isadjtosink(:)
        REAL (wp) :: imb(:)
        INTEGER :: COUNT(:), mark1(:), mark2(:), list(:)
        LOGICAL :: augcap

        augcap = .TRUE.

        lp = 0

        IF (msglvl>0) THEN
          WRITE (lp,'(/A)') '### inside mknetwork()'
          WRITE (lp,'(A,I8,A,I8)') 'nvtx', nvtx
        END IF


        ! write(0,*) 'a_n,a_n1,a_n2,a_ne',   &
        ! a_n,a_n1,a_n2,a_ne
        ! write(0,*) 'a_ptr',a_ptr(1:a_n)
        ! write(0,*) 'a_row',a_row(1:a_ne)
        ! write(0,*) 'partition',partition(1:a_n)
        ! write(0,*) 'map',map(1:a_n)
        ! write(0,*) 'vwts',vwts(1:nvtx)
        ! write(0,*) 'wtb',wtb
        ! write(0,*) 'wtw',wtw

        isadjtosource = 0
        isadjtosink = 0

        ! Determine mapping of global variables of matrix to separator set
        DO k = 1, nvtx
          i = partition(a_n1+a_n2+k)
          sep_map(i) = k
        END DO

        ! We could use a single array although the logic for generating it is
        ! marginally more complicated.
        ! For nodes in separator, set isadj as
        ! 1 if only connected to source
        ! 2 if only connected to sink
        ! 3 if connected to source and sink
        ! isadj = 0


        ! Run through nodes in separator S and generate edges
        nedge = 0
        DO k = 1, nvtx
          i = partition(a_n1+a_n2+k)
          j1 = a_ptr(i)
          IF (i==a_n) THEN
            j2 = a_ne
          ELSE
            j2 = a_ptr(i+1) - 1
          END IF

          ! Run through vertices connected to vertex i
          DO jj = j1, j2
            j = a_row(jj)
            ! Find out in which partition node j lies using map array
            IF (map(j)==1) THEN
              ! If in partition B add vertex k to AdjToSource
              isadjtosource(k) = 1
            END IF
            IF (map(j)==2) THEN
              ! If in partition W add vertex k to AdjToSink
              isadjtosink(k) = 1
            END IF
            IF (map(j)==0) THEN
              ! If in separator add edge accumulating number of edges
              nedge = nedge + 1
              ! Edge has this orientation to emulate matlab code
              sedge(nedge,2) = sep_map(i)
              sedge(nedge,1) = sep_map(j)
            END IF
          END DO

        END DO


        ! narc1 is number of vertices in separator and is equal to the number
        ! of
        ! added edges u- to u+
        narc1 = nvtx
        ! narc2 is number of vertices in separator connected to source
        narc2 = 0
        ! narc3 is number of vertices in separator connected to sink
        narc3 = 0


        ! count the number of arcs
        ! Can't do easily in above loop because of multiple edges from
        ! source/sink
        ! to vertex in the separator set

        DO k = 1, nvtx
          IF (isadjtosource(k)==1) narc2 = narc2 + 1
          IF (isadjtosink(k)==1) narc3 = narc3 + 1
        END DO

        narc4 = 0
        DO ii = 1, nedge
          u = sedge(ii,1)
          v = sedge(ii,2)
          ! --- ignore self edges ---
          IF (u==v) CYCLE
          ! --- omit edges with essential vertices ---
          IF (isadjtosink(u)==1 .AND. isadjtosource(u)==1) CYCLE
          IF (isadjtosink(v)==1 .AND. isadjtosource(v)==1) CYCLE
          ! --- omit pairs both adjacent to source ---
          IF ((isadjtosource(u)==1) .AND. (isadjtosource(v)==1)) CYCLE
          ! --- omit pairs both adjacent to sink ---
          IF ((isadjtosink(u)==1) .AND. (isadjtosink(v)==1)) CYCLE
          ! Qualifying arc found
          narc4 = narc4 + 1
        END DO

        nnode = 2*nvtx + 2
        netw%nnode = nnode
        narc = narc1 + narc2 + narc3 + narc4
        netw%narc = narc
        netw%source = 1
        netw%sink = 2*nvtx + 2

        IF (msglvl>0) THEN
          WRITE (lp,'(I8,A)') narc1, ' internal arcs'
          WRITE (lp,'(I8,A)') narc2, ' arcs from source'
          WRITE (lp,'(I8,A)') narc3, ' arcs from sink'
          WRITE (lp,'(I8,A)') narc4, ' edge arcs'
          WRITE (lp,'(I8,A)') narc, ' total arcs'
        END IF


        ! create the arc arrays



        ! Allocations done here but could be easily moved up the path although
        ! values very dependent on separator size.
        ALLOCATE (netw%firsts(narc),netw%seconds(narc),netw%capacities(narc))
        ALLOCATE (netw%flows(narc))
        ALLOCATE (netw%inheads(nnode),netw%outheads(nnode),netw%nextin(narc), &
          netw%nextout(narc))

        netw%firsts = -1
        netw%seconds = -1


        ! (u-,u+) arcs first

        iarc = 0
        DO u = 1, nvtx
          iarc = iarc + 1
          netw%firsts(iarc) = 2*u
          netw%seconds(iarc) = 2*u + 1
          ! We set capacities after computing imbalance penalty
          ! netw%capacities(iarc) = vwts(u)
        END DO

        IF (msglvl>0) WRITE (lp,'(A,I8)') 'after (u-,u+) arcs, iarc = ', iarc


        ! (source,u) arcs

        DO u = 1, nvtx
          IF (isadjtosource(u)==1) THEN
            iarc = iarc + 1
            netw%firsts(iarc) = netw%source
            netw%seconds(iarc) = 2*u
            netw%capacities(iarc) = HUGE(1)/2
          END IF
        END DO

        IF (msglvl>0) WRITE (lp,'(A,I8)') 'after (source,u-) arcs, iarc = ', &
          iarc


        ! (u,sink) arcs

        DO u = 1, nvtx
          IF (msglvl>5) WRITE (lp,'(A,I4,A,I8)') 'isAdjToSink(', u, ')= ', &
            isadjtosink(u)
          IF (isadjtosink(u)==1) THEN
            iarc = iarc + 1
            netw%firsts(iarc) = 2*u + 1
            netw%seconds(iarc) = netw%sink
            netw%capacities(iarc) = HUGE(1)/2
          END IF
        END DO

        IF (msglvl>0) WRITE (lp,'(A,I8)') 'after (u,sink) arcs, iarc = ', iarc


        ! (u+,v-) arcs

        DO ii = 1, nedge
          u = sedge(ii,1)
          v = sedge(ii,2)
          IF ((u/=v) .AND. (isadjtosource(u)/=1 .OR. isadjtosource( &
              v)/=1) .AND. (isadjtosink(u)/=1 .OR. isadjtosink( &
              v)/=1) .AND. (isadjtosource(u)/=1 .OR. isadjtosink( &
              u)/=1) .AND. (isadjtosource(v)/=1 .OR. isadjtosink(v)/=1)) THEN
            iarc = iarc + 1
            netw%firsts(iarc) = 2*u + 1
            netw%seconds(iarc) = 2*v
            netw%capacities(iarc) = HUGE(1)/2
          END IF
        END DO

        IF (msglvl>0) WRITE (lp,'(A,I8)') 'after (u+,v-) arcs, iarc = ', iarc


        ! Generate the head vectors for in/out edges
        ! and the in/out link vectors for the arcs


        netw%inheads = -1
        netw%outheads = -1
        netw%nextin = -1
        netw%nextout = -1
        DO ii = narc, 1, -1
          u = netw%firsts(ii)
          v = netw%seconds(ii)
          IF (msglvl>1) WRITE (lp,'(A,I8,A,I8,A,I8,A)') 'ii', ii, 'arc (', u, &
            ',', v, ')'
          netw%nextin(ii) = netw%inheads(v)
          netw%inheads(v) = ii
          netw%nextout(ii) = netw%outheads(u)
          netw%outheads(u) = ii
        END DO

        IF (augcap) THEN
          ! Generate wtS
          wts = 0
          DO i = 1, nvtx
            wts = wts + vwts(i)
          END DO
          IF (msglvl>0) WRITE (lp,*) 'Calling findpenalty'
          ! Compute network capacities for separator arcs.
          CALL findpenalty(msglvl,a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
            map,vwts,wtb,wtw,wts,isadjtosource,isadjtosink,mark1,mark2,list, &
            imb)
          IF (msglvl>0) WRITE (lp,*) 'Exiting findpenalty'
          DO i = 1, nvtx
            netw%capacities(i) = mark1(i)
          END DO

        ELSE
          DO i = 1, nvtx
            netw%capacities(i) = vwts(i)
          END DO
        END IF

        ! ISD Would this not be initialized in routine that computes flows?
        netw%flows = 0

        IF (msglvl>0) WRITE (lp,'(A/)') '### leaving mknetwork()'

      END SUBROUTINE mk_network


      SUBROUTINE findmaxflow(netw,pred,list,tags,deltas)
        ! Find a maximum flow through the network


        IMPLICIT NONE
        TYPE (network), INTENT (INOUT) :: netw

        INTEGER avail, iarc, lp, nnode, sink, source, stats(2), tag

        ! Work arrays ... all of length nnode

        INTEGER :: pred(:), list(:), tags(:), deltas(:)

        lp = 6

        nnode = netw%nnode
        source = netw%source
        sink = netw%sink

        ! Network flows initialized to zero
        netw%flows = 0

        ! tag is just used to count which arc from source is being used to
        ! start augmenting path
        tag = 0
        iarc = netw%outheads(source)
        ! Run through all nodes connected to source
        DO
          IF (iarc==-1) EXIT

          DO
            avail = netw%capacities(iarc) - netw%flows(iarc)
            IF (avail>0) THEN
              ! Set in findaugpath

              ! use BFS to find path

              ! avail is an output from findaugpath giving available flow on
              ! augmenting
              ! path.  Not to be confused with dummy avail above. The
              ! augmenting
              ! path is given through the array pred.

              CALL findaugpath(netw,iarc,avail,pred,stats,list,tags,deltas)


              ! Go to next arc from source node if no augmenting path has been
              ! found
              IF ((avail==0) .OR. (pred(sink)==0)) EXIT

              ! Update flows
              CALL augmentpath(netw,avail,pred)

            ELSE
              EXIT
            END IF
          END DO

          iarc = netw%nextout(iarc)
          tag = tag + 1

        END DO

      END SUBROUTINE findmaxflow

      SUBROUTINE findmincut(netw,mark1,mark2,list)

        ! Finds one or two mincuts, one nearest the source, one nearest the
        ! sink

        ! Input parameters

        ! netw    -- network object
        ! msglvl  -- message level

        ! Output parameters
        ! mark1 --- to identify cut set nearest source
        ! mark2 --- to identify cut set nearest sink


        ! Workspace
        ! list --- to hold list of nodes being searched

        IMPLICIT NONE
        INTEGER, INTENT (OUT) :: mark1(:), mark2(:), list(:)
        TYPE (network), INTENT (INOUT) :: netw

        ! Local variables
        INTEGER iarc, last, lp, nnode, now, sink, source, x, z

        lp = 6

        nnode = netw%nnode
        source = netw%source
        sink = netw%sink


        ! breadth first traversal from source


        mark1 = 2
        mark1(source) = 1

        list = 0
        now = 1
        last = 1
        list(now) = source
        ! while now <= last
        DO
          IF (now>last) EXIT
          x = list(now)
          now = now + 1
          iarc = netw%outheads(x)
          ! while iarc ~= -1
          ! Run through all arcs starting at node x and putting node at end of
          ! arc
          ! on list if there is spare capacity on arc
          DO
            IF (iarc==-1) EXIT
            z = netw%seconds(iarc)
            IF (mark1(z)==1) THEN
            ELSE
              IF (netw%flows(iarc)<netw%capacities(iarc)) THEN
                last = last + 1
                list(last) = z
                mark1(z) = 1
              END IF
            END IF
            iarc = netw%nextout(iarc)
          END DO
          iarc = netw%inheads(x)
          ! while iarc ~= -1
          ! Run through all arcs terminating at node x and putting node at
          ! start of arc
          ! on list if there is spare capacity on arc
          DO
            IF (iarc==-1) EXIT
            z = netw%firsts(iarc)
            IF (mark1(z)==1) THEN
            ELSE
              IF (netw%flows(iarc)>0) THEN
                last = last + 1
                list(last) = z
                mark1(z) = 1
              END IF
            END IF
            iarc = netw%nextin(iarc)
          END DO
        END DO

        ! breadth first traversal from sink


        mark2 = 1
        mark2(sink) = 2

        list = 0
        now = 1
        last = 1
        list(now) = sink
        ! while now <= last
        DO
          IF (now>last) EXIT
          x = list(now)
          now = now + 1
          iarc = netw%outheads(x)
          ! while iarc ~= -1
          ! Run through all arcs starting at node x and putting node at end of
          ! arc
          ! on list if there is spare capacity on arc
          DO
            IF (iarc==-1) EXIT
            z = netw%seconds(iarc)
            IF (mark2(z)==2) THEN
            ELSE
              IF (netw%flows(iarc)>0) THEN
                last = last + 1
                list(last) = z
                mark2(z) = 2
              END IF
            END IF
            iarc = netw%nextout(iarc)
          END DO
          iarc = netw%inheads(x)
          ! while iarc ~= -1
          ! Run through all arcs terminating at node x and putting node at
          ! start of arc
          ! on list if there is spare capacity on arc
          DO
            IF (iarc==-1) EXIT
            z = netw%firsts(iarc)
            IF (mark2(z)==2) THEN
            ELSE
              IF (netw%flows(iarc)<netw%capacities(iarc)) THEN
                last = last + 1
                list(last) = z
                mark2(z) = 2
              END IF
            END IF
            iarc = netw%nextin(iarc)
          END DO
        END DO


      END SUBROUTINE findmincut

      SUBROUTINE findpenalty(msglvl,a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
          map,vwts,wtb,wtw,wts,count,head,mark,mark1,list,imb)

        ! count is isAdjToSource in call
        ! head is isAdjToSink in call

        ! Computes augmented capacities using weighting based on balance

        ! Input parameters

        ! msglvl  -- message level

        ! Output parameters
        ! mark --- records level set of vertex in cutset and two-level set
        ! also
        ! mark1 --- records level set of vertex in cutset starting from sink


        IMPLICIT NONE
        INTEGER, INTENT (IN) :: msglvl, a_n, a_ne, a_n1, a_n2
        INTEGER, INTENT (IN) :: a_ptr(:), a_row(:)
        INTEGER, INTENT (IN) :: vwts(:), partition(:), map(:)
        INTEGER, INTENT (INOUT) :: COUNT(:), head(:)
        INTEGER, INTENT (INOUT) :: wtb, wtw, wts
        INTEGER, INTENT (OUT) :: mark(:), mark1(:)

        ! Workspace
        ! list --- to hold list of nodes being searched
        ! length bounded by a_n
        INTEGER :: list(:)
        REAL (wp) :: imb(:)

        ! Local variables
        INTEGER inode, last, lp, maxl, minl, x, z
        INTEGER i, j1, j2, jj, k, nvtx
        INTEGER begin_lev, end_lev, ilev, penp

        lp = 0

        ! Set imbalance weighting
        penp = 100

        ! nvtx is number of vertices in separator
        nvtx = a_n - a_n1 - a_n2

        IF (msglvl>0) THEN
          WRITE (lp,'(A)') ''
          WRITE (lp,'(A)') '### inside findpenalty()'
        END IF


        ! Breadth first traversal from source
        ! Source is defined as all vertices in black partition

        IF (msglvl>0) WRITE (lp,'(A)') 'breadth first traversal from source'

        mark = 0
        ! Black vertices at level 0
        ! Run through separator vertices
        last = 0
        DO k = a_n1 + a_n2 + 1, a_n
          z = partition(k)
          ! Check if adjacent to source
          IF (COUNT(k-a_n1-a_n2)==1) THEN
            last = last + 1
            list(last) = z
            mark(z) = 1
          END IF
        END DO
        end_lev = 0

        ! Each pass through this loop determines all nodes in level set k
        ! a_n is just a dummy end
        DO k = 2, a_n
          ! Run through all nodes in the previous level set
          begin_lev = end_lev + 1
          end_lev = last
          DO ilev = begin_lev, end_lev
            x = list(ilev)
            IF (msglvl>1) WRITE (lp,'(A,I10)') 'Processing vertex', x
            ! Run through vertices connected to vertex x
            j1 = a_ptr(x)
            IF (x==a_n) THEN
              j2 = a_ne
            ELSE
              j2 = a_ptr(x+1) - 1
            END IF
            DO jj = j1, j2
              z = a_row(jj)
              ! Jump if vertex not in separator
              IF (map(z)/=0) CYCLE
              ! Jump if vertex visited already
              IF (mark(z)/=0) CYCLE
              mark(z) = k
              ! write(0,*) 'z,mark(z)',z,mark(z)
              ! Add node z to list
              last = last + 1
              list(last) = z
            END DO
          END DO
          IF (last==end_lev) EXIT
        END DO ! end of processing of nodes on level k


        ! breadth first traversal from sink
        ! Sink is defined as all vertices in white partition

        IF (msglvl>0) WRITE (lp,'(A)') 'breadth first traversal from the sink'

        mark1 = 0
        ! Keep sink at level 0
        maxl = -a_n
        minl = a_n

        ! White nodes at level 0
        ! Put all separator vertices connected to sink in list
        last = 0
        DO k = a_n1 + a_n2 + 1, a_n
          z = partition(k)
          ! Check if adjacent to source
          IF (head(k-a_n1-a_n2)==1) THEN
            last = last + 1
            list(last) = z
            mark1(z) = 1
            ! k = 1
            mark(z) = mark(z) - 1
            minl = MIN(minl,mark(z))
            maxl = MAX(maxl,mark(z))
          END IF
        END DO
        end_lev = 0

        ! Each pass through this loop determines all nodes in level set k
        ! a_n is just a dummy end
        DO k = 2, a_n
          ! Run through all nodes in the previous level set
          begin_lev = end_lev + 1
          end_lev = last
          DO ilev = begin_lev, end_lev
            x = list(ilev)
            IF (msglvl>1) WRITE (lp,'(A,I10)') 'Processing vertex', x
            ! Run through vertices connected to vertex x
            j1 = a_ptr(x)
            IF (x==a_n) THEN
              j2 = a_ne
            ELSE
              j2 = a_ptr(x+1) - 1
            END IF
            DO jj = j1, j2
              z = a_row(jj)
              ! Jump if vertex not in separator
              IF (map(z)/=0) CYCLE
              ! Jump if vertex visited already
              IF (mark1(z)/=0) CYCLE
              mark1(z) = k
              ! write(0,*) 'z,mark1(z)',z,mark1(z)
              ! write(0,*) 'mark(z)',mark(z)
              mark(z) = mark(z) - k
              ! write(0,*) 'z,mark(z)',z,mark(z)
              minl = MIN(minl,mark(z))
              maxl = MAX(maxl,mark(z))
              ! Add node z to list
              last = last + 1
              list(last) = z
            END DO
          END DO
          IF (last==end_lev) EXIT
        END DO ! end of processing of nodes on level k

        ! Compute half-level sets
        IF (msglvl>1) WRITE (lp,'(A,2I4)') 'minl, maxl ', minl, maxl

        ! We will number levels from 1 to maxl-minl+1
        ! count will hold total weight of all vertices in half-level set
        ! Nodes in level set are accessed by linked list in list with headers
        ! in head
        COUNT(1:maxl-minl+1) = 0
        head(1:maxl-minl+1) = -1

        ! Map the mark values into the local coordinates where separator
        ! vertices numbered from 1 to nvtx
        DO i = 1, nvtx
          k = partition(a_n1+a_n2+i)
          mark1(i) = mark(k)
        END DO

        ! Run through all notes in cutset, resetting numbering of level set
        ! putting
        ! them in level set for level and accumulating weight of half-level
        ! set
        DO i = 1, nvtx
          mark1(i) = mark1(i) - minl + 1
          list(i) = head(mark1(i))
          head(mark1(i)) = i
          COUNT(mark1(i)) = COUNT(mark1(i)) + vwts(i)
        END DO

        IF (msglvl>1) THEN
          WRITE (lp,'(A)') 'Number of vertices in each half-level set'
          DO i = 1, maxl - minl + 1
            WRITE (lp,'(2I10)') i, COUNT(i)
          END DO
        END IF

        ! Run through half-level sets computing imbalances
        ! wtB is weight of B
        ! wtW is set to total weight of rest of network
        wtw = wtw + wts
        IF (msglvl>1) WRITE (lp,('(A,3I10)')) 'wtB,wtW,wtS', wtb, wtw, wts
        IF (maxl-minl==0) THEN
          ! Only one level set
          imb(1) = MAX(REAL(wtb)/REAL(wtw),REAL(wtw)/REAL(wtb))
        ELSE
          wtw = wtw - COUNT(1)
         ! IF (msglvl>0) WRITE (14,'(A)') &
         !   'Half-level set   width   |B|,|W|, imbalance'
          DO k = 1, maxl - minl
            wtw = wtw - COUNT(k+1)
            imb(k) = MAX(REAL(wtb)/REAL(wtw),REAL(wtw)/REAL(wtb))
         !   IF (msglvl>0) WRITE (14,'(I10,4G12.2)') k, COUNT(k), wtb, wtw, &
         !     imb(k)
            wtb = wtb + COUNT(k)
          END DO
        END IF

        IF (msglvl>1) THEN
          WRITE (lp,'(A)') 'Imbalances'
          DO i = 1, maxl - minl
            WRITE (lp,'(I10,G12.2)') i, imb(i)
          END DO
        END IF

        ! Run through nodes in level set assigning penalty to them
        DO k = 1, maxl - minl + 1
          inode = head(k)
          DO
            IF (inode==-1) EXIT
            IF (msglvl>1) WRITE (lp,('(A,2I4)')) 'level and node', k, inode
            IF (k==1) THEN
              mark(inode) = FLOOR(penp*imb(1)*vwts(inode))
            ELSE
              IF (k==maxl-minl+1) mark(inode) = FLOOR(penp*imb(maxl-minl)*vwts &
                (inode))
              IF (k>1 .AND. k<maxl-minl+1) mark(inode) = FLOOR(penp*MIN(imb( &
                k),imb(k-1))*vwts(inode))
            END IF
            inode = list(inode)
          END DO
        END DO

        IF (msglvl>1) THEN
          WRITE (lp,'(A)') 'Computed penalties'
          DO i = 1, nvtx
            WRITE (lp,'(2I10)') i, mark(i)
          END DO
        END IF

        IF (msglvl>0) WRITE (lp,'(A/)') '### leaving findpenalty()'

      END SUBROUTINE findpenalty

      SUBROUTINE augmentpath(netw,delta,pred)

        ! Reset flows on augmenting path


        ! Input
        ! delta   -- increment flow
        ! pred    -- tree predecessor vector, size nnode
        ! msglvl  -- message level

        ! Input/output
        ! network -- network object

        IMPLICIT NONE
        INTEGER, INTENT (IN) :: delta
        TYPE (network), INTENT (INOUT) :: netw

        INTEGER iarc, lp, sink, source, v, w

        INTEGER, INTENT (IN) :: pred(:)


        lp = 6

        source = netw%source
        sink = netw%sink

        ! Should set an error flag
        IF (delta<=0 .OR. pred(sink)<=0) THEN
          WRITE (lp,'(A,I4,A,I4)') 'ERROR : delta', delta, ', pred(sink) = ', &
            pred(sink)
          RETURN
        END IF


        ! work back from the sink resetting network flows

        w = sink
        ! while w ~= source
        DO
          IF (w==source) EXIT
          iarc = pred(w)
          IF (netw%firsts(iarc)==w) THEN
            v = netw%seconds(iarc)
            netw%flows(iarc) = netw%flows(iarc) - delta
          ELSE IF (netw%seconds(iarc)==w) THEN
            v = netw%firsts(iarc)
            netw%flows(iarc) = netw%flows(iarc) + delta
          END IF
          w = v
        END DO

      END SUBROUTINE augmentpath

      SUBROUTINE findaugpath(netw,iarc_m,avail,pred,stats,list,tags,deltas)

        ! Find an augmenting path starting from arc iarc_m
        ! here we use a breadth first search (BFS),
        ! in findaugpath2() we use a depth first search (DFS)
        ! in findaugpath3() we will use a max distance
        ! from source to grow the tree to find a path


        ! input --
        ! netw -- network object
        ! iarc_m -- label for starting arc (u,v)


        ! output --
        ! avail -- if nonzero, available flow on augmenting path
        ! pred -- tree predecessor vector, size nnode
        ! source <-- pred^m(sink) <-- pred(pred(sink)) ...
        ! <-- pred(sink) <-- sink
        ! stats -- statistics
        ! stats(1) = # nodes visited in search
        ! stats(2) = # arcs visited


        ! working --
        ! list -- stack vector used for depth first search, size nnode
        ! tags -- mark vector, size nnode
        ! deltas -- increment flow vector, size nnode

        IMPLICIT NONE
        INTEGER, INTENT (IN) :: iarc_m
        INTEGER, INTENT (OUT) :: avail, stats(2)
        TYPE (network), INTENT (IN) :: netw
        INTEGER, INTENT (OUT) :: pred(:)

        INTEGER iarc, last, lp, nnode, now, root, sink, source, v, w

        INTEGER :: list(:), tags(:), deltas(:)

        INTEGER narc, u, n_nodevisit, n_arcvisit

        lp = 6

        ! As input variable is intent(in) we set a local value that we will
        ! update
        iarc = iarc_m

        nnode = netw%nnode
        narc = netw%narc
        source = netw%source
        sink = netw%sink

        ! write(0,*) 'sink',sink

        stats(1) = 0
        stats(2) = 0

        list = 0

        ! Initial working storage and array pred

        tags = 0
        deltas = 0
        pred = 0


        ! check that (source,u) is an edge
        ! that is, that iarc as input is an edge from the source node

        u = netw%seconds(iarc)

        ! This will never be the case because all arcs iarc_m come from source
        ! node
        IF (netw%firsts(iarc)/=source) THEN
          WRITE (lp,'(A,I4,A)') 'u', u, 'is not adjacent to source'
          RETURN
        END IF


        ! check for available capacity

        avail = netw%capacities(iarc) - netw%flows(iarc)
        IF (avail==0) RETURN


        ! Find augmenting path using an alternating tree

        root = u
        now = 1
        last = 1
        list(1) = root
        ! tags is used to tell whether node has been visited on this attempt
        ! to find an augmenting path
        tags(root) = root
        tags(source) = root
        tags(sink) = -1
        pred(sink) = -1
        deltas(root) = avail
        pred(root) = iarc
        n_nodevisit = 0
        n_arcvisit = 0

        ! while now <= last
        DO
          IF (now>last) EXIT
          v = list(now)
          now = now + 1
          n_nodevisit = n_nodevisit + 1

          iarc = netw%outheads(v)
          ! while iarc ~= -1
          DO
            ! Run through all edges emanating from v
            ! First is v^- to v^+
            IF (iarc==-1) EXIT
            w = netw%seconds(iarc)
            n_arcvisit = n_arcvisit + 1

            IF (tags(w)/=root) THEN
              ! Node w has not yet been visited

              IF (netw%capacities(iarc)>netw%flows(iarc)) THEN
                avail = netw%capacities(iarc) - netw%flows(iarc)

                IF (avail>deltas(v)) THEN
                  avail = deltas(v)
                END IF
                deltas(w) = avail
                pred(w) = iarc
                ! Flag w as being visited
                tags(w) = root

                IF (w==sink) EXIT

                last = last + 1
                list(last) = w

              END IF
            END IF
            ! Go to next arc from v
            iarc = netw%nextout(iarc)
          END DO
          IF (w==sink) EXIT


          iarc = netw%inheads(v)
          ! while iarc ~= -1
          DO
            ! Run through all edges coming in to v
            IF (iarc==-1) EXIT
            w = netw%firsts(iarc)
            n_arcvisit = n_arcvisit + 1
            IF (tags(w)/=root) THEN
              IF (netw%flows(iarc)>0) THEN
                IF (avail>netw%flows(iarc)) THEN
                  avail = netw%flows(iarc)
                END IF
                deltas(w) = avail
                pred(w) = iarc
                tags(w) = root
                last = last + 1
                list(last) = w
              END IF
            END IF
            iarc = netw%nextin(iarc)
          END DO
          ! Don't think you can reach this statement
          IF (w==sink) EXIT
        END DO

        ! Flag to show augmenting path not found
        IF (w/=sink) avail = 0

        stats(1) = n_nodevisit
        stats(2) = n_arcvisit


      END SUBROUTINE findaugpath

      SUBROUTINE evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,map,alpha,costf,stats, &
          stats10)
        ! Matlab call
        ! function stats = evalBSW ( A, map, alpha, beta, msglvl )

        ! stats = EVALBSW ( A, map, alpha, beta )

        ! input ---
        ! map[nvtx] -- map from vertices to region
        ! map[u] == 0 --> u in S
        ! map[u] == 1 --> u in B
        ! map[u] == 2 --> u in W
        ! alpha --- acceptability parameter
        ! beta  --- imbalance penalty parameter

        ! output --
        ! stats[1] -- weight of vertices in S
        ! stats[2] -- weight of vertices in B
        ! stats[3] -- weight of vertices in W
        ! stats[4] -- weight of edges in A_{S,S}
        ! stats[5] -- weight of edges in A_{S,B}
        ! stats[6] -- weight of edges in A_{S,W}
        ! stats[7] -- weight of edges in A_{B,B}
        ! stats[8] -- weight of edges in A_{B,W}
        ! stats[9] -- 1 if acceptable, 0 if not
        ! acceptable --> alpha*min(|B|,|W|) >= max(|B|,|W|)
        ! stats10 -- cost of partition
        ! cost = |S|*(1 + (beta*| |B| - |W| |)/(|B|+|S|+|W|)) ;

        ! created -- 12jan12, cca

        IMPLICIT NONE
        INTEGER, INTENT (IN) :: a_n
        INTEGER, INTENT (IN) :: a_ne
        INTEGER, INTENT (IN) :: map(:), a_ptr(:), a_row(:), a_weight(:)
        REAL (wp), INTENT (IN) :: alpha
        INTEGER, INTENT (IN) :: costf
        INTEGER, INTENT (OUT) :: stats(9)
        REAL (wp), INTENT (OUT) :: stats10
        INTEGER minbw, maxbw, nss, nsb, nsw, nbb, nww, nvtx, ns, nb, nw
        INTEGER j, j1, j2, jj, u, v
        REAL (wp) diffbw, beta
        LOGICAL :: imbal

        beta = 0.5_wp
        nvtx = a_n
        stats(1:9) = -1
        ns = 0
        nb = 0
        nw = 0
        DO u = 1, nvtx
          IF (map(u)==0) THEN
            ns = ns + a_weight(u)
          ELSE IF (map(u)==1) THEN
            nb = nb + a_weight(u)
          ELSE IF (map(u)==2) THEN
            nw = nw + a_weight(u)
          END IF
        END DO
        stats(1) = ns
        stats(2) = nb
        stats(3) = nw
        minbw = MIN(nb,nw)
        maxbw = MAX(nb,nw)
        diffbw = REAL(ABS(nb-nw))/REAL(ns+nb+nw)
        IF (.FALSE.) THEN
          nss = 0
          nsb = 0
          nsw = 0
          nbb = 0
          nww = 0
          ! [rows, cols, ents] = find(A) ;
          ! nzA = length(rows) ;
          DO j = 1, a_n
            j1 = a_ptr(j)
            IF (j==a_n) THEN
              j2 = a_ne
            ELSE
              j2 = a_ptr(j+1) - 1
            END IF
            v = j
            DO jj = j1, j2
              u = a_row(jj)
              ! v = cols(ii) ;
              IF (map(u)==0) THEN
                IF (map(v)==0) THEN
                  nss = nss + 1
                ELSE IF (map(v)==1) THEN
                  nsb = nsb + 1
                ELSE IF (map(v)==2) THEN
                  nsw = nsw + 1
                END IF
              ELSE IF (map(u)==1) THEN
                IF (map(v)==1) THEN
                  nbb = nbb + 1
                END IF
              ELSE IF (map(u)==2) THEN
                IF (map(v)==2) THEN
                  nww = nww + 1
                END IF
              END IF
            END DO
          END DO
          stats(4) = nss
          stats(5) = nsb
          stats(6) = nsw
          stats(7) = nbb
          stats(8) = nww
        END IF
        ! stats[9] -- 1 if acceptable, 0 if not
        ! acceptable --> alpha*min(|B|,|W|) >= max(|B|,|W|)
        ! stats10 -- cost of partition
        IF (alpha*minbw>=maxbw) THEN
          stats(9) = 1
        ELSE
          stats(9) = 0
        END IF

        IF (alpha > nb+nw+ns-2) THEN
            imbal = .FALSE.
        ELSE
            imbal = .TRUE.
        END IF

        CALL cost_function(nb,nw,ns,nb+nw+ns, &
          alpha,imbal,costf,stats10)

      END SUBROUTINE evalbsw

      SUBROUTINE nd_supervars(n,ptr,row,perm,invp,nsvar,svar,info,st)
        ! Detects supervariables - modified version of subroutine from hsl_mc78
        INTEGER, INTENT (INOUT) :: n ! Dimension of system
        INTEGER, DIMENSION (n+1), INTENT (IN) :: ptr ! Column pointers
        INTEGER, DIMENSION (ptr(n+1)-1), INTENT (IN) :: row ! Row indices
        INTEGER, DIMENSION (n), INTENT (INOUT) :: perm
        ! perm(i) must hold position of i in the pivot sequence.
        ! On exit, holds the pivot order to be used by factorization.
        INTEGER, DIMENSION (n), INTENT (INOUT) :: invp ! inverse of perm
        INTEGER, INTENT (OUT) :: nsvar ! number of supervariables
        INTEGER, DIMENSION (n), INTENT (OUT) :: svar ! number of vars in each
        ! svar
        INTEGER, INTENT (OUT) :: info
        INTEGER, INTENT (OUT) :: st

        LOGICAL :: full_rank ! flags if supervariable 1 has ever become
        ! empty.
        ! If it has not, then the varaibles in s.v. 1 are those that never
        ! occur
        INTEGER :: i
        INTEGER (long) :: ii
        INTEGER :: j
        INTEGER :: idx ! current index
        INTEGER :: next_sv ! head of free sv linked list
        INTEGER :: nsv ! new supervariable to move j to
        INTEGER :: piv ! current pivot
        INTEGER :: col ! current column of A
        INTEGER :: sv ! current supervariable
        INTEGER :: svc ! temporary holding supervariable count
        INTEGER, DIMENSION (:), ALLOCATABLE :: sv_new ! Maps each
        ! supervariable to
        ! a new supervariable with which it is associated.
        INTEGER, DIMENSION (:), ALLOCATABLE :: sv_seen ! Flags whether
        ! svariables have
        ! been seen in the current column. sv_seen(j) is set to col when svar
        ! j
        ! has been encountered.
        INTEGER, DIMENSION (:), ALLOCATABLE :: sv_count ! number of variables
        ! in sv.

        info = 0 ! by default completed succefully

        ALLOCATE (sv_new(n+1),sv_seen(n+1),sv_count(n+1),STAT=st)
        IF (st/=0) THEN
          info = ND_ERR_MEMORY_ALLOC
          RETURN
        END IF

        svar(:) = 1
        sv_count(1) = n
        sv_seen(1) = 0

        ! Setup linked list of free super variables
        next_sv = 2
        DO i = 2, n
          sv_seen(i) = i + 1
        END DO
        sv_seen(n+1) = -1

        ! Determine supervariables using modified Duff and Reid algorithm
        full_rank = .FALSE.
        DO col = 1, n
          IF (ptr(col+1)/=ptr(col)) THEN
            ! If column is not empty, add implicit diagonal entry
            j = col
            sv = svar(j)
            IF (sv_count(sv)==1) THEN ! Are we only (remaining) var in sv
              full_rank = full_rank .OR. (sv==1)
              ! MUST BE the first time that sv has been seen for this
              ! column, so just leave j in sv, and go to next variable.
              ! (Also there can be no other vars in this block pivot)
            ELSE
              ! There is at least one other variable remaining in sv
              ! MUST BE first occurence of sv in the current row/column,
              ! so define a new supervariable and associate it with sv.
              sv_seen(sv) = col
              sv_new(sv) = next_sv
              nsv = next_sv
              next_sv = sv_seen(next_sv)
              sv_new(nsv) = nsv ! avoids problems with duplicates
              sv_seen(nsv) = col
              ! Now move j from sv to nsv
              nsv = sv_new(sv)
              svar(j) = nsv
              sv_count(sv) = sv_count(sv) - 1
              sv_count(nsv) = 1
              ! This sv cannot be empty as initial sv_count was > 1
            END IF
          END IF
          DO ii = ptr(col), ptr(col+1) - 1
            j = row(ii)
            sv = svar(j)
            IF (sv_count(sv)==1) THEN ! Are we only (remaining) var in sv
              full_rank = full_rank .OR. (sv==1)
              ! If so, and this is first time that sv has been seen for this
              ! column, then we can just leave j in sv, and go to next
              ! variable.
              IF (sv_seen(sv)<col) CYCLE
              ! Otherwise, we have already defined a new supervariable
              ! associated
              ! with sv. Move j to this variable, then retire (now empty) sv.
              nsv = sv_new(sv)
              IF (sv==nsv) CYCLE
              svar(j) = nsv
              sv_count(nsv) = sv_count(nsv) + 1
              ! Old sv is now empty, add it to top of free stack
              sv_seen(sv) = next_sv
              next_sv = sv
            ELSE
              ! There is at least one other variable remaining in sv
              IF (sv_seen(sv)<col) THEN
                ! this is the first occurence of sv in the current row/column,
                ! so define a new supervariable and associate it with sv.
                sv_seen(sv) = col
                sv_new(sv) = next_sv
                sv_new(next_sv) = next_sv ! avoids problems with duplicates
                next_sv = sv_seen(next_sv)
                sv_count(sv_new(sv)) = 0
                sv_seen(sv_new(sv)) = col
              END IF
              ! Now move j from sv to nsv
              nsv = sv_new(sv)
              svar(j) = nsv
              sv_count(sv) = sv_count(sv) - 1
              sv_count(nsv) = sv_count(nsv) + 1
              ! This sv cannot be empty as sv_count was > 1
            END IF
          END DO
        END DO


        ! Now modify pivot order such that all variables in each supervariable
        ! are
        ! consecutive. Do so by iterating over pivots in elimination order. If
        ! a
        ! pivot has not already been listed, then order that pivot followed by
        ! any other pivots in that supervariable.

        ! We will build a new inverse permutation in invp, and then find perm
        ! afterwards. First copy invp to perm:
        perm(:) = invp(:)
        ! Next we iterate over the pivots that have not been ordered already
        ! Note: as we begin, all entries of sv_seen are less than or equal to
        ! n+1
        ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable
        ! has been
        ! ordered.
        idx = 1
        nsvar = 0
        DO piv = 1, n
          IF (sv_seen(piv)>n+1) CYCLE ! already ordered
          ! Record information for supervariable
          sv = svar(perm(piv))
          IF ( .NOT. full_rank .AND. sv==1) CYCLE ! Don't touch unused vars
          nsvar = nsvar + 1
          svc = sv_count(sv)
          sv_new(nsvar) = svc ! store # vars in s.v. to copy to svar
          ! later
          j = piv
          ! Find all variables that are members of sv and order them.
          DO WHILE (svc>0)
            DO j = j, n
              IF (svar(perm(j))==sv) EXIT
            END DO
            sv_seen(j) = n + 2 ! flag as ordered
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
          END DO
        END DO
        ! Push unused variables to end - these are those vars still in s.v. 1
        IF ( .NOT. full_rank) THEN
          svc = sv_count(1)
          ! Find all variables that are members of sv and order them.
          j = 1
          DO WHILE (svc>0)
            DO j = j, n
              IF (svar(perm(j))==1) EXIT
            END DO
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
          END DO
          n = n - sv_count(1)
        END IF

        ! Recover perm as inverse of invp
        DO piv = 1, n
          perm(invp(piv)) = piv
        END DO
        ! sv_new has been used to store number of variables in each svar, copy
        ! into
        ! svar where it is returned.
        svar(1:nsvar) = sv_new(1:nsvar)

        DEALLOCATE (sv_new,sv_seen,sv_count,STAT=st)
        IF (st/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF
      END SUBROUTINE nd_supervars


      ! This subroutine takes a set of supervariables and compresses the
      ! supplied
      ! matrix using them.

      SUBROUTINE nd_compress_by_svar(n,ne,ptr,row,invp,nsvar,svar,ptr2, &
          row2,info,st)
        INTEGER, INTENT (IN) :: n ! Dimension of system
        INTEGER, INTENT (IN) :: ne ! Number off-diagonal zeros in system
        INTEGER, DIMENSION (n+1), INTENT (IN) :: ptr ! Column pointers
        INTEGER, DIMENSION (ptr(n+1)-1), INTENT (IN) :: row ! Row indices
        INTEGER, DIMENSION (n), INTENT (IN) :: invp ! inverse of perm
        INTEGER, INTENT (IN) :: nsvar
        INTEGER, DIMENSION (nsvar), INTENT (IN) :: svar ! super variables of A
        INTEGER, DIMENSION (nsvar+1), INTENT (OUT) :: ptr2
        INTEGER, DIMENSION (ne), INTENT (OUT) :: row2
        INTEGER, INTENT (OUT) :: info
        INTEGER, INTENT (OUT) :: st

        INTEGER :: piv, svc, sv, col
        INTEGER :: j, idx
        INTEGER, DIMENSION (:), ALLOCATABLE :: flag, sv_map

        info = 0 ! by default completed succefully

        ALLOCATE (flag(nsvar),sv_map(n),STAT=st)
        IF (st/=0) THEN
          info = ND_ERR_MEMORY_ALLOC
          RETURN
        END IF
        flag(:) = 0

        ! Setup sv_map
        piv = 1
        DO svc = 1, nsvar
          DO piv = piv, piv + svar(svc) - 1
            sv_map(invp(piv)) = svc
          END DO
        END DO

        piv = 1
        idx = 1
        DO svc = 1, nsvar
          col = invp(piv)
          ptr2(svc) = idx
          DO j = ptr(col), ptr(col+1) - 1
            sv = sv_map(row(j))
            IF (flag(sv)==piv) CYCLE ! Already dealt with this supervariable
            ! Add row entry for this sv
            row2(idx) = sv
            flag(sv) = piv
            idx = idx + 1
          END DO
          piv = piv + svar(svc)
        END DO
        ptr2(svc) = idx

        DEALLOCATE (flag,sv_map,STAT=st)
        IF (st/=0) THEN
          info = ND_ERR_MEMORY_DEALLOC
          RETURN
        END IF
      END SUBROUTINE nd_compress_by_svar

    END MODULE spral_nd

