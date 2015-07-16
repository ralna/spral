    PROGRAM example
      USE spral_nestd
      IMPLICIT NONE

      ! Local variables
      INTEGER :: n, ne
      INTEGER, DIMENSION (:), ALLOCATABLE :: row, ptr, perm

      TYPE (mc70_options) :: options
      TYPE (mc70_inform) :: inform

      ! Read in the order n of the matrix and the number
      ! of non-zeros in its lower triangular part.
      READ (5,*) mtx, n, ne

      ! Allocate arrays
      ALLOCATE (row(ne),ptr(n+1),perm(n),STAT=st)
      IF (st/=0) THEN
        WRITE (6,*) ' Allocation error'
        STOP
      END IF

      ! Read in pointers
      READ (5,*) ptr(1:n+1)

      ! Read in row indices for lower triangular part
      READ (5,*) row(1:ne)

      ! Call nested dissection and switch to approximate minimum degree when
      ! sub matrix has order less than or equal to 4
      options%amd_switch1 = 4
      options%amd_call = 3
      mtx = 0
      CALL nestd_order(mtx,n,ptr,row,perm,options,inform)

      ! Print out nested dissection ordering
      WRITE (6,'(a)') ' Permutation : '
      WRITE (6,'(8i8)') perm

      ! Deallocate all arrays
      DEALLOCATE (row,ptr,perm)

    END PROGRAM example
