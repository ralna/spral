    PROGRAM main
      USE spral_nd
      IMPLICIT NONE

      LOGICAL :: ok

      CALL test_errors(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_input(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_dense(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_ashcraft(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_levelset(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_shift(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_multigrid(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_amd(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_misc(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_refine(ok)
      IF (.NOT. ok) GO TO 80
      CALL test_fm(ok)
      IF (.NOT. ok) GO TO 80


80    IF (ok) THEN
        WRITE (*,'(a)') 'All tests successfully completed'
      ELSE
        WRITE (*,'(a)') 'WARNING: Tests NOT successfully completed'
      END IF

    CONTAINS

      SUBROUTINE reset_control(control_orig,control_reset)


        TYPE (nd_options), INTENT (IN) :: control_orig

        TYPE (nd_options), INTENT (OUT) :: control_reset

        control_reset%print_level = control_orig%print_level
        control_reset%unit_diagnostics = control_orig%unit_diagnostics
        control_reset%unit_error = control_orig%unit_error
        control_reset%amd_call = control_orig%amd_call
        control_reset%amd_switch1 = control_orig%amd_switch1
        control_reset%amd_switch2 = control_orig%amd_switch2
        control_reset%cost_function = control_orig%cost_function
        control_reset%partition_method = control_orig%partition_method
        control_reset%matching = control_orig%matching
        control_reset%coarse_partition_method = control_orig% &
          coarse_partition_method
        control_reset%refinement = control_orig%refinement
        control_reset%refinement_band = control_orig%refinement_band
        control_reset%remove_dense_rows = control_orig%remove_dense_rows
        control_reset%stop_coarsening2 = control_orig%stop_coarsening2
        control_reset%stop_coarsening1 = control_orig%stop_coarsening1
        control_reset%ml_call = control_orig%ml_call
        control_reset%min_reduction = control_orig%min_reduction
        control_reset%max_reduction = control_orig%max_reduction
        control_reset%balance = control_orig%balance
        control_reset%max_improve_cycles = control_orig%max_improve_cycles
        control_reset%find_supervariables = control_orig%find_supervariables


      END SUBROUTINE reset_control

      SUBROUTINE testpr(n,perm)
        ! .. Scalar Arguments ..
        INTEGER n
        ! ..
        ! .. Array Arguments ..
        INTEGER perm(n)
        ! ..
        ! .. Local Scalars ..
        INTEGER i, ip
        ! .. Local Arrays
        INTEGER, ALLOCATABLE, DIMENSION (:) :: iw
        ! ..
        ! Initialize array used to test for valid permutation.
        ALLOCATE (iw(n))
        iw(1:n) = 1
        DO i = 1, n
          ip = abs(perm(i))
          IF (iw(ip)==0) THEN
            WRITE (68,*) i, ip
            GO TO 10
          ELSE
            iw(ip) = 0
          END IF
        END DO
        DEALLOCATE (iw)
        RETURN
10      WRITE (6,'(A)') '**** Error in permutation'
        DEALLOCATE (iw)
      END SUBROUTINE testpr

      SUBROUTINE ndef(n,nz,a,ip,ind,w,sym)
        ! This subroutine augments the diagonal to make the matrix pos def
        ! Then large entries put in off-diagonal to make it a little not so.
        ! .. Parameters ..
        DOUBLE PRECISION zero, one
        PARAMETER (zero=0.0D0,one=1.0D0)
        ! ..
        ! .. Scalar Arguments ..
        INTEGER n, nz
        LOGICAL sym
        ! ..
        ! .. Array Arguments ..
        DOUBLE PRECISION a(nz), w(n)
        INTEGER ind(nz), ip(*)
        ! ..
        ! .. Local Scalars ..
        DOUBLE PRECISION rmax
        INTEGER i, i1, i2, idiag, ii, ioff, j, maxoff, mprint, numoff
        ! ..
        ! .. Intrinsic Functions ..
        INTRINSIC abs, max
        ! ..
        numoff = 0
        ! !! MAXOFF was 1 .. now 10 to try to have more 2 x 2 pivots
        maxoff = 10
        mprint = 6
        DO i = 1, n
          w(i) = zero
        END DO
        DO j = 1, n
          rmax = zero
          IF (sym) rmax = w(j)
          idiag = 0
          ioff = 0
          i1 = ip(j)
          i2 = ip(j+1) - 1
          IF (i2<i1) THEN
            GO TO 30
          ELSE
            DO ii = i1, i2
              i = ind(ii)
              rmax = rmax + abs(a(ii))
              w(i) = w(i) + abs(a(ii))
              IF (i==j) idiag = ii
              IF (i/=j) ioff = ii
            END DO
            IF (idiag==0) THEN
              GO TO 30
            ELSE
              a(idiag) = rmax + one
              IF (ioff/=0 .AND. numoff<maxoff) THEN
                a(ioff) = 1.1*a(idiag)
                ! !! added
                a(idiag) = zero
                numoff = numoff + 1
              END IF
            END IF
          END IF
        END DO
        IF ( .NOT. sym) THEN
          DO j = 1, n
            i1 = ip(j)
            i2 = ip(j+1) - 1
            DO ii = i1, i2
              i = ind(ii)
              IF (i==j) GO TO 10
            END DO
            GO TO 20
10          a(ii) = max(a(ii),w(i)+one)
20          CONTINUE
          END DO
        END IF
        GO TO 40
30      WRITE (mprint,fmt=90000) j
        nz = -1
40      RETURN
90000   FORMAT (' NO DIAGONAL ENTRY IN COLUMN ',I8,/, &
          ' SO CANNOT MAKE MATRIX DIAGONALLY DOMINANT')
      END SUBROUTINE ndef

      SUBROUTINE pdef(n,nz,a,ip,ind,w,sym)
        ! THIS SUBROUTINE AUGMENTS THE DIAGONAL TO MAKE THE MATRIX POS DEF.
        ! .. Parameters ..
        DOUBLE PRECISION zero, one
        PARAMETER (zero=0.0D0,one=1.0D0)
        ! ..
        ! .. Scalar Arguments ..
        INTEGER n, nz
        LOGICAL sym
        ! ..
        ! .. Array Arguments ..
        DOUBLE PRECISION a(nz), w(n)
        INTEGER ind(nz), ip(*)
        ! ..
        ! .. Local Scalars ..
        DOUBLE PRECISION rmax
        INTEGER i, i1, i2, idiag, ii, j, mprint
        ! ..
        ! .. Intrinsic Functions ..
        INTRINSIC abs, max
        ! ..
        mprint = 6
        DO i = 1, n
          w(i) = zero
        END DO
        DO j = 1, n
          rmax = zero
          IF (sym) rmax = w(j)
          idiag = 0
          i1 = ip(j)
          i2 = ip(j+1) - 1
          IF (i2<i1) THEN
            GO TO 30
          ELSE
            DO ii = i1, i2
              i = ind(ii)
              IF (i/=j) THEN
                rmax = rmax + abs(a(ii))
                w(i) = w(i) + abs(a(ii))
              END IF
              IF (i==j) idiag = ii
            END DO
            IF (idiag==0) THEN
              GO TO 30
            ELSE
              a(idiag) = rmax + one
            END IF
          END IF
        END DO
        IF ( .NOT. sym) THEN
          DO j = 1, n
            i1 = ip(j)
            i2 = ip(j+1) - 1
            DO ii = i1, i2
              i = ind(ii)
              IF (i==j) GO TO 10
            END DO
            GO TO 20
10          a(ii) = max(a(ii),w(i)+one)
20          CONTINUE
          END DO
        END IF
        GO TO 40
30      WRITE (mprint,fmt=90000) j
        nz = -1
40      RETURN
90000   FORMAT (' NO DIAGONAL ENTRY IN COLUMN ',I8,/, &
          ' SO CANNOT MAKE MATRIX DIAGONALLY DOMINANT')
      END SUBROUTINE pdef


     SUBROUTINE test_errors(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: n, ne, st, err_count, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control
      TYPE (nd_inform) :: info

      ok = .TRUE.

      ! --------------------------------------
      ! Test error flag n<1
      ! --------------------------------------
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a)') '****Test n<1 error ****'
      n = 0
      ne = 0
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1

      control%amd_switch1 = 2
      err_count = 0
      DO i = 0, 1
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag==-3) err_count = err_count + 1
        CALL nd_order(1,n,ptr,row,perm,control,info,seps)
        IF (info%flag==-3) err_count = err_count + 1
      END DO

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20

      IF (err_count/=4) THEN
        WRITE (6,'(a)') 'Code failure in test section'
        ok = .FALSE.
      END IF
      RETURN

10    WRITE (*,'(a)') 'Allocation error during test'
      RETURN

20    WRITE (*,'(a)') 'Deallocation error during test'


     END SUBROUTINE test_errors

! *****************************************************************

     SUBROUTINE test_input(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i, j
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      LOGICAL :: corr
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info


      ok = .TRUE.

      test_count = 0
      ! --------------------------------------
      ! Test diagonal matrix
      ! --------------------------------------
      testi_count = 0
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a)') '**** Test inputs: diagonal matrix ****'
      n = 10
      ne = 0
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = 1

      control%amd_switch1 = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        corr = .TRUE.
        DO j = 1, n
          corr = (corr .AND. (perm(j)==j))
        END DO
        IF (info%flag>=0 .AND. corr) testi_count = testi_count + 1
        CALL nd_order(1,n,ptr,row,perm,control,info,seps)
        corr = .TRUE.
        DO j = 1, n
          corr = (corr .AND. (perm(j)==j))
        END DO
        IF (info%flag>=0 .AND. corr) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20

      IF (testi_count/=2) THEN
        WRITE (6,'(a29)') 'Code failure in test section '
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expansion of matrix from lower triangle input
      ! --------------------------------------
      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
      n = 10
      ne = 9
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n-1) = (/ (i,i=1,n-1) /)
      ptr(n:n+1) = n
      row(1:ne) = (/ (i,i=2,n) /)

      control%amd_switch1 = 2
      control%amd_call = 2
      DO i = 0, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expansion of matrix from lower triangle format and that diags are
      ! removed
      ! --------------------------------------

      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
      n = 10
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n) = (/ (2*i-1,i=1,n) /)
      ptr(n+1) = ne + 1
      row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
        10 /)

      control%amd_switch1 = 2
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expansion of matrix from lower and upper triangle input with no
      ! diags
      ! --------------------------------------
      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
      n = 10
      ne = 18
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1
      ptr(2:n) = (/ (2*(i-1),i=2,n) /)
      ptr(n+1) = ne + 1
      row(1) = 2
      DO i = 2, n - 1
        row(ptr(i)) = i - 1
        row(ptr(i)+1) = i + 1
      END DO
      row(ptr(n)) = n - 1

      control%amd_switch1 = 2
      control%amd_call = 2
      DO i = 2, 2
        control%print_level = i
        CALL nd_order(1,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expansion of matrix from lower and upper triangle input with diag
      ! entry
      ! --------------------------------------
      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a34,i5,a4)') '****Test inputs: matrix expansion ', test, ' ***'
      n = 10
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1
      ptr(2:n) = (/ (2*(i-1),i=2,n) /)
      ptr(n+1) = ne + 1
      row(1) = 2
      DO i = 2, n - 1
        row(ptr(i)) = i - 1
        row(ptr(i)+1) = i + 1
      END DO
      row(ptr(n)) = n - 1
      row(ne) = n

      control%amd_switch1 = 2
      control%amd_call = 2
      DO i = 0, 2
        control%print_level = i
        CALL nd_order(1,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF





      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

     END SUBROUTINE test_input


! *****************************************************************


     SUBROUTINE test_dense(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0
      ! --------------------------------------
      ! Test dense row removal - row 1 dense
      ! --------------------------------------
      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
      n = 800
      ne = n
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1
      ptr(2:n+1) = ne + 1
      row(1:ne) = (/ (i,i=1,n) /)

      control%amd_switch1 = 2
      control%amd_call = 20
      DO i = 0, 2
        control%print_level = i
        control%remove_dense_rows = .TRUE.
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==799 .AND. info%nzsuper==0) &
          testi_count = testi_count + 1
        control%remove_dense_rows = .FALSE.
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==800 .AND. info%nzsuper==1598) &
          testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=6) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test dense row removal - first row dense
      ! --------------------------------------
      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
      n = 800
      ne = 2*n - 3
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1
      ptr(2:n) = (/ (n+i-2,i=2,n) /)
      ptr(n+1) = ne + 1
      row(1:ptr(2)-1) = (/ (i,i=2,n) /)
      row(ptr(2):ptr(n)-1) = (/ (i+1,i=2,n-1) /)

      control%amd_switch1 = 2
      control%amd_call = 20
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==799 .AND. info%nzsuper==1596) &
          testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test dense row removal - last row dense
      ! --------------------------------------
      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
      n = 800
      ne = n - 1
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n-1) = (/ (i,i=1,n-1) /)
      ptr(n:n+1) = ne + 1
      row(1:ne) = n

      control%amd_switch1 = 2
      control%amd_call = 20
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==799 .AND. info%nzsuper==0) &
          testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test dense row removal - first two rows dense but row 2 has max degree
      ! --------------------------------------
      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Dense row removal test', test, ' ***'
      n = 800
      ne = n - 3 + n - 2
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1
      ptr(2) = n - 2
      ptr(3:n+1) = ne + 1
      row(ptr(1):ptr(2)-1) = (/ (i+3,i=1,n-3) /)
      row(ptr(2):ne) = (/ (i+2,i=1,n-2) /)

      control%amd_switch1 = 2
      control%amd_call = 20
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==798 .AND. info%nzsuper==0) &
          testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_dense


! *****************************************************************


     SUBROUTINE test_ashcraft(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0


      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
      n = 10
      ne = 14
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
      row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      DO i = 0, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
      n = 10
      ne = 15
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 4, 6, 8, 9, 10, 11, 13, 15, 16, 16 /)
      row(1:ne) = (/ 2, 3, 4, 3, 5, 4, 8, 6, 7, 9, 8, 10, 9, 10, 10 /)

      control%amd_switch1 = 2
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
      n = 13
      ne = 16
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 4, 6, 8, 10, 12, 13, 13, 14, 14, 16, 17, 17 /)
      row(1:ne) = (/ 2, 3, 4, 4, 11, 5, 6, 7, 8, 9, 10, 8, 10, 12, 13, 13 /)

      control%amd_switch1 = 2
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Ashcraft test ', test, ' ***'
      n = 13
      ne = 16
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 4, 6, 8, 10, 12, 13, 13, 14, 14, 16, 17, 17 /)
      row(1:ne) = (/ 2, 3, 4, 4, 11, 5, 6, 7, 8, 9, 10, 8, 10, 12, 13, 13 /)

      control%amd_switch1 = 2
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%balance = 20.0
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_ashcraft


! *****************************************************************


     SUBROUTINE test_levelset(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0


      ! --------------------------------------
      ! Test one-sided levelset method
      ! --------------------------------------

      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Level set test ', test, ' ***'
      n = 10
      ne = 14
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
      row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 2
      control%amd_call = 2
      control%balance = 20.0
      DO i = 0, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test one-sided levelset method
      ! --------------------------------------

      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Level set test ', test, ' ***'
      n = 10
      ne = 14
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 12, 13, 14, 15, 15 /)
      row(1:ne) = (/ 2, 7, 7, 3, 4, 8, 5, 9, 6, 10, 10, 8, 9, 10 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 2
      control%balance = 1.5
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_levelset

! *****************************************************************


     SUBROUTINE test_shift(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0

     ! --------------------------------------
      ! Test DM-style refinement with 2-sided partition
      ! --------------------------------------

      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 8
      ne = 11
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 4, 5, 6, 8, 9, 11, 12, 12 /)
      row(1:ne) = (/ 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 8 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      DO i = 2, 2
        control%print_level = i
        control%refinement = 2
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
        control%refinement = 5
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
        control%print_level = i
        control%refinement = 2
        CALL nd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1
        control%refinement = 5
        CALL nd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=4) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition
      ! --------------------------------------

      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 21
      ne = 37
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 9, 11, 13, 15, 16, 20, 22, 24, 26, 28, 29, &
        29, 32, 33, 35, 37, 38, 38 /)
      row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 16, &
        17, 11, 13, 12, 13, 13, 15, 14, 15, 15, 17, 18, 19, 18, 19, 20, 20, &
        21, 21 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%refinement = 2
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 1-sided partition
      ! --------------------------------------

      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
        27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
        13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
        24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 2
      control%refinement = 2
      control%amd_call = 2
      control%balance = 33.0
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition
      ! --------------------------------------

      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 7
      ne = 6
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 6, 7, 7 /)
      row(1:ne) = (/ 7, 7, 7, 7, 7, 7 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 2
      control%refinement = 2
      control%amd_call = 2
      control%balance = 12.0
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 1-sided partition
      ! --------------------------------------

      test = 5
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 7
      ne = 6
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 6, 7, 7 /)
      row(1:ne) = (/ 7, 7, 7, 7, 7, 7 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 2
      control%refinement = 5
      control%amd_call = 2
      control%balance = 12.0
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test DM refinement with 1-sided partition with balanced partition
      ! --------------------------------------

      test = 6
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
        27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
        13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
        24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 2
      control%refinement = 5
      control%balance = 2.0
      control%amd_call = 2
      DO i = 0, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 7
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
        27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
        13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
        24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%refinement = 5
      control%balance = 2.0
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 8
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%refinement = 2
      control%balance = 1.05
      control%refinement_band = n
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 9
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 13
      ne = 15
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 3, 6, 7, 8, 9, 11, 11, 13, 14, 15, 16, 16 /)
      row(1:ne) = (/ 2, 3, 4, 5, 6, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%refinement = 2
      control%balance = 1.30
      control%refinement_band = n
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 10
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a26,i5,a4)') '****Shift refinement test ', test, ' ***'
      n = 13
      ne = 15
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 16 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 7, 7, 8, 9, 10, 11, 11, 11, 12, 13 /)

      control%amd_switch1 = 5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%refinement = 2
      control%balance = 1.30
      control%refinement_band = n
      control%amd_call = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF





      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_shift

! *****************************************************************


     SUBROUTINE test_multigrid(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i, j
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0

      ! --------------------------------------
      ! Test multigrid, 2-sided partitioning
      ! --------------------------------------

      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 5
      control%balance = 1.05
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 2
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test multigrid, 1-sided partitioning
      ! --------------------------------------

      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 5
      control%balance = 1.05
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        DO j = 0, 2
          control%print_level = i
          control%matching = j
          CALL nd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=9) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid coarsening
      ! --------------------------------------

      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 5
      control%balance = 1.05
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 2
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test multigrid coarsening
      ! --------------------------------------

      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 5
      control%balance = 1.05
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        DO j = 0, 2
          control%print_level = i
          control%matching = j
          CALL nd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=9) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test multigrid coarsening
      ! --------------------------------------

      test = 5
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 12
      ne = 18
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 4, 6, 7, 9, 10, 12, 13, 16, 18, 19, 19 /)
      row(1:ne) = (/ 2, 9, 9, 4, 10, 10, 6, 11, 11, 8, 12, 12, 10, 11, 12, 11, &
        12, 12 /)

      control%amd_switch1 = 5
      control%balance = 1.05
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      control%min_reduction = 0.4
      control%matching = 1
      DO i = 0, 2
        control%print_level = i
        control%matching = j
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid refinement
      ! --------------------------------------

      test = 6
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 5
      control%balance = 1.05
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 7
        DO j = 0, 2
          control%refinement = i
          control%matching = j
          CALL nd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=24) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test no. multigrid - automatic choice refinement
      ! --------------------------------------

      test = 7
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 9
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
      row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
        9 /)

      control%amd_switch1 = 4
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 4
      control%print_level = 2
      DO i = 0, 7
        DO j = 0, 2
          control%refinement = 3
          control%matching = j
          CALL nd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=24) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid - automatic choice refinement
      ! --------------------------------------

      test = 8
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 9
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
      row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
        9 /)

      control%amd_switch1 = 4
      control%balance = 1.0
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 2
      DO i = 0, 7
        DO j = 0, 2
          control%refinement = i
          control%matching = j
          CALL nd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=24) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF




      ! --------------------------------------
      ! Test automatic multilevel-no multilevel check
      ! --------------------------------------

      test = 9
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 200
      ne = 199
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10


      ptr(1:n) = (/ (i,i=1,n) /)
      ptr(n+1) = ne + 1
      row(1:ne) = (/ (i,i=2,n) /)

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%balance = 8.0
      control%partition_method = 2
      control%coarse_partition_method = 2
      control%ml_call = 10
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%refinement = 7
      DO i = 0, 0
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test automatic multilevel-no multilevel check
      ! --------------------------------------

      test = 10
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 200
      ne = 7*(n-7)+21
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10


      ptr(1:n-6) = (/ (1+7*(i-1),i=1,n-6) /)

      ptr(n-5) = 1 + 6*(n-6) + 6
      ptr(n-4) = 1 + 6*(n-6) + 11
      ptr(n-3) = 1 + 6*(n-6) + 15
      ptr(n-2) = 1 + 6*(n-6) + 18
      ptr(n-1) = 1 + 6*(n-6) + 20
      ptr(n) = 1 + 6*(n-6) + 21
      ptr(n+1) = ne + 1
      j = 1
      DO i = 1,n-7
       row(j) = (i+1)
       row(j+1) = (i+2)
       row(j+2) = (i+3)
       row(j+3) = (i+4)
       row(j+4) = (i+5)
       row(j+5) = (i+6)
       row(j+6) = (i+7)
       j= j+7
      END DO
      
      i = n-6
      row(j) = (i+1)
      row(j+1) = (i+2)
      row(j+2) = (i+3)
      row(j+3) = (i+4)
      row(j+4) = (i+5)
      row(j+4) = (i+6)
      j = j+6
      i = n-5
      row(j) = (i+1)
      row(j+1) = (i+2)
      row(j+2) = (i+3)
      row(j+3) = (i+4)
      row(j+4) = (i+5)
      j = j+5
      i = n-4
      row(j) = (i+1)
      row(j+1) = (i+2)
      row(j+2) = (i+3)
      row(j+3) = (i+4)
      j = j+4
      i = n-3
      row(j) = (i+1)
      row(j+1) = (i+2)
      row(j+2) = (i+3)
      j = j+3
      i = n-2
      row(j) = (i+1)
      row(j+1) = (i+2)
      j = j+2
      i = n-1
      row(j) = (i+1)

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%balance = 8.0
      control%partition_method = 2
      control%ml_call = 10
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%refinement = 7
      DO i = 1, 2
        control%coarse_partition_method = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=2) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid reduction ratio
      ! --------------------------------------

      test = 11
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 20
      ne = n-1
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10


      ptr(1) = 1

      ptr(2:n+1) = ne + 1

      row(1:ne) = (/ (i,i=2,n) /)

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%balance = 2.0
      control%partition_method = 1
      control%ml_call = 10
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%refinement = 7
      control%find_supervariables = .false.
      control%print_level = 2
      DO i = 1, 2
        control%coarse_partition_method = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=2) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for full galerkin_graph_rap coverage
      ! --------------------------------------

      test = 12
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Multigrid  test ', test, '****'
      n = 16
      ne = 42
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      j = 1
      DO i = 1, 16
        ptr(i) = j
        IF (mod(i,4)/=0) THEN
          row(j) = i + 1
          j = j + 1
        END IF
        IF (i<13) THEN
          IF (mod(i,4)==1) THEN
            row(j) = i + 4
            j = j + 1
            row(j) = i + 5
            j = j + 1
          ELSE
            IF (mod(i,4)==0) THEN
              row(j) = i + 3
              j = j + 1
              row(j) = i + 4
              j = j + 1
            ELSE
              row(j) = i + 3
              j = j + 1
              row(j) = i + 4
              j = j + 1
              row(j) = i + 5
              j = j + 1
            END IF
          END IF
        END IF
      END DO
      ptr(n+1) = j

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%min_reduction = 0.5
      control%balance = 1.0
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 1
      control%print_level = 0
      DO i = 0, 0
        control%refinement = 1
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_multigrid

! *****************************************************************


     SUBROUTINE test_amd(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0

      ! --------------------------------------
      ! Call AMD with no nested
      ! --------------------------------------

      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****AMD test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 4, 6, 8, 10, 12, 13, 15, 16, 18, 20, 21, 22, 25, 26, &
        27, 29, 31, 32, 34, 36, 38, 39, 40, 42, 44, 45, 47, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 11, 12, &
        13, 14, 14, 14, 15, 17, 18, 18, 16, 17, 19, 20, 21, 21, 20, 22, 23, &
        24, 24, 25, 23, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 31 /)


      control%amd_call = 40
      DO i = 2, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Call AMD with no nested
      ! --------------------------------------

      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****AMD test ', test, '****'
      n = 31
      ne = 1
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1
      ptr(2:n+1) = 2
      row(1:ne) = (/ 2 /)

      control%amd_call = 40
      DO i = 2, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_amd

! *****************************************************************


     SUBROUTINE test_refine(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i, j
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0

      ! --------------------------------------
      ! Test for refinement controls
      ! --------------------------------------

      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 14
      ne = 16
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
      row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

      control%amd_switch1 = 6
      control%amd_call = 2
      control%find_supervariables = .FALSE.
      control%partition_method = 0
      control%refinement = 0
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0 .AND. info%nzsuper==32) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for refinement controls
      ! --------------------------------------

      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 14
      ne = 16
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
      row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

      control%amd_switch1 = 6
      control%amd_call = 2
      control%find_supervariables = .FALSE.
      control%partition_method = 0
      control%refinement = 0
      DO i = 0, 0
        control%print_level = i
        DO j = 1, 7
          control%refinement = j
          CALL nd_order(0,n,ptr,row,perm,control,info)
          IF (info%flag>=0 .AND. info%nzsuper==32) testi_count = testi_count + &
            1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=7) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF





      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 4
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 2
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%refinement_band = -1
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 7
      ne = 12
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
      row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

      control%amd_switch1 = 3
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 5
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 13
      ne = 27
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 7, 12, 16, 19, 21, 22, 23, 24, 25, 26, 27, 28, 28 /)
      row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, &
        7, 7, 8, 9, 10, 11, 12, 13 /)

      control%amd_switch1 = 3
      control%stop_coarsening1 = 2
      control%min_reduction = 0.5
      control%balance = 1.0
      control%partition_method = 1
      control%ml_call = 4
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%print_level = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for full refine_trim coverage
      ! --------------------------------------

      test = 6
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 16
      ne = 22
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 6, 7, 8, 9, 12, 14, 16, 17, 18, 20, 21, 22, 23, &
        23 /)
      row(1:ne) = (/ 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 12, 13, 13, 14, &
        14, 15, 15, 16, 16 /)

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%min_reduction = 0.5
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 1
      control%print_level = 0
      DO i = 0, 0
        control%refinement = 1
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for full refine_trim coverage
      ! --------------------------------------

      test = 7
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 16
      ne = 22
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 14, 15, 16, 18, 19, 20, 22, &
        23, 23 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 10, 11, 12, 13, &
        14, 15, 15, 16, 16 /)

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%min_reduction = 0.5
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 1
      control%print_level = 2
      DO i = 0, 0
        control%refinement = 1
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for full refine_trim coverage
      ! --------------------------------------

      test = 8
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 25
      ne = 44
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 5, 8, 10, 12, 14, 16, 18, 20, 21, 22, 23, 29, 30, &
        30, 32, 34, 36, 38, 39, 41, 43, 44, 45, 45 /)
      row(1:ne) = (/ 2, 3, 4, 5, 4, 6, 7, 5, 8, 9, 10, 7, 13, 8, 13, 9, 13, &
        10, 13, 13, 12, 13, 14, 16, 17, 18, 19, 20, 15, 17, 21, 18, 21, 19, &
        22, 20, 23, 23, 22, 24, 23, 24, 24, 25 /)

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%min_reduction = 0.5
      control%balance = 1.5
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 1
      control%print_level = 2
      DO i = 0, 0
        control%refinement = 1
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for not checking balance
      ! --------------------------------------

      test = 9
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 25
      ne = 44
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 5, 8, 10, 12, 14, 16, 18, 20, 21, 22, 23, 29, 30, &
        30, 32, 34, 36, 38, 39, 41, 43, 44, 45, 45 /)
      row(1:ne) = (/ 2, 3, 4, 5, 4, 6, 7, 5, 8, 9, 10, 7, 13, 8, 13, 9, 13, &
        10, 13, 13, 12, 13, 14, 16, 17, 18, 19, 20, 15, 17, 21, 18, 21, 19, &
        22, 20, 23, 23, 22, 24, 23, 24, 24, 25 /)

      control%amd_switch1 = 4
      control%stop_coarsening1 = 3
      control%min_reduction = 0.5
      control%balance = real(n+1)
      control%partition_method = 1
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 1
      control%print_level = 0
      DO i = 0, 0
        control%refinement = 7
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 10
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 9
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
      row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
        9 /)

      control%amd_switch1 = 4
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 2
      control%max_improve_cycles = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 11
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 9
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 6, 10, 14, 17, 19, 20, 20 /)
      row(1:ne) = (/ 2, 3, 3, 4, 4, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 8, 9, &
        9 /)

      control%amd_switch1 = 4
      control%balance = 10.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 2
      control%max_improve_cycles = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 12
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 7
      ne = 12
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 5, 8, 10, 11, 12, 13, 13 /)
      row(1:ne) = (/ 2, 3, 4, 5, 3, 4, 5, 4, 5, 5, 6, 7 /)

      control%amd_switch1 = 3
      control%balance = 8.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 13
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 31
      ne = 49
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 5, 7, 8, 10, 12, 13, 15, 17, 18, 20, 22, 23, 25, &
        27, 28, 29, 32, 33, 34, 37, 39, 40, 42, 44, 45, 46, 48, 49, 50, 50 /)
      row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 12, 10, 12, 13, &
        14, 15, 13, 15, 16, 18, 19, 16, 19, 17, 19, 20, 21, 22, 22, 23, 23, &
        24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 31 /)

      control%amd_switch1 = 4
      control%balance = 8.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test automatic choice of shift within expand-refine cycle
      ! --------------------------------------

      test = 14
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 29
      ne = 67
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 7, 9,11, 14,15,16,19,22,26,28,31,35,37,40,44,46,&
            51, 55,59,62,64,66,67,67,68,69,69,69 /)
      row(1:ne) = (/ 2,3,4,5,6,7,5,6,5,7,6,7,8,8,8,9,10,11,10,13,14,11,12,13,&
                  14,12,13,13,16,17,14,15,16,17,15,16,16,19,20,17,18,19,20,18,&
               19,19,21,23,24,25,21,22,26,27,22,23,28,26,27,24,25,25,27,29 /)

      control%amd_switch1 = 4
      control%balance = 1.01
      control%partition_method = 1
      control%coarse_partition_method = 2
      control%refinement_band = 0
      control%stop_coarsening1 = 3
      control%stop_coarsening2 = 3
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%print_level = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test automatic choice of shift within expand-refine cycle
      ! --------------------------------------

      test = 15
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Refinement test ', test, '****'
      n = 21
      ne = 34
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1,2,3,5,7,8,10,12,13,15,18,19,21,22,25,27,29,30,32,34,&
                       35,35 /)
      row(1:ne) = (/ 2,4,4,6,5,7,8,7,9,8,10,11,10,12,11,12,13,13,13,14,14,15,&
           16,17,16,18,17,18,18,19,20,20,21,21 /)

      control%amd_switch1 = 4
      control%balance = 1.01
      control%partition_method = 1
      control%coarse_partition_method = 2
      control%refinement_band = 0
      control%stop_coarsening1 = 3
      control%stop_coarsening2 = 3
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%print_level = 2
      DO i = 0, 7
        control%refinement = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_refine




! *****************************************************************


     SUBROUTINE test_misc(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i, j, k
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0
      ! --------------------------------------
      ! Test diagonal submatrix - not diag!!!!
      ! --------------------------------------
      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 10
      ne = 10
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n-2) = (/ (i,i=1,n-2) /)
      ptr(n-1) = ptr(n-2) + 2
      ptr(n:n+1) = ne + 1
      row(1:ne-2) = n
      row(ne-1) = n - 1
      row(ne) = n

      control%amd_switch1 = 2
      control%amd_call = 2
      DO i = 1, 1
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==16) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF




      ! --------------------------------------
      ! Test independent component detection
      ! --------------------------------------
      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 5
      ne = 6
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 3, 4, 6, 7 /)
      row(1:ne) = (/ 2, 1, 4, 3, 5, 4 /)

      control%amd_switch1 = 2
      control%amd_call = 2
      DO i = 0, 2
        control%print_level = i
        control%coarse_partition_method = 2
        CALL nd_order(1,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1

        control%coarse_partition_method = 1
        CALL nd_order(1,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1

      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=6) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF






      ! --------------------------------------
      ! Test for supervariable detection turned off
      ! --------------------------------------

      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 10
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n) = (/ (2*i-1,i=1,n) /)
      ptr(n+1) = ne + 1
      row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
        10 /)

      control%amd_switch1 = 2
      control%amd_call = 2
      control%find_supervariables = .FALSE.
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for full matrix
      ! --------------------------------------

      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 10
      ne = 45
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1
      DO i = 2, n
        ptr(i) = ptr(i-1) + n - i + 1
      END DO
      ptr(n+1) = ne + 1

      j = 1
      DO i = 1, n - 1
        DO k = i + 1, n
          row(j) = k
          j = j + 1
        END DO
      END DO

      control%amd_switch1 = 2
      control%amd_call = 2
      control%find_supervariables = .FALSE.
      DO i = 2, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==90) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test cost function cost2
      ! --------------------------------------

      test = 5
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 10
      ne = 19
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n) = (/ (2*i-1,i=1,n) /)
      ptr(n+1) = ne + 1
      row(1:ne) = (/ 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, &
        10 /)

      control%amd_switch1 = 2
      control%amd_call = 2
      control%find_supervariables = .FALSE.
      control%cost_function = 2
      DO i = 0, 0
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for no separator return, partition 2 larger than amd_switch1 but
      ! partition
      ! 1 equal to amd_switch1
      ! --------------------------------------

      test = 6
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 14
      ne = 16
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 14, 15, 16, 17, 17 /)
      row(1:ne) = (/ 2, 3, 4, 5, 6, 7, 8, 8, 8, 9, 10, 11, 12, 13, 13, 14 /)

      control%amd_switch1 = 6
      control%amd_call = 2
      control%find_supervariables = .FALSE.
      control%partition_method = 0
      DO i = 2, 2
        control%print_level = i
        CALL nd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0 .AND. info%nzsuper==32) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF






      ! --------------------------------------
      ! Test matrix extraction involving end row
      ! --------------------------------------

      test = 7
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 14
      ne = 26
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 6, 8, 10, 13, 17, 19, 23, 24, 26, 27, 27, 27, &
        27 /)
      row(1:ne) = (/ 2, 3, 3, 4, 5, 5, 14, 5, 13, 6, 13, 14, 8, 9, 13, 14, 10, &
        13, 9, 10, 11, 13, 11, 11, 12, 12 /)

      control%amd_switch1 = 7
      control%stop_coarsening1 = 2
      control%min_reduction = 0.5
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 0
      control%print_level = 2
      DO i = 0, 0
        control%refinement = 7
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF



      ! --------------------------------------
      ! Test matrix extraction involving end row
      ! --------------------------------------

      test = 8
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 16
      ne = 31
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1, 3, 6, 8, 10, 13, 14, 16, 21, 25, 27, 29, 31, 31, 31, &
        32, 32 /)
      row(1:ne) = (/ 2, 3, 3, 4, 5, 5, 6, 5, 15, 6, 15, 16, 16, 11, 15, 9, 11, &
        12, 15, 16, 10, 12, 13, 16, 16, 13, 12, 14, 13, 14, 16 /)

      control%amd_switch1 = 7
      control%stop_coarsening1 = 2
      control%min_reduction = 0.5
      control%balance = 1.0
      control%partition_method = 0
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 0
      control%print_level = 2
      DO i = 0, 0
        control%refinement = 7
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test no balanced partition found so switch to multilevel
      ! --------------------------------------

      test = 9
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Misc test ', test, '****'
      n = 18
      ne = 41
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = (/ 1,2,3,8,10,12,15,19,21,24,28,30,33,37,39,40,42,42,42 /)
      row(1:ne) = (/ 2,3,4,5,6,7,8,6,7,7,8,7,9,10,8,9,10,11,10,11,10,12,13,11,&
           12,13,14,13,14,13,15,16,14,15,16,17,16,17,16,17,18 /)

      control%amd_switch1 = 7
      control%stop_coarsening1 = 2
      control%min_reduction = 0.5
      control%balance = 1.0
      control%partition_method = 2
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 0
      control%print_level = 2
      DO i = 0, 0
        control%refinement = 7
        CALL nd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_misc




! *****************************************************************


     SUBROUTINE test_fm(ok)
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: test_count, testi_count, test
      INTEGER :: n, ne, st, i, j, k, an1,an2,swgt,aw1,aw2,aws
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr,row,wgt,work,part
      TYPE (nd_options) :: control, control_orig
      TYPE (nd_inform) :: info

      ok = .TRUE.
      testi_count = 0
      test_count = 0

      ! --------------------------------------
      ! Test fm with empty partition 2
      ! --------------------------------------
      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****FM test ', test, '****'
      n = 11
      ne = 14
      swgt = n
      ALLOCATE (ptr(n),row(ne),wgt(n),work(8*n+swgt),part(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n) = (/ 1,2,4,6,8,9,11,12,13,14,15  /)
      row(1:ne) = (/ 2,3,4,5,6,6,7,8,8,9,9,10,10,11  /)
      wgt(:)=1
      an1 = 3
      an2 = 0
      aw1 = an1
      aw2 = an2
      aws = swgt - aw1-aw2
      part(1:n) = (/ 1,2,4,8,10,11,3,5,6,7,9 /)

      DO i = 1, 1
        control%print_level = i
        CALL nd_refine_fm(n,ne,ptr,row,wgt,swgt,an1,an2,aw1,aw2,aws,part,work,&
          control)
        testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,wgt,work,part,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test fm with empty partition 1
      ! --------------------------------------
      test = 2
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****FM test ', test, '****'
      n = 11
      ne = 14
      swgt = n
      ALLOCATE (ptr(n),row(ne),wgt(n),work(8*n+swgt),part(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n) = (/ 1,2,4,6,8,9,11,12,13,14,15  /)
      row(1:ne) = (/ 2,3,4,5,6,6,7,8,8,9,9,10,10,11  /)
      wgt(:)=1
      an1 = 0
      an2 = 3
      aw1 = an1
      aw2 = an2
      aws = swgt - aw1-aw2
      part(1:n) = (/ 1,2,4,8,10,11,3,5,6,7,9 /)

      DO i = 1, 1
        control%print_level = i
        CALL nd_refine_fm(n,ne,ptr,row,wgt,swgt,an1,an2,aw1,aw2,aws,part,work,&
          control)
        testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,wgt,work,part,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        ok = .FALSE.
        RETURN
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      RETURN

10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      RETURN

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test

    END SUBROUTINE test_fm

    END PROGRAM main


