    PROGRAM main
      USE spral_nestd
      USE hsl_mc69_double
      IMPLICIT NONE


      ! --------------------------------------
      ! Parameters
      ! --------------------------------------
      INTEGER, PARAMETER :: myint1 = kind(1), wp = kind(1.0D+0), &
        myint = kind(1), leniw = 1000000, lena = 1000000, maxn = 10000, &
        maxnz = 1000000, lenirn = 1000000, maxrhs = 5
      REAL (kind=wp), PARAMETER :: zero = 0.0_wp
      ! --------------------------------------
      ! Local variables
      ! --------------------------------------
      INTEGER :: i, n, ne, st, test, jj, k, netot, nf, mc69flag
      INTEGER :: j, kase, lcase, lfact, liw, lkeep, lp, lrhs, ncase, ne1, &
        nres, numed0, numed1, tour, ix, iy
      INTEGER :: err_count, test_count, testi_count
      INTEGER :: irn1(lenirn), iw(leniw), jcn1(lenirn), jcolst(maxn+1), &
        keep(5*maxn+2*maxnz+42), ym11_icntl(10)

      INTEGER, ALLOCATABLE, DIMENSION (:) :: ptr, row, perm, seps
      INTEGER, DIMENSION (:), ALLOCATABLE :: irn, row_out, ptr_out, jcn
      REAL (kind=wp) :: eps, resmax
      REAL (kind=wp) :: aa1(maxnz), w(maxn,maxrhs)
      LOGICAL :: def, llcase, ltest, undef, corr
      CHARACTER (len=8) :: key1
      TYPE (nestd_options) :: control
      TYPE (nestd_options) :: control_orig
      TYPE (nestd_inform) :: info

      ! .. External Functions ..
      REAL (kind=wp) fa14ad, fd15ad
      EXTERNAL fa14ad, fd15ad


      ! --------------------------------------
      ! Test error flag n<1
      ! --------------------------------------
      test = 1
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
      n = 0
      ne = 0
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1) = 1

      control%amd_switch1 = 2
      err_count = 0
      DO i = 0, 1
        control%print_level = i
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag==-3) err_count = err_count + 1
        CALL nestd_order(1,n,ptr,row,perm,control,info,seps)
        IF (info%flag==-3) err_count = err_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20

      IF (err_count/=4) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      END IF

      ! ******** Basic tests ********

      test_count = 0
      ! --------------------------------------
      ! Test diagonal matrix
      ! --------------------------------------
      test = 2
      testi_count = 0
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
      n = 10
      ne = 0
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10

      ptr(1:n+1) = 1

      control%amd_switch1 = 2
      DO i = 0, 0
        control%print_level = i
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        corr = .TRUE.
        DO j = 1, n
          corr = corr .AND. (perm(j)==j)
        END DO
        IF (info%flag>=0 .AND. corr) testi_count = testi_count + 1
        CALL nestd_order(1,n,ptr,row,perm,control,info,seps)
        corr = .TRUE.
        DO j = 1, n
          corr = corr .AND. (perm(j)==j)
        END DO
        IF (info%flag>=0 .AND. corr) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20

      IF (testi_count/=2) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test diagonal submatrix
      ! --------------------------------------
      test = 3
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      DO i = 0, 0
        control%print_level = i
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==16) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expansion of matrix from lower triangle input
      ! --------------------------------------
      test = 4
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expansion of matrix from lower triangle format and that diags are
      ! removed
      ! --------------------------------------

      test = 5
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expansion of matrix from lower and upper triangle input with no
      ! diags
      ! --------------------------------------
      test = 6
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(1,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expansion of matrix from lower and upper triangle input with diag
      ! entry
      ! --------------------------------------
      test = 7
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(1,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test dense row removal - row 1 dense
      ! --------------------------------------
      test = 8
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==799 .AND. info%nzsuper==0) &
          testi_count = testi_count + 1
        control%remove_dense_rows = .FALSE.
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==800 .AND. info%nzsuper==1598) &
          testi_count = testi_count + 1
        WRITE (*,*) 'i', i
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=6) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test dense row removal - first row dense
      ! --------------------------------------
      test = 9
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, ' ****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==799 .AND. info%nzsuper==1596) &
          testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test dense row removal - last row dense
      ! --------------------------------------
      test = 10
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==799 .AND. info%nzsuper==0) &
          testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test dense row removal - first two rows dense but row 2 has max degree
      ! --------------------------------------
      test = 11
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nsuper==798 .AND. info%nzsuper==0) &
          testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test independent component detection
      ! --------------------------------------
      test = 12
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(1,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1

        control%coarse_partition_method = 1
        CALL nestd_order(1,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1

        control%coarse_partition_method = 3
        CALL nestd_order(1,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=9) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 13
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 14
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 15
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test Ashcraft method
      ! --------------------------------------

      test = 16
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test one-sided levelset method
      ! --------------------------------------

      test = 17
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test one-sided levelset method
      ! --------------------------------------

      test = 18
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM-style refinement with 2-sided partition
      ! --------------------------------------

      test = 19
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
        control%refinement = 5
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
        control%print_level = i
        control%refinement = 2
        CALL nestd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1
        control%refinement = 5
        CALL nestd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=4) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition
      ! --------------------------------------

      test = 20
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 1-sided partition
      ! --------------------------------------

      test = 21
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition
      ! --------------------------------------

      test = 22
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 1-sided partition
      ! --------------------------------------

      test = 23
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test DM refinement with 1-sided partition with balanced partition
      ! --------------------------------------

      test = 24
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 25
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 26
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 27
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test DM refinement with 2-sided partition with balanced partition
      ! --------------------------------------

      test = 28
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid, 2-sided partitioning
      ! --------------------------------------

      test = 29
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%partition_method = 2
      control%coarse_partition_method = 2
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        control%print_level = i
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test multigrid, 1-sided partitioning
      ! --------------------------------------

      test = 30
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%partition_method = 2
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        DO j = 0, 2
          control%print_level = i
          control%matching = j
          CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=9) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Call AMD with no nested
      ! --------------------------------------

      test = 29
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Call AMD with no nested
      ! --------------------------------------

      test = 30
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF



      ! --------------------------------------
      ! Tests for AMD coverage
      ! --------------------------------------

      ALLOCATE (irn(500000),jcn(500000),ptr(50000))
      DO n = 10, 50
        ! N is the matrix dimension
        DO nf = 1, n/5
          ! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
            ! Add full rows
            DO k = jj + 1, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
          DO jj = nf + 1, n
            netot = netot + 1
            iw(netot) = jj
            jcn(netot) = jj - 1
          END DO

          ! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF

          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))
          control%amd_switch1 = n - 1
          CALL nestd_order(0,n,ptr_out,row_out(1:ptr_out(n+ &
            1)-1),perm,control,info,seps)
          CALL reset_control(control_orig,control)
          CALL testpr(n,perm)
        END DO
      END DO

      DO n = 10, 100
        ! N is the matrix dimension
        DO nf = 1, n/5
          ! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
            ! Add quasi-dense rows
            DO k = jj + 7, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
          ! Add a full row
          jj = nf + 1
          DO k = jj + 1, n
            netot = netot + 1
            iw(netot) = k
            jcn(netot) = jj
          END DO
          DO k = 1, jj - 1
            netot = netot + 1
            jcn(netot) = k
            iw(netot) = jj
          END DO
          ! Add sparse rows
          DO jj = nf + 2, n
            netot = netot + 1
            iw(netot) = jj
            jcn(netot) = jj - 1
          END DO

          ! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF
          ! call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))
          control%amd_switch1 = n - 1
          CALL nestd_order(0,n,ptr_out,row_out(1:ptr_out(n+ &
            1)-1),perm,control,info,seps)
          CALL reset_control(control_orig,control)
          CALL testpr(n,perm)
        END DO
      END DO

      DO n = 10, 100
        ! N is the matrix dimension
        DO nf = 1, n/5
          ! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
            ! Add quasi-dense rows
            DO k = jj + 1, n - n/4
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO

          ! Add a full rows
          DO jj = nf + 1, nf + 1
            DO k = jj + 1, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
            DO k = 1, jj - 1
              netot = netot + 1
              jcn(netot) = k
              iw(netot) = jj
            END DO
          END DO
          DO jj = nf + 2, n
            netot = netot + 1
            iw(netot) = jj
            jcn(netot) = jj - 1
          END DO

          ! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF
          ! call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))
          control%amd_switch1 = n - 1
          CALL nestd_order(0,n,ptr_out,row_out(1:ptr_out(n+ &
            1)-1),perm,control,info,seps)
          CALL reset_control(control_orig,control)
          CALL testpr(n,perm)
        END DO
      END DO


      n = 1000
      ! N is the matrix dimension
      DO nf = 1, n/4 + 1, n/8
        ! NF is the number of full rows
        netot = 0
        DO jj = 1, nf
          ! Add quasi-dense rows
          DO k = jj + 1, n - 2*jj
            netot = netot + 1
            iw(netot) = k
            jcn(netot) = jj
          END DO
        END DO
        DO jj = nf + 2, n
          netot = netot + 1
          iw(netot) = jj
          jcn(netot) = jj - 2
        END DO

        ! Copy column indices into IW
        DO jj = 1, netot
          iw(netot+jj) = jcn(jj)
        END DO
        ne = netot
        irn(1:ne) = iw(1:ne)
        IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
        ALLOCATE (ptr_out(n+1))

        CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
        IF (mc69flag<0) THEN
          GO TO 70
        END IF
        ! call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
        IF (allocated(perm)) DEALLOCATE (perm)
        IF (allocated(seps)) DEALLOCATE (seps)
        ALLOCATE (perm(n),seps(n))
        control%amd_switch1 = n - 1
        CALL nestd_order(0,n,ptr_out,row_out(1:ptr_out(n+ &
          1)-1),perm,control,info,seps)
        CALL reset_control(control_orig,control)
        CALL testpr(n,perm)
        WRITE (6,'(a10)') 'perm='
        WRITE (6,'(10I6)') (perm(1:min(20,n)))
      END DO

      DO n = 1000, 1000, 1000
        ! N is the matrix dimension
        DO nf = 1, 1
          ! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
            ! Add quasi-dense rows
            DO k = jj + 10, jj + 100
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
            DO k = jj + 110, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
          ! Add sparse rows
          DO jj = 2, 40
            DO k = 98 + jj, 98 + jj
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
          ! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF
          ! call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))


          control%amd_switch1 = n - 1
          CALL nestd_order(0,n,ptr_out,row_out(1:ptr_out(n+ &
            1)-1),perm,control,info,seps)
          CALL reset_control(control_orig,control)
          CALL testpr(n,perm)
        END DO
      END DO

      ltest = .FALSE.

      resmax = zero
      nres = 0

      numed0 = 0
      numed1 = 0

      lrhs = maxn
      lkeep = 5*maxn + 2*maxnz + 42
      keep(:) = 0

      lp = 6

      ! Initialize LLCASE.   Again only used for serious debug.
      llcase = .FALSE.

      ! Generate value of eps
      eps = fd15ad('E')

      CALL fa14id(ix)
      CALL ym11id(ym11_icntl,iy)

      tour = 1

      ! DO ncase = -5, 24
      DO ncase = -5, 24
        lcase = 20
        IF (ncase<0) lcase = 10
        IF (ncase==0) lcase = 11
        IF (ncase==-5) lcase = 6
        IF (ncase==21) lcase = 10
        IF (ncase==22) lcase = 2
        IF (ncase==23) lcase = 1
        IF (ncase==24) lcase = 1
        DO kase = 1, lcase
          control%print_level = 2
          WRITE (lp,'(//A,2I8)') '***** Case *****', ncase, kase

          IF (llcase) THEN
            GO TO 40
          ELSE
            ! Set LLCASE
            ! LTEST  = TOUR.EQ.1 .AND. NCASE.EQ.-4.AND. KASE.EQ.4
            llcase = ltest

            ! Set LIW and LFACT
            liw = leniw
            lfact = lena

            ! jump to reading special cases

            ! Set N, IN
            n = 50*ncase
            IF (ncase<0) THEN
              n = 100
            END IF
            IF (ncase==0) n = 12

            ! Set NE
            CALL fa14bd(ix,n,ne)
            CALL fa14bd(ix,n,ne1)
            ne1 = ne*ne1
            IF (fa14ad(ix,1)>=0.01) ne1 = max(ne1,1+n/2)
            ! IF (fa04ad(1)>=0.01) ne1 = max(ne1,1+n/2)
            ! Generate test matrix
            ym11_icntl(2) = 0
            ym11_icntl(3) = 0
            key1 = 'file'
            CALL ym11ad(n,n,ne1,ne,irn1,aa1,jcolst,iw,ym11_icntl,key1,iy)
            DO j = 1, ne
              aa1(j) = float(j)
            END DO

            def = .FALSE.
            undef = .FALSE.
            IF (ncase>0 .AND. (kase==3 .OR. kase==5)) undef = .TRUE.
            IF (ncase>0 .AND. (kase==6 .OR. kase==8)) def = .TRUE.
            IF (undef) CALL ndef(n,ne,aa1,jcolst,irn1,w,.TRUE.)
            IF (def) CALL pdef(n,ne,aa1,jcolst,irn1,w,.TRUE.)
            ! Input negative definite matrix
            IF (ncase>0 .AND. kase==8) THEN
              DO k = 1, ne
                aa1(k) = -aa1(k)
              END DO
            END IF
            IF (ne==-1) THEN
              GO TO 100
            END IF

            IF (llcase) THEN
              control%print_level = 2
            END IF

            IF (ncase==-3) THEN
              IF (kase==1) control%print_level = 0
            END IF
            IF (ncase==-2) THEN
              jcn1(1) = irn1(2)
              IF (kase==1) control%print_level = 0
            END IF
            IF (ncase==-1) THEN
              jcn1(1) = irn1(2)
            END IF

            IF (ncase==-5) THEN
              IF (kase==1) n = 1
              IF (kase==1 .OR. kase==2) control%print_level = 0
            END IF
            IF (ncase==22 .AND. kase==1) lkeep = 2
            IF (ltest) THEN
              control%print_level = 1
              WRITE (lp,'(A)') 'Testing test ... MC68 (first)'
            END IF

            IF (allocated(perm)) DEALLOCATE (perm)
            IF (allocated(seps)) DEALLOCATE (seps)
            ALLOCATE (perm(n),seps(n))
            IF (ncase>10) control%amd_switch1 = max(50,n/4)
            control%print_level = 0
            CALL nestd_order(0,n,jcolst(1:n+1),irn1(1:jcolst(n+ &
              1)-1),perm,control,info,seps)
            control%print_level = 0
            IF (ncase>20) control%amd_switch1 = 50
            CALL testpr(n,perm)

            IF (ncase==-5 .AND. (kase==1 .OR. kase==2)) &
              control%print_level = -1
            IF (ncase==22 .AND. kase==1) lkeep = 5*maxn + 2*maxnz + 42
            IF (ncase/=22 .OR. kase/=1) THEN
              IF (.TRUE.) THEN
                IF (ncase<=20 .AND. kase<=4) THEN
                  IF (ncase==-5) THEN
                    ! n = 10
                    IF (kase==3) keep(1) = keep(2)
                    IF (kase==4) keep(1) = 0
                  END IF
                  IF (ncase==-4) THEN
                    DO i = 1, n
                      keep(i) = i
                    END DO
                  END IF
                  IF (ltest) THEN
                    control%print_level = 1
                    WRITE (lp,'(A)') 'Testing test ... MC68 (second)'
                  END IF

                  IF (allocated(perm)) DEALLOCATE (perm)
                  IF (allocated(seps)) DEALLOCATE (seps)
                  ALLOCATE (perm(n),seps(n))
                  control%print_level = 0
                  CALL nestd_order(0,n,jcolst(1:n+1),irn1(1:jcolst(n+ &
                    1)-1),perm,control,info,seps)
                  control%print_level = 0
                  CALL testpr(n,perm)

                END IF


              ELSE IF (ncase>0) THEN
                GO TO 30
              END IF
            END IF
          END IF
100       CONTINUE
        END DO
      END DO
      ! End of extra loop

      GO TO 50

30    CONTINUE
      STOP
40    STOP
50    CONTINUE


70    IF (allocated(irn)) DEALLOCATE (irn)
      IF (allocated(ptr)) DEALLOCATE (ptr)
      IF (allocated(perm)) DEALLOCATE (perm)
      IF (allocated(seps)) DEALLOCATE (seps)
      IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
      IF (allocated(row_out)) DEALLOCATE (row_out)


      ! --------------------------------------
      ! Test for supervariable detection turned off
      ! --------------------------------------

      test = 32
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==18) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for full matrix
      ! --------------------------------------

      test = 33
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0 .AND. info%nzsuper==90) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for no separator return, partition 2 larger than amd_switch1 but
      ! partition
      ! 1 equal to amd_switch1
      ! --------------------------------------

      test = 34
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0 .AND. info%nzsuper==32) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for refinement controls
      ! --------------------------------------

      test = 35
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info)
        IF (info%flag>=0 .AND. info%nzsuper==32) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for refinement controls
      ! --------------------------------------

      test = 36
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
          CALL nestd_order(0,n,ptr,row,perm,control,info)
          IF (info%flag>=0 .AND. info%nzsuper==32) testi_count = testi_count + &
            1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=7) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid coarsening
      ! --------------------------------------

      test = 37
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%coarse_partition_method = 2
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        control%print_level = i
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test multigrid coarsening
      ! --------------------------------------

      test = 38
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        DO j = 0, 2
          control%print_level = i
          control%matching = j
          CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=9) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test multigrid coarsening
      ! --------------------------------------

      test = 39
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      control%min_reduction = 0.4
      control%matching = 1
      DO i = 0, 2
        control%print_level = i
        control%matching = j
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=3) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid refinement
      ! --------------------------------------

      test = 40
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 7
        DO j = 0, 2
          control%refinement = i
          control%matching = j
          CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=24) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test no. multigrid - automatic choice refinement
      ! --------------------------------------

      test = 41
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
          CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=24) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid - automatic choice refinement
      ! --------------------------------------

      test = 42
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%stop_coarsening1 = 2
      DO i = 0, 7
        DO j = 0, 2
          control%refinement = i
          control%matching = j
          CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=24) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF



      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 43
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 44
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 45
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 46
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test automatic multilevel-no multilevel check
      ! --------------------------------------

      test = 47
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
      n = 200
      ne = 199
      ALLOCATE (ptr(n+1),row(ne),perm(n),seps(n),STAT=st)
      IF (st/=0) GO TO 10


      ptr(1:n) = (/ (i,i=1,n) /)
      ptr(n+1) = ne + 1
      row(1:ne) = (/ (i,i=2,n) /)

      control%amd_switch1 = 4
      control%balance = 8.0
      control%partition_method = 2
      control%coarse_partition_method = 2
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%refinement = 7
      DO i = 0, 0
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 48
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 49
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test expand-refine cycles
      ! --------------------------------------

      test = 50
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 2
      control%print_level = 2
      DO i = 0, 7
        control%refinement = i
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=8) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF




      ! --------------------------------------
      ! Test matrix extraction involving end row
      ! --------------------------------------

      test = 51
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF



      ! --------------------------------------
      ! Test matrix extraction involving end row
      ! --------------------------------------

      test = 52
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for full galerkin_graph_rap coverage
      ! --------------------------------------

      test = 52
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%coarse_partition_method = 1
      control%amd_call = 2
      control%max_improve_cycles = 1
      control%print_level = 0
      DO i = 0, 0
        control%refinement = 1
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test multigrid, growth partitioning
      ! --------------------------------------

      test = 53
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
      control%partition_method = 2
      control%coarse_partition_method = 3
      control%amd_call = 2
      control%stop_coarsening1 = 5
      DO i = 0, 2
        DO j = 0, 2
          control%print_level = i
          control%matching = j
          CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
          IF (info%flag>=0) testi_count = testi_count + 1
        END DO
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=9) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      ! --------------------------------------
      ! Test for full refine_trim coverage
      ! --------------------------------------

      test = 54
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for full refine_trim coverage
      ! --------------------------------------

      test = 55
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF


      ! --------------------------------------
      ! Test for full refine_trim coverage
      ! --------------------------------------

      test = 56
      WRITE (6,'(a1)') ' '
      WRITE (6,'(a20,i5,a4)') '****Test section = ', test, '****'
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
        CALL nestd_order(0,n,ptr,row,perm,control,info,seps)
        IF (info%flag>=0) testi_count = testi_count + 1
      END DO
      CALL reset_control(control_orig,control)

      DEALLOCATE (ptr,row,perm,seps,STAT=st)
      IF (st/=0) GO TO 20
      IF (testi_count/=1) THEN
        WRITE (6,'(a29,i5)') 'Code failure in test section ', test
        GO TO 80
      ELSE
        test_count = test_count + 1
        testi_count = 0
      END IF

      GO TO 90


10    WRITE (*,'(a,i4)') 'Allocation error during test section ', test
      GO TO 80

20    WRITE (*,'(a,i4)') 'Deallocation error during test section ', test
      GO TO 80

90    WRITE (*,'(a)') 'All tests successfully completed'

80    CONTINUE

    CONTAINS

      SUBROUTINE reset_control(control_orig,control_reset)


        TYPE (nestd_options), INTENT (IN) :: control_orig

        TYPE (nestd_options), INTENT (OUT) :: control_reset

        control_reset%print_level = control_orig%print_level
        control_reset%unit_diagnostics = control_orig%unit_diagnostics
        control_reset%unit_error = control_orig%unit_error
        control_reset%amd_call = control_orig%amd_call
        control_reset%amd_switch1 = control_orig%amd_switch1
        control_reset%amd_switch2 = control_orig%amd_switch2
        control_reset%coarse_partition_method = control_orig% &
          coarse_partition_method
        control_reset%matching = control_orig%matching
        control_reset%coarse_partition_method = control_orig% &
          coarse_partition_method
        control_reset%refinement = control_orig%refinement
        control_reset%refinement_band = control_orig%refinement_band
        control_reset%remove_dense_rows = control_orig%remove_dense_rows
        control_reset%stop_coarsening2 = control_orig%stop_coarsening2
        control_reset%stop_coarsening1 = control_orig%stop_coarsening1
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

    END PROGRAM main


