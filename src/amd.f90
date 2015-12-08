MODULE spral_amd
      IMPLICIT NONE

      PRIVATE
      PUBLIC amd_order

      ! ---------------------------------------------------
      ! Precision
      ! ---------------------------------------------------
      INTEGER, PARAMETER :: wp = KIND(1.0D0)

      ! ---------------------------------------------------
      ! Partition flags
      ! ---------------------------------------------------
      INTEGER, PARAMETER :: nd_sep_flag = 1 ! node in separator

      ! ---------------------------------------------------
      ! The main code
      ! ---------------------------------------------------

    CONTAINS


      ! -------------------------------------------------------------------
      ! hamd is an implementation of the halo_amd method by Pellegrini, Roman
      ! and Amestoy. This is a modified version of AMD given to us by Tim Davis
      ! under the terms of the BSD licence.

      ! We use the term Le to denote the set of all supervariables in element
      ! E.
      ! -------------------------------------------------------------------
      SUBROUTINE amd_order(n,ne,lirn,ip,irn,sep,perm,work)
        ! ..
        ! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: lirn, n, ne
        ! ..
        ! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: irn(lirn), ip(n)
        INTEGER, OPTIONAL, INTENT (IN) :: sep(n)
        INTEGER, INTENT (OUT) :: perm(n)
        INTEGER, INTENT (OUT) :: work(7*n)


        ! N must be set to the matrix order.
        ! Restriction:  N .ge. 1
        ! NE is set to number of non-zeros stored

        ! lirn must be set to the length of irn. It is not altered. On input,
        ! the matrix is stored in irn (1..ne).
        ! *** We do not recommend running this algorithm with ***
        ! ***      lirn .LT. NE+1 + N.                      ***
        ! *** Better performance will be obtained if          ***
        ! ***      lirn .GE. NE+1 + N                       ***
        ! *** or better yet                                   ***
        ! ***      lirn .GT. 1.2 *(NE+1)                    ***
        ! Restriction: lirn .GE. NE

        ! irn(1..NE) must be set to  hold the patterns of the rows of
        ! the matrix.  The matrix must be symmetric, and both upper and
        ! lower triangular parts must be present.  The diagonal must not be
        ! present.  Row I is held as follows:
        ! irn(ip(I)...ip(I) + work(len+I) - 1) must hold the list of
        ! column indices for entries in row I (simple
        ! supervariables), excluding the diagonal.  All
        ! supervariables start with one row/column each
        ! (supervariable I is just row I). If work(len+I) is zero on
        ! input, then ip(I) is ignored on input. Note that the
        ! rows need not be in any particular order, and there may
        ! be empty space between the rows.
        ! During execution, the supervariable I experiences fill-in. This
        ! is represented by constructing a list of the elements that cause
        ! fill-in in supervariable I:
        ! IE(ip(i)...ip(I) + perm(I) - 1) is the list of elements
        ! that contain I. This list is kept short by removing
        ! absorbed elements. irn(ip(I)+perm(I)...ip(I)+work(len+I)-1)
        ! is the list of supervariables in I. This list is kept
        ! short by removing nonprincipal variables, and any entry
        ! J that is also contained in at least one of the
        ! elements in the list for I.
        ! When supervariable I is selected as pivot, we create an element E
        ! of the same name (E=I):
        ! IE(ip(E)..ip(E)+work(len+E)-1) is the list of supervariables
        ! in element E.
        ! An element represents the fill-in that occurs when supervariable
        ! I is selected as pivot.
        ! CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
        ! The contents of irn are undefined on output.

        ! ip(i) must be set to the the index in irn of the start of row I, or
        ! be
        ! zero if row I has no off-diagonal entries. During execution,
        ! it is used for both supervariables and elements:
        ! * Principal supervariable I:  index into irn of the
        ! list of supervariable I.  A supervariable
        ! represents one or more rows of the matrix
        ! with identical pattern.
        ! * Non-principal supervariable I:  if I has been absorbed
        ! into another supervariable J, then ip(I) = -J.
        ! That is, J has the same pattern as I.
        ! Note that J might later be absorbed into another
        ! supervariable J2, in which case ip(I) is still -J,
        ! and ip(J) = -J2.
        ! * Unabsorbed element E:  the index into irn of the list
        ! of element E.  Element E is created when
        ! the supervariable of the same name is selected as
        ! the pivot.
        ! * Absorbed element E:  if element E is absorbed into element
        ! E2, then ip(E) = -E2.  This occurs when one of its
        ! variables is eliminated and when the pattern of
        ! E (that is, Le) is found to be a subset of the pattern
        ! of E2 (that is, Le2).  If element E is "null" (it has
        ! no entries outside its pivot block), then ip(E) = 0.

        ! On output, ip holds the assembly tree/forest, which implicitly
        ! represents a pivot order with identical fill-in as the actual
        ! order (via a depth-first search of the tree). If work(nv+I) .GT. 0,
        ! then I represents a node in the assembly tree, and the parent of
        ! I is -ip(I), or zero if I is a root. If work(nv+I)=0, then
        ! (I,-ip(I))
        ! represents an edge in a subtree, the root of which is a node in
        ! the assembly tree.

        ! sep is an array that is used indicate whether a row/column is in the
        ! seperator
        ! or not. If row i is in the seperator, then sep(i) must equal
        ! nd_sep_flag;
        ! otherwise, sep(i) most be set to a value that is not equal to
        ! nd_sep_flag

        ! perm(I) need not be set. See the description of irn above. At the
        ! start of execution, perm(I) is set to zero. For a supervariable,
        ! perm(I) is the number of elements in the list for supervariable
        ! I. For an element, perm(E) is the negation of the position in the
        ! pivot sequence of the supervariable that generated it. perm(I)=0
        ! if I is nonprincipal.
        ! On output perm(1..N) holds the permutation. That is, if K = perm(I),
        ! then row K is the ith pivot row.  Row K of A appears as the
        ! I-th row in the permuted matrix, PAP^T.

        ! CONTROL is of type nd_options and contains control
        ! parameters and must be set by the user.


        ! Local arrays:
        ! ---------------

        ! WORK(NV+1:NV+N) During execution, ABS(WORK(NV+I)) is equal to the
        ! number of rows represented by the principal supervariable I. If I
        ! is a nonprincipal variable, then WORK(NV+I) = 0. Initially,
        ! WORK(NV+I) = 1
        ! for all I.  WORK(NV+I) .LT. 0 signifies that I is a principal
        ! variable
        ! in the pattern Lme of the current pivot element ME. On termination,
        ! WORK(NV+E) holds the true degree of element E at the time it was
        ! created (including the diagonal part).

        ! WORK(LAST+1:LAST+N) In a degree list, work(last+I) is the
        ! supervariable preceding I, or zero if I is the head of the list.
        ! In a hash bucket, work(last+I) is the hash key for I.
        ! work(last+work(head+HASH))
        ! is also used as the head of a hash bucket if work(head+HASH)
        ! contains
        ! a degree list (see HEAD, below).
        ! On output, work(last+1..last+N) holds the permutation (the same as
        ! the
        ! 'PERM' argument in Sparspak). That is, if I = work(last+K), then row
        ! I
        ! is the Kth pivot row.  Row work(last+K) of A is the K-th row in the
        ! permuted matrix, PAP^T.

        ! work(len+I) is initialised to hold the number of entries in row I of
        ! the
        ! matrix, excluding the diagonal.  The contents of work(len+1..N) are
        ! undefined on output.

        ! DEGREE If I is a supervariable and sparse,
        ! then work(degree+I) holds the current approximation of the external
        ! degree of row I (an upper bound). The external degree is the
        ! number of entries in row I, minus ABS(work(nv+I)) (the diagonal
        ! part). The bound is equal to the external degree if perm(I) is
        ! less than or equal to two. We also use the term "external degree"
        ! for elements E to refer to |Le \ Lme|. If I is full in the reduced
        ! matrix, then work(degree+I)=N+1. If I is dense in the reduced
        ! matrix,
        ! then work(degree+I)=N+1+last_approximate_external_deg of I.

        ! work(head+DEG) is used for degree lists.
        ! work(head+DEG) is the first supervariable in a degree list (all
        ! supervariables I in a degree list DEG have the same approximate
        ! degree, namely, DEG = work(degree+I)). If the list DEG is empty then
        ! work(head+DEG) = 0.
        ! During supervariable detection work(head+HASH) also serves as a
        ! pointer to a hash bucket.
        ! If work(head+HASH) .GT. 0, there is a degree list of degree HASH.
        ! The
        ! hash bucket head pointer is work(last+work(head+HASH)).
        ! If work(head+HASH) = 0, then the degree list and hash bucket are
        ! both empty.
        ! If work(head+HASH) .LT. 0, then the degree list is empty, and
        ! -work(head+HASH) is the head of the hash bucket.
        ! After supervariable detection is complete, all hash buckets are
        ! empty, and the (work(last+work(head+HASH)) = 0) condition is
        ! restored for
        ! the non-empty degree lists.

        ! work(denxt+I)  For supervariable I, work(denxt+I) is
        ! the supervariable following I in a link list, or zero if I is
        ! the last in the list. Used for two kinds of lists: degree lists
        ! and hash buckets (a supervariable can be in only one kind of
        ! list at a time). For element E, work(denxt+E) is the number of
        ! variables with dense or full rows in the element E.

        ! work(w+I) The flag array W determines the status
        ! of elements and variables, and the external degree of elements.
        ! For elements:
        ! if work(w+E) = 0, then the element E is absorbed.
        ! if work(w+E) .GE. WFLG, then work(w+E)-WFLG is the size of the set
        ! |Le \ Lme|, in terms of nonzeros (the sum of ABS(work(nv+I))
        ! for each principal variable I that is both in the
        ! pattern of element E and NOT in the pattern of the
        ! current pivot element, ME).
        ! if WFLG .GT. WE(E) .GT. 0, then E is not absorbed and has
        ! not yet been seen in the scan of the element lists in
        ! the computation of |Le\Lme| in loop 150 below.
        ! ***SD: change comment to remove reference to label***
        ! For variables:
        ! during supervariable detection, if work(w+J) .NE. WFLG then J is
        ! not in the pattern of variable I.
        ! The W array is initialized by setting work(w+I) = 1 for all I, and
        ! by
        ! setting WFLG = 2. It is reinitialized if WFLG becomes too large
        ! (to ensure that WFLG+N does not cause integer overflow).

        ! PFREE must be set to the position in irn of the first free variable.
        ! During execution, additional data is placed in irn, and PFREE is
        ! modified so that components  of irn from PFREE are free.
        ! On output, PFREE is set equal to the size of irn that would have
        ! caused no compressions to occur.  If NCMPA is zero, then
        ! PFREE (on output) is less than or equal to lirn, and the space
        ! irn(PFREE+1 ... lirn) was not used. Otherwise, PFREE (on output)
        ! is greater than lirn, and all the memory in irn was used.

        ! Local variables:
        ! ---------------

        ! DEG:        the degree of a variable or element
        ! DEGME:      size (no. of variables), |Lme|, of the current element,
        ! ME (= work(degree+ME))
        ! DEXT:       external degree, |Le \ Lme|, of some element E
        ! DMAX:       largest |Le| seen so far
        ! E:          an element
        ! permME:     the length, perm(ME), of element list of pivotal var.
        ! ELN:        the length, perm(...), of an element list
        ! HASH:       the computed value of the hash function
        ! HMOD:       the hash function is computed modulo HMOD = MAX(1,N-1)
        ! I:          a supervariable
        ! IDUMMY:     loop counter
        ! ILAST:      the entry in a link list preceding I
        ! INEXT:      the entry in a link list following I
        ! IOVFLO:     local copy of ICNTL(5)
        ! J:          a supervariable
        ! JDUMMY:     loop counter
        ! JLAST:      the entry in a link list preceding J
        ! JNEXT:      the entry in a link list, or path, following J
        ! K:          the pivot order of an element or variable
        ! KNT1:       loop counter used during element construction
        ! KNT2:       loop counter used during element construction
        ! KNT3:       loop counter used during element construction
        ! LENJ:       work(len+J)
        ! LN:         length of a supervariable list
        ! MAXMEM:     amount of memory needed for no compressions
        ! ME:         current supervariable being eliminated, and the
        ! current element created by eliminating that
        ! supervariable
        ! MEM:        memory in use assuming no compressions have occurred
        ! MINDEG:     current approximate minimum degree
        ! NCMPA:      counter for the number of times irn was compressed
        ! NEL:        number of pivots selected so far
        ! NEWMEM:     amount of new memory needed for current pivot element
        ! NLEFT:      N-NEL, the number of nonpivotal rows/columns remaining
        ! NRLADU:     counter for the forecast number of reals in matrix
        ! factor
        ! NVI:        the number of variables in a supervariable I (=
        ! work(nv+I))
        ! NVJ:        the number of variables in a supervariable J (=
        ! work(nv+J))
        ! NVPIV:      number of pivots in current element
        ! P:          pointer into lots of things
        ! P1:         ip (i) for some variable i (start of element list)
        ! P2:         ip (i) + perm (i) -  1 for some var. i (end of el. list)
        ! P3:         index of first supervariable in clean list
        ! PJ:         pointer into an element or variable
        ! PDST:       destination pointer, for compression
        ! PEND:       end of memory to compress
        ! PFREE is set to the position in irn of the first free variable.
        ! At the end, PFREE is set equal to the size of irn that would have
        ! caused no compressions to occur.  If NCMPA is zero, then
        ! PFREE (on output) is less than or equal to lirn, and the space
        ! irn(PFREE+1 ... lirn) was not used. Otherwise, PFREE (on output)
        ! is greater than lirn, and all the memory in irn was used.
        ! PME:        pointer into the current element (PME1...PME2)
        ! PME1:       the current element, ME, is stored in irn(PME1...PME2)
        ! PME2:       the end of the current element
        ! PN:         pointer into a "clean" variable, also used to compress
        ! PSRC:       source pointer, for compression
        ! SLENME:     number of variables in variable list of pivotal variable
        ! WE:         work(w+E)
        ! WFLG:       used for flagging the W array.  See description of W.
        ! WNVI:       WFLG-work(nv+I)
        ! X:          either a supervariable or an element

        ! OPS:        counter for forecast number of flops

        ! IDENSE is true if supervariable I is dense

        ! -------------------------------------------------------------------
        ! FUNCTIONS CALLED:
        ! -------------------------------------------------------------------

        ! ====================================================================
        ! INITIALIZATIONS
        ! ====================================================================

        ! ..
        ! .. Local Arrays ..


        ! ..
        ! .. Local Scalars ..
        INTEGER :: nv, last, degree, head, denxt, w, len
        INTEGER deg, degme, dext, dmax, e, permme, eln, hash, hmod, i, idummy, &
          ilast, inext, iovflo, j, jdummy, jlast, jnext, k, knt1, knt2, knt3, &
          lenj, ln, maxmem, me, mem, mindeg, ncmpa, nel, newmem, nleft, nvi, &
          nvj, nvpiv, p, p1, p2, p3, pdst, pend, pj, pme, pme1, pme2, pn, &
          psrc, slenme, we, wflg, wnvi, x, pfree, nosep, l, temp
        ! ..
        ! .. Intrinsic Functions ..
        INTRINSIC ABS, MAX, MIN, MOD
        ! ..
        me = 0
        nosep = 0
        DO i = 1, n
          IF (sep(i)==nd_sep_flag) nosep = nosep + 1
        END DO
        nv = 0
        last = nv + n
        degree = last + n
        head = degree + n
        denxt = head + n
        w = denxt + n
        len = w + n

        dmax = 0
        hmod = MAX(1,n-1)
        iovflo = HUGE(0)
        pfree = ne + 1
        mem = pfree - 1
        maxmem = mem
        mindeg = 1
        ncmpa = 0
        nel = 0
        wflg = 2

        ! Assign len and pfree
        work(len+1:len+n-1) = ip(2:n) - ip(1:n-1)
        work(len+n) = ne + 1 - ip(n)

        ! ----------------------------------------------------------
        ! initialize arrays and eliminate rows with no off-diag. nz.
        ! ----------------------------------------------------------
        work(last+1:last+n) = 0
        work(head+1:head+n) = 0
        work(nv+1:nv+n) = 1
        work(degree+1:degree+n) = work(len+1:len+n)
        DO i = 1, n
          IF (work(degree+i)==0 .AND. sep(i)/=nd_sep_flag) THEN
            nel = nel + 1
            perm(i) = -nel
            ip(i) = 0
            work(w+i) = 0
          ELSE
            work(w+i) = 1
            perm(i) = 0
          END IF
        END DO
        ! ----------------------------------------------------------------
        ! initialize degree lists
        ! ----------------------------------------------------------------
        DO i = 1, n
          deg = work(degree+i)
          IF (deg>0 .AND. sep(i)/=nd_sep_flag) THEN
            ! ----------------------------------------------------------
            ! place i in the degree list corresponding to its degree
            ! or in the dense row list if i is dense
            ! ----------------------------------------------------------
            ! place i in the degree list corresponding to its degree
            inext = work(head+deg)
            IF (inext/=0) work(last+inext) = i
            work(denxt+i) = inext
            work(head+deg) = i
          END IF
        END DO

        DO WHILE (nel<n-nosep)

          ! ==================================================================
          ! GET PIVOT OF MINIMUM APPROXIMATE DEGREE
          ! ==================================================================
          ! -------------------------------------------------------------
          ! find next supervariable for elimination
          ! -------------------------------------------------------------
          DO deg = mindeg, n
            me = work(head+deg)
            IF (me>0) GO TO 10
          END DO
10        mindeg = deg

          ! -------------------------------------------------------------
          ! remove chosen variable from linked list
          ! -------------------------------------------------------------
          inext = work(denxt+me)
          IF (inext/=0) work(last+inext) = 0
          work(head+deg) = inext
          ! -------------------------------------------------------------
          ! me represents the elimination of pivots nel+1 to nel+work(nv+me).
          ! place me itself as the first in this set.  It will be moved
          ! to the nel+work(nv+me) position when the permutation vectors are
          ! computed.
          ! -------------------------------------------------------------
          permme = perm(me)
          perm(me) = -(nel+1)
          nvpiv = work(nv+me)
          nel = nel + nvpiv
          work(denxt+me) = 0

          ! ==================================================================
          ! ==
          ! CONSTRUCT NEW ELEMENT
          ! ==================================================================
          ! ==

          ! -------------------------------------------------------------
          ! At this point, me is the pivotal supervariable.  It will be
          ! converted into the current element.  Scan list of the
          ! pivotal supervariable, me, setting tree pointers and
          ! constructing new list of supervariables for the new element,
          ! me.  p is a pointer to the current position in the old list.
          ! -------------------------------------------------------------

          ! flag the variable "me" as being in the front by negating
          ! work(nv+me)
          work(nv+me) = -nvpiv
          degme = 0
          IF (permme==0) THEN
            ! ----------------------------------------------------------
            ! There are no elements involved.
            ! Construct the new element in place.
            ! ----------------------------------------------------------
            pme1 = ip(me)
            pme2 = pme1 - 1
            DO p = pme1, pme1 + work(len+me) - 1
              i = irn(p)
              nvi = work(nv+i)
              IF (nvi>0) THEN
                ! ----------------------------------------------------
                ! i is a principal variable not yet placed in the
                ! generated element. Store i in new list
                ! ----------------------------------------------------
                degme = degme + nvi
                ! flag i as being in Lme by negating nv (i)
                work(nv+i) = -nvi
                pme2 = pme2 + 1
                irn(pme2) = i

                ! ----------------------------------------------------
                ! remove variable i from degree list.
                ! ----------------------------------------------------
                IF (sep(i)/=nd_sep_flag) THEN
                  ilast = work(last+i)
                  inext = work(denxt+i)
                  IF (inext/=0) work(last+inext) = ilast
                  IF (ilast/=0) THEN
                    work(denxt+ilast) = inext
                  ELSE
                    ! i is at the head of the degree list
                    temp = work(degree+i)
                    work(head+temp) = inext
                  END IF
                END IF
              END IF
            END DO
            ! this element takes no new memory in irn:
            newmem = 0
          ELSE
            ! ----------------------------------------------------------
            ! construct the new element in empty space, irn (pfree ...)
            ! ----------------------------------------------------------
            p = ip(me)
            pme1 = pfree
            slenme = work(len+me) - permme
            DO knt1 = 1, permme
              ! search the elements in me.
              e = irn(p)
              p = p + 1
              pj = ip(e)
              ln = work(len+e)
              ! -------------------------------------------------------
              ! search for different supervariables and add them to the
              ! new list, compressing when necessary.
              ! -------------------------------------------------------
              DO knt2 = 1, ln
                i = irn(pj)
                pj = pj + 1
                nvi = work(nv+i)
                IF (nvi>0) THEN
                  ! -------------------------------------------------
                  ! compress irn, if necessary
                  ! -------------------------------------------------
                  IF (pfree>lirn) THEN
                    ! prepare for compressing irn by adjusting
                    ! pointers and lengths so that the lists being
                    ! searched in the inner and outer loops contain
                    ! only the remaining entries.
                    ! ***** SD: Seperate compression subroutine tried
                    ! but found to be inefficient in comparison ****
                    ip(me) = p
                    work(len+me) = work(len+me) - knt1
                    ! Check if anything left in supervariable ME
                    IF (work(len+me)==0) ip(me) = 0
                    ip(e) = pj
                    work(len+e) = ln - knt2
                    ! Check if anything left in element E
                    IF (work(len+e)==0) ip(e) = 0
                    ncmpa = ncmpa + 1
                    ! store first item in ip
                    ! set first entry to -item
                    DO j = 1, n
                      pn = ip(j)
                      IF (pn>0) THEN
                        ip(j) = irn(pn)
                        irn(pn) = -j
                      END IF
                    END DO

                    ! psrc/pdst point to source/destination
                    pdst = 1
                    psrc = 1
                    pend = pme1 - 1

                    ! while loop:
                    DO idummy = 1, lirn
                      IF (psrc>pend) THEN
                        GO TO 20
                      ELSE
                        ! search for next negative entry
                        j = -irn(psrc)
                        psrc = psrc + 1
                        IF (j>0) THEN
                          irn(pdst) = ip(j)
                          ip(j) = pdst
                          pdst = pdst + 1
                          ! copy from source to destination
                          lenj = work(len+j)
                          DO knt3 = 0, lenj - 2
                            irn(pdst+knt3) = irn(psrc+knt3)
                          END DO
                          pdst = pdst + lenj - 1
                          psrc = psrc + lenj - 1
                        END IF
                      END IF
                    END DO

                    ! move the new partially-constructed element
20                  p1 = pdst
                    DO psrc = pme1, pfree - 1
                      irn(pdst) = irn(psrc)
                      pdst = pdst + 1
                    END DO
                    pme1 = p1
                    pfree = pdst
                    pj = ip(e)
                    p = ip(me)
                  END IF

                  ! -------------------------------------------------
                  ! i is a principal variable not yet placed in Lme
                  ! store i in new list
                  ! -------------------------------------------------
                  degme = degme + nvi
                  ! flag i as being in Lme by negating nv (i)
                  work(nv+i) = -nvi
                  irn(pfree) = i
                  pfree = pfree + 1

                  ! -------------------------------------------------
                  ! remove variable i from degree link list
                  ! -------------------------------------------------
                  IF (sep(i)/=nd_sep_flag) THEN
                    ilast = work(last+i)
                    inext = work(denxt+i)
                    IF (inext/=0) work(last+inext) = ilast
                    IF (ilast/=0) THEN
                      work(denxt+ilast) = inext
                    ELSE
                      ! i is at the head of the degree list
                      temp = work(degree+i)
                      work(head+temp) = inext
                    END IF
                  END IF
                END IF
              END DO

              ! set tree pointer and flag to indicate element e is
              ! absorbed into new element me (the parent of e is me)
              IF (e/=me) THEN
                IF (sep(e)/=nd_sep_flag) THEN
                  ip(e) = -me
                  work(w+e) = 0
                ELSE
                  ip(e) = 0
                  work(w+e) = 0
                END IF
              END IF
            END DO

            ! search the supervariables in me.
            knt1 = permme + 1
            e = me
            pj = p
            ln = slenme

            ! -------------------------------------------------------
            ! search for different supervariables and add them to the
            ! new list, compressing when necessary.
            ! -------------------------------------------------------
            DO knt2 = 1, ln
              i = irn(pj)
              pj = pj + 1
              nvi = work(nv+i)
              IF (nvi>0) THEN
                ! -------------------------------------------------
                ! compress irn, if necessary
                ! -------------------------------------------------
                IF (pfree>lirn) THEN
                  ! prepare for compressing irn by adjusting
                  ! pointers and lengths so that the lists being
                  ! searched in the inner and outer loops contain
                  ! only the remaining entries.
                  ip(me) = p
                  work(len+me) = work(len+me) - knt1
                  ! Check if anything left in supervariable ME
                  IF (work(len+me)==0) ip(me) = 0
                  ip(e) = pj
                  work(len+e) = ln - knt2
                  ! Check if anything left in element E
                  IF (work(len+e)==0) ip(e) = 0
                  ncmpa = ncmpa + 1
                  ! store first item in ip
                  ! set first entry to -item
                  DO j = 1, n
                    pn = ip(j)
                    IF (pn>0) THEN
                      ip(j) = irn(pn)
                      irn(pn) = -j
                    END IF
                  END DO

                  ! psrc/pdst point to source/destination
                  pdst = 1
                  psrc = 1
                  pend = pme1 - 1

                  ! while loop:
                  ! 122              CONTINUE
                  DO idummy = 1, lirn
                    IF (psrc>pend) THEN
                      GO TO 30
                    ELSE
                      ! search for next negative entry
                      j = -irn(psrc)
                      psrc = psrc + 1
                      IF (j>0) THEN
                        irn(pdst) = ip(j)
                        ip(j) = pdst
                        pdst = pdst + 1
                        ! copy from source to destination
                        lenj = work(len+j)
                        DO knt3 = 0, lenj - 2
                          irn(pdst+knt3) = irn(psrc+knt3)
                        END DO
                        pdst = pdst + lenj - 1
                        psrc = psrc + lenj - 1
                      END IF
                    END IF
                  END DO

                  ! move the new partially-constructed element
30                p1 = pdst
                  DO psrc = pme1, pfree - 1
                    irn(pdst) = irn(psrc)
                    pdst = pdst + 1
                  END DO
                  pme1 = p1
                  pfree = pdst
                  pj = ip(e)
                  p = ip(me)
                END IF

                ! -------------------------------------------------
                ! i is a principal variable not yet placed in Lme
                ! store i in new list
                ! -------------------------------------------------
                degme = degme + nvi
                ! flag i as being in Lme by negating nv (i)
                work(nv+i) = -nvi
                irn(pfree) = i
                pfree = pfree + 1

                ! -------------------------------------------------
                ! remove variable i from degree link list
                ! -------------------------------------------------
                IF (sep(i)/=nd_sep_flag) THEN
                  ilast = work(last+i)
                  inext = work(denxt+i)
                  IF (inext/=0) work(last+inext) = ilast
                  IF (ilast/=0) THEN
                    work(denxt+ilast) = inext
                  ELSE
                    ! i is at the head of the degree list
                    temp = work(degree+i)
                    work(head+temp) = inext
                  END IF
                END IF
              END IF
            END DO

            pme2 = pfree - 1
            ! this element takes newmem new memory in irn (possibly zero)
            newmem = pfree - pme1
            mem = mem + newmem
            maxmem = MAX(maxmem,mem)
          END IF

          ! -------------------------------------------------------------
          ! me has now been converted into an element in irn (pme1..pme2)
          ! -------------------------------------------------------------
          ! degme holds the external degree of new element
          work(degree+me) = degme
          ip(me) = pme1
          work(len+me) = pme2 - pme1 + 1

          ! -------------------------------------------------------------
          ! make sure that wflg is not too large.  With the current
          ! value of wflg, wflg+n must not cause integer overflow
          ! -------------------------------------------------------------
          IF (wflg>iovflo-n) THEN
            DO x = 1, n
              IF (work(w+x)/=0) work(w+x) = 1
            END DO
            wflg = 2
          END IF

          ! ==================================================================
          ! ==
          ! COMPUTE (work(w+e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
          ! where G' is the subgraph of G containing just the sparse rows)
          ! ==================================================================
          ! ==
          ! -------------------------------------------------------------
          ! Scan 1:  compute the external degrees of elements touched
          ! with respect to the current element.  That is:
          ! (w (e) - wflg) = |Le \ Lme|
          ! for each element e involving a supervariable in Lme.
          ! The notation Le refers to the pattern (list of
          ! supervariables) of a previous element e, where e is not yet
          ! absorbed, stored in irn (ip (e) + 1 ... ip (e) + irn (ip (e))).
          ! The notation Lme refers to the pattern of the current element
          ! (stored in irn (pme1..pme2)).
          ! -------------------------------------------------------------
          DO pme = pme1, pme2
            i = irn(pme)
            eln = perm(i)
            IF (eln>0) THEN
              ! note that nv (i) has been negated to denote i in Lme:
              nvi = -work(nv+i)
              wnvi = wflg - nvi
              DO p = ip(i), ip(i) + eln - 1
                e = irn(p)
                we = work(w+e)
                IF (we>=wflg) THEN
                  ! unabsorbed element e has been seen in this loop
                  we = we - nvi
                ELSE IF (we/=0) THEN
                  ! e is an unabsorbed element - this is
                  ! the first we have seen e in all of Scan 1
                  we = work(degree+e) + wnvi
                END IF
                work(w+e) = we
              END DO
            END IF
          END DO

          ! ==================================================================
          ! ==
          ! DEGREE UPDATE AND ELEMENT ABSORPTION
          ! ==================================================================
          ! ==

          ! -------------------------------------------------------------
          ! Scan 2:  for each sparse i in Lme, sum up the external degrees
          ! of each Le for the elements e appearing within i, plus the
          ! supervariables in i.  Place i in hash list.
          ! -------------------------------------------------------------

          DO pme = pme1, pme2
            i = irn(pme)
            ! remove absorbed elements from the list for i
            p1 = ip(i)
            p2 = p1 + perm(i) - 1
            pn = p1
            hash = 0
            deg = 0

            ! -------------------------------------------------------
            ! scan the element list associated with supervariable i
            ! -------------------------------------------------------
            DO p = p1, p2
              e = irn(p)
              ! dext = | Le | - | (Le \cap Lme)\D | - work(denxt+e)
              IF (work(w+e)/=0) THEN
                dext = work(w+e) - wflg
                IF (dext>0) THEN
                  deg = deg + dext
                  irn(pn) = e
                  pn = pn + 1
                  hash = hash + e
                ELSE IF (dext==0) THEN
                  ! aggressive absorption: e is not adjacent to me, but
                  ! |Le(G') \ Lme(G')| is 0, so absorb it into me
                  IF (sep(e)/=nd_sep_flag) THEN
                    ip(e) = -me
                    work(w+e) = 0
                  ELSE
                    ip(e) = 0
                    work(w+e) = 0
                  END IF
                END IF
              END IF
            END DO

            ! count the number of elements in i (including me):
            perm(i) = pn - p1 + 1

            ! ----------------------------------------------------------
            ! scan the supervariables in the list associated with i
            ! ----------------------------------------------------------
            p3 = pn
            DO p = p2 + 1, p1 + work(len+i) - 1
              j = irn(p)
              nvj = work(nv+j)
              IF (nvj>0) THEN
                ! j is unabsorbed, and not in Lme.
                ! add to degree and add to new list
                deg = deg + nvj
                irn(pn) = j
                pn = pn + 1
                hash = hash + j
              END IF
            END DO

            ! ----------------------------------------------------------
            ! update the degree and check for mass elimination
            ! ----------------------------------------------------------
            IF (deg==0 .AND. sep(i)/=nd_sep_flag) THEN
              ! -------------------------------------------------------
              ! mass elimination - supervariable i can be eliminated
              ! -------------------------------------------------------
              ip(i) = -me
              nvi = -work(nv+i)
              degme = degme - nvi
              nvpiv = nvpiv + nvi
              nel = nel + nvi
              work(nv+i) = 0
              perm(i) = 0
            ELSE
              ! -------------------------------------------------------
              ! update the upper-bound degree of i
              ! A bound for the new external degree is the old bound plus
              ! the size of the generated element
              ! -------------------------------------------------------

              ! the following degree does not yet include the size
              ! of the current element, which is added later:
              work(degree+i) = MIN(deg,work(degree+i))

              ! -------------------------------------------------------
              ! add me to the list for i
              ! -------------------------------------------------------
              ! move first supervariable to end of list
              irn(pn) = irn(p3)
              ! move first element to end of element part of list
              irn(p3) = irn(p1)
              ! add new element to front of list.
              irn(p1) = me
              ! store the new length of the list in len (i)
              work(len+i) = pn - p1 + 1

              ! -------------------------------------------------------
              ! place in hash bucket.  Save hash key of i in last (i).
              ! -------------------------------------------------------
              hash = ABS(MOD(hash,hmod)) + 1
              j = work(head+hash)
              IF (j<=0) THEN
                ! the degree list is empty, hash head is -j
                work(denxt+i) = -j
                work(head+hash) = -i
              ELSE
                ! degree list is not empty - has j as its head
                ! last is hash head
                work(denxt+i) = work(last+j)
                work(last+j) = i
              END IF
              work(last+i) = hash
            END IF
          END DO
          work(degree+me) = degme

          ! -------------------------------------------------------------
          ! Clear the counter array, w (...), by incrementing wflg.
          ! -------------------------------------------------------------
          dmax = MAX(dmax,degme)
          wflg = wflg + dmax

          ! make sure that wflg+n does not cause integer overflow
          IF (wflg>=iovflo-n) THEN
            DO x = 1, n
              IF (work(w+x)/=0) work(w+x) = 1
            END DO
            wflg = 2
          END IF
          ! at this point, w (1..n) .lt. wflg holds

          ! ==================================================================
          ! ==
          ! SUPERVARIABLE DETECTION
          ! ==================================================================
          ! ==
          DO pme = pme1, pme2
            i = irn(pme)
            IF ((work(nv+i)<0)) THEN
              ! replace i by head of its hash bucket, and set the hash
              ! bucket header to zero

              ! -------------------------------------------------------
              ! examine all hash buckets with 2 or more variables.  We
              ! do this by examing all unique hash keys for super-
              ! variables in the pattern Lme of the current element, me
              ! -------------------------------------------------------
              hash = work(last+i)
              ! let i = head of hash bucket, and empty the hash bucket
              j = work(head+hash)
              IF (j/=0) THEN
                IF (j<0) THEN
                  ! degree list is empty
                  i = -j
                  work(head+hash) = 0
                ELSE
                  ! degree list is not empty, restore last () of head
                  i = work(last+j)
                  work(last+j) = 0
                END IF
                IF (i/=0) THEN

                  ! while loop:
                  DO jdummy = 1, n
                    IF (work(denxt+i)==0) THEN
                      GO TO 70
                    ELSE
                      ! ----------------------------------------------------
                      ! this bucket has one or more variables following i.
                      ! scan all of them to see if i can absorb any entries
                      ! that follow i in hash bucket.  Scatter i into w.
                      ! ----------------------------------------------------
                      ln = work(len+i)
                      eln = perm(i)
                      ! do not flag the first element in the list (me)
                      DO p = ip(i) + 1, ip(i) + ln - 1
                        work(w+irn(p)) = wflg
                      END DO

                      ! ----------------------------------------------------
                      ! scan every other entry j following i in bucket
                      ! ----------------------------------------------------
                      jlast = i
                      j = work(denxt+i)

                      ! while loop:
                      DO idummy = 1, n
                        IF (j==0) THEN
                          GO TO 60
                        ELSE

                          ! -------------------------------------------------
                          ! check if j and i have identical nonzero pattern
                          ! -------------------------------------------------
                          ! jump if i and j do not have same size data
                          ! structure
                          ! jump if i and j do not have same number adj elts
                          IF (work(len+j)==ln .AND. perm(j)==eln .AND. &
                              sep(i)==sep(j)) THEN
                            ! do not flag the first element in the list (me)

                            DO p = ip(j) + 1, ip(j) + ln - 1
                              ! jump if an entry (irn(p)) is in j but not in i
                              IF (work(w+irn(p))/=wflg) GO TO 40
                            END DO

                            ! ------------------------------------------------
                            ! -
                            ! found it!  j can be absorbed into i
                            ! ------------------------------------------------
                            ! -
                            ip(j) = -i
                            ! both nv (i) and nv (j) are negated since they
                            ! are in Lme, and the absolute values of each
                            ! are the number of variables in i and j:
                            work(nv+i) = work(nv+i) + work(nv+j)
                            work(nv+j) = 0
                            perm(j) = 0
                            ! delete j from hash bucket
                            j = work(denxt+j)
                            work(denxt+jlast) = j
                            GO TO 50
                          END IF

                          ! -------------------------------------------------
40                        CONTINUE
                          ! j cannot be absorbed into i
                          ! -------------------------------------------------
                          jlast = j
                          j = work(denxt+j)
                        END IF
50                      CONTINUE
                      END DO

                      ! ----------------------------------------------------
                      ! no more variables can be absorbed into i
                      ! go to next i in bucket and clear flag array
                      ! ----------------------------------------------------
60                    wflg = wflg + 1
                      i = work(denxt+i)
                      IF (i==0) GO TO 70
                    END IF
                  END DO
                END IF
              END IF
            END IF
70          CONTINUE
          END DO

          ! ==================================================================
          ! ==
          ! RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM
          ! ELEMENT
          ! Squeeze out absorbed variables
          ! ==================================================================
          ! ==
          p = pme1
          nleft = n - nel
          DO pme = pme1, pme2
            i = irn(pme)
            nvi = -work(nv+i)
            IF (nvi>0) THEN
              ! i is a principal variable in Lme
              ! restore nv (i) to signify that i is principal
              work(nv+i) = nvi
              ! -------------------------------------------------------
              ! compute the external degree (add size of current elem)
              ! -------------------------------------------------------
              deg = MIN(work(degree+i)+degme-nvi,nleft-nvi)
              work(degree+i) = deg
              ! -------------------------------------------------------
              ! place the supervariable at the head of the degree list
              ! -------------------------------------------------------
              IF (sep(i)/=nd_sep_flag) THEN
                inext = work(head+deg)
                IF (inext/=0) work(last+inext) = i
                work(denxt+i) = inext
                work(last+i) = 0
                work(head+deg) = i
                ! -------------------------------------------------------
                ! save the new degree, and find the minimum degree
                ! -------------------------------------------------------
                mindeg = MIN(mindeg,deg)
              END IF
              ! -------------------------------------------------------
              ! place the supervariable in the element pattern
              ! -------------------------------------------------------
              irn(p) = i
              p = p + 1
            END IF
          END DO

          ! ==================================================================
          ! ===
          ! FINALIZE THE NEW ELEMENT
          ! ==================================================================
          ! ===
          work(nv+me) = nvpiv + degme
          ! nv (me) is now the degree of pivot (including diagonal part)
          ! save the length of the list for the new element me
          work(len+me) = p - pme1
          IF (work(len+me)==0) THEN
            ! there is nothing left of the current pivot element
            ip(me) = 0
            work(w+me) = 0
          END IF
          IF (newmem/=0) THEN
            ! element was not constructed in place: deallocate part
            ! of it (final size is less than or equal to newmem,
            ! since newly nonprincipal variables have been removed).
            pfree = p
            mem = mem - newmem + work(len+me)
          END IF

          ! ==================================================================
          ! ===
          ! END WHILE (selecting pivots)

        END DO
        ! ===================================================================
        ! COMPUTE THE PERMUTATION VECTORS
        ! ===================================================================

        ! ----------------------------------------------------------------
        ! The time taken by the following code is O(n).  At this
        ! point, perm (e) = -k has been done for all elements e,
        ! and perm (i) = 0 has been done for all nonprincipal
        ! variables i.  At this point, there are no principal
        ! supervariables left, and all elements are absorbed.
        ! ----------------------------------------------------------------

        ! ----------------------------------------------------------------
        ! compute the ordering of unordered nonprincipal variables
        ! ----------------------------------------------------------------
        l = n
        DO i = 1, n
          IF (perm(i)==0 .AND. sep(i)/=nd_sep_flag) THEN
            ! ----------------------------------------------------------
            ! i is an un-ordered row.  Traverse the tree from i until
            ! reaching an element, e.  The element, e, was the
            ! principal supervariable of i and all nodes in the path
            ! from i to when e was selected as pivot.
            ! ----------------------------------------------------------
            j = -ip(i)
            ! while (j is a variable) do:
            DO jdummy = 1, n
              IF (perm(j)<0) THEN
                GO TO 80
              ELSE
                j = -ip(j)
              END IF
            END DO
80          e = j
            ! ----------------------------------------------------------
            ! get the current pivot ordering of e
            ! ----------------------------------------------------------
            k = -perm(e)

            ! ----------------------------------------------------------
            ! traverse the path again from i to e, and compress the
            ! path (all nodes point to e).  Path compression allows
            ! this code to compute in O(n) time.  Order the unordered
            ! nodes in the path, and place the element e at the end.
            ! ----------------------------------------------------------
            j = i
            ! while (j is a variable) do:
            DO idummy = 1, n
              IF (perm(j)<0) THEN
                GO TO 90
              ELSE
                jnext = -ip(j)
                ip(j) = -e
                IF (perm(j)==0) THEN
                  ! j is an unordered row
                  perm(j) = k
                  k = k + 1
                END IF
                j = jnext
              END IF
            END DO
            ! leave perm (e) negative, so we know it is an element
90          perm(e) = -k
          ELSE
            IF (sep(i)==nd_sep_flag) THEN
              perm(i) = l
              l = l - 1
            END IF
          END IF
        END DO

        ! ----------------------------------------------------------------
        ! reset the permutation (perm (1..n)) to be positive and ignore the
        ! halo
        ! ----------------------------------------------------------------
        DO i = 1, n
          k = ABS(perm(i))
          ip(k) = i
        END DO
        perm(1:n) = ip(1:n)

        ! ====================================================================
        ! RETURN THE MEMORY USAGE IN irn AND SET INFORMATION ARRAYS
        ! ====================================================================
        ! If maxmem is less than or equal to lirn, then no compressions
        ! occurred, and irn (maxmem+1 ... lirn) was unused.  Otherwise
        ! compressions did occur, and lirn would have had to have been
        ! greater than or equal to maxmem for no compressions to occur.
        ! Return the value of maxmem in the pfree argument.


        pfree = maxmem
      END SUBROUTINE amd_order

END MODULE spral_amd
