module spral_nd_maxflow
   use spral_nd_types
   use spral_nd_util
   implicit none

   private
   public :: nd_refine_max_flow

   ! *****************************************************************

   ! This is a version of maxflow with no bandgraph and assumption that
   ! there are no duplicated edges.
   type network
      integer :: nnode ! number of nodes in the network
      integer :: narc ! number of arcs in the network
      integer :: source ! source node = 1
      integer :: sink ! sink node = nnode
      integer, allocatable :: inheads(:) ! for each node u first incoming arc
         ! for u
      integer, allocatable :: outheads(:) ! for each node u first outgoing arc
         ! for u
      integer, allocatable :: firsts(:) ! for each arc e first node for arc e
      integer, allocatable :: seconds(:) ! for each arc e second node for arc e
      integer, allocatable :: capacities(:) ! for each arc e capacity of arc e
      integer, allocatable :: flows(:) ! for each arc e flow through arc e
      integer, allocatable :: nextin(:) ! for each arc e next incoming arc
      integer, allocatable :: nextout(:) ! for each arc e next outgoing arc
   end type network
contains
! ---------------------------------------------------
! nd_refine_max_flow
! ---------------------------------------------------
! Given a partition, trim the partition using blocks to make it minimal
subroutine nd_refine_max_flow(a_n,a_ne,a_ptr,a_row,a_weight,a_n1, &
    a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,options)
  integer, intent(in) :: a_n ! order of matrix
  integer, intent(in) :: a_ne ! number of entries in matrix
  integer, intent(in) :: a_ptr(a_n) ! On input a_ptr(i) contains
  ! position in a_row that entries for column i start.
  integer, intent(in) :: a_row(a_ne) ! On input a_row contains row
  ! indices of the non-zero rows. Diagonal entries have been removed
  ! and the matrix expanded.
  integer, intent(in) :: a_weight(a_n) ! On input a_weight(i) contains
  ! the weight of column i
  integer, intent(inout) :: a_n1 ! Size of partition 1
  integer, intent(inout) :: a_n2 ! Size of partition 2
  integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
  ! hted
  ! size of partitions and separator
  integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end. This is updated to the new
  ! partition
  integer, intent(out) :: work(8) ! Work array
  type (nd_options), intent(in) :: options

  ! ---------------------------------------------
  ! Local variables
  integer :: msglvl
  real(wp) :: cost, balance_tol

  msglvl = 0
  if (options%print_level.eq.1 .and. options%unit_diagnostics.ge.0) &
    msglvl = 1
  if (options%print_level.ge.2 .and. options%unit_diagnostics.ge.0) &
    msglvl = 3


  if (a_n-a_n1-a_n2.gt.1) then
    balance_tol = max(1.0_wp,options%balance)
    call nd_maxflow(a_n,a_ne,a_ptr,a_row,a_weight,options%cost_function,&
      a_n1,a_n2, &
      a_weight_1,a_weight_2,a_weight_sep,partition,balance_tol,msglvl, &
      options%unit_diagnostics, &
      work(1:8),cost)
  end if


end subroutine nd_refine_max_flow


! ---------------------------------------------------
! nd_maxflow
! ---------------------------------------------------
! Given a partition, get better partition using maxflow algorithm
subroutine nd_maxflow(a_n,a_ne,a_ptr,a_row,a_weight,costf,a_n1,a_n2, &
    a_weight_1,a_weight_2,a_weight_sep,partition,balance_tol,msglvl,lp,stats, &
    cost)

  ! Input matrix: a_n, a_ne, a_ptr, a_row
  integer, intent(in) :: a_n ! order of matrix
  integer, intent(in) :: a_ne ! number of entries in matrix (lower and
  ! upper triangle)
  integer, intent(in) :: a_ptr(:) ! On input, a_ptr(i) contains
  ! position in a_row that entries for column i start.
  integer, intent(in) :: a_row(:) ! On input, a_row contains row
  ! indices of the nonzero entries. Diagonal entries have been
  ! removed and the matrix expanded.
  ! At the moment weights are not used at all
  integer, intent(in) :: a_weight(a_n) ! On input, a_weight(i) contains
  ! the weight of column i
  integer, intent(in) :: costf ! Determines which cost function is used
  ! Data on partition a_n1, a_n2, partition ... will be updated
  integer, intent(inout) :: a_n1 ! Size of partition 1 (ie B)
  integer, intent(inout) :: a_n2 ! Size of partition 2 (ie W)
  integer, intent(inout) :: a_weight_1, a_weight_2, a_weight_sep ! Weig
  ! hted
  ! size of partitions and separator
  integer, intent(inout) :: partition(a_n) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end. This is updated to the new
  ! partition.

  ! Parameters balance_tol (for balance) for cost function
  real(wp), intent(in) :: balance_tol
  integer, intent(in) :: msglvl
  integer, intent(in) :: lp
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
  integer, intent(out) :: stats(8)
  real(wp), intent(out) :: cost

  type (network) :: netw

  ! Work arrays
  ! map,mapL,mapR of length a_n
  ! dmapL,dmapR,vwts of length a_ns
  ! mark1,mark2,pred,list of length number of nodes in network that
  ! is bounded by 2*a_ns+2

  ! Eventually these arrays will be allocated higher up the chain and
  ! will be passed as parameters
  integer, allocatable :: map(:), mapl(:), mapr(:)
  integer, allocatable :: dmapl(:), dmapr(:)
  integer, allocatable :: vwts(:)
  integer, allocatable :: sedge(:,:)
  integer, allocatable :: mark1(:), mark2(:), pred(:), list(:)
  real(wp), allocatable :: imb(:)

  ! Local variables
  integer :: a_ns, i, istart_s, j1, k, wtw, wtb, statsr(9), &
    statsl(9)
  integer nedge, matsiz
  real(wp) :: costr, costl

  ! write(9,*) 'Entering maxflow'
  ! write(0,*) 'Entering maxflow'
  allocate (map(a_n),mapl(a_n),mapr(a_n))


  ! Number vertices in separator
  a_ns = a_n - a_n1 - a_n2

  allocate (dmapl(2*a_ns),dmapr(2*a_ns),vwts(a_ns))

  ! Allocate network work arrays.  Length is upper bound
  ! allocate (mark1(2*a_ns+2),mark2(2*a_ns+2),pred(2*a_ns+2), &
  ! list(2*a_ns+2))
  matsiz = max(a_n,2*a_ns+2)
  allocate (mark1(matsiz),mark2(matsiz),pred(2*a_ns+2),list(matsiz))
  allocate (imb(2*a_ns))

  ! Set up map array to define in what partition each vertex lies
  ! At same time set weights for partition (can check with Sue's input)

  call nd_convert_partition_flags(a_n,a_n1,a_n2,partition,1,2,0, &
    map(1:a_n))

  wtb = a_weight_1
  wtw = a_weight_2

  do i = a_n1 + a_n2 + 1, a_n
    k = partition(i)
    vwts(i-a_n1-a_n2) = a_weight(k)
  end do


  ! Count edges to get upper bound on size of sedge array
  nedge = 0
  do k = 1, a_ns
    i = partition(a_n1+a_n2+k)
    j1 = a_ptr(i)
    if (i.eq.a_n) then
      nedge = nedge + a_ne + 1 - a_ptr(i)
    else
      nedge = nedge + a_ptr(i+1) - a_ptr(i)
    end if
  end do
  allocate (sedge(nedge,2))


  ! Generate network graph.  The structure for our maxflow algorithm

  ! Work arrays dmapL and dmapR used to hold isadjsource and isadjsink
  ! Work array  mapL used to hold sep_map
  ! Source is associated with partition B (size a_n1)
  ! Sink is associated with partition W   (size a_n2)
  ! write(9,*) 'Calling mk_network'
  ! write(0,*) 'Calling mk_network'
  call mk_network(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,map,a_ns, &
    msglvl,lp,netw,vwts,wtb,wtw,sedge,mapl,dmapl,dmapr,list,mark1, &
    mark2,imb)
  ! write(9,*) 'Leaving mk_network'
  ! write(0,*) 'Leaving mk_network'


  ! solve a max flow problem to find the two new maps dmapL and dmapR

  call solvemaxflow(netw,dmapl,dmapr,mark1,mark2,pred,list)


  if (msglvl.gt.2) then
    write (lp,*) 'dmapL ...'
    write (lp,'(10I4)') dmapl
    write (lp,*) 'dmapR ...'
    write (lp,'(10I4)') dmapr
  end if

  mapl = map
  mapr = map
  istart_s = a_n1 + a_n2
  do i = 1, a_ns
    mapl(partition(istart_s+i)) = dmapl(i)
    mapr(partition(istart_s+i)) = dmapr(i)
  end do

  if (msglvl.gt.2) then
    write (lp,*) 'mapL ...'
    write (lp,'(10I4)') mapl
    write (lp,*) 'mapR ...'
    write (lp,'(10I4)') mapr
  end if

  ! Use evaluation function to choose best partition from among these
  ! two
  ! Use Sue's weighted code
  call evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,mapl,balance_tol,costf,statsl,costl)
  call evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,mapr,balance_tol,costf,statsr,costr)


  ! Find the better of the two partitions


  if (statsl(9).eq.1 .and. statsr(9).eq.1) then
    if (msglvl.gt.0) write (lp,'(A)') 'both maps are acceptable'
    if (costl.le.costr) then
      map = mapl
      stats = statsl(1:8)
      cost = costl
      if (msglvl.gt.0) write (lp,'(A)') 'left map accepted'
    else
      map = mapr
      stats = statsr(1:8)
      cost = costr
      if (msglvl.gt.0) write (lp,'(A)') 'right map accepted'
    end if
  else if (statsl(9).eq.1) then
    map = mapl
    stats = statsl(1:8)
    cost = costl
    if (msglvl.gt.0) write (lp,'(A)') &
      'right map not acceptable, left map accepted'
  else if (statsr(9).eq.1) then
    map = mapr
    stats = statsr(1:8)
    cost = costr
    if (msglvl.gt.0) write (lp,'(A)') &
      'left map not acceptable, right map accepted'
  else
    if (msglvl.gt.0) write (lp,'(A)') 'NEITHER map acceptable'
    if (costl.le.costr) then
      map = mapl
      stats = statsl(1:8)
      cost = costl
      if (msglvl.gt.0) write (lp,'(A)') 'left map accepted'
    else
      map = mapr
      stats = statsr(1:8)
      cost = costr
      if (msglvl.gt.0) write (lp,'(A)') 'right map accepted'
    end if
  end if
  a_weight_1 = stats(2)
  a_weight_2 = stats(3)
  a_weight_sep = stats(1)

  ! Now update partition
  ! First count number of vertices in each part
  a_n1 = 0
  a_n2 = 0
  a_ns = 0
  do i = 1, a_n
    if (map(i).eq.1) a_n1 = a_n1 + 1
    if (map(i).eq.2) a_n2 = a_n2 + 1
    if (map(i).eq.0) a_ns = a_ns + 1
  end do


  call nd_convert_flags_partition(a_n,a_n1,a_n2,map(1:a_n),1,2, &
    partition(1:a_n))

  deallocate (map)
  deallocate (dmapl,dmapr)
  deallocate (mapl,mapr)
  deallocate (mark1,mark2,pred,list,vwts)
  deallocate (sedge)


end subroutine nd_maxflow

subroutine solvemaxflow(netw,maps1,maps2,mark1,mark2,pred,list)
  ! Find two partitions of a wide separator graph by solving a max flow
  ! problem




  ! output --

  ! mapS1[n_S] -- first map from wide separator to {0,1,2} = {S,B,W}
  ! mapS2[n_S] -- second map from wide separator to {0,1,2} = {S,B,W}


  ! input/output --

  ! network graph

  ! work arrays

  ! mark1, mark2, pred, list of length nnode


  type (network), intent(inout) :: netw
  integer, intent(out) :: maps1(:), maps2(:)

  integer :: mark1(:), mark2(:), pred(:), list(:)

  ! Local variables
  integer ii, lp, narc, nnode, u

  lp = 6


  nnode = netw%nnode
  narc = netw%narc


  ! Find maxflow through the network using the Ford-Fulkerson algorithm

  call findmaxflow(netw,pred,list,mark1,mark2)



  ! Find the two mincuts

  call findmincut(netw,mark1,mark2,list)


  ! Use mark1 and mark2 to generate maps1 and maps2


  maps1 = -1
  maps2 = -1

  do ii = 2, nnode - 1, 2
    u = ii/2
    if (mark1(ii).eq.1) then
      if (mark1(ii+1).eq.1) then
        maps1(u) = 1
      else
        maps1(u) = 0
      end if
    end if
    if (mark1(ii).eq.2) then
      if (mark1(ii+1).eq.2) then
        maps1(u) = 2
      else
        maps1(u) = 0
      end if
    end if
  end do

  do ii = 2, nnode - 1, 2
    u = ii/2
    if (mark2(ii).eq.1) then
      if (mark2(ii+1).eq.1) then
        maps2(u) = 1
      else
        maps2(u) = 0
      end if
    end if
    if (mark2(ii).eq.2) then
      if (mark2(ii+1).eq.2) then
        maps2(u) = 2
      else
        maps2(u) = 0
      end if
    end if
  end do

end subroutine solvemaxflow


subroutine mk_network(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,map,nvtx, &
    msglvl,lp,netw,vwts,wtb,wtw,sedge,sep_map,isadjtosource,isadjtosink, &
    list,mark1,mark2,imb)
  ! Create and return a network structure


  ! Input matrix: a_n, a_ne, a_ptr, a_row
  integer, intent(in) :: a_n ! order of matrix
  integer, intent(in) :: a_ne ! number of entries in matrix (lower and
  ! upper triangle)
  integer, intent(in) :: a_ptr(:) ! On input, a_ptr(i) contains
  ! position in a_row that entries for column i start.
  integer, intent(in) :: a_row(:) ! On input, a_row contains row
  ! indices of the nonzero entries. Diagonal entries have been
  ! removed and the matrix expanded.

  integer, intent(in) :: partition(a_n) ! First a_n1 entries contain
  ! list of (local) indices in partition 1; next a_n2 entries
  ! contain list of (local) entries in partition 2; entries in
  ! separator are listed at the end.
  integer, intent(in) :: map(a_n) ! First a_n1 entries contain
  integer, intent(in) :: msglvl, lp, nvtx
  integer, intent(in) :: vwts(:)
  integer, intent(inout) :: wtb, wtw
  integer, intent(inout) :: a_n1 ! Size of partition 1 (ie B)
  integer, intent(inout) :: a_n2 ! Size of partition 2 (ie W)

  ! Note that we still do the allocations here.  Doing it further up the
  ! call tree would need to allocate much more storage.
  type (network), intent(out) :: netw

  integer i, iarc, ii, jj, narc1, narc2, narc3, narc4, narc, nnode, &
    nedge, u, v, wts
  integer j1, j2, j, k

  ! Work arrays sedge(nedge,2),sep_map(a_n),isAdjToSource(nvtx),
  ! isAdjToSink(nvtx)
  integer :: sedge(:,:)
  integer :: sep_map(:)
  integer :: isadjtosource(:), isadjtosink(:)
  real(wp) :: imb(:)
  integer :: mark1(:), mark2(:), list(:)
  logical :: augcap

  augcap = .true.

  if (msglvl.gt.0) then
    write (lp,'(/A)') '### inside mknetwork()'
    write (lp,'(A,I8,A,I8)') 'nvtx', nvtx
  end if


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
  do k = 1, nvtx
    i = partition(a_n1+a_n2+k)
    sep_map(i) = k
  end do

  ! We could use a single array although the logic for generating it is
  ! marginally more complicated.
  ! For nodes in separator, set isadj as
  ! 1 if only connected to source
  ! 2 if only connected to sink
  ! 3 if connected to source and sink
  ! isadj = 0


  ! Run through nodes in separator S and generate edges
  nedge = 0
  do k = 1, nvtx
    i = partition(a_n1+a_n2+k)
    j1 = a_ptr(i)
    if (i.eq.a_n) then
      j2 = a_ne
    else
      j2 = a_ptr(i+1) - 1
    end if

    ! Run through vertices connected to vertex i
    do jj = j1, j2
      j = a_row(jj)
      ! Find out in which partition node j lies using map array
      if (map(j).eq.1) then
        ! If in partition B add vertex k to AdjToSource
        isadjtosource(k) = 1
      end if
      if (map(j).eq.2) then
        ! If in partition W add vertex k to AdjToSink
        isadjtosink(k) = 1
      end if
      if (map(j).eq.0) then
        ! If in separator add edge accumulating number of edges
        nedge = nedge + 1
        ! Edge has this orientation to emulate matlab code
        sedge(nedge,2) = sep_map(i)
        sedge(nedge,1) = sep_map(j)
      end if
    end do

  end do


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

  do k = 1, nvtx
    if (isadjtosource(k).eq.1) narc2 = narc2 + 1
    if (isadjtosink(k).eq.1) narc3 = narc3 + 1
  end do

  narc4 = 0
  do ii = 1, nedge
    u = sedge(ii,1)
    v = sedge(ii,2)
    ! --- ignore self edges ---
    if (u.eq.v) cycle
    ! --- omit edges with essential vertices ---
    if (isadjtosink(u).eq.1 .and. isadjtosource(u).eq.1) cycle
    if (isadjtosink(v).eq.1 .and. isadjtosource(v).eq.1) cycle
    ! --- omit pairs both adjacent to source ---
    if ((isadjtosource(u).eq.1) .and. (isadjtosource(v).eq.1)) cycle
    ! --- omit pairs both adjacent to sink ---
    if ((isadjtosink(u).eq.1) .and. (isadjtosink(v).eq.1)) cycle
    ! Qualifying arc found
    narc4 = narc4 + 1
  end do

  nnode = 2*nvtx + 2
  netw%nnode = nnode
  narc = narc1 + narc2 + narc3 + narc4
  netw%narc = narc
  netw%source = 1
  netw%sink = 2*nvtx + 2

  if (msglvl.gt.0) then
    write (lp,'(I8,A)') narc1, ' internal arcs'
    write (lp,'(I8,A)') narc2, ' arcs from source'
    write (lp,'(I8,A)') narc3, ' arcs from sink'
    write (lp,'(I8,A)') narc4, ' edge arcs'
    write (lp,'(I8,A)') narc, ' total arcs'
  end if


  ! create the arc arrays



  ! Allocations done here but could be easily moved up the path although
  ! values very dependent on separator size.
  allocate (netw%firsts(narc),netw%seconds(narc),netw%capacities(narc))
  allocate (netw%flows(narc))
  allocate (netw%inheads(nnode),netw%outheads(nnode),netw%nextin(narc), &
    netw%nextout(narc))

  netw%firsts = -1
  netw%seconds = -1


  ! (u-,u+) arcs first

  iarc = 0
  do u = 1, nvtx
    iarc = iarc + 1
    netw%firsts(iarc) = 2*u
    netw%seconds(iarc) = 2*u + 1
    ! We set capacities after computing imbalance penalty
    ! netw%capacities(iarc) = vwts(u)
  end do

  if (msglvl.gt.0) write (lp,'(A,I8)') 'after (u-,u+) arcs, iarc = ', iarc


  ! (source,u) arcs

  do u = 1, nvtx
    if (isadjtosource(u).eq.1) then
      iarc = iarc + 1
      netw%firsts(iarc) = netw%source
      netw%seconds(iarc) = 2*u
      netw%capacities(iarc) = huge(1)/2
    end if
  end do

  if (msglvl.gt.0) write (lp,'(A,I8)') 'after (source,u-) arcs, iarc = ', &
    iarc


  ! (u,sink) arcs

  do u = 1, nvtx
    if (msglvl.gt.5) write (lp,'(A,I4,A,I8)') 'isAdjToSink(', u, ')= ', &
      isadjtosink(u)
    if (isadjtosink(u).eq.1) then
      iarc = iarc + 1
      netw%firsts(iarc) = 2*u + 1
      netw%seconds(iarc) = netw%sink
      netw%capacities(iarc) = huge(1)/2
    end if
  end do

  if (msglvl.gt.0) write (lp,'(A,I8)') 'after (u,sink) arcs, iarc = ', iarc


  ! (u+,v-) arcs

  do ii = 1, nedge
    u = sedge(ii,1)
    v = sedge(ii,2)
    if ((u.ne.v) .and. (isadjtosource(u).ne.1 .or. isadjtosource( &
        v).ne.1) .and. (isadjtosink(u).ne.1 .or. isadjtosink( &
        v).ne.1) .and. (isadjtosource(u).ne.1 .or. isadjtosink( &
        u).ne.1) .and. (isadjtosource(v).ne.1 .or. isadjtosink(v).ne.1)) then
      iarc = iarc + 1
      netw%firsts(iarc) = 2*u + 1
      netw%seconds(iarc) = 2*v
      netw%capacities(iarc) = huge(1)/2
    end if
  end do

  if (msglvl.gt.0) write (lp,'(A,I8)') 'after (u+,v-) arcs, iarc = ', iarc


  ! Generate the head vectors for in/out edges
  ! and the in/out link vectors for the arcs


  netw%inheads = -1
  netw%outheads = -1
  netw%nextin = -1
  netw%nextout = -1
  do ii = narc, 1, -1
    u = netw%firsts(ii)
    v = netw%seconds(ii)
    if (msglvl.gt.1) write (lp,'(A,I8,A,I8,A,I8,A)') 'ii', ii, 'arc (', u, &
      ',', v, ')'
    netw%nextin(ii) = netw%inheads(v)
    netw%inheads(v) = ii
    netw%nextout(ii) = netw%outheads(u)
    netw%outheads(u) = ii
  end do

  if (augcap) then
    ! Generate wtS
    wts = 0
    do i = 1, nvtx
      wts = wts + vwts(i)
    end do
    if (msglvl.gt.0) write (lp,*) 'Calling findpenalty'
    ! Compute network capacities for separator arcs.
    call findpenalty(msglvl,lp,a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
      map,vwts,wtb,wtw,wts,isadjtosource,isadjtosink,mark1,mark2,list, &
      imb)
    if (msglvl.gt.0) write (lp,*) 'Exiting findpenalty'
    do i = 1, nvtx
      netw%capacities(i) = mark1(i)
    end do

  else
    do i = 1, nvtx
      netw%capacities(i) = vwts(i)
    end do
  end if

  ! ISD Would this not be initialized in routine that computes flows?
  netw%flows = 0

  if (msglvl.gt.0) write (lp,'(A/)') '### leaving mknetwork()'

end subroutine mk_network


subroutine findmaxflow(netw,pred,list,tags,deltas)
  ! Find a maximum flow through the network


  type (network), intent(inout) :: netw

  integer avail, iarc, lp, nnode, sink, source, stats(2), tag

  ! Work arrays ... all of length nnode

  integer :: pred(:), list(:), tags(:), deltas(:)

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
  do
    if (iarc.eq.-1) exit

    do
      avail = netw%capacities(iarc) - netw%flows(iarc)
      if (avail.gt.0) then
        ! Set in findaugpath

        ! use BFS to find path

        ! avail is an output from findaugpath giving available flow on
        ! augmenting
        ! path.  Not to be confused with dummy avail above. The
        ! augmenting
        ! path is given through the array pred.

        call findaugpath(netw,iarc,avail,pred,stats,list,tags,deltas)


        ! Go to next arc from source node if no augmenting path has been
        ! found
        if ((avail.eq.0) .or. (pred(sink).eq.0)) exit

        ! Update flows
        call augmentpath(netw,avail,pred)

      else
        exit
      end if
    end do

    iarc = netw%nextout(iarc)
    tag = tag + 1

  end do

end subroutine findmaxflow

subroutine findmincut(netw,mark1,mark2,list)

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

  integer, intent(out) :: mark1(:), mark2(:), list(:)
  type (network), intent(inout) :: netw

  ! Local variables
  integer iarc, last, lp, nnode, now, sink, source, x, z

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
  do
    if (now.gt.last) exit
    x = list(now)
    now = now + 1
    iarc = netw%outheads(x)
    ! while iarc ~= -1
    ! Run through all arcs starting at node x and putting node at end of
    ! arc
    ! on list if there is spare capacity on arc
    do
      if (iarc.eq.-1) exit
      z = netw%seconds(iarc)
      if (mark1(z).eq.1) then
      else
        if (netw%flows(iarc).lt.netw%capacities(iarc)) then
          last = last + 1
          list(last) = z
          mark1(z) = 1
        end if
      end if
      iarc = netw%nextout(iarc)
    end do
    iarc = netw%inheads(x)
    ! while iarc ~= -1
    ! Run through all arcs terminating at node x and putting node at
    ! start of arc
    ! on list if there is spare capacity on arc
    do
      if (iarc.eq.-1) exit
      z = netw%firsts(iarc)
      if (mark1(z).eq.1) then
      else
        if (netw%flows(iarc).gt.0) then
          last = last + 1
          list(last) = z
          mark1(z) = 1
        end if
      end if
      iarc = netw%nextin(iarc)
    end do
  end do

  ! breadth first traversal from sink


  mark2 = 1
  mark2(sink) = 2

  list = 0
  now = 1
  last = 1
  list(now) = sink
  ! while now <= last
  do
    if (now.gt.last) exit
    x = list(now)
    now = now + 1
    iarc = netw%outheads(x)
    ! while iarc ~= -1
    ! Run through all arcs starting at node x and putting node at end of
    ! arc
    ! on list if there is spare capacity on arc
    do
      if (iarc.eq.-1) exit
      z = netw%seconds(iarc)
      if (mark2(z).eq.2) then
      else
        if (netw%flows(iarc).gt.0) then
          last = last + 1
          list(last) = z
          mark2(z) = 2
        end if
      end if
      iarc = netw%nextout(iarc)
    end do
    iarc = netw%inheads(x)
    ! while iarc ~= -1
    ! Run through all arcs terminating at node x and putting node at
    ! start of arc
    ! on list if there is spare capacity on arc
    do
      if (iarc.eq.-1) exit
      z = netw%firsts(iarc)
      if (mark2(z).eq.2) then
      else
        if (netw%flows(iarc).lt.netw%capacities(iarc)) then
          last = last + 1
          list(last) = z
          mark2(z) = 2
        end if
      end if
      iarc = netw%nextin(iarc)
    end do
  end do


end subroutine findmincut

subroutine findpenalty(msglvl,lp,a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition, &
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


  integer, intent(in) :: msglvl, lp, a_n, a_ne, a_n1, a_n2
  integer, intent(in) :: a_ptr(:), a_row(:)
  integer, intent(in) :: vwts(:), partition(:), map(:)
  integer, intent(inout) :: count(:), head(:)
  integer, intent(inout) :: wtb, wtw, wts
  integer, intent(out) :: mark(:), mark1(:)

  ! Workspace
  ! list --- to hold list of nodes being searched
  ! length bounded by a_n
  integer :: list(:)
  real(wp) :: imb(:)

  ! Local variables
  integer inode, last, maxl, minl, x, z
  integer i, j1, j2, jj, k, nvtx
  integer begin_lev, end_lev, ilev, penp

  ! Set imbalance weighting
  penp = 100

  ! nvtx is number of vertices in separator
  nvtx = a_n - a_n1 - a_n2

  if (msglvl.gt.0) then
    write (lp,'(A)') ''
    write (lp,'(A)') '### inside findpenalty()'
  end if


  ! Breadth first traversal from source
  ! Source is defined as all vertices in black partition

  if (msglvl.gt.0) write (lp,'(A)') 'breadth first traversal from source'

  mark = 0
  ! Black vertices at level 0
  ! Run through separator vertices
  last = 0
  do k = a_n1 + a_n2 + 1, a_n
    z = partition(k)
    ! Check if adjacent to source
    if (count(k-a_n1-a_n2).eq.1) then
      last = last + 1
      list(last) = z
      mark(z) = 1
    end if
  end do
  end_lev = 0

  ! Each pass through this loop determines all nodes in level set k
  ! a_n is just a dummy end
  do k = 2, a_n
    ! Run through all nodes in the previous level set
    begin_lev = end_lev + 1
    end_lev = last
    do ilev = begin_lev, end_lev
      x = list(ilev)
      if (msglvl.gt.1) write (lp,'(A,I10)') 'Processing vertex', x
      ! Run through vertices connected to vertex x
      j1 = a_ptr(x)
      if (x.eq.a_n) then
        j2 = a_ne
      else
        j2 = a_ptr(x+1) - 1
      end if
      do jj = j1, j2
        z = a_row(jj)
        ! Jump if vertex not in separator
        if (map(z).ne.0) cycle
        ! Jump if vertex visited already
        if (mark(z).ne.0) cycle
        mark(z) = k
        ! write(0,*) 'z,mark(z)',z,mark(z)
        ! Add node z to list
        last = last + 1
        list(last) = z
      end do
    end do
    if (last.eq.end_lev) exit
  end do ! end of processing of nodes on level k


  ! breadth first traversal from sink
  ! Sink is defined as all vertices in white partition

  if (msglvl.gt.0) write (lp,'(A)') 'breadth first traversal from the sink'

  mark1 = 0
  ! Keep sink at level 0
  maxl = -a_n
  minl = a_n

  ! White nodes at level 0
  ! Put all separator vertices connected to sink in list
  last = 0
  do k = a_n1 + a_n2 + 1, a_n
    z = partition(k)
    ! Check if adjacent to source
    if (head(k-a_n1-a_n2).eq.1) then
      last = last + 1
      list(last) = z
      mark1(z) = 1
      ! k = 1
      mark(z) = mark(z) - 1
      minl = min(minl,mark(z))
      maxl = max(maxl,mark(z))
    end if
  end do
  end_lev = 0

  ! Each pass through this loop determines all nodes in level set k
  ! a_n is just a dummy end
  do k = 2, a_n
    ! Run through all nodes in the previous level set
    begin_lev = end_lev + 1
    end_lev = last
    do ilev = begin_lev, end_lev
      x = list(ilev)
      if (msglvl.gt.1) write (lp,'(A,I10)') 'Processing vertex', x
      ! Run through vertices connected to vertex x
      j1 = a_ptr(x)
      if (x.eq.a_n) then
        j2 = a_ne
      else
        j2 = a_ptr(x+1) - 1
      end if
      do jj = j1, j2
        z = a_row(jj)
        ! Jump if vertex not in separator
        if (map(z).ne.0) cycle
        ! Jump if vertex visited already
        if (mark1(z).ne.0) cycle
        mark1(z) = k
        ! write(0,*) 'z,mark1(z)',z,mark1(z)
        ! write(0,*) 'mark(z)',mark(z)
        mark(z) = mark(z) - k
        ! write(0,*) 'z,mark(z)',z,mark(z)
        minl = min(minl,mark(z))
        maxl = max(maxl,mark(z))
        ! Add node z to list
        last = last + 1
        list(last) = z
      end do
    end do
    if (last.eq.end_lev) exit
  end do ! end of processing of nodes on level k

  ! Compute half-level sets
  if (msglvl.gt.1) write (lp,'(A,2I4)') 'minl, maxl ', minl, maxl

  ! We will number levels from 1 to maxl-minl+1
  ! count will hold total weight of all vertices in half-level set
  ! Nodes in level set are accessed by linked list in list with headers
  ! in head
  count(1:maxl-minl+1) = 0
  head(1:maxl-minl+1) = -1

  ! Map the mark values into the local coordinates where separator
  ! vertices numbered from 1 to nvtx
  do i = 1, nvtx
    k = partition(a_n1+a_n2+i)
    mark1(i) = mark(k)
  end do

  ! Run through all notes in cutset, resetting numbering of level set
  ! putting
  ! them in level set for level and accumulating weight of half-level
  ! set
  do i = 1, nvtx
    mark1(i) = mark1(i) - minl + 1
    list(i) = head(mark1(i))
    head(mark1(i)) = i
    count(mark1(i)) = count(mark1(i)) + vwts(i)
  end do

  if (msglvl.gt.1) then
    write (lp,'(A)') 'Number of vertices in each half-level set'
    do i = 1, maxl - minl + 1
      write (lp,'(2I10)') i, count(i)
    end do
  end if

  ! Run through half-level sets computing imbalances
  ! wtB is weight of B
  ! wtW is set to total weight of rest of network
  wtw = wtw + wts
  if (msglvl.gt.1) write (lp,('(A,3I10)')) 'wtB,wtW,wtS', wtb, wtw, wts
  if (maxl-minl.eq.0) then
    ! Only one level set
    imb(1) = max(real(wtb)/real(wtw),real(wtw)/real(wtb))
  else
    wtw = wtw - count(1)
   ! if (msglvl.gt.0) write (14,'(A)') &
   !   'Half-level set   width   |B|,|W|, imbalance'
    do k = 1, maxl - minl
      wtw = wtw - count(k+1)
      imb(k) = max(real(wtb)/real(wtw),real(wtw)/real(wtb))
   !   if (msglvl.gt.0) write (14,'(I10,4G12.2)') k, count(k), wtb, wtw, &
   !     imb(k)
      wtb = wtb + count(k)
    end do
  end if

  if (msglvl.gt.1) then
    write (lp,'(A)') 'Imbalances'
    do i = 1, maxl - minl
      write (lp,'(I10,G12.2)') i, imb(i)
    end do
  end if

  ! Run through nodes in level set assigning penalty to them
  do k = 1, maxl - minl + 1
    inode = head(k)
    do
      if (inode.eq.-1) exit
      if (msglvl.gt.1) write (lp,('(A,2I4)')) 'level and node', k, inode
      if (k.eq.1) then
        mark(inode) = floor(penp*imb(1)*vwts(inode))
      else
        if (k.eq.maxl-minl+1) mark(inode) = floor(penp*imb(maxl-minl)*vwts &
          (inode))
        if (k.gt.1 .and. k.lt.maxl-minl+1) mark(inode) = floor(penp*min(imb( &
          k),imb(k-1))*vwts(inode))
      end if
      inode = list(inode)
    end do
  end do

  if (msglvl.gt.1) then
    write (lp,'(A)') 'Computed penalties'
    do i = 1, nvtx
      write (lp,'(2I10)') i, mark(i)
    end do
  end if

  if (msglvl.gt.0) write (lp,'(A/)') '### leaving findpenalty()'

end subroutine findpenalty

subroutine augmentpath(netw,delta,pred)

  ! Reset flows on augmenting path


  ! Input
  ! delta   -- increment flow
  ! pred    -- tree predecessor vector, size nnode
  ! msglvl  -- message level

  ! Input/output
  ! network -- network object

  integer, intent(in) :: delta
  type (network), intent(inout) :: netw

  integer iarc, lp, sink, source, v, w

  integer, intent(in) :: pred(:)


  lp = 6

  source = netw%source
  sink = netw%sink

  ! Should set an error flag
  if (delta.le.0 .or. pred(sink).le.0) then
    write (lp,'(A,I4,A,I4)') 'ERROR : delta', delta, ', pred(sink) = ', &
      pred(sink)
    return
  end if


  ! work back from the sink resetting network flows

  w = sink
  ! while w ~= source
  do
    if (w.eq.source) exit
    iarc = pred(w)
    if (netw%firsts(iarc).eq.w) then
      v = netw%seconds(iarc)
      netw%flows(iarc) = netw%flows(iarc) - delta
    else if (netw%seconds(iarc).eq.w) then
      v = netw%firsts(iarc)
      netw%flows(iarc) = netw%flows(iarc) + delta
    end if
    w = v
  end do

end subroutine augmentpath

subroutine findaugpath(netw,iarc_m,avail,pred,stats,list,tags,deltas)

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

  integer, intent(in) :: iarc_m
  integer, intent(out) :: avail, stats(2)
  type (network), intent(in) :: netw
  integer, intent(out) :: pred(:)

  integer iarc, last, lp, nnode, now, root, sink, source, v, w

  integer :: list(:), tags(:), deltas(:)

  integer narc, u, n_nodevisit, n_arcvisit

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
  if (netw%firsts(iarc).ne.source) then
    write (lp,'(A,I4,A)') 'u', u, 'is not adjacent to source'
    return
  end if


  ! check for available capacity

  avail = netw%capacities(iarc) - netw%flows(iarc)
  if (avail.eq.0) return


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
  do
    if (now.gt.last) exit
    v = list(now)
    now = now + 1
    n_nodevisit = n_nodevisit + 1

    iarc = netw%outheads(v)
    ! while iarc ~= -1
    do
      ! Run through all edges emanating from v
      ! First is v^- to v^+
      if (iarc.eq.-1) exit
      w = netw%seconds(iarc)
      n_arcvisit = n_arcvisit + 1

      if (tags(w).ne.root) then
        ! Node w has not yet been visited

        if (netw%capacities(iarc).gt.netw%flows(iarc)) then
          avail = netw%capacities(iarc) - netw%flows(iarc)

          if (avail.gt.deltas(v)) then
            avail = deltas(v)
          end if
          deltas(w) = avail
          pred(w) = iarc
          ! Flag w as being visited
          tags(w) = root

          if (w.eq.sink) exit

          last = last + 1
          list(last) = w

        end if
      end if
      ! Go to next arc from v
      iarc = netw%nextout(iarc)
    end do
    if (w.eq.sink) exit


    iarc = netw%inheads(v)
    ! while iarc ~= -1
    do
      ! Run through all edges coming in to v
      if (iarc.eq.-1) exit
      w = netw%firsts(iarc)
      n_arcvisit = n_arcvisit + 1
      if (tags(w).ne.root) then
        if (netw%flows(iarc).gt.0) then
          if (avail.gt.netw%flows(iarc)) then
            avail = netw%flows(iarc)
          end if
          deltas(w) = avail
          pred(w) = iarc
          tags(w) = root
          last = last + 1
          list(last) = w
        end if
      end if
      iarc = netw%nextin(iarc)
    end do
    ! Don't think you can reach this statement
    if (w.eq.sink) exit
  end do

  ! Flag to show augmenting path not found
  if (w.ne.sink) avail = 0

  stats(1) = n_nodevisit
  stats(2) = n_arcvisit


end subroutine findaugpath

subroutine evalbsw(a_n,a_ne,a_ptr,a_row,a_weight,map,balance_tol,costf,stats, &
    stats10)
  ! Matlab call
  ! function stats = evalBSW ( A, map, balance_tol, beta, msglvl )

  ! stats = EVALBSW ( A, map, balance_tol, beta )

  ! input ---
  ! map[nvtx] -- map from vertices to region
  ! map[u] .eq. 0 --> u in S
  ! map[u] .eq. 1 --> u in B
  ! map[u] .eq. 2 --> u in W
  ! balance_tol --- acceptability parameter
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
  ! acceptable --> balance_tol*min(|B|,|W|) >= max(|B|,|W|)
  ! stats10 -- cost of partition
  ! cost = |S|*(1 + (beta*| |B| - |W| |)/(|B|+|S|+|W|)) ;

  ! created -- 12jan12, cca

  integer, intent(in) :: a_n
  integer, intent(in) :: a_ne
  integer, intent(in) :: map(:), a_ptr(:), a_row(:), a_weight(:)
  real(wp), intent(in) :: balance_tol
  integer, intent(in) :: costf
  integer, intent(out) :: stats(9)
  real(wp), intent(out) :: stats10
  integer minbw, maxbw, nss, nsb, nsw, nbb, nww, nvtx, ns, nb, nw
  integer j, j1, j2, jj, u, v
  real(wp) diffbw, beta
  logical :: imbal

  beta = 0.5_wp
  nvtx = a_n
  stats(1:9) = -1
  ns = 0
  nb = 0
  nw = 0
  do u = 1, nvtx
    if (map(u).eq.0) then
      ns = ns + a_weight(u)
    else if (map(u).eq.1) then
      nb = nb + a_weight(u)
    else if (map(u).eq.2) then
      nw = nw + a_weight(u)
    end if
  end do
  stats(1) = ns
  stats(2) = nb
  stats(3) = nw
  minbw = min(nb,nw)
  maxbw = max(nb,nw)
  diffbw = real(abs(nb-nw))/real(ns+nb+nw)
  if (.false.) then
    nss = 0
    nsb = 0
    nsw = 0
    nbb = 0
    nww = 0
    ! [rows, cols, ents] = find(A) ;
    ! nzA = length(rows) ;
    do j = 1, a_n
      j1 = a_ptr(j)
      if (j.eq.a_n) then
        j2 = a_ne
      else
        j2 = a_ptr(j+1) - 1
      end if
      v = j
      do jj = j1, j2
        u = a_row(jj)
        ! v = cols(ii) ;
        if (map(u).eq.0) then
          if (map(v).eq.0) then
            nss = nss + 1
          else if (map(v).eq.1) then
            nsb = nsb + 1
          else if (map(v).eq.2) then
            nsw = nsw + 1
          end if
        else if (map(u).eq.1) then
          if (map(v).eq.1) then
            nbb = nbb + 1
          end if
        else if (map(u).eq.2) then
          if (map(v).eq.2) then
            nww = nww + 1
          end if
        end if
      end do
    end do
    stats(4) = nss
    stats(5) = nsb
    stats(6) = nsw
    stats(7) = nbb
    stats(8) = nww
  end if
  ! stats[9] -- 1 if acceptable, 0 if not
  ! acceptable --> balance_tol*min(|B|,|W|) >= max(|B|,|W|)
  ! stats10 -- cost of partition
  if (balance_tol*minbw.ge.maxbw) then
    stats(9) = 1
  else
    stats(9) = 0
  end if

  imbal = (balance_tol .le. nb+nw+ns-2)

  call cost_function(nb,nw,ns,nb+nw+ns, &
    balance_tol,imbal,costf,stats10)

end subroutine evalbsw

end module spral_nd_maxflow
