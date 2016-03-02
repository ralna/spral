module spral_nd_multilevel
   use spral_nd_maxflow
   use spral_nd_partition
   use spral_nd_refine
   use spral_nd_types
   use spral_nd_util
   implicit none

   private
   public :: multilevel_partition

contains

subroutine multilevel_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, &
      partition, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, options,   &
      info, lwork, work, grids)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, dimension(a_n), intent(in) :: a_ptr
   integer, dimension(a_ne), intent(in) :: a_row
   integer, dimension(a_n), intent(in) :: a_weight
   integer, intent(in) :: sumweight ! sum of entries in a_weight
   integer, dimension(a_n), intent(out) :: partition
   integer, intent(out) :: a_n1 ! number of entries in partition 1
   integer, intent(out) :: a_n2 ! number of entries in partition 2
   integer, intent(out) :: a_weight_1, a_weight_2, a_weight_sep ! Weighted
      ! size of partitions and separator
   type (nd_options), intent(in) :: options
   integer, intent(inout) :: info
   integer, intent(in) :: lwork ! length of work array: must be atleast
      ! 9a_n + sumweight
   integer, intent(out) :: work(lwork) ! work array
   type (nd_multigrid), dimension(max(1,options%stop_coarsening2)), target, &
      intent(inout) :: grids ! the multilevel of graphs (matrices)

   integer :: level
   type(nd_multigrid), pointer :: grid
   integer :: i, j, k, inv1, inv2, ins, st, grid_ne, cexit
   integer :: stop_coarsening1

   info = 0

   call nd_print_diagnostic(1, options, 'Start multilevel_partition:')

   !
   ! construct the grid at this level
   !
   st = 0
   call nd_matrix_construct(grids(1)%graph,a_n,a_n,a_ne,st)
   if (st.lt.0) then
      info = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info, options, ' multilevel_partition')
      return
   end if

   grids(1)%graph%ptr(1:a_n) = a_ptr(1:a_n)
   grids(1)%graph%ptr(a_n+1) = a_ne + 1
   grids(1)%graph%col(1:a_ne) = a_row(1:a_ne)

   do i = 1, a_n
      do j = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr) - 1
         k = a_row(j)
         grids(1)%graph%val(j) = a_weight(i)*a_weight(k)
      end do
   end do

   grids(1)%size = a_n
   grids(1)%level = 1

   call nd_alloc(grids(1)%where, a_n, st)
   if (st.lt.0) then
      info = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info, options, ' multilevel_partition')
      return
   end if

   call nd_alloc(grids(1)%row_wgt, a_n, info)
   if (info.lt.0) then
      info = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info, options, ' multilevel_partition')
      return
   end if

   ! Initialise row weights
   grids(1)%row_wgt(1:a_n) = a_weight(1:a_n)

   !
   ! Build coarse grid hierarchy
   !
   stop_coarsening1 = max(2,options%stop_coarsening1)
   !stop_coarsening1 = max(stop_coarsening1, a_n/30)
   do level = 1, options%stop_coarsening2-1 ! NB we are are creating level+1
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
         call level_print(options%unit_diagnostics, 'size of grid on level ', &
            level, ' is ', real(grids(level)%size,wp))

      ! Test to see if the matrix size too small
      if (grids(level)%size.le.stop_coarsening1) then
         if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
            call level_print(options%unit_diagnostics, &
               'Terminating coarsening as grid%size <= &
               &options%stop_coarsening1 at level ', level)
         exit ! Stop, we're at the coarsest level
      end if

      ! Allocate coarser level
      grids(level+1)%level = level+1

      call coarsen(grids(level), grids(level+1), work(1:grids(level)%size), &
         options, cexit, st)
      if (st.lt.0) then
         info = ND_ERR_MEMORY_ALLOC
         return
      endif
      if(cexit.ne.0) exit ! Stop coarsening and partition
   end do

   ! Perform coarse patitioning (if it fails, keep tring on finer grids)
   do
      grid => grids(level)
      grid_ne = grid%graph%ptr(grid%graph%n+1) - 1
      call nd_coarse_partition(grid%graph%n, grid_ne, grid%graph%ptr, &
         grid%graph%col, grid%row_wgt, sumweight, grid%part_div(1), &
         grid%part_div(2), grid%where, lwork, work, options, info)
      if (info.lt.0) return

      ! check if partition was succesful returned
      if (grid%part_div(1).ne.0 .and. grid%part_div(2).ne.0) exit

      ! Even if not, give up if we're at the coarsest level
      if ( level.eq.1 ) exit

      ! Unlikely to get here because 99.999% of cases caught in full
      ! matrix check above. Follows same procedure as when full matrix found
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
         write (options%unit_diagnostics,'(a,i10,a)') &
            'at level ', grid%level, ' no partition found'

      ! Try partitioning previous grid level
      level = level - 1
      if(level.lt.1) then
         info = ND_ERR_INTERNAL
         return
      endif
   end do

   ! Prolongation
   do level = level, 2, -1
      grid => grids(level)
      call prolong(grids(level-1), grids(level), sumweight, a_weight_1, &
         a_weight_2, a_weight_sep,work(1:9*grids(level-1)%graph%n+sumweight), &
         options)
      if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) &
         call level_print(options%unit_diagnostics, ' after post smoothing ', &
            level-1)
   end do


   !
   ! Convert from flags to partition
   ! 

   inv1 = 1
   inv2 = grids(1)%part_div(1) + 1
   ins = grids(1)%part_div(1) + grids(1)%part_div(2) + 1

   a_weight_1 = 0
   a_weight_2 = 0
   a_weight_sep = 0
   do i = 1, a_n
      select case (grids(1)%where(i))
      case (ND_PART1_FLAG)
         partition(inv1) = i
         inv1 = inv1 + 1
         a_weight_1 = a_weight_1 + a_weight(i)
      case (ND_PART2_FLAG)
         partition(inv2) = i
         inv2 = inv2 + 1
         a_weight_2 = a_weight_2 + a_weight(i)
      case default
         partition(ins) = i
         ins = ins + 1
         a_weight_sep = a_weight_sep + a_weight(i)
      end select
   end do

   a_n1 = grids(1)%part_div(1)
   a_n2 = grids(1)%part_div(2)

   !write (*,'(a)') ' '
   !write (*,'(a)') 'Multilevel partition found'
   !write (*,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, ', a_n2=', a_n2, &
   !   ', a_n_sep=', a_n - a_n1 - a_n2
   !write (*,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', a_weight_1, &
   !   ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
   !   sumweight - a_weight_1 - a_weight_2

   call nd_print_diagnostic(2, options, &
      'multilevel_partition: successful completion' &
      )
end subroutine multilevel_partition

subroutine coarsen(grid, cgrid, work, options, cexit, st)
   type(nd_multigrid), target, intent(inout) :: grid
   type(nd_multigrid), target, intent(inout) :: cgrid
   integer, dimension(3*grid%size), intent(out) :: work
   type(nd_options), intent(in) :: options
   integer, intent(out) :: cexit
   integer, intent(out) :: st

   integer :: cnvtx ! number of vertices (rows) in the coarse matrix
   integer :: cnedge ! number of edge (entries) in the coarse matrix
   type (nd_matrix), pointer :: cgraph ! the coarse graph
   type (nd_matrix), pointer :: graph ! the fine graph
   real(wp) :: reduction ! actual grid reduction factor
   real(wp) :: min_reduction ! min grid reduction factor
   real(wp) :: max_reduction ! max grid reduction factor

   real(wp), parameter :: dense_test = 1.0

   ! Initialize return values
   cexit = 0   ! Keep coarsening
   st = 0      ! No error

   call nd_alloc(grid%match, grid%graph%n, st)
   if(st.ne.0) return
   select case(options%matching)
   case(:ND_MATCH_COMMON_NEIGHBOURS)
      call match_common_neighbours(grid%graph%n, grid%graph%ptr, &
         grid%graph%col, cgrid%size, grid%match, work(1:grid%size))
   case(ND_MATCH_SHEM:)
      call match_sorted_heavy_edge(grid%graph%n, grid%graph%ptr, &
         grid%graph%col, grid%graph%val, cgrid%size, grid%match, &
         work(1:2*grid%size))
   end select

   ! ensure coarse grid quantities are allocated to sufficient size
   cnvtx = cgrid%size
   call nd_alloc(cgrid%where, cnvtx, st)
   if (st.ne.0) return

   call nd_alloc(cgrid%row_wgt, cnvtx, st)
   if (st.ne.0) return

   ! see if the grid reduction is achieved, if not, set the allowed
   ! maximum level to current level and partition this level
   reduction = real(cgrid%size) / grid%size
   min_reduction = max( 0.01_wp,options%min_reduction )
   max_reduction = min( 1.0_wp, max(0.5_wp,options%max_reduction) )
   if (  reduction.gt.max_reduction .or. reduction.lt.min_reduction .or. &
         cgrid%size.lt.4 ) then
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) then
         write (options%unit_diagnostics, '(a,i10,a,f12.4,i4)') &
            'at level ', grid%level, &
            ' further coarsening gives reduction factor', &
            cgrid%size/real(grid%size)
         write (options%unit_diagnostics,'(a,i10)') &
            'current size = ', grid%size
      endif

      ! stop coarsening
      cexit = 2
      return
   end if

   ! restriction ================

   ! form the coarse grid graph and matrix
   ! cmatrix = P^T*matrix = R*matrix
   graph => grid%graph
   cgraph => cgrid%graph

   ! Construct compressed graph
   cgraph%m = cgrid%size
   cgraph%n = cgrid%size
   call nd_alloc(cgraph%ptr, cgraph%n+1, st)
   if (st.ne.0) return
   call compress_matrix(graph%n, graph%ptr, graph%col, graph%val, grid%match, &
      cgraph%ptr, work) ! NB: Only fills cgraph%ptr, to count #entries
   cgraph%ne = cgraph%ptr(cgraph%n+1)-1
   call nd_alloc(cgraph%col, cgraph%ne, st)
   if (st.ne.0) return
   call nd_alloc(cgraph%val, cgraph%ne, st)
   if (st.ne.0) return
   call compress_matrix(graph%n, graph%ptr, graph%col, graph%val, grid%match, &
      cgraph%ptr, work, row_out=cgraph%col, val_out=cgraph%val)

   ! check if matrix is full
   cnedge = cgrid%graph%ptr(cgrid%graph%n+1)-1
   if (real(cnedge)/cgrid%graph%n .ge. dense_test*(cgrid%graph%n-1)) then
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
         write (options%unit_diagnostics,'(a,i10,a)') &
            'at level ', grid%level, ' further coarsening gives full matrix'

      cexit = 3
      return
   end if

   ! row weight cw = R*w
   call compress_vector(grid%size, grid%match, grid%row_wgt, cgrid%row_wgt)

end subroutine coarsen

subroutine compress_vector(n, match, v_in, v_out)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: match
   integer, dimension(n), intent(in) :: v_in
   integer, dimension(*), intent(out) :: v_out

   integer :: i, j, k

   k = 1
   do i = 1, n
      j = match(i)
      if(j.lt.i) cycle
      v_out(k) = v_in(i)
      if(i.ne.j) v_out(k) = v_out(k) + v_in(j)
      k = k + 1
   end do
end subroutine compress_vector

subroutine prolong(grid, cgrid, sumweight, a_weight_1, a_weight_2, &
      a_weight_sep, work, options)
   type(nd_multigrid), target, intent(inout) :: grid
   type(nd_multigrid), target, intent(inout) :: cgrid
   integer, intent(in) :: sumweight
   integer, intent(out) :: a_weight_1
   integer, intent(out) :: a_weight_2
   integer, intent(out) :: a_weight_sep
   integer, dimension(9*grid%graph%n+sumweight), intent(out) :: work
   type(nd_options), intent(in) :: options

   integer, dimension(:), pointer :: fwhere ! partition on fine grid

   integer :: i, j, a_ne, p1, p2, psep
   integer :: work_ptr, partition_ptr

   fwhere => grid%where(1:grid%size)

   ! injection of the order from coarse grid to the
   ! fine grid, since cwhere(i) is the index of the
   ! i-th vertex in the new ordering, the order
   ! of this vertex should be where(i)
   ! grid%where = P*order_on_coarse_grid
   ! here P is a special matrix with only one non-zero entry per row
   grid%part_div(1:2) = 0
   call prolong_vector(grid%size, grid%match, cgrid%where, grid%where)

   do i = 1, grid%size
      if (fwhere(i).eq.ND_PART1_FLAG) then
         grid%part_div(1) = grid%part_div(1) + 1
      else
         if (fwhere(i).eq.ND_PART2_FLAG) then
            grid%part_div(2) = grid%part_div(2) + 1
         end if
      end if
   end do
   a_weight_1 = 0
   a_weight_2 = 0
   a_weight_sep = 0

   ! Set partition
   partition_ptr = 0
   work_ptr = partition_ptr + grid%graph%n
   p1 = 1
   p2 = grid%part_div(1) + 1
   psep = grid%part_div(1) + grid%part_div(2) + 1
   do i = 1, grid%size
      select case (grid%where(i))
      case (ND_PART1_FLAG)
         work(partition_ptr+p1) = i
         a_weight_1 = a_weight_1 + grid%row_wgt(i)
         p1 = p1 + 1
      case (ND_PART2_FLAG)
         work(partition_ptr+p2) = i
         a_weight_2 = a_weight_2 + grid%row_wgt(i)
         p2 = p2 + 1
      case (ND_SEP_FLAG)
         work(partition_ptr+psep) = i
         a_weight_sep = a_weight_sep + grid%row_wgt(i)
         psep = psep + 1
      end select
   end do

   a_ne = grid%graph%ptr(grid%graph%n+1) - 1
   if (a_weight_sep.gt.0) then ! Do not refine if separable graph
      call refine_partition(grid%graph%n, a_ne, grid%graph%ptr, &
         grid%graph%col, grid%row_wgt, sumweight, grid%part_div(1), &
         grid%part_div(2), work(partition_ptr+1:partition_ptr+grid%graph%n), &
         a_weight_1, a_weight_2, a_weight_sep, options, &
         work(work_ptr+1:work_ptr+8*grid%graph%n+sumweight) )
   end if

   do i = 1, grid%part_div(1)
      j = work(partition_ptr+i)
      grid%where(j) = ND_PART1_FLAG
   end do
   do i = grid%part_div(1) + 1, grid%part_div(1) + grid%part_div(2)
      j = work(partition_ptr+i)
      grid%where(j) = ND_PART2_FLAG
   end do
   do i = grid%part_div(1) + grid%part_div(2) + 1, grid%graph%n
      j = work(partition_ptr+i)
      grid%where(j) = ND_SEP_FLAG
   end do
end subroutine prolong

subroutine prolong_vector(n, match, vec_in, vec_out)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: match
   integer, dimension(*), intent(in) :: vec_in
   integer, dimension(n), intent(out) :: vec_out

   integer :: i, j, k

   k = 1
   do i = 1, n
      j = match(i)
      if(j.lt.i) cycle
      vec_out(i) = vec_in(k)
      vec_out(j) = vec_in(k)
      k = k + 1
   end do
end subroutine prolong_vector

! ***************************************************************
! ---------------------------------------------------
! nd_partition matrix
! ---------------------------------------------------
! Partition the matrix and if one (or more) of the generated submatrices
! is
! small enough, apply halo amd

subroutine nd_coarse_partition(a_n,a_ne,a_ptr,a_row,a_weight, &
    sumweight,a_n1,a_n2,where1,lwork,work,options,info)

  integer, intent(in) :: a_n ! dimension of subproblem ND is applied to
  integer, intent(in) :: a_ne ! no. nonzeros of subproblem
  integer, intent(inout) :: a_ptr(a_n) ! On input a_ptr(i) contains
  ! position in a_row that entries for column i start. This is then
  ! used to hold positions for submatrices after partitioning
  integer, intent(inout) :: a_row(a_ne) ! On input a_row contains row
  ! indices of the non-zero rows. Diagonal entries have been removed
  ! and the matrix expanded.This is then used to hold row indices for
  ! submatrices after partitioning
  integer, intent(inout) :: a_weight(a_n) ! On input a_weight(i)
  ! contains
  ! the weight of column i. This is then used to hold the weights for
  ! the submatrices after partitioning.
  integer, intent(in) :: sumweight ! Sum entries in a_weight.
  ! Unchanged.
  integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
  integer, intent(out) :: where1(a_n) ! Computed partition
  integer, intent(in) :: lwork ! .ge. 9*a_n+sumweight
  integer, intent(out) :: work(lwork)
  type (nd_options), intent(in) :: options
  integer, intent(inout) :: info
  ! real(wp), optional, intent(out) :: real_work(a_n)

  ! ---------------------------------------------
  ! Local variables
  integer :: unit_diagnostics ! unit on which to print diagnostics
  logical :: printi, printd, use_multilevel
  integer :: partition_ptr ! pointer into work array
  integer :: work_ptr ! pointer into work array
  integer :: st
  integer :: a_weight_1, a_weight_2, a_weight_sep
  integer, allocatable :: work1(:)
  real(wp) :: dummy, dummy1

  ! ---------------------------------------------
  ! Printing levels
  unit_diagnostics = options%unit_diagnostics
  printi = (options%print_level.eq.1 .and. unit_diagnostics.ge.0)
  printd = (options%print_level.ge.2 .and. unit_diagnostics.ge.0)

  if (printi .or. printd) then
    write (unit_diagnostics,'(a)') ' '
    write (unit_diagnostics,'(a)') 'Start finding a coarse partition'
    write (unit_diagnostics,'(a,i10,a,i10)') 'a_n=', a_n, ', a_ne=', &
      a_ne
  end if

  ! Find the partition
  partition_ptr = 0 ! length a_n
  work_ptr = partition_ptr + a_n ! max length needed 9*a_n+a_ne

  allocate (work1(a_n),stat=st)
  if (st.ne.0) then
    info = ND_ERR_MEMORY_ALLOC
    return
  end if


  select case (options%coarse_partition_method)
  case (:1)
    ! Half level set method
    use_multilevel = .false.
    call nd_half_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
      a_n2,a_weight_1,a_weight_2,a_weight_sep, &
      work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
      options,dummy,dummy1,use_multilevel,info)
    if(info.ne.0) return ! it's all gone horribly wrong

    if (printi .or. printd) then
      write (unit_diagnostics,'(a)') ' '
      write (unit_diagnostics,'(a)') 'Initial half-level set partition'
      write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
        ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
      write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
        a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
        sumweight - a_weight_1 - a_weight_2
    end if

  case (2)
    ! Level set method
    use_multilevel = .false.
    call nd_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
      a_n2,a_weight_1,a_weight_2,a_weight_sep, &
      work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
      options,dummy,dummy1,use_multilevel, info)
    if(info.ne.0) return ! it's all gone horribly wrong

    if (printi .or. printd) then
      write (unit_diagnostics,'(a)') ' '
      write (unit_diagnostics,'(a)') 'Initial level set partition'
      write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
        ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
      write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
        a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
        sumweight - a_weight_1 - a_weight_2
    end if
  case(3:)
    ! Region-growing edge seperator vertex cover
    use_multilevel = .false.
    call region_grow_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1, &
      a_n2,a_weight_1,a_weight_2,a_weight_sep, &
      work1(partition_ptr+1:partition_ptr+a_n),work(1:9*a_n+sumweight), &
      options,dummy,dummy1,use_multilevel, info)
    if(info.ne.0) return ! it's all gone horribly wrong

    if (printi .or. printd) then
      write (unit_diagnostics,'(a)') ' '
      write (unit_diagnostics,'(a)') 'Initial region-based partition'
      write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
        ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
      write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
        a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
        sumweight - a_weight_1 - a_weight_2
    end if

  end select



  if (a_n1.ne.0 .and. a_n2.ne.0 .and. a_n.ge.3) then
    if (a_n1+a_n2.lt.a_n) then
      ! Call standard refinement, but WITHOUT any improvement cycles
      call refine_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, &
         a_n1, a_n2, work1(partition_ptr+1:partition_ptr+a_n), a_weight_1, &
         a_weight_2, a_weight_sep, options, work(1:8*a_n+sumweight), &
         max_improve_cycles=0)
    end if
  else
    go to 10
  end if


  call nd_convert_partition_flags(a_n,a_n1,a_n2, &
    work1(partition_ptr+1:partition_ptr+a_n),ND_PART1_FLAG, &
    ND_PART2_FLAG,ND_SEP_FLAG,where1(1:a_n))

  deallocate (work1,stat=st)
  if (st.ne.0) then
    info = ND_ERR_MEMORY_ALLOC
    return
  end if



  if (printi .or. printd .or. .false.) then
    write (unit_diagnostics,'(a)') ' '
    write (unit_diagnostics,'(a)') 'Initial coarse partition found'
    write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_n1 =', a_n1, &
      ', a_n2=', a_n2, ', a_n_sep=', a_n - a_n1 - a_n2
    write (unit_diagnostics,'(a,i10,a,i10,a,i10)') 'a_weight_1 =', &
      a_weight_1, ', a_weight_2=', a_weight_2, ', a_weight_sep=', &
      sumweight - a_weight_1 - a_weight_2
  end if
  go to 20

10      if (printi .or. printd) then
    write (unit_diagnostics,'(a)') ' '
    write (unit_diagnostics,'(a)') 'No partition found'
  end if

20      info = 0
  if (printi .or. printd) then
    call nd_print_message(info,unit_diagnostics, &
      'nd_coarse_partition')
  end if
  return

end subroutine nd_coarse_partition

!
! Construct data structure for storing sparse matrix.
!
! Arrays in nd_matrix will only be (re)allocated if they are not long
! enough. On exit,
! size(p%val) <-  max(ne, size(p%val)
! size(p%col) <-  max(ne, size(p%col)
! size(p%ptr) <-  max(m+1, size(p%ptr)
!
subroutine nd_matrix_construct(p, m, n, ne, st)
   type (nd_matrix), intent(inout) :: p
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, intent(in) :: ne ! number entries
   integer, intent(out) :: st

   p%m = m
   p%n = n
   p%ne = ne

   call nd_alloc(p%ptr, m+1, st)
   if(st.ne.0) return
   call nd_alloc(p%col, ne, st)
   if(st.ne.0) return
   call nd_alloc(p%val, ne, st)
   if(st.ne.0) return
end subroutine nd_matrix_construct

! ***************************************************************

subroutine nd_alloc(v, n, st)
   integer, intent(inout), allocatable :: v(:)
   integer, intent(in) :: n
   integer, intent(out) :: st

   st = 0 ! By default no allocation error

   ! Check for immediate return if array already correct size
   if (allocated(v)) then
      if(size(v).ge.n) return
   endif

   ! Otherwise, ensure deallocated (ignore any error)
   deallocate (v, stat=st)
   
   ! Allocate to correct size (caller should handle any error)
   allocate (v(n), stat=st)
end subroutine nd_alloc


! *******************************************************************

!
! Find matching using heavy-edge collapsing starting from smallest degree
!
subroutine match_sorted_heavy_edge(n, ptr, row, val, cn, match, work)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(ptr(n+1)-1), intent(in) :: val
   integer, intent(out) :: cn ! number of vertices after matching [ie coarse]
   integer, dimension(n), intent(out) :: match
   integer, dimension(2*n), intent(out) :: work

   ! working variables
   integer :: v, u, j, vv

   ! maximum weight and index of edges connected to the current vertex
   integer :: bestv, best

   ! Sort nodes based on degree
   call construct_degree(n, ptr, work(n+1:2*n))
   call sort_by_degree(n, work(1:n), work(n+1:2*n))

   ! Initialise everything to unmatched
   match(:) = -1

   ! loop over each vertex and match along the heaviest edge
   cn = 0
   do vv = 1, n
      v = work(vv)
      if (match(v).ne.-1) cycle ! v is already matched
      bestv = -huge(0)
      best = v ! Default to unmatched [ie with itself]
      ! Loop over unmatched entries in col v, find one with largest val
      do j = ptr(v), ptr(v+1) - 1
         u = row(j)
         if(match(u).ne.-1) cycle ! u already matched
         ! heavy edge matching
         ! if u is unmatched and value of the entry in col. u is greater
         ! than maxwgt, select u as the matching.
         if (val(j).gt.bestv) then
            bestv = val(j)
            best = u
         end if
      end do
      ! Match v and best
      match(v) = best
      match(best) = v
      ! increase number of vertices in coarse graph by 1
      cn = cn + 1
   end do
end subroutine match_sorted_heavy_edge

! *******************************************************************

subroutine construct_degree(n, ptr, degree)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(n), intent(out) :: degree

   integer :: i, j

   do i = 1, n
      degree(i) = ptr(i+1) - ptr(i)
   end do
end subroutine construct_degree

!
! Returns an order with vertices in order of increasing degree
!
subroutine sort_by_degree(n, order, val)
   integer, intent(in) :: n
   integer, dimension(n), intent(out) :: order
   integer, dimension(n), intent(in) :: val 

   integer :: start, last, u, temp

   !
   ! Heapsort
   !
   ! Our heap has the following properties:
   ! An element must be greater than either of its children (if they exist)

   do u = 1, n
      order(u) = u
   end do

   ! Construct heap with smallest element at top
   start = n / 2 ! parent of n
   do while(start.ge.1)
      ! Sift start down to correct location in heap
      call heap_sift_down(start, n, order, val)
      
      ! Go down to next parent
      start = start - 1
   end do

   ! Keep popping elements off heap and grow ordered list from rear
   last = n
   do while (last.ge.1)
      temp = order(last)
      order(last) = order(1)
      order(1) = temp
      last = last - 1
      call heap_sift_down(1, last, order, val)
   end do
end subroutine sort_by_degree

subroutine heap_sift_down(first, last, order, val)
   integer, intent(in) :: first
   integer, intent(in) :: last
   integer, dimension(*), intent(inout) :: order
   integer, dimension(*), intent(in) :: val

   integer :: root, swapidx, swapval, leftidx, rightidx, leftval, rightval
   integer :: temp

   root = first
   do while(2*root.le.last) ! i.e. while root has any children
      ! Initialise swap to be current root
      swapidx = root
      swapval = val(order(root))

      ! Check if left child is greater than current best
      leftidx = 2*root
      leftval = val(order(leftidx))
      if(leftval.gt.swapval) then
         swapidx = leftidx
         swapval = leftval
      endif

      ! Check if right child is larger than current best
      rightidx = leftidx + 1
      if(rightidx.le.last) then ! right child exists
         rightval = val(order(rightidx))
         if(rightval.gt.swapval) then
            swapidx = rightidx
            swapval = rightval
         endif
      endif

      ! If root is largest element, heap property now holds
      if(swapidx.eq.root) return

      ! Otherwise, swap root and largest child and descend to child
      temp = order(root)
      order(root) = order(swapidx)
      order(swapidx) = temp
      root = swapidx
   end do
end subroutine heap_sift_down

! *******************************************************************

!
! Find matching using common neighbours criterion
!
subroutine match_common_neighbours(n, ptr, row, cn, match, work)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, intent(out) :: cn ! number of vertices in coarse grid
   integer, dimension(n), intent(out) :: match
   integer, dimension(n), intent(out) :: work

   ! working variables
   integer :: v, u, j

   ! maximum no. neighbours and index of edges connected to the current vertex
   integer :: bestv, best, num

   ! initialise the matching status
   match(1:n) = -1
   work(:) = 0

   ! loop over each vertex and match based on number of neighbours in common
   cn = 0
   do v = 1, n
      if (match(v).ne.-1) cycle ! Already matched

      best = v ! Default to unmatched [i.e. with self]

      ! Flag neighbours of v
      work(v) = v
      do j = ptr(v), ptr(v+1) - 1
         u = row(j)
         work(u) = v
      end do

      ! For each unmatched neighbour of v, count the number of
      ! neighbours it has in common with v
      bestv = 0
      do j = ptr(v), ptr(v+1) - 1
         u = row(j)
         if (match(u).ne.-1) cycle ! u already matched
         ! Count number of flagged neighbours
         num = count( work(row(ptr(u) : ptr(u+1)-1)) .eq. v)
         ! Compare to best
         if (num.gt.bestv) then
            bestv = num
            best = u
         end if
      end do

      ! match v and best
      match(v) = best
      match(best) = v

      ! increase number of vertices in coarse graph by 1
      cn = cn + 1
   end do
end subroutine match_common_neighbours

! *******************************************************************
subroutine level_print(mp,title1,level,title2,res)

  character (len=*), intent(in) :: title1
  integer, intent(in) :: mp, level
  real(wp), optional, intent(in) :: res
  character (len=*), optional, intent(in) :: title2
  integer :: char_len1, char_len2

  char_len1 = len_trim(title1)

  if (present(res) .and. present(title2)) then
    char_len2 = len_trim(title2)
    write (mp,'(a,i4,a,g14.3)') title1, level, title2, res
  else
    write (mp,'(a,i4)') title1, level
  end if

end subroutine level_print

!
! Given a matching, produce the compressed matrix obtained by merging matched
! rows and columns.
!
subroutine compress_matrix(n, ptr_in, row_in, val_in, match, ptr_out, work, &
      row_out, val_out)
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr_in
   integer, dimension(ptr_in(n+1)-1), intent(in) :: row_in
   integer, dimension(ptr_in(n+1)-1), intent(in) :: val_in
   integer, dimension(n), intent(in) :: match
   integer, dimension(*), intent(out) :: ptr_out
   integer, dimension(3*n), target, intent(out) :: work
   integer, dimension(*), optional, intent(out) :: row_out
   integer, dimension(*), optional, intent(out) :: val_out

   integer :: col, idx, j, k, u, v
   integer, dimension(:), pointer :: seen
   integer, dimension(:), pointer :: map
   integer, dimension(:), pointer :: loc

   map  => work(    1 :   n)
   seen => work(  n+1 : 2*n)
   loc  => work(2*n+1 : 3*n)

   ! Build map array
   col = 1 ! Variable to map to
   do u = 1, n
      v = match(u)
      if(v.lt.u) cycle ! Already handled column
      map(u) = col
      map(v) = col
      col = col + 1
   end do
   
   ! Compress matrix
   col = 1 ! Insert column
   idx = 1 ! Insert location
   seen(:) = 0
   do u = 1, n
      v = match(u)
      if(v.lt.u) cycle ! Already handled column
      ptr_out(col) = idx
      seen(col) = u ! Avoid diagonals
      ! Loop over entries in column u
      do j = ptr_in(u), ptr_in(u+1)-1
         k = map(row_in(j))
         if(seen(k).ge.u) then
            ! Entry already present - just add value
            if(k.ne.col .and. present(val_out)) &
               val_out(loc(k)) = val_out(loc(k)) + val_in(j)
         else
            ! New entry - insert
            if(present(row_out)) row_out(idx) = k
               if(present(val_out)) then
                  loc(k) = idx
                  val_out(loc(k)) = val_in(j)
               endif
            seen(k) = u
            idx = idx + 1
         end if
      end do
      if(u.ne.v) then
         ! Loop over entries in column v
         do j = ptr_in(v), ptr_in(v+1)-1
            k = map(row_in(j))
            if(seen(k).ge.u) then
               ! Entry already present - just add value
               if(k.ne.col .and. present(val_out)) &
                  val_out(loc(k)) = val_out(loc(k)) + val_in(j)
            else
               ! New entry - insert
               if(present(row_out)) row_out(idx) = k
               if(present(val_out)) then
                  loc(k) = idx
                  val_out(loc(k)) = val_in(j)
               endif
               seen(k) = u
               idx = idx + 1
            end if
         end do
      end if
      col = col + 1
   end do
   ptr_out(col) = idx ! Set entry n+1
end subroutine compress_matrix

end module spral_nd_multilevel
