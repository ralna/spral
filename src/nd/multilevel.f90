module spral_nd_multilevel
   use spral_nd_maxflow
   use spral_nd_partition
   use spral_nd_refine
   use spral_nd_types
   use spral_nd_util
   implicit none

   private
   public :: multilevel_partition, mg_grid_destroy

contains

subroutine multilevel_partition(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, &
      partition, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, options,   &
      info1, lwork, work, basegrid)
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
   integer, intent(inout) :: info1
   integer, intent(in) :: lwork ! length of work array: must be atleast
      ! 9a_n + sumweight
   integer, intent(out) :: work(lwork) ! work array
   type (nd_multigrid), target, intent(inout) :: basegrid ! the multilevel of graphs 
      ! (matrices)

   type(nd_multigrid), pointer :: grid
   integer :: i, j, k, inv1, inv2, ins, info, st, grid_ne, cexit

   info1 = 0

   call nd_print_diagnostic(1, options, 'Start multilevel_partition:')

   !
   ! construct the grid at this level
   !
   st = 0
   if (.not.allocated(basegrid%graph)) allocate (basegrid%graph, stat=st)
   if (st.ne.0) then
      info1 = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info1, options, ' multilevel_partition')
      return
   end if
   call nd_matrix_construct(basegrid%graph,a_n,a_n,a_ne,st)
   if (st.lt.0) then
      info1 = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info1, options, ' multilevel_partition')
      return
   end if

   basegrid%graph%ptr(1:a_n) = a_ptr(1:a_n)
   basegrid%graph%ptr(a_n+1) = a_ne + 1
   basegrid%graph%col(1:a_ne) = a_row(1:a_ne)

   do i = 1, a_n
      do j = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr) - 1
         k = a_row(j)
         basegrid%graph%val(j) = a_weight(i)*a_weight(k)
      end do
   end do

   basegrid%size = a_n
   basegrid%level = 1

   call nd_alloc(basegrid%where, a_n, st)
   if (st.lt.0) then
      info1 = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info1, options, ' multilevel_partition')
      return
   end if

   call nd_alloc(basegrid%row_wgt, a_n, info1)
   if (info1.lt.0) then
      info1 = ND_ERR_MEMORY_ALLOC
      call nd_print_error(info1, options, ' multilevel_partition')
      return
   end if

   ! Initialise row weights
   basegrid%row_wgt(1:a_n) = a_weight(1:a_n)

   !
   ! Build coarse grid hierarchy
   !
   grid => basegrid
   do
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
         call level_print(options%unit_diagnostics, 'size of grid on level ', &
            grid%level, ' is ', real(grid%size,wp))

      call coarsen(grid, work(1:grid%size), options, cexit, st)
      if (st.lt.0) then
         info = ND_ERR_MEMORY_ALLOC
         return
      endif
      if(cexit.ne.0) exit ! Stop coarsening and partition
      grid => grid%coarse
   end do

   ! Perform coarse patitioning (if it fails, keep tring on finer grids)
   do
      grid_ne = grid%graph%ptr(grid%graph%n+1) - 1
      call nd_coarse_partition(grid%graph%n, grid_ne, grid%graph%ptr, &
         grid%graph%col, grid%row_wgt, sumweight, grid%part_div(1), &
         grid%part_div(2), grid%where, lwork, work, options, info)
      if (info.lt.0) return

      ! check if partition was succesful returned
      if (grid%part_div(1).ne.0 .and. grid%part_div(2).ne.0) exit

      ! Unlikely to get here because 99.999% of cases caught in full
      ! matrix check above. Follows same procedure as when full matrix found
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
         write (options%unit_diagnostics,'(a,i10,a)') &
            'at level ', grid%level, ' no partition found'

      ! Try partitioning previous grid level
      grid => grid%fine
   end do

   ! Prolongation
   do while(.not.associated(grid, basegrid))
      call prolong(grid%fine, sumweight, a_weight_1, a_weight_2, a_weight_sep,&
         work(1:9*grid%fine%graph%n+sumweight), options)
      if (options%print_level.ge.2 .and. options%unit_diagnostics.gt.0) &
         call level_print(options%unit_diagnostics, ' after post smoothing ', &
            grid%fine%level)
      grid => grid%fine
   end do


   !
   ! Convert from flags to partition
   ! 

   inv1 = 1
   inv2 = basegrid%part_div(1) + 1
   ins = basegrid%part_div(1) + basegrid%part_div(2) + 1

   a_weight_1 = 0
   a_weight_2 = 0
   a_weight_sep = 0
   do i = 1, a_n
      select case (basegrid%where(i))
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

   a_n1 = basegrid%part_div(1)
   a_n2 = basegrid%part_div(2)

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

subroutine coarsen(grid, work, options, cexit, st)
   type(nd_multigrid), target, intent(inout) :: grid
   integer, dimension(3*grid%size), intent(out) :: work
   type(nd_options), intent(in) :: options
   integer, intent(out) :: cexit
   integer, intent(out) :: st

   integer :: lwk
   type (nd_multigrid), pointer :: cgrid ! the coarse level grid
   integer :: cnvtx ! number of vertices (rows) in the coarse matrix
   integer :: cnedge ! number of edge (entries) in the coarse matrix
   type (nd_matrix), pointer :: p ! the coarse grid prolongator
   type (nd_matrix), pointer :: r ! the coarse grid restrictor (= p')
   type (nd_matrix), pointer :: cgraph ! the coarse graph
   type (nd_matrix), pointer :: graph ! the fine graph
   integer, dimension(:), pointer :: row_wgt ! fine graph vertex weights
   integer, dimension(:), pointer :: crow_wgt ! coarse graph vertex weights
   integer :: stop_coarsening1
   real(wp) :: reduction ! actual grid reduction factor
   real(wp) :: min_reduction ! min grid reduction factor
   real(wp) :: max_reduction ! max grid reduction factor

   ! Initialize return values
   cexit = 0   ! Keep coarsening
   st = 0      ! No error

   ! Test to see if this is either the last level or
   ! if the matrix size too small
   stop_coarsening1 = max(2,options%stop_coarsening1)
   if (grid%level.ge.options%stop_coarsening2 .or. &
         grid%size.le.stop_coarsening1) then
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
         call level_print(options%unit_diagnostics, 'end of level ', grid%level)

      cexit = 1 ! Stop, we're at the coarsest level
      return
   end if

   ! Coarsest level not yet reached so carry on coarsening
   select case(options%matching)
   case(:ND_MATCH_COMMON_NEIGHBOURS)
      lwk = 2*grid%size
      call coarsen_cn(grid, lwk, work(1:lwk), st)
   case(ND_MATCH_HEAVY)
      lwk = grid%size
      call coarsen_hec(grid, lwk, work(1:lwk), st)
   case(ND_MATCH_BEST:)
      lwk = 3*grid%size
      call coarsen_best(grid, lwk, work(1:lwk), st)
   end select
   if (st.ne.0) return

   ! ensure coarse grid quantities are allocated to sufficient size
   cgrid => grid%coarse
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

      ! recurse
      cexit = 2
      return
   end if

   ! restriction ================

   ! form the coarse grid graph and matrix
   ! cmatrix = P^T*matrix = R*matrix
   p => cgrid%p
   r => cgrid%r
   graph => grid%graph
   cgraph => cgrid%graph

   ! get the coarse matrix
   lwk = 3*grid%size
   call galerkin_graph(graph,p,r,cgraph,st,lwk,work(1:lwk))
   if (st.ne.0) return

   ! check if matrix is full
   cnedge = cgrid%graph%ptr(cgrid%graph%n+1)-1
   if (real(cnedge)/cgrid%graph%n .ge. cgrid%graph%n-1) then
      if (options%print_level.ge.1 .and. options%unit_diagnostics.gt.0) &
         write (options%unit_diagnostics,'(a,i10,a)') &
            'at level ', grid%level, ' further coarsening gives full matrix'

      cexit = 3
      return
   end if

   ! row weight cw = R*w
   row_wgt => grid%row_wgt(1:grid%size)
   crow_wgt => cgrid%row_wgt(1:cgrid%size)
   call nd_matrix_multiply_vec(r, row_wgt, crow_wgt)

end subroutine coarsen

subroutine prolong(grid, sumweight, a_weight_1, a_weight_2, a_weight_sep, work, options)
   type(nd_multigrid), target, intent(inout) :: grid
   integer, intent(in) :: sumweight
   integer, intent(out) :: a_weight_1
   integer, intent(out) :: a_weight_2
   integer, intent(out) :: a_weight_sep
   integer, dimension(9*grid%graph%n+sumweight), intent(out) :: work
   type(nd_options), intent(in) :: options

   type (nd_multigrid), pointer :: cgrid ! the coarse level grid
   type (nd_matrix), pointer :: p ! the coarse grid prolongator
   integer, dimension(:), pointer :: fwhere ! partition on fine grid
   integer, dimension(:), pointer :: cwhere ! partition on coarse grid

   integer :: i, j, a_ne, p1, p2, psep
   integer :: work_ptr, partition_ptr

   cgrid => grid%coarse
   p => cgrid%p
   fwhere => grid%where(1:grid%size)
   cwhere => cgrid%where(1:cgrid%size)

   ! injection of the order from coarse grid to the
   ! fine grid, since cwhere(i) is the index of the
   ! i-th vertex in the new ordering, the order
   ! of this vertex should be where(i)
   ! grid%where = P*order_on_coarse_grid
   ! here P is a special matrix with only one non-zero entry per row
   grid%part_div(1:2) = 0
   call nd_matrix_multiply_vec(p, cwhere, fwhere)

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
  integer :: partition_method
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
  if (options%coarse_partition_method.le.1) then
    partition_method = 1
  else
    partition_method = 2
  end if

  partition_ptr = 0 ! length a_n
  work_ptr = partition_ptr + a_n ! max length needed 9*a_n+a_ne

  allocate (work1(a_n),stat=st)
  if (st.ne.0) then
    info = ND_ERR_MEMORY_ALLOC
    return
  end if


  select case (partition_method)
  case (1)
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

! *****************************************************************

recursive subroutine mg_grid_destroy(grid)
  ! deallocate a grid structure
  type (nd_multigrid) :: grid

  if (associated(grid%coarse)) then

    call mg_grid_destroy(grid%coarse)

    if (grid%level.ne.1) then

      call multigrid_deallocate(grid)

    else

      call multigrid_deallocate_first(grid)

    end if

  else

    if (grid%level.ne.1) then

      call multigrid_deallocate_last(grid)

    else

      call multigrid_deallocate_first(grid)

    end if

  end if

end subroutine mg_grid_destroy


!
! Deallocate a grid (at given level between last and first)
!
subroutine multigrid_deallocate(grid)
  type (nd_multigrid) :: grid

  call nd_matrix_destruct(grid%graph)
  call nd_matrix_destruct(grid%p)
  call nd_matrix_destruct(grid%r)

  if(associated(grid%coarse)) deallocate(grid%coarse)
  deallocate (grid%graph,grid%p,grid%r,grid%where,grid%row_wgt)
  nullify (grid%coarse)
end subroutine multigrid_deallocate

!
! Deallocate a grid (at the last level).
! In this case the matrix grid%graph has not been formed yet
!
subroutine multigrid_deallocate_last(grid)
   type (nd_multigrid) :: grid

   call nd_matrix_destruct(grid%p)
   call nd_matrix_destruct(grid%r)

   if(associated(grid%coarse)) deallocate(grid%coarse)
   deallocate (grid%graph,grid%p,grid%r,grid%where,grid%row_wgt)
   nullify (grid%coarse)
end subroutine multigrid_deallocate_last

!
! Deallocate a grid (at the first level).
! In this case the matrix grid%p does not exist
!
subroutine multigrid_deallocate_first(grid)
  type (nd_multigrid) :: grid

  if (allocated(grid%graph)) &
    call nd_matrix_destruct(grid%graph)

  deallocate (grid%where, grid%row_wgt)
end subroutine multigrid_deallocate_first

! ***************************************************************
subroutine coarsen_hec(grid,lwork,work,st)
  ! coarsen the grid using heavy-edge collapsing and set up the
  ! coarse grid equation, the prolongator and restrictor

  type (nd_multigrid), intent(inout), TARGET :: grid
  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)
  integer, intent(out) :: st

  st = 0
  if (.not.associated(grid%coarse)) allocate (grid%coarse, stat=st)
  if(st.ne.0) return

  grid%coarse%fine => grid

  ! find the prolongator
  call prolng_heavy_edge(grid,lwork,work,st)


  grid%coarse%level = grid%level + 1

end subroutine coarsen_hec


! ***************************************************************
subroutine coarsen_cn(grid,lwork,work,st)
  ! coarsen the grid using common neighbours collapsing and set up the
  ! coarse grid equation, the prolongator and restrictor

  type (nd_multigrid), intent(inout), TARGET :: grid

  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)
  integer, intent(out) :: st

  st = 0
  if (.not.associated(grid%coarse)) allocate (grid%coarse, stat=st)
  if(st.ne.0) return

  grid%coarse%fine => grid

  ! find the prolongator

  call prolng_common_neigh(grid,lwork,work,st)

  grid%coarse%level = grid%level + 1

end subroutine coarsen_cn

! ***************************************************************
subroutine coarsen_best(grid,lwork,work,st)
  ! coarsen the grid using common neighbours collapsing and set up the
  ! coarse grid equation, the prolongator and restrictor

  integer, intent(inout) :: st
  type (nd_multigrid), intent(inout), TARGET :: grid

  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)

  st = 0
  if ( .not. associated(grid%coarse)) allocate (grid%coarse, stat=st)
  if(st.ne.0) return

  grid%coarse%fine => grid

  ! find the prolongator

  call prolng_best(grid,lwork,work,st)

  grid%coarse%level = grid%level + 1

end subroutine coarsen_best

! ***************************************************************
subroutine nd_matrix_multiply_vec(matrix,x,y)
  ! subroutine nd_matrix_multiply_vec(matrix,x,y)

  ! y = matrix*x where x and y are integer vectors. Entries of
  ! matrix is assumed to be one. Dimension of y
  ! is checked and returned if it is smaller than the row dimension
  ! of x    !

  ! matrix: of the derived type nd_matrix, intent(in),
  ! the sparse matrix in compressed sparse row format
  type (nd_matrix), intent(in) :: matrix

  ! x: integer array of intent(in), a vector to be
  ! multiplied with the matrix
  integer, intent(in), dimension(*) :: x

  ! y: integer array of intent(out), the result of
  ! matrix*x or matrix^T*x
  integer, intent(out), dimension(*) :: y

  ! local ==========
  integer :: m, n, i, l1, l2

  m = matrix%m
  n = matrix%n

  do i = 1, m
    l1 = matrix%ptr(i)
    l2 = matrix%ptr(i+1) - 1
    y(i) = sum(x(matrix%col(l1:l2)))
  end do
end subroutine nd_matrix_multiply_vec


!
! Ensures all allocatable components of matrix are deallocated
!
! NB: Any deallocation errors are ignored
subroutine nd_matrix_destruct(matrix)
  type (nd_matrix), intent(inout) :: matrix

  integer :: st

  deallocate (matrix%ptr, stat=st)
  deallocate (matrix%col, stat=st)
  deallocate (matrix%val, stat=st)
end subroutine nd_matrix_destruct

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


! ********************************************************

subroutine prolng_heavy_edge(grid,lwork,work,st)

  ! calculate the prolongator for heavy-edge collapsing:
  ! match the vertices of the heaviest edges

  ! input fine grid
  type (nd_multigrid), target, intent(inout) :: grid
  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)
  integer, intent(out) :: st

  ! coarse grid based on the fine grid
  type (nd_multigrid), pointer :: cgrid

  ! the fine grid row connectivity graph
  type (nd_matrix), pointer :: graph

  ! the coarse grid prolongator
  type (nd_matrix), pointer :: p

  ! the coarse grid restrictor
  type (nd_matrix), pointer :: r

  ! the number of fine and coarse grid vertices
  integer :: nvtx, cnvtx

  ! working variables
  integer :: v, u, j, i, k
  integer :: nz

  ! whether a vertex is matched already
  integer, parameter :: unmatched = -1

  ! matching status of each vertex
  integer :: ptr_match

  ! maximum weight and index of edges connected to the current vertex
  integer :: maxwgt
  integer :: maxind

  ! allocate the prolongation matrix pointers
  cgrid => grid%coarse
  graph => grid%graph

  ! allocate the graph and matrix pointer and the mincut pointer
  ! so that everything is defined 
  st = 0
  if (.not.allocated(cgrid%graph)) allocate (cgrid%graph, stat=st)
  if (st.ne.0) return

  nvtx = graph%n

  ! prolongator start here ================================

  ! initialise the matching status and randomly permute the vertex order
  ptr_match = 0

  work(ptr_match+1:ptr_match+nvtx) = unmatched

  ! loop over each vertex and match along the heaviest edge
  cnvtx = 0
  nz = 0
  do i = 1, nvtx
    v = i
    ! If already matched, next vertex please
    if (work(ptr_match+v).ne.unmatched) cycle
    maxwgt = -huge(0)
    ! in the case no match is found then match itself
    maxind = v
    ! Loop over entries in row v
    do j = graph%ptr(v), graph%ptr(v+1) - 1
      ! u is col index of en entry in row v (so u is neighbor of v)
      u = graph%col(j)
      ! heavy edge matching
      ! if u is unmatched and value of the entry in col. u is greater
      ! than maxwgt, select u as the matching.
      if (work(ptr_match+u).eq.unmatched .and. maxwgt.lt.abs(graph%val(j))) &
          then
        maxwgt = abs(graph%val(j))
        maxind = u
      end if
    end do
    ! NOTE: maxind .ge. v
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
    if (maxind.ne.v) then
      nz = nz + 1
    end if
  end do

  ! storage allocation for col. indices and values of prolongation
  ! matrix P (nvtx * cnvtx)
  if (.not.allocated(cgrid%p)) allocate (cgrid%p, stat=st)
  if (st.ne.0) return
  p => cgrid%p
  call nd_matrix_construct(p,nvtx,cnvtx,nz,st)
  if(st.ne.0) return


  ! storage allocation for col. indices and values of restiction
  ! matrix R (cnvtx * nvtx)
  if (.not.allocated(cgrid%r)) allocate (cgrid%r, stat=st)
  if (st.ne.0) return
  r => cgrid%r
  call nd_matrix_construct(r,cnvtx,nvtx,nz,st)
  if(st.ne.0) return

  r%val(1:nz) = 1

  ! store restriction matrix
  r%ptr(cnvtx+1) = nz + 1

  j = 1
  k = 1
  do i = 1, nvtx
    if (work(ptr_match+i).eq.i) then
      r%ptr(k) = j
      r%col(j) = i
      j = j + 1
      k = k + 1
    else
      if (work(ptr_match+i).gt.i) then
        r%ptr(k) = j
        r%col(j) = i
        r%col(j+1) = work(ptr_match+i)
        j = j + 2
        k = k + 1
      end if
    end if
  end do

  ! store prolongation matrix

  p%ptr(1) = 1
  do i = 1, nvtx
    p%ptr(i+1) = p%ptr(i) + 1
  end do

  p%val(1:nz) = 1

  j = 1
  do i = 1, nvtx
    k = work(ptr_match+i)
    if (k.eq.i) then
      p%col(p%ptr(i)) = j
      j = j + 1
    else
      if (k.gt.i) then
        p%col(p%ptr(i)) = j
        p%col(p%ptr(k)) = j
        j = j + 1
      end if
    end if
  end do

  ! size of coarse grid
  cgrid%size = cnvtx

end subroutine prolng_heavy_edge

! *******************************************************************

subroutine prolng_common_neigh(grid,lwork,work,st)

  ! calculate the prolongator:
  ! match the vertices of with most neighbours in common

  ! input fine grid
  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)
  type (nd_multigrid), target, intent(inout) :: grid
  integer, intent(out) :: st

  ! coarse grid based on the fine grid
  type (nd_multigrid), pointer :: cgrid

  ! the fine grid row connectivity graph
  type (nd_matrix), pointer :: graph

  ! the coarse grid prolongator
  type (nd_matrix), pointer :: p

  ! the fine grid restrictor
  type (nd_matrix), pointer :: r

  ! the number of fine and coarse grid vertices
  integer :: nvtx, cnvtx

  ! working variables
  integer :: v, u, w, j, i, k

  ! whether a vertex is matched already
  integer, parameter :: unmatched = -1

  ! matching status of each vertex
  integer :: ptr_match
  ! flag array to flag up neighbours of a  node
  integer :: ptr_flag

  ! maximum no. neighbours and index of edges connected to the current
  ! vertex
  integer :: max_neigh, maxind, num

  integer :: nz

  ! allocate the prolongation matrix pointers
  cgrid => grid%coarse
  graph => grid%graph

  ! allocate the graph and matrix pointer and the mincut pointer
  ! so that everything is defined
  st = 0
  if (.not.allocated(cgrid%graph)) allocate (cgrid%graph, stat=st)
  if (st.ne.0) return

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
  do i = 1, nvtx
    v = i
    ! If already matched, next vertex please
    if (work(ptr_match+v).ne.unmatched) cycle
    ! access the col. indices of row v

    ! in the case no match is found then match itself
    maxind = v
    ! Loop over entries in row v and set flag for each entry
    work(ptr_flag+v) = i
    do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
      ! u is col index of en entry in row v (so u is neighbor of v)
      u = grid%graph%col(j)
      work(ptr_flag+u) = i
    end do
    ! For each unmatched neighbour of v, count the number of
    ! neighbours it has in common with v
    max_neigh = 0
    do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
      u = grid%graph%col(j)
      ! cycle is u is already matched
      if (work(ptr_match+u).ne.unmatched) cycle
      num = 0
      do k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
        w = grid%graph%col(k)
        if (work(ptr_flag+w).eq.i) num = num + 1
      end do
      if (num.gt.max_neigh) then
        max_neigh = num
        maxind = u
      end if
    end do

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
    if (maxind.ne.v) then
      nz = nz + 1
    end if
  end do

  ! storage allocation for col. indices and values of prolongation
  ! matrix P (order nvtx * cnvtx)
  if ( .not. allocated(cgrid%p)) allocate (cgrid%p, stat=st)
  if (st.ne.0) return
  p => cgrid%p
  call nd_matrix_construct(p,nvtx,cnvtx,nz,st)
  if(st.ne.0) return
  p%val(1:nz) = 0

  ! storage allocation for col. indices and values of restiction
  ! matrix R (cnvtx * nvtx)
  if ( .not. allocated(cgrid%r)) allocate (cgrid%r, stat=st)
  if (st.ne.0) return
  r => cgrid%r
  call nd_matrix_construct(r,cnvtx,nvtx,nz,st)
  if(st.ne.0) return

  r%val(1:nz) = 1

  ! store restriction matrix
  r%ptr(cnvtx+1) = nz + 1

  j = 1
  k = 1
  do i = 1, nvtx
    if (work(ptr_match+i).eq.i) then
      r%ptr(k) = j
      r%col(j) = i
      j = j + 1
      k = k + 1
    else
      if (work(ptr_match+i).gt.i) then
        r%ptr(k) = j
        r%col(j) = i
        r%col(j+1) = work(ptr_match+i)
        j = j + 2
        k = k + 1
      end if
    end if
  end do


  ! store prolongation matrix

  p%ptr(1) = 1
  do i = 1, nvtx
    p%ptr(i+1) = p%ptr(i) + 1
  end do

  p%val(1:nz) = 1

  j = 1
  do i = 1, nvtx
    k = work(ptr_match+i)
    if (k.eq.i) then
      p%col(p%ptr(i)) = j
      j = j + 1
    else
      if (k.gt.i) then
        p%col(p%ptr(i)) = j
        p%col(p%ptr(k)) = j
        j = j + 1
      end if
    end if
  end do


  ! size of coarse grid
  cgrid%size = cnvtx

end subroutine prolng_common_neigh

! ********************************************
subroutine prolng_best(grid,lwork,work,st)

  ! calculate the prolongator for heavy-edge collapsing:
  ! match the vertices of the heaviest edges

  ! input fine grid
  type (nd_multigrid), target, intent(inout) :: grid
  integer, intent(in) :: lwork
  integer, intent(out), TARGET :: work(lwork)
  integer, intent(out) :: st

  ! coarse grid based on the fine grid
  type (nd_multigrid), pointer :: cgrid

  ! the fine grid row connectivity graph
  type (nd_matrix), pointer :: graph

  ! the coarse grid prolongator
  type (nd_matrix), pointer :: p

  ! the coarse grid restrictor
  type (nd_matrix), pointer :: r

  ! the number of fine and coarse grid vertices
  integer :: nvtx, cnvtx, cnvtx1

  ! working variables
  integer :: v, u, j, i, k
  integer :: nz

  ! whether a vertex is matched already
  integer, parameter :: unmatched = -1

  ! matching status of each vertex
  integer :: ptr_match, ptr_match1, ptr_flag, max_neigh, num, w

  ! maximum weight and index of edges connected to the current vertex
  integer :: maxwgt
  integer :: maxind

  integer, pointer, dimension(:) :: matching

  ! allocate the prolongation matrix pointers
  cgrid => grid%coarse
  graph => grid%graph

  ! allocate the graph and matrix pointer and the mincut pointer
  ! so that everything is defined
  st = 0
  if (.not.allocated(cgrid%graph)) allocate (cgrid%graph, stat=st)
  if (st.ne.0) return

  nvtx = graph%n

  ptr_match = 0
  ptr_match1 = ptr_match + nvtx

  ! -----------------------------------------------------------------
  ! Find heavy-edge matching

  ! initialise the matching status and randomly permute the vertex order

  work(ptr_match+1:ptr_match+nvtx) = unmatched

  ! loop over each vertex and match along the heaviest edge
  cnvtx = 0
  do i = 1, nvtx
    v = i
    ! If already matched, next vertex please
    if (work(ptr_match+v).ne.unmatched) cycle
    maxwgt = -huge(0)
    ! in the case no match is found then match itself
    maxind = v
    ! Loop over entries in row v
    do j = graph%ptr(v), graph%ptr(v+1) - 1
      ! u is col index of en entry in row v (so u is neighbor of v)
      u = graph%col(j)
      ! heavy edge matching
      ! if u is unmatched and value of the entry in col. u is greater
      ! than maxwgt, select u as the matching.
      if (work(ptr_match+u).eq.unmatched .and. maxwgt.lt.abs(graph%val(j))) &
          then
        maxwgt = abs(graph%val(j))
        maxind = u
      end if
    end do
    ! NOTE: maxind .ge. v
    ! the neighbor with heaviest weight
    work(ptr_match+v) = maxind
    ! mark maxind as having been matched
    work(ptr_match+maxind) = v
    ! increase number of vertices in coarse graph by 1
    cnvtx = cnvtx + 1
    ! construct the prolongation matrix: find vertex v and maxind is
    ! linked
    ! with the coarse grid vertex cnvtx
  end do
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
  do i = 1, nvtx
    v = i
    ! If already matched, next vertex please
    if (work(ptr_match1+v).ne.unmatched) cycle
    ! access the col. indices of row v

    ! in the case no match is found then match itself
    maxind = v
    ! Loop over entries in row v and set flag for each entry
    work(ptr_flag+v) = i
    do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
      ! u is col index of en entry in row v (so u is neighbor of v)
      u = grid%graph%col(j)
      work(ptr_flag+u) = i
    end do
    ! For each unmatched neighbour of v, count the number of
    ! neighbours it has in common with v
    max_neigh = 0
    do j = grid%graph%ptr(v), grid%graph%ptr(v+1) - 1
      u = grid%graph%col(j)
      ! cycle is u is already matched
      if (work(ptr_match1+u).ne.unmatched) cycle
      num = 0
      do k = grid%graph%ptr(u), grid%graph%ptr(u+1) - 1
        w = grid%graph%col(k)
        if (work(ptr_flag+w).eq.i) num = num + 1
      end do
      if (num.gt.max_neigh) then
        max_neigh = num
        maxind = u
      end if
    end do

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
    if (maxind.ne.v) then
      nz = nz + 1
    end if
  end do

  ! --------------------------------------------------------------------
  ! -
  if (cnvtx.le.cnvtx1) then
    ! use heavy-edge matching
    matching => work(ptr_match+1:ptr_match+nvtx)
  else
    ! use common neighbours matching
    matching => work(ptr_match1+1:ptr_match1+nvtx)
    cnvtx = cnvtx1
  end if


  ! storage allocation for col. indices and values of prolongation
  ! matrix P (nvtx * cnvtx)
  if ( .not. allocated(cgrid%p)) allocate (cgrid%p, stat=st)
  if (st.ne.0) return
  p => cgrid%p
  call nd_matrix_construct(p,nvtx,cnvtx,nz,st)
  if(st.ne.0) return

  ! storage allocation for col. indices and values of restiction
  ! matrix R (cnvtx * nvtx)
  if ( .not. allocated(cgrid%r)) allocate (cgrid%r, stat=st)
  if (st.ne.0) return
  r => cgrid%r
  call nd_matrix_construct(r,cnvtx,nvtx,nz,st)
  if(st.ne.0) return

  r%val(1:nz) = 1

  ! store restriction matrix
  r%ptr(cnvtx+1) = nz + 1

  j = 1
  k = 1
  do i = 1, nvtx
    if (matching(i).eq.i) then
      r%ptr(k) = j
      r%col(j) = i
      j = j + 1
      k = k + 1
    else
      if (matching(i).gt.i) then
        r%ptr(k) = j
        r%col(j) = i
        r%col(j+1) = matching(i)
        j = j + 2
        k = k + 1
      end if
    end if
  end do

  ! store prolongation matrix

  p%ptr(1) = 1
  do i = 1, nvtx
    p%ptr(i+1) = p%ptr(i) + 1
  end do

  p%val(1:nz) = 1

  j = 1
  do i = 1, nvtx
    k = matching(i)
    if (k.eq.i) then
      p%col(p%ptr(i)) = j
      j = j + 1
    else
      if (k.gt.i) then
        p%col(p%ptr(i)) = j
        p%col(p%ptr(k)) = j
        j = j + 1
      end if
    end if
  end do

  ! size of coarse grid
  cgrid%size = cnvtx
  nullify (matching)

end subroutine prolng_best


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



! *************************************************
subroutine galerkin_graph(matrix,p,r,cmatrix,st,lwork,work)

  ! Given matrix on fine grid and a prolongation operator p,
  ! find the coarse matrix R*A*P

  ! matrix: fine grid matrix
  type (nd_matrix), intent(in) :: matrix
  ! p: prolongation operator
  type (nd_matrix), intent(in) :: p
  ! r: restriction operator
  type (nd_matrix), intent(in) :: r
  ! cmatrix: coarse grid matrix
  type (nd_matrix), intent(inout) :: cmatrix
  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)

  ! nvtx,cnvtx: size of fine and coarse grid
  integer :: nvtx, cnvtx
  integer :: nz

  integer, intent(inout) :: st

  ! call mc65_matrix_transpose(p,r,info65)
  ! if (info65.lt.0) then
  ! info = info65
  ! return
  ! end if
  nvtx = matrix%n
  cnvtx = p%n

  ! get the size of the coarse matrix first
  call galerkin_graph_rap_size(nvtx,cnvtx,nz,p%ptr(nvtx+1)-1,p%col, &
    p%ptr,matrix%ptr(nvtx+1)-1,matrix%col,matrix%ptr,r%ptr(cnvtx+1)-1, &
    r%col,r%ptr,lwork,work(1:lwork))

  call nd_matrix_construct(cmatrix,cnvtx,cnvtx,nz,st)
  if(st.ne.0) return

  call galerkin_graph_rap(nvtx,cnvtx,p%ptr(nvtx+1)-1,p%val,p%col,p%ptr, &
    matrix%ptr(nvtx+1)-1,matrix%val,matrix%col,matrix%ptr, &
    r%ptr(cnvtx+1)-1,r%val,r%col,r%ptr,nz,cmatrix%val,cmatrix%col, &
    cmatrix%ptr,lwork,work(1:lwork))

end subroutine galerkin_graph

! *************************************************

subroutine galerkin_graph_rap_size(nvtx,cnvtx,nz,nzp,pcol,pptr,nzaa, &
    acol,aptr,nzr,rcol,rptr,lwork,work)
  ! get the number of nonzeros in R*A*P
  ! nvtx: size of aa matrix
  ! cnvtx: size of ca matrix
  integer, intent(in) :: nvtx, cnvtx
  ! nz: number of nonzeros in R*A*P
  integer, intent(out) :: nz

  ! P: matrix
  integer, intent(in) :: nzp
  integer, intent(in), dimension(nzp) :: pcol
  integer, intent(in), dimension(nvtx+1) :: pptr
  ! aa: matrix
  integer, intent(in) :: nzaa
  integer, intent(in), dimension(nzaa) :: acol
  integer, intent(in), dimension(nvtx+1) :: aptr
  ! R: matrix
  integer, intent(in) :: nzr
  integer, intent(in), dimension(nzr) :: rcol
  integer, intent(in), dimension(cnvtx+1) :: rptr

  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)

  ! mask: masking array to see if an entry has been seen before
  integer :: ptr_mask
  ! i,j,k: loop index
  integer :: i, j, k
  ! nz: number of nonzeros so far in ca
  integer :: nz1
  ! various neighbors
  integer :: neigh, neighneigh

  ! col: column index of a row of r*matrix
  integer :: ptr_col

  ptr_mask = 0
  ptr_col = ptr_mask + nvtx
  work(ptr_mask+1:ptr_mask+nvtx) = 0
  nz = 0
  ! loop over coarse grid points
  do i = 1, cnvtx
    ! first form row i of (r*matrix)
    nz1 = 0
    ! for each vertex D that restricts to C (including itself).
    do j = rptr(i), rptr(i+1) - 1
      neigh = rcol(j)
      ! find D's neighbor
      do k = aptr(neigh), aptr(neigh+1) - 1
        neighneigh = acol(k)
        if (work(ptr_mask+neighneigh).ne.i) then
          nz1 = nz1 + 1
          work(ptr_col+nz1) = neighneigh
          work(ptr_mask+neighneigh) = i
        end if
      end do
    end do
    ! form row i of (r*matrix)*p
    do j = 1, nz1
      neigh = work(ptr_col+j)
      do k = pptr(neigh), pptr(neigh+1) - 1
        neighneigh = pcol(k)
        if (work(ptr_mask+neighneigh).ne.-i .and. neighneigh.ne.i) then
          nz = nz + 1
          work(ptr_mask+neighneigh) = -i
        end if
      end do
    end do
  end do

end subroutine galerkin_graph_rap_size
! ******************************************************
subroutine galerkin_graph_rap(nvtx,cnvtx,nzp,pa,pcol,pptr,nzaa,aa,acol, &
    aptr,nzr,ra,rcol,rptr,nzca,ca,ccol,cptr,lwork,work)
  ! multiply R*A*P to get CA
  ! nvtx: size of aa matrix
  ! cnvtx: size of ca matrix
  integer, intent(in) :: nvtx, cnvtx
  ! p: matrix
  integer, intent(in) :: nzp
  integer, intent(in), dimension(nzp) :: pa
  integer, intent(in), dimension(nzp) :: pcol
  integer, intent(in), dimension(nvtx+1) :: pptr
  ! aa: matrix
  integer, intent(in) :: nzaa
  integer, intent(in), dimension(:) :: aa
  integer, intent(in), dimension(nzaa) :: acol
  integer, intent(in), dimension(:) :: aptr
  ! r: matrix
  integer, intent(in) :: nzr
  integer, intent(in), dimension(nzr) :: ra
  integer, intent(in), dimension(nzr) :: rcol
  integer, intent(in), dimension(cnvtx+1) :: rptr
  ! ca: matrix
  integer, intent(in) :: nzca
  integer, intent(inout), dimension(nzca) :: ca
  integer, intent(inout), dimension(nzca) :: ccol
  integer, intent(inout), dimension(cnvtx+1) :: cptr

  integer, intent(in) :: lwork
  integer, intent(out) :: work(lwork)


  ! mask: masking array to see if an entry has been seen before
  integer :: ptr_mask
  ! i,j,k,l: loop index
  integer :: i, j, k
  ! nz: number of nonzeros so far in ca
  integer :: nz, nzz, nz1
  ! various neighbors
  integer :: neigh, neighneigh
  ! r_ij: (i,j) element of r
  integer :: r_ij
  ! col: column index of a row of r*matrix
  ! a: values of a row of r*matrix
  integer :: ptr_col, ptr_a

  ptr_mask = 0
  ptr_col = ptr_mask + nvtx
  ptr_a = ptr_col + nvtx
  ! now get the entries of the coarse matrix
  cptr(1) = 1
  work(ptr_mask+1:ptr_mask+nvtx) = 0
  nz = 0
  ! loop over every coarse grid point
  do i = 1, cnvtx
    ! first form row i of (r*matrix)
    nz1 = 0
    ! foreach each vertex D that restricts to C (including itself).
    do j = rptr(i), rptr(i+1) - 1
      neigh = rcol(j)
      r_ij = ra(j)
      ! find D's neighbor
      do k = aptr(neigh), aptr(neigh+1) - 1
        neighneigh = acol(k)
        nzz = work(ptr_mask+neighneigh)
        if (nzz.eq.0) then
          nz1 = nz1 + 1
          work(ptr_col+nz1) = neighneigh
          work(ptr_a+nz1) = r_ij*aa(k)
          work(ptr_mask+neighneigh) = nz1
        else
          work(ptr_a+nzz) = work(ptr_a+nzz) + r_ij*aa(k)
        end if
      end do
    end do
    do j = 1, nz1
      work(ptr_mask+work(ptr_col+j)) = 0
    end do

    ! form row i of (r*matrix)*p
    do j = 1, nz1
      neigh = work(ptr_col+j)
      r_ij = work(ptr_a+j)
      do k = pptr(neigh), pptr(neigh+1) - 1
        neighneigh = pcol(k)
        if (neighneigh.eq.i) cycle
        nzz = work(ptr_mask+neighneigh)
        if (nzz.eq.0) then
          nz = nz + 1
          work(ptr_mask+neighneigh) = nz
          ca(nz) = r_ij*pa(k)
          ccol(nz) = neighneigh
        else
          ca(nzz) = ca(nzz) + r_ij*pa(k)
        end if
      end do
    end do

    do j = cptr(i), nz
      work(ptr_mask+ccol(j)) = 0
    end do
    cptr(i+1) = nz + 1
  end do


end subroutine galerkin_graph_rap


end module spral_nd_multilevel
