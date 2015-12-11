module spral_nd_partition
   use spral_nd_types
   use spral_nd_util
   implicit none

   private
   public :: nd_half_level_set, nd_level_set

contains

!
! Partition the matrix using the half level set (Ashcraft) method
!
subroutine nd_half_level_set(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, &
      level, a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, partition, &
      work, options, info, band, depth, use_multilevel)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(in) :: a_ptr(a_n)
   integer, intent(in) :: a_row(a_ne)
   integer, intent(in) :: a_weight(a_n)
   integer, intent(in) :: sumweight ! sum of entries in a_weight
   integer, intent(in) :: level ! current level of nested dissection
   integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
   integer, intent(out) :: a_weight_1, a_weight_2, a_weight_sep ! Weighted
      ! size of partitions and separator
   integer, intent(out) :: partition(a_n) ! First a_n1 entries will contain
      ! list of (local) indices in partition 1; next a_n2 entries will
      ! contain list of (local) entries in partition 2; entries in
      ! separator are listed at the end
   integer, intent(out) :: work(9*a_n+sumweight) ! used as work array
   type (nd_options), intent(in) :: options
   integer, intent(inout) :: info
   real(wp), intent(out) :: band ! If level = 0, then on output
      ! band = 100*L/a_n, where L is the size of the largest levelset
   real(wp), intent(out) :: depth ! If level = 0, then on output
      ! depth = num_levels_nend
   logical, intent(inout) :: use_multilevel ! are we allowed to use a
      ! multilevel partitioning strategy

   integer :: nstrt, nend
   integer :: i, j, dptr, p1sz, p2sz, sepsz, k
   integer :: unit_diagnostics ! unit on which to print diagnostics
   integer :: mask_p, level_p, level_ptr_p, level2_p, level2_ptr_p, &
      work_p
   integer :: num_levels_nend ! no. levels in structure rooted at nend
   integer :: num_levels_nstrt ! no. levels in structure rooted at nstrt
   integer :: num_entries ! no. entries in level set structure
   integer :: best_sep_start
   integer :: distance
   integer :: distance_ptr
   integer :: lwidth, mindeg, degree
   integer :: ww
   real(wp) :: bestval
   real(wp) :: val
   real(wp) :: balance_tol
   logical :: printi, printd
   logical :: imbal, use_multilevel_copy

   integer, parameter :: max_search = 5

   ! ---------------------------------------------
   ! Printing levels
   unit_diagnostics = options%unit_diagnostics
   printi = (options%print_level.eq.1 .and. unit_diagnostics.ge.0)
   printd = (options%print_level.ge.2 .and. unit_diagnostics.ge.0)

   call nd_print_diagnostic(1, options, ' ')
   call nd_print_diagnostic(1, options, 'Use two-sided level set method')

   ! Initialize return vars
   band = -1
   depth = -1

   ! If we're going to use multilevel regardless, immediate return
   if (options%partition_method.eq.1 .and. use_multilevel) return
   if (options%partition_method.gt.1 .and. level.gt.0 .and. use_multilevel) &
      return

   ! Initialize various internal variables
   balance_tol = max(1.0_wp, options%balance)
   imbal = (balance_tol .le. sumweight-2.0_wp)
   p2sz = 0
   sepsz = 0
   use_multilevel_copy = use_multilevel

   ! Find pseudoperipheral nodes nstart and nend, and the level structure
   ! rooted at nend
   mask_p = 0 ! size a_n
   level_ptr_p = mask_p + a_n ! size a_n
   level_p = level_ptr_p + a_n ! size a_n
   level2_ptr_p = level_p + a_n ! size a_n
   level2_p = level2_ptr_p + a_n ! size a_n
   work_p = level2_p + a_n ! size 2*a_n

   nend = -1

   ! Choose nstrt as node with minimum degree
   mindeg = sumweight + 1
   do i = 1, a_n
      degree = 0
      do j = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
         degree = degree + a_weight(a_row(j))
      end do
      if (degree.lt.mindeg) then
         mindeg = degree
         nstrt = i
      end if
   end do

   call nd_find_pseudo(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,          &
      work(level_ptr_p+1:level_ptr_p+a_n), work(level_p+1:level_p+a_n),       &
      nstrt, nend, max_search, work(work_p+1:work_p+2*a_n), num_levels_nend,  &
      num_entries, lwidth)

   if (num_entries.lt.a_n) then
      ! matrix is separable
      a_n1 = num_entries
      a_n2 = a_n - a_n1
      a_weight_sep = 0
      a_weight_1 = 0
      work(work_p+1:work_p+a_n) = 0
      do i = 1, num_entries
         j = work(level_p+i)
         partition(i) = j
         work(work_p+j) = 1
         a_weight_1 = a_weight_1 + a_weight(j)
      end do
      a_weight_2 = sumweight - a_weight_1
      j = num_entries + 1
      do i = 1, a_n
         if (work(work_p+i).eq.0) then
            partition(j) = i
            j = j + 1
         end if
      end do
      if (level.eq.0) band = -real(lwidth,wp)

      use_multilevel = .false. ! Will be reset for each component anyway
      return
   end if

   ! ********************************************************************
   if (level.eq.0) then
      band = (100.0_wp * lwidth) / a_n
      depth = (100.0_wp * num_levels_nend) / a_n
   end if
   if (options%stop_coarsening2.le.0 .or. options%partition_method.lt.1) &
      use_multilevel = .false.
   if (options%partition_method.ge.2 .and. use_multilevel) then
      if ( (100.0_wp*lwidth)/sumweight.le.3.0 .or. a_n.lt.options%ml_call ) &
         use_multilevel = .false.
   end if

   if (use_multilevel) return


   ! Find level structure rooted at nstrt
   work(mask_p+1:mask_p+a_n) = 1
   call nd_level_struct(nstrt, a_n, a_ne, a_ptr, a_row,                    &
      work(mask_p+1:mask_p+a_n), work(level2_ptr_p+1:level2_ptr_p+a_n),    &
      work(level2_p+1:level2_p+a_n), num_levels_nstrt, lwidth, num_entries)

   ! Calculate difference in distances from nstart and nend, and set up
   ! lists D_i of nodes with same distance
   distance_ptr = work_p
   distance = mask_p
   call nd_distance(a_n, num_levels_nend, work(level_ptr_p+1:level_ptr_p+a_n), &
      work(level_p+1:level_p+a_n), num_levels_nstrt,                           &
      work(level2_ptr_p+1:level2_ptr_p+a_n), work(level2_p+1:level2_p+a_n),    &
      work(distance_ptr+1:distance_ptr+2*a_n-1), work(distance+1:distance+a_n) )

   ! Do not need the information in work(level_ptr_p+1:level2_ptr_p)
   ! Calculate total weight in each distance level
   dptr = level_ptr_p + a_n
   work(dptr+1-a_n:dptr+a_n-1) = 0
   do i = 1 - num_levels_nstrt, num_levels_nend - 2
      do j = work(distance_ptr+a_n+i), work(distance_ptr+a_n+i+1) - 1
         work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
      end do
   end do
   i = num_levels_nend - 1
   do j = work(distance_ptr+a_n+i), a_n
      work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
   end do

   ! Find first possible separator
   ww = 2
   p1sz = 0
   do i = 1 - num_levels_nstrt, num_levels_nend - ww - 2
      p1sz = p1sz + work(dptr+i)
      sepsz = work(dptr+i+1)
      do j = 1, ww - 1
         sepsz = sepsz + work(dptr+i+1+j)
      end do
      p2sz = sumweight - p1sz - sepsz
      if (p1sz.gt.0 .and. sepsz.gt.0) exit
   end do

   if (i+ww.ge.num_levels_nend-1 .or. p2sz.eq.0) then
      ! Not possible to find separator
      ! This can only be invoked for a fully connected graph. The partition
      ! subroutine checks for this case so it should never be called
      a_n1 = a_n
      a_n2 = 0
      partition(1:a_n) = (/ (i,i=1,a_n) /)
      return
   end if

   ! Entries in level i will form first possible partition
   best_sep_start = i + 1
   a_n1 = work(distance_ptr+a_n+i+1) - 1
   a_n2 = a_n - work(distance_ptr+a_n+i+1+ww) + 1
   a_weight_1 = p1sz
   a_weight_2 = p2sz
   a_weight_sep = sepsz
   call cost_function(p1sz,p2sz,sepsz,sumweight,balance_tol,imbal,&
      options%cost_function,bestval)

   ! Search for best separator using tau

   do j = i + 1, num_levels_nend - ww - 2
      p1sz = p1sz + work(dptr+j)
      sepsz = work(dptr+j+1)
      do k = 1, ww - 1
         sepsz = sepsz + work(dptr+j+1+k)
      end do
      p2sz = sumweight - sepsz - p1sz
      if (p2sz.eq.0) exit
      call cost_function(p1sz,p2sz,sepsz,sumweight,balance_tol,imbal,&
         options%cost_function,val)
      if (val.lt.bestval) then
         bestval = val
         best_sep_start = j + 1
         a_n1 = work(distance_ptr+a_n+j+1) - 1
         a_n2 = a_n - work(distance_ptr+a_n+j+1+ww) + 1
         a_weight_1 = p1sz
         a_weight_2 = p2sz
         a_weight_sep = sepsz
      end if
   end do

   if (imbal .and. use_multilevel_copy .and. options%partition_method.ge.2) &
         then
      if (real(max(a_weight_1,a_weight_2)) / min(a_weight_1,a_weight_2) .gt. &
            balance_tol) then
         use_multilevel = .true.
         return
      end if
   end if

   ! Rearrange partition
   ! Entries in partition 1
   j = 1
   do i = 1, work(distance_ptr+a_n+best_sep_start) - 1
      partition(j) = work(distance+i)
      j = j + 1
   end do

   ! Entries in partition 2
   do i = work(distance_ptr+a_n+best_sep_start+ww), a_n
      partition(j) = work(distance+i)
      j = j + 1
   end do

   ! Entries in separator
   do i = work(distance_ptr+a_n+best_sep_start), &
         work(distance_ptr+a_n+best_sep_start+ww) - 1
      partition(j) = work(distance+i)
      j = j + 1
   end do

   info = 0
   if (printi .or. printd) &
      call nd_print_message(info,unit_diagnostics,'nd_half_level_set')

end subroutine nd_half_level_set

!
! Partition the matrix using the level set method
!
subroutine nd_level_set(a_n, a_ne, a_ptr, a_row, a_weight, sumweight, level,  &
      a_n1, a_n2, a_weight_1, a_weight_2, a_weight_sep, partition, work,      &
      options, info, band, depth, use_multilevel)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(in) :: a_ptr(a_n)
   integer, intent(in) :: a_row(a_ne)
   integer, intent(in) :: a_weight(a_n)
   integer, intent(in) :: sumweight ! sum of entries in a_weight
   integer, intent(in) :: level ! current nested dissection level
   integer, intent(out) :: a_n1, a_n2 ! size of the two submatrices
   integer, intent(out) :: a_weight_1, a_weight_2, a_weight_sep ! Weighted
      ! size of partitions and separator
   integer, intent(out) :: partition(a_n) ! First a_n1 entries will contain
      ! list of (local) indices in partition 1; next a_n2 entries will
      ! contain list of (local) entries in partition 2; entries in
      ! separator are listed at the end
   integer, intent(out) :: work(9*a_n+sumweight)
   type (nd_options), intent(in) :: options
   integer, intent(inout) :: info
   real(wp), intent(out) :: band ! If level = 0, then on output
      ! band = 100*L/a_n, where L is the size of the largest levelset
   real(wp), intent(out) :: depth ! If level = 0, then on output
      ! band = num_levels_nend
   logical, intent(inout) :: use_multilevel ! are we allowed to use a multilevel
      ! partitioning strategy

   integer :: unit_diagnostics ! unit on which to print diagnostics
   integer :: nstrt, nend
   logical :: printi, printd
   integer :: level_p, level_ptr_p, work_p
   integer :: num_levels_nend ! no. levels in structure rooted at nend
   integer :: num_entries ! no. entries in level structure rooted at nend
   integer :: best_sep_start
   integer :: i, j, p1sz, p2sz, sepsz, lwidth
   integer :: mindeg, degree, max_search
   real(wp) :: bestval
   real(wp) :: val
   real(wp) :: balance_tol
   logical :: imbal

   unit_diagnostics = options%unit_diagnostics
   printi = (options%print_level.eq.1 .and. unit_diagnostics.ge.0)
   printd = (options%print_level.ge.2 .and. unit_diagnostics.ge.0)

   if (printi .or. printd) then
      write (unit_diagnostics,'(a)') ' '
      write (unit_diagnostics,'(a)') 'Use one-sided level set method'
   end if

   ! Initialize return values
   band = -1
   depth = -1

   ! If we're going to use multilevel regardless, immediate return
   if (options%partition_method.eq.1 .and. use_multilevel) return
   if (options%partition_method.gt.1 .and. level.gt.0 .and. use_multilevel) &
      return

   ! Initialize internal variables
   balance_tol = max(real(1.0,wp),options%balance)
   imbal = (balance_tol .le. (sumweight-2))

   ! Find pseudoperipheral nodes nstart and nend, and the level structure
   ! rooted at nend
   level_ptr_p = 0 ! size a_n
   level_p = level_ptr_p + a_n ! size a_n
   work_p = level_p + a_n ! size 2*a_n
   mindeg = sumweight + 1
   do i = 1, a_n
      degree = 0
      do j = a_ptr(i), nd_get_ptr(i+1, a_n, a_ne, a_ptr)-1
         degree = degree + a_weight(a_row(j))
      end do
      if (degree.lt.mindeg) then
         mindeg = degree
         nstrt = i
      end if
   end do
   max_search = 5

   call nd_find_pseudo(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,          &
      work(level_ptr_p+1:level_ptr_p+a_n), work(level_p+1:level_p+a_n),       &
      nstrt, nend, max_search, work(work_p+1:work_p+2*a_n), num_levels_nend,  &
      num_entries, lwidth)

   if (num_entries.lt.a_n) then
      ! matrix is separable
      a_n1 = num_entries
      a_n2 = a_n - a_n1
      work(work_p+1:work_p+a_n) = 0
      a_weight_1 = 0
      do i = 1, num_entries
         j = work(level_p+i)
         partition(i) = j
         work(work_p+j) = 1
         a_weight_1 = a_weight_1 + a_weight(j)
      end do
      j = num_entries + 1
      a_weight_2 = 0
      do i = 1, a_n
         if (work(work_p+i).eq.0) then
            partition(j) = i
            a_weight_2 = a_weight_2 + a_weight(j)
            j = j + 1
         end if
      end do
      a_weight_sep = 0
      if (level.eq.0) band = -real(lwidth,wp)

      use_multilevel = .false. ! Will be reset for each component
      return
   end if

   if (level.eq.0) then
      band = (100.0_wp * lwidth) / sumweight
      depth = (100.0_wp * num_levels_nend) / sumweight
   end if

   if (options%partition_method.le.0 .or. options%stop_coarsening2.le.0) &
      use_multilevel = .false.
   if (options%partition_method.ge.2 .and. use_multilevel) then
      if ( (100.0_wp*lwidth)/sumweight.le.3.0 .or. a_n.lt. options%ml_call ) &
         use_multilevel = .false.
   end if

   if (use_multilevel) return

   if (num_levels_nend.le.2) then
      ! Not possible to find separator
      ! This can only be invoked for a full connected graph. The partition
      ! subroutine checks for this case so it should never be called
      a_n1 = a_n
      a_n2 = 0
      partition(1:a_n) = (/ (i,i=1,a_n) /)
      return
   end if

   ! Calculate total weight in each level set
   work(work_p+1:work_p+num_levels_nend) = 0
   do i = 1, num_levels_nend - 1
      do j = work(level_ptr_p+i), work(level_ptr_p+i+1) - 1
         work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
      end do
   end do
   i = num_levels_nend
   do j = work(level_ptr_p+i), a_n
      work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
   end do


   ! First possible separator contains all of the nodes in level 2
   p1sz = work(work_p+1)
   sepsz = work(work_p+2)
   p2sz = sum(work(work_p+3:work_p+num_levels_nend))
   a_weight_1 = p1sz
   a_weight_2 = p2sz
   a_weight_sep = sepsz
   best_sep_start = 2
   a_n1 = work(level_ptr_p+2) - 1
   a_n2 = a_n - work(level_ptr_p+3) + 1
   call cost_function(p1sz,p2sz,sepsz,sumweight,balance_tol,imbal,&
      options%cost_function,bestval)

   ! Search for best separator using tau
   do j = 2, num_levels_nend - 4
      p1sz = p1sz + work(work_p+j)
      sepsz = work(work_p+j+1)
      p2sz = p2sz - work(work_p+j+1)
      call cost_function(p1sz,p2sz,sepsz,sumweight,balance_tol,imbal,&
         options%cost_function,val)
      if (val.lt.bestval) then
         bestval = val
         best_sep_start = j + 1
         a_n1 = work(level_ptr_p+j+1) - 1
         a_n2 = a_n - work(level_ptr_p+j+2) + 1
         a_weight_1 = p1sz
         a_weight_2 = p2sz
         a_weight_sep = sepsz
      end if
   end do


   ! Rearrange partition
   ! Entries in partition 1
   j = 1
   do i = 1, work(level_ptr_p+best_sep_start) - 1
      partition(j) = work(level_p+i)
      j = j + 1
   end do

   ! Entries in partition 2
   do i = work(level_ptr_p+best_sep_start+1), a_n
      partition(j) = work(level_p+i)
      j = j + 1
   end do

   ! Entries in separator
   do i = work(level_ptr_p+best_sep_start), work(level_ptr_p+ &
         best_sep_start+1) - 1
      partition(j) = work(level_p+i)
      j = j + 1
   end do

   info = 0
   if (printi .or. printd) &
      call nd_print_message(info,unit_diagnostics,'nd_level_set')

end subroutine nd_level_set

!
! Find pseudoperipheral pairs of nodes for irreducible graph
!
! Based on MC60HD
subroutine nd_find_pseudo(a_n, a_ne, a_ptr, a_row, a_weight, sumweight,        &
      level_ptr, level, nstrt, nend, max_search, work, num_levels, num_entries,&
      lwidth)
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(in) :: a_ptr(a_n)
   integer, intent(in) :: a_row(a_ne)
   integer, intent(in) :: a_weight(a_n)
   integer, intent(in) :: sumweight ! sum of entries in a_weight
   integer, intent(out) :: level_ptr(a_n) ! On output level_ptr(i) contains
      ! position in level that entries for level i start.
   integer, intent(out) :: level(a_n) ! On output level contains lists of
      ! rows according to the level set that they are in
   integer, intent(inout) :: nstrt ! Starting pseudoperipheral node
   integer, intent(out) :: nend ! End pseudoperipheral node
   integer, intent(in) :: max_search
   integer, target, intent(out) :: work(2*a_n)
   integer, intent(out) :: num_levels
   integer, intent(out) :: num_entries ! number of entries in level structure
   integer, intent(out) :: lwidth

   integer :: i, j, k, l, ll
   integer :: mindeg, maxdep, main, lwidth1
   integer, dimension(:), pointer :: mask, list
   integer :: nstop ! ending pseudoperipheral node
   integer :: node ! Index of graph node
   integer :: lsize ! size of final levelset
   integer :: nlsize ! no. nodes of differing degrees in final level set
   integer :: minwid ! minimum levelset width
   integer :: nlvl ! number of levels in level structure
   integer :: minwid1

   j = 0
   mask => work(1:a_n)
   list => work(a_n+1:2*a_n)
   ! Initialise mask(:) and list(:)
   mask(:) = 1
   list(:) = 0
   level_ptr(:) = 0


   ! First guess for starting node is input nstrt

   ! Generate level structure for node nstrt
   call nd_level_struct(nstrt, a_n, a_ne, a_ptr, a_row, &
      mask, level_ptr, level, maxdep, lwidth, num_entries)
   if (num_entries.lt.a_n) then
      ! matrix is separable
      num_levels = maxdep
      return
   end if

   nstop = level(a_n)
   main_loop: &
   do main = 1, min(a_n,10) ! Provides max # itr to try and improve with

      ! Store nodes in the final level set and their (weighted) degrees
      lsize = 0
      j = level_ptr(maxdep)
      do i = j, a_n
         node = level(i)
         lsize = lsize + 1
         list(lsize) = node
         level_ptr(node) = 0
         do l = a_ptr(node), nd_get_ptr(node+1, a_n, a_ne, a_ptr)-1
            level_ptr(node) = level_ptr(node) + a_weight(a_row(l))
         end do
      end do

      ! Choose at most max_search nodes
      do nlsize = 1, max_search
         ! Look for candiate with least degree
         mindeg = sumweight + 1
         ! mindeg = -1
         do i = nlsize, lsize
            if (level_ptr(list(i)).lt.mindeg) then
               j = i
               mindeg = level_ptr(list(i))
            end if
         end do
         ! Jump out of loop if no candidates left
         if (mindeg.eq.sumweight+1) exit
         ! Swap chose candidate to next position
         node = list(j)
         list(j) = list(nlsize)
         list(nlsize) = node
         ! Rule out the neighbours of the chosen node
         do i = a_ptr(node), nd_get_ptr(node+1, a_n, a_ne, a_ptr)-1
            level_ptr(a_row(i)) = sumweight + 1
         end do
      end do
      nlsize = nlsize - 1

      ! Loop over nodes in list
      minwid = huge(a_n)
      minwid1 = huge(a_n)

      do i = 1, nlsize
         node = list(i)

         ! Form rooted level structures for node
         call nd_level_struct(node, a_n, a_ne, a_ptr, a_row, &
            mask, level_ptr, level, nlvl, lwidth, num_entries)
         lwidth1 = 0
         do k = 1, nlvl - 1
            do l = level_ptr(k), level_ptr(k+1) - 1
               ll = level(l)
               lwidth1 = lwidth1 + k*a_weight(ll)
            end do
         end do
         do l = level_ptr(nlvl), a_n
            ll = level(l)
            lwidth1 = lwidth1 + nlvl*a_weight(ll)
         end do

         if (nlvl.gt.maxdep) then
            ! Level structure of greater depth. Begin a new iteration.
            nstrt = node
            maxdep = nlvl
            cycle main_loop
         else
            if (lwidth1.lt.minwid1) then
               nstop = node
               minwid = lwidth
               minwid1 = lwidth1
            end if
         end if
      end do
      ! If we reach this point, there has been no improvement, stop trying
      exit main_loop
   end do main_loop

   if (nstop.ne.node) then
      call nd_level_struct(node, a_n, a_ne, a_ptr, a_row, &
         mask, level_ptr, level, nlvl, lwidth, num_entries)
   end if
   num_levels = maxdep
   nend = nstop
end subroutine nd_find_pseudo

!
! Given a root, calculate the level structure of a given graph from that root
!
! Based on MC60LD
subroutine nd_level_struct(root, a_n, a_ne, a_ptr, a_row, mask, level_ptr, &
      level, num_levels, lwidth, num_entries)
   integer, intent(in) :: root ! Root node for level structure
   integer, intent(in) :: a_n
   integer, intent(in) :: a_ne
   integer, intent(in) :: a_ptr(a_n)
   integer, intent(in) :: a_row(a_ne)
   integer, intent(inout) :: mask(a_n) ! Must be all 1s on entry, returned to
      ! all 1s on exit. mask(node) > 0 for all visible nodes.
   integer, intent(out) :: level_ptr(a_n) ! On output level_ptr(i) contains
      ! position in level that entries for level i start.
   integer, intent(out) :: level(a_n) ! On output level contains lists of
      ! rows according to the level set that they are in
   integer, intent(out) :: num_levels ! On output num_levels contains the
      ! number of levels
   integer, intent(out) :: lwidth ! On output, contains the width of the
      ! structure
   integer, intent(out) :: num_entries ! On output, contains number of
      ! entries in the tree structure containing root

   integer :: lvlend ! End of previous level in level
   integer :: lnbr ! Next position in level
   integer :: lbegin ! Beginning of present level in level
   integer :: lw ! Level width
   integer :: nbr ! neighbour node
   integer :: node ! current node
   integer :: nvars
   integer :: i, j

   ! Establish level 1
   level_ptr(:) = 0
   mask(root) = -1
   level(1) = root
   lvlend = 0
   nvars = 0
   lnbr = 1
   lwidth = 1
   do num_levels = 1, a_n
      ! Generate next level by finding all unmasked neighbours of nodes in
      ! present level
      lbegin = lvlend + 1
      lvlend = lnbr
      level_ptr(num_levels) = lbegin
      lw = 0
      do i = lbegin, lvlend
         node = level(i)
         do j = a_ptr(node), nd_get_ptr(node+1, a_n, a_ne, a_ptr)-1
            nbr = a_row(j)
            if (mask(nbr).eq.1) then
               lnbr = lnbr + 1
               level(lnbr) = nbr
               mask(nbr) = -1
               lw = lw + 1
            end if
         end do
      end do
      lwidth = max(lwidth,lw)
      nvars = nvars + lw
      ! If no neighbours found, we are done
      if (lnbr.eq.lvlend) exit
   end do
   ! Reset mask
   do i = 1, lnbr
      mask(level(i)) = 1
   end do
   num_entries = lnbr
end subroutine nd_level_struct


! ---------------------------------------------------
! nd_distance
! ---------------------------------------------------
! Given two level structures, calculate the difference in each nodes
! distance
! from the start and end node
subroutine nd_distance(a_n,num_levels_nend,level_ptr_nend,level_nend, &
    num_levels_nstrt,level_ptr_nstrt,level_nstrt,distance_ptr,distance)
  integer, intent(in) :: a_n ! order of matrix
  integer, intent(in) :: num_levels_nend ! number of levels with root
  ! nend
  integer, intent(in) :: level_ptr_nend(a_n) ! level_ptr(i) contains
  ! position in level that entries for level i start (root = nend)
  integer, intent(in) :: level_nend(a_n) ! Contains lists of rows
  ! according to the level set that they are in (root = nend)
  integer, intent(in) :: num_levels_nstrt ! no. of levels with root
  ! nstrt
  integer, intent(inout) :: level_ptr_nstrt(a_n) ! level_ptr(i)
  ! contains
  ! position in level that entries for level i start (root = nstrt)
  ! Reused during subroutine
  integer, intent(in) :: level_nstrt(a_n) ! Contains lists of rows
  ! according to the level set that they are in (root = nstrt)
  integer, intent(out) :: distance_ptr(2*a_n-1) ! distance(i) contains
  ! position in distance where entries with distance i-a_n
  integer, intent(out) :: distance(a_n) ! Contains lists of rows
  ! ordered
  ! according to their distance

  ! ---------------------------------------------
  ! Local variables
  integer :: j, k
  integer :: dptr ! used as offset into distance_ptr
  integer :: lev ! stores current level

  dptr = a_n
  distance_ptr(dptr+1-a_n:dptr+a_n-1) = 0
  distance(1:a_n) = 0

  ! set distance(i) to hold the level that i belongs to in the structure
  ! rooted at nend (= distance from end node + 1)
  do lev = 1, num_levels_nend - 1
    do j = level_ptr_nend(lev), level_ptr_nend(lev+1) - 1
      distance(level_nend(j)) = lev
    end do
  end do
  do j = level_ptr_nend(num_levels_nend), a_n
    distance(level_nend(j)) = num_levels_nend
  end do

  ! now consider level structure rooted at start node
  do lev = 1, num_levels_nstrt - 1
    do j = level_ptr_nstrt(lev), level_ptr_nstrt(lev+1) - 1
      distance(level_nstrt(j)) = distance(level_nstrt(j)) - lev
      k = distance(level_nstrt(j))
      distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
    end do
  end do
  do j = level_ptr_nstrt(num_levels_nstrt), a_n
    distance(level_nstrt(j)) = distance(level_nstrt(j)) - lev
    k = distance(level_nstrt(j))
    distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
  end do

  ! Copy distance into level_ptr_nstrt to save memory
  level_ptr_nstrt(1:a_n) = distance(1:a_n)

  ! Set distance_ptr(i) to point one place to the right of where entries
  ! with distance i will be
  do j = 1 - a_n, a_n - 1
    if (distance_ptr(dptr+j).gt.0) then
      exit
    else
      distance_ptr(dptr+j) = 1
    end if
  end do
  distance_ptr(dptr+j) = distance_ptr(dptr+j) + 1
  do k = j + 1, a_n - 1
    distance_ptr(dptr+k) = distance_ptr(dptr+k) + distance_ptr(dptr+k-1)
  end do

  ! Set up lists of rows with same value of distance
  do j = 1, a_n
    k = level_ptr_nstrt(j)
    distance_ptr(dptr+k) = distance_ptr(dptr+k) - 1
    distance(distance_ptr(dptr+k)) = j
  end do


end subroutine nd_distance

end module spral_nd_partition
