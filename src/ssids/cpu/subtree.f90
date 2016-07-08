module spral_ssids_cpu_subtree
   use, intrinsic :: iso_c_binding
   use spral_ssids_alloc, only : smalloc_setup, smfreeall
   use spral_ssids_contrib, only : contrib_type
   use spral_ssids_cpu_iface ! fixme only
   use spral_ssids_cpu_solve ! fixme only
   use spral_ssids_datatypes
   use spral_ssids_inform, only : ssids_inform_base
   use spral_ssids_subtree, only : symbolic_subtree_base, numeric_subtree_base
   implicit none

   private
   public :: cpu_symbolic_subtree, construct_cpu_symbolic_subtree
   public :: cpu_numeric_subtree

   type, extends(symbolic_subtree_base) :: cpu_symbolic_subtree
      integer :: n
      integer :: nnodes
      integer(long) :: nfactor ! Total number of entries in L (expected)
      integer, dimension(:), pointer :: sptr
      integer(long), dimension(:), pointer :: rptr
      integer, dimension(:), pointer :: rlist
      integer, dimension(:), pointer :: nptr
      integer, dimension(:,:), pointer :: nlist
      type(C_PTR) :: csubtree
   contains
      procedure :: factor
      final :: symbolic_final
   end type cpu_symbolic_subtree

   type, extends(numeric_subtree_base) :: cpu_numeric_subtree
      logical(C_BOOL) :: posdef
      type(cpu_symbolic_subtree), pointer :: symbolic
      type(smalloc_type), pointer :: alloc=>null() ! Linked list of memory pages
         ! pointed to by nodes variable
      type(node_type), dimension(:), allocatable :: nodes
      type(cpu_node_data), dimension(:), allocatable :: cnodes
      type(contrib_type), pointer :: contrib
      type(C_PTR) :: csubtree
   contains
      procedure :: get_contrib
      procedure :: solve_fwd
      procedure :: solve_diag
      procedure :: solve_diag_bwd
      procedure :: solve_bwd
      procedure :: enquire_posdef
      procedure :: enquire_indef
      procedure :: alter
      final :: numeric_final
   end type cpu_numeric_subtree

   interface
      type(C_PTR) function c_create_symbolic_subtree(nnodes, sptr, sparent, &
            rptr, rlist)&
            bind(C, name="spral_ssids_cpu_create_symbolic_subtree")
         use, intrinsic :: iso_c_binding
         implicit none
         integer(C_INT), value :: nnodes
         integer(C_INT), dimension(nnodes+1), intent(in) :: sptr
         integer(C_INT), dimension(nnodes+1), intent(in) :: sparent
         integer(C_LONG), dimension(nnodes+1), intent(in) :: rptr
         integer(C_INT), dimension(*), intent(in) :: rlist
      end function c_create_symbolic_subtree

      subroutine c_destroy_symbolic_subtree(subtree) &
            bind(C, name="spral_ssids_cpu_destroy_symbolic_subtree")
         use, intrinsic :: iso_c_binding
         implicit none
         type(C_PTR), value :: subtree
      end subroutine c_destroy_symbolic_subtree
      subroutine c_factor_cpu(pos_def, subtree, n, nnodes, nodes, aval, &
            scaling, alloc, options, stats) &
            bind(C, name="spral_ssids_factor_cpu_dbl")
         use, intrinsic :: iso_c_binding
         import :: cpu_node_data, cpu_factor_options, cpu_factor_stats
         implicit none
         logical(C_BOOL), value :: pos_def
         type(C_PTR), value :: subtree
         integer(C_INT), value :: n
         integer(C_INT), value :: nnodes
         type(cpu_node_data), dimension(nnodes), intent(inout) :: nodes
         real(C_DOUBLE), dimension(*), intent(in) :: aval
         type(C_PTR), value :: scaling
         type(C_PTR), value :: alloc
         type(cpu_factor_options), intent(in) :: options
         type(cpu_factor_stats), intent(out) :: stats
      end subroutine c_factor_cpu

      type(C_PTR) function c_create_cpu_subtree(posdef, symbolic_subtree, &
            nnodes, nodes) &
            bind(C, name="spral_ssids_create_cpu_subtree_dbl")
         use, intrinsic :: iso_c_binding
         import :: cpu_node_data
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: symbolic_subtree
         integer(C_INT), value :: nnodes
         type(cpu_node_data), dimension(nnodes), intent(inout) :: nodes
      end function c_create_cpu_subtree

      subroutine c_destroy_cpu_subtree(posdef, subtree) &
            bind(C, name="spral_ssids_destroy_cpu_subtree_dbl")
         use, intrinsic :: iso_c_binding
         implicit none
         logical(C_BOOL), value :: posdef
         type(C_PTR), value :: subtree
      end subroutine c_destroy_cpu_subtree
   end interface

contains

function construct_cpu_symbolic_subtree(n, nnodes, sptr, sparent, rptr, &
      rlist, nptr, nlist) result(this)
   class(cpu_symbolic_subtree), pointer :: this
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), target, intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), target, intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), target, intent(in) :: rlist
   integer, dimension(nnodes+1), target, intent(in) :: nptr
   integer, dimension(2,nptr(nnodes+1)-1), target, intent(in) :: nlist

   integer :: node
   integer(long) :: blkm, blkn
   integer :: st

   ! Allocate output
   allocate(this, stat=st)
   if(st.ne.0) then
      nullify(this)
      return
   endif

   ! Store basic details
   this%n = n
   this%nnodes = nnodes
   this%nptr => nptr
   this%nlist => nlist
   this%sptr => sptr
   this%rptr => rptr
   this%rlist => rlist

   ! Count size of factors
   this%nfactor = 0
   do node = 1, nnodes
      blkm = rptr(node+1) - rptr(node)
      blkn = sptr(node+1) - sptr(node)
      this%nfactor = this%nfactor + blkm*blkn
   end do

   ! Call C++ subtree analyse
   this%csubtree = &
      c_create_symbolic_subtree(nnodes, sptr, sparent, rptr, rlist)

end function construct_cpu_symbolic_subtree

subroutine symbolic_final(this)
   type(cpu_symbolic_subtree) :: this

   call c_destroy_symbolic_subtree(this%csubtree)
end subroutine symbolic_final

function factor(this, posdef, aval, options, inform, scaling)
   class(numeric_subtree_base), pointer :: factor
   class(cpu_symbolic_subtree), target, intent(inout) :: this
   logical(C_BOOL), intent(in) :: posdef
   real(wp), dimension(*), intent(in) :: aval
   class(ssids_options), intent(in) :: options
   class(ssids_inform_base), intent(inout) :: inform
   real(wp), dimension(*), target, optional, intent(in) :: scaling

   type(cpu_numeric_subtree), pointer :: cpu_factor
   type(cpu_factor_options) :: coptions
   type(cpu_factor_stats) :: cstats
   type(C_PTR) :: cscaling
   integer :: st

   ! Leave output as null until successful exit
   nullify(factor)

   ! Allocate cpu_factor for output
   allocate(cpu_factor, stat=st)
   if(st.ne.0) goto 10
   cpu_factor%symbolic => this

   ! Setup memory allocator
   ! * options%multiplier * n             integers (for nodes(:)%perm)
   ! * options%multiplier * (nfactor+2*n) reals    (for nodes(:)%lcol)
   allocate(cpu_factor%alloc, stat=st)
   if(st.ne.0) goto 10
   call smalloc_setup(cpu_factor%alloc, &
      max(this%n+0_long, int(options%multiplier*this%n, kind=long)), &
      max(this%nfactor+2*this%n, &
         int(options%multiplier*real(this%nfactor,wp)+2*this%n,kind=long)), st)
   if(st.ne.0) goto 10

   ! Allocate Fortran nodes structure
   allocate(cpu_factor%nodes(this%nnodes+1), stat=st)
   if(st.ne.0) goto 10
   cpu_factor%nodes(:)%ndelay = 0

   ! Allocate cnodes and setup for main call
   allocate(cpu_factor%cnodes(this%nnodes+1), stat=st)
   if(st.ne.0) goto 10
   call setup_cpu_data(this%nnodes, cpu_factor%cnodes, this%nptr, this%nlist)

   ! Setup C++ data structure
   cpu_factor%posdef = posdef
   cpu_factor%csubtree = &
      c_create_cpu_subtree(posdef, this%csubtree, this%nnodes, cpu_factor%cnodes)

   ! Call C++ factor routine
   cscaling = C_NULL_PTR
   if(present(scaling)) cscaling = C_LOC(scaling)
   call cpu_copy_options_in(options, coptions)
   call c_factor_cpu(posdef, cpu_factor%csubtree, this%n, &
      this%nnodes, cpu_factor%cnodes, aval, cscaling, C_LOC(cpu_factor%alloc), &
      coptions, cstats)

   ! Extract to Fortran data structures
   call extract_cpu_data(this%nnodes, cpu_factor%cnodes, cpu_factor%nodes, &
      cstats, inform)

   ! Succes, set result and return
   factor => cpu_factor
   return

   ! Allocation error handler
   10 continue
   inform%flag = SSIDS_ERROR_ALLOCATION
   inform%stat = st
   deallocate(cpu_factor, stat=st)
   return
end function factor

subroutine numeric_final(this)
   type(cpu_numeric_subtree) :: this

   call c_destroy_cpu_subtree(this%posdef, this%csubtree)

   call smfreeall(this%alloc)
   deallocate(this%alloc)
   nullify(this%alloc)
end subroutine numeric_final

function get_contrib(this)
   class(contrib_type), pointer :: get_contrib
   class(cpu_numeric_subtree), intent(in) :: this

   get_contrib => this%contrib
end function get_contrib

subroutine solve_fwd(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform
 
   logical :: fposdef
   integer :: blkm
   integer :: blkn
   integer :: nd
   integer :: nelim
   integer :: node
   real(wp), dimension(:), allocatable :: xlocal
   integer, dimension(:), allocatable :: map

   fposdef = this%posdef
   associate(symbolic=>this%symbolic, nodes=>this%nodes)
      allocate(xlocal(nrhs*symbolic%n), map(symbolic%n), stat=inform%stat)
      if(inform%stat.ne.0) return

      ! Forwards solve Ly = b
      do node = 1, symbolic%nnodes
         nelim = nodes(node)%nelim
         if (nelim.eq.0) cycle
         nd = nodes(node)%ndelay
         blkn = symbolic%sptr(node+1) - symbolic%sptr(node) + nd
         blkm = int(symbolic%rptr(node+1) - symbolic%rptr(node)) + nd
         
         if(nrhs.eq.1) then
            call solve_fwd_one(fposdef, symbolic%rlist(symbolic%rptr(node):), &
               x, blkm, blkn, nelim, nd, &
               nodes(node)%lcol, nodes(node)%perm, xlocal, map)
         else
            call solve_fwd_mult(fposdef, symbolic%rlist(symbolic%rptr(node):), &
               nrhs, x, ldx, blkm, blkn, nelim, nd, &
               nodes(node)%lcol, nodes(node)%perm, xlocal, map)
         end if
      end do
   end associate
end subroutine solve_fwd

subroutine solve_diag(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform

   logical :: fposdef
   integer :: blkm
   integer :: blkn
   integer :: nd
   integer :: nelim
   integer :: node

   fposdef = this%posdef
   associate(symbolic=>this%symbolic, nodes=>this%nodes)
      do node = symbolic%nnodes, 1, -1
         nelim = nodes(node)%nelim
         if (nelim.eq.0) cycle
         nd = nodes(node)%ndelay
         blkn = symbolic%sptr(node+1) - symbolic%sptr(node) + nd
         blkm = int(symbolic%rptr(node+1) - symbolic%rptr(node)) + nd
         
         if(nrhs.eq.1) then
            call solve_diag_one(x, nelim, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
               nodes(node)%ismptr%imem(nodes(node)%ismsa)) ! node%perm
         else
            call solve_diag_mult(nrhs, x, ldx, nelim, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
               nodes(node)%ismptr%imem(nodes(node)%ismsa)) ! node%perm
         end if
      end do
   end associate
end subroutine solve_diag

subroutine solve_diag_bwd(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform

   logical :: fposdef

   fposdef = this%posdef
   call subtree_bwd_solve(this%symbolic%nnodes, 1, SSIDS_SOLVE_JOB_DIAG_BWD, &
      fposdef, this%symbolic%nnodes, this%nodes, this%symbolic%sptr, &
      this%symbolic%rptr, this%symbolic%rlist, nrhs, x, ldx, inform%stat)
end subroutine solve_diag_bwd

subroutine solve_bwd(this, nrhs, x, ldx, inform)
   class(cpu_numeric_subtree), intent(inout) :: this
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   class(ssids_inform_base), intent(inout) :: inform

   logical :: fposdef

   fposdef = this%posdef
   call subtree_bwd_solve(this%symbolic%nnodes, 1, SSIDS_SOLVE_JOB_BWD, &
      fposdef, this%symbolic%nnodes, this%nodes, this%symbolic%sptr, &
      this%symbolic%rptr, this%symbolic%rlist, nrhs, x, ldx, inform%stat)
end subroutine solve_bwd

subroutine enquire_posdef(this, d)
   class(cpu_numeric_subtree), target, intent(in) :: this
   real(wp), dimension(*), intent(out) :: d

   integer :: blkn, blkm
   integer(long) :: i
   integer :: j
   integer :: node
   integer :: piv

   type(node_type), pointer :: nptr

   associate(symbolic => this%symbolic)
      piv = 1
      do node = 1, symbolic%nnodes
         nptr => this%nodes(node)
         blkn = symbolic%sptr(node+1) - symbolic%sptr(node)
         blkm = int(symbolic%rptr(node+1) - symbolic%rptr(node))
         i = 1
         do j = 1, blkn
            d(piv) = nptr%lcol(i)
            i = i + blkm + 1
            piv = piv + 1
         end do
      end do 
   end associate
end subroutine enquire_posdef

subroutine enquire_indef(this, invp, piv_order, d)
   class(cpu_numeric_subtree), target, intent(in) :: this
   integer, dimension(*), intent(in) :: invp
   integer, dimension(*), optional, intent(out) :: piv_order
   real(wp), dimension(2,*), optional, intent(out) :: d

   integer :: blkn, blkm
   integer :: j, k
   integer :: nd
   integer :: node
   integer(long) :: offset
   integer :: piv

   type(node_type), pointer :: nptr

   associate(symbolic => this%symbolic)
      piv = 1
      do node = 1, symbolic%nnodes
         nptr => this%nodes(node)
         j = 1
         nd = nptr%ndelay
         blkn = symbolic%sptr(node+1) - symbolic%sptr(node) + nd
         blkm = int(symbolic%rptr(node+1) - symbolic%rptr(node)) + nd
         offset = blkm*(blkn+0_long)
         do while(j .le. nptr%nelim)
            if (nptr%lcol(offset+2*j).ne.0) then
               ! 2x2 pivot
               if(present(piv_order))  then
                  k = invp( nptr%perm(j) )
                  piv_order(k) = -piv
                  k = invp( nptr%perm(j+1) )
                  piv_order(k) = -(piv+1)
               end if
               if(present(d)) then
                  d(1,piv) = nptr%lcol(offset+2*j-1)
                  d(2,piv) = nptr%lcol(offset+2*j)
                  d(1,piv+1) = nptr%lcol(offset+2*j+1)
                  d(2,piv+1) = 0
               end if
               piv = piv + 2
               j = j + 2
            else
               ! 1x1 pivot
               if(present(piv_order)) then
                  k = invp( nptr%perm(j) )
                  piv_order(k) = piv
               end if
               if(present(d)) then
                  d(1,piv) = nptr%lcol(offset+2*j-1)
                  d(2,piv) = 0
               end if
               piv = piv + 1
               j = j + 1
            end if
         end do
      end do
   end associate
end subroutine enquire_indef

subroutine alter(this, d)
   class(cpu_numeric_subtree), target, intent(inout) :: this
   real(wp), dimension(2,*), intent(in) :: d

   integer :: blkm, blkn
   integer(long) :: ip
   integer :: j
   integer :: nd
   integer :: node
   integer :: piv

   type(node_type), pointer :: nptr

   associate(symbolic => this%symbolic)
      piv = 1
      do node = 1, symbolic%nnodes
         nptr => this%nodes(node)
         nd = nptr%ndelay
         blkn = symbolic%sptr(node+1) - symbolic%sptr(node) + nd
         blkm = int(symbolic%rptr(node+1) - symbolic%rptr(node)) + nd
         ip = blkm*(blkn+0_long) + 1
         do j = 1, nptr%nelim
            nptr%lcol(ip)   = d(1,piv)
            nptr%lcol(ip+1) = d(2,piv)
            ip = ip + 2
            piv = piv + 1
         end do
      end do
   end associate
end subroutine alter

end module spral_ssids_cpu_subtree
