****************************************************************
:f:mod:`spral_ssids` - Sparse Symmetric Indefinite Direct Solver
****************************************************************
.. f:module:: spral_ssids
   :synopsis: Sparse Symmetric Indefinite Direct Solver

=======
Purpose
=======

This package solves one or more sets of :math:`n\times n`
sparse **symmetric** equations  :math:`AX = B` using a multifrontal method.
The following cases are covered:

1. :math:`A` is **indefinite**.
SSIDS computes the sparse factorization

.. math::

   A =  PLD(PL)^T

where :math:`P` is a permutation matrix, :math:`L` is unit lower triangular,
and :math:`D` is block diagonal with blocks of size :math:`1 \times 1`
and :math:`2 \times 2`.

2. :math:`A` is **positive definite**.
SSIDS computes the **sparse Cholesky factorization**

.. math::

   A =  PL(PL)^T

where :math:`P` is a permutation matrix and :math:`L` is lower triangular.
*However, as SSIDS is designed primarily for indefinite
systems, this may be slower than a dedicated Cholesky solver.*

The code optionally supports hybrid computation using one or more NVIDIA GPUs.

SSIDS returns bit-compatible results.

An option exists to scale the matrix. In this case, the factorization of
the scaled matrix  :math:`\overline{A} = S A S` is computed,
where :math:`S` is a diagonal scaling matrix.

Version history
---------------

[For detail, see ChangeLog]

2016-09-XX Version 2.0.0
   * Add support for CPU and hybrid execution.
   * Add support for 64-bit ptr(:) on input matrix :math:`A`.
   * options%presolve removed.

2014-03-17 Version 1.0.0
   * Initial release: GPU only

==============
Usage overview
==============

Solving :math:`AX=B` using SSIDS is a four stage process.

1. Call :f:subr:`ssids_analyse()` to perform a symbolic factorization, stored
   in `akeep`.
2. Call :f:subr:`ssids_factor()` to perform a numeric factorization, stored in
   `fkeep`. More than one numeric factorization can refer to the same `akeep`.
3. Call :f:subr:`ssids_solve()` to perform a solve with the factors. More than
   one solve can be performed with the same `fkeep`.
4. Once all desired solutions have been performed, free memory with
   :f:subr:`ssids_free()`.

In addition, advanced users may use the following routines:

* :f:subr:`ssids_enquire_posdef()` and :f:subr:`ssids_enquire_indef()` return
  the diagonal entries of the factors and the pivot sequence.
* :f:subr:`ssids_alter()` allows altering the diagonal entries of the factors.


.. note::

   **Bit-compatibility:**
   If used with bit-compatible BLAS and compiler options, this routine will
   return bit compatible results. That is, consecutive runs with the same data
   on the same machine produces exactly the same solution.

=================
Basic Subroutines
=================

Note: For the most efficient use of the routine, CSC format should be used
without checking.

.. f:subroutine:: ssids_analyse(check,n,ptr,row,akeep,options,inform[,order,val,topology])

   Perform the analyse (symbolic) phase of the factorization for a matrix
   supplied in :doc:`CSC format<csc_format>`. The resulting symbolic factors
   stored in :f:type:`akeep` should be passed unaltered in the following call to
   ssids_factor().

   :p logical check [in]: if true, matrix data is checked. Out-of-range entries
      are dropped and duplicate entries are summed.
   :p integer n [in]: number of columns in :math:`A`.
   :p integer(long) ptr(n+1) [in]: column pointers for :math:`A` (see :doc:`CSC format<csc_format>`).
   :p integer row(ptr(n+1)-1) [in]: row indices for :math:`A` (see :doc:`CSC format<csc_format>`).
   :p ssids_akeep akeep [out]: returns symbolic factorization, to be passed unchanged to subsequent routines.
   :p ssids_options options [in]: specifies algorithm options to be used (see [ssids:type:option]).
   :p ssids_inform inform [out]: returns information about the execution of the routine (see [ssids:type:inform]).
   :o integer order(n) [inout]: on entry a user-supplied ordering (options%ordering=0). On return, the actual ordering used (if present).
   :o real val(ptr(n+1)-1) [in]: non-zero values for :math:`A` (see :doc:`CSC format<csc_format>`). Only used if a matching-based ordering is requested.
   :o numa_region topology(*) [in]: If present, specifies the machine topology
      to be exploited. The size of the machine is the number of independent
      NUMA regions. Region `i` will use `topology(i)%nproc` threads and is
      associated with the GPUs in the array `topology(i)%gpus`, which may have
      length 0. If not present, these parameters are auto-detected using
      the hwloc library (if detected at compiled time) or the environment
      variable `OMP_NUM_THREADS` if the hwloc library is not available. See
      the :ref:`method section <ssids_method>` for details of how work is
      divided.

   .. note::

      If a user-supplied ordering is used, it may be altered by this routine,
      with the altered version returned in order(:). This version will be
      equivalent to the original ordering, except that some supernodes may have
      been amalgamated, a topographic ordering may have been applied to the tree
      and the order of columns within a supernode may have been adjusted to
      improve cache locality.

   .. note::

      A version where `ptr` is of kind default integer is also provided for
      backwards compatibility.

.. f:subroutine:: ssids_analyse_coord(n,ne,row,col,akeep,options,inform[,order,val, topology])

   As :f:subr:`ssids_analyse()`, but for coordinate data. The variant parameters
   are:

   :p integer(long) ne [in]: number of non-zero entries in :math:`A`.
   :p integer row(ne) [in]: row indices for :math:`A` (see :doc:`Coordinate format<coord_format>`).
   :p integer col(ne) [in]: column indices for :math:`A` (see :doc:`Coordinate format<coord_format>`).

.. f:subroutine::  ssids_factor(posdef,val,akeep,fkeep,options,inform[,scale,ptr,row])

   :p logical posdef [in]: true if matrix is matrix is positive-definite
   :p real val(*) [in]: non-zero values for :math:`A` in same format as for
      the call to :f:subr:`ssids_analyse()` or :f:subr:`ssids_analyse_coord()`.
   :p ssids_akeep akeep [in]: symbolic factorization returned by preceding call
      to :f:subr:`ssids_analyse()` or :f:subr:`ssids_analyse_coord()`.
   :p ssids_fkeep fkeep [inout]: returns numeric factorization, to be passed
      unchanged to subsequent routines.
   :p ssids_options options [in]: specifies algorithm options to be used
      (see [ssids:type:option]).
   :p ssids_inform inform [out]: returns information about the execution of the
      routine (see [ssids:type:inform]).
   :o real scale(n) [inout]: diagonal scaling. scale(i) contains entry
      :math:`S_ii` of :math:`S`. Must be supplied by user if
      options%scaling=0 (user-supplied scaling). On exit, return scaling used.
   :o integer(long) ptr(n+1) [in]: column pointers for :math:`A`, only required
      if :f:type:`akeep` was obtained by running :f:subr:`ssids_analyse()` with
      :f:type:`check` =.true., in which case it must be unchanged since that
      call.
   :o integer row(ptr(n+1)-1) [in]: row indices for :math:`A`, only required if
      :f:type:`akeep` was obtained by running :f:subr:`ssids_analyse()` with
      :f:type:`check` =.true., in which case it must be unchanged since that
      call.

   .. note::

      A version where `ptr` is of kind default integer is also provided for
      backwards compatibility.

.. f:subroutine:: ssids_solve(x,akeep,fkeep,options,inform[,job])
   
   Solve (for a single right-hand side) one of the following equations:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0 (or absent) | :math:`Ax=b`             |
   +---------------+--------------------------+
   | 1             | :math:`PLx=Sb`           |
   +---------------+--------------------------+
   | 2             | :math:`Dx=b`             |
   +---------------+--------------------------+
   | 3             | :math:`(PL)^TS^{-1}x=b`  |
   +---------------+--------------------------+
   | 4             | :math:`D(PL)^TS^{-1}x=b` |
   +---------------+--------------------------+

   Recall :math:`A` has been factorized as either:
   
   * :math:`SAS = (PL)(PL)^T~` (positive-definite case); or
   * :math:`SAS = (PL)D(PL)^T` (indefinite case).

   :p real x(n) [in]: right-hand side :math:`b` on entry, solution :math:`x`
      on exit.
   :p ssids_akeep akeep [in]: symbolic factorization returned by preceding
      call to :f:subr:`ssids_analyse()` or :f:subr:`ssids_analyse_coord()`.
   :p ssids_fkeep fkeep [in]: numeric factorization returned by preceding
      call to :f:subr:`ssids_factor()`.
   :p ssids_options options [in]: specifies algorithm options to be used
      (see [ssids:type:option]).
   :p ssids_inform inform [out]: returns information about the execution of the
      routine (see [ssids:type:inform]).
   :o integer job [in]: specifies equation to solve, as per above table.

.. f:subroutine:: ssids_solve(nrhs,x,ldx,akeep,fkeep,options,inform[,job])
   
   Solve (for multiple right-hand sides) one of the following equations:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0 (or absent) | :math:`AX=B`             |
   +---------------+--------------------------+
   | 1             | :math:`PLX=SB`           |
   +---------------+--------------------------+
   | 2             | :math:`DX=B`             |
   +---------------+--------------------------+
   | 3             | :math:`(PL)^TS^{-1}X=B`  |
   +---------------+--------------------------+
   | 4             | :math:`D(PL)^TS^{-1}X=B` |
   +---------------+--------------------------+

   Recall :math:`A` has been factorized as either:
   
   * :math:`SAS = (PL)(PL)^T~` (positive-definite case); or
   * :math:`SAS = (PL)D(PL)^T` (indefinite case).

   :p integer nrhs [in]: number of right-hand sides.
   :p real x(ldx,nrhs) [inout]: right-hand sides :math:`B` on entry,
      solutions :math:`X` on exit.
   :p integer ldx [in]: leading dimension of :f:type:`x`.
   :p ssids_akeep akeep [in]: symbolic factorization returned by preceding
      call to :f:subr:`ssids_analyse()` or :f:subr:`ssids_analyse_coord()`.
   :p ssids_fkeep fkeep [in]: numeric factorization returned by preceding
      call to :f:subr:`ssids_factor()`.
   :p ssids_options options [in]: specifies algorithm options to be used
      (see [ssids:type:option]).
   :p ssids_inform inform [out]: returns information about the execution of the
      routine (see [ssids:type:inform]).
   :o integer job [in]: specifies equation to solve, as per above table.

.. f:subroutine:: ssids_free([akeep,fkeep,]cuda_error)

   Frees memory and resources associated with :f:type:`akeep` and/or
   :f:type:`fkeep`, at least one of which must be present.

   :p ssids_akeep akeep [inout]: symbolic factors to be freed.
   :p ssids_fkeep fkeep [inout]: numeric factors to be freed.
   :p integer cuda_error [out]: 0 on success, or a CUDA error code on failure.

   .. warning::

      :f:subr:`ssids_free` must be called by the user. Merely deallocating
      :f:type:`akeep` or :f:type:`fkeep`, or allowing them to go out of scope
      will result in memory leaks due to underlying memory allocation performed
      using C++ and CUDA mechanisms. :f:type:`akeep` should only be deallocated
      after all associated numeric factorizations :f:type:`fkeep` have been
      freed.

====================
Advanced subroutines
====================

.. f:subroutine:: ssids_enquire_posdef(akeep,fkeep,options,inform,d)

   Return the diagonal entries of the Cholesky factor.

   :p ssids_akeep akeep [in]: symbolic factorization returned by preceding
      call to :f:subr:`ssids_analyse()` or :f:subr:`ssids_analyse_coord()`.
   :p ssids_fkeep fkeep [in]: numeric factorization returned by preceding
      call to :f:subr:`ssids_factor()`.
   :p ssids_options options [in]: specifies algorithm options to be used
      (see [ssids:type:option]).
   :p ssids_inform inform [out]: returns information about the execution of the
      routine (see [ssids:type:inform]).
   :p real d (n) [out]: returns the diagonal of :math:`L`. d(i) stores the
      entry :math:`L_{ii}`.

.. f:subroutine:: ssids_enquire_indef(akeep,fkeep,options,inform[,piv_order,d])

   Return the pivot order and/or values of :math:`D` of the Symmetric Indefinite
   Factorization.

   :p ssids_akeep akeep [in]: symbolic factorization returned by preceding
      call to :f:subr:`ssids_analyse()` or :f:subr:`ssids_analyse_coord()`.
   :p ssids_fkeep fkeep [in]: numeric factorization returned by preceding
      call to :f:subr:`ssids_factor()`.
   :p ssids_options options [in]: specifies algorithm options to be used
      (see [ssids:type:option]).
   :p ssids_inform inform [out]: returns information about the execution of the
      routine (see [ssids:type:inform]).
   :o integer piv_order (n) [out]: returns the pivot order.
      :math:`|\,\texttt{piv\_order(i)}|` gives the position of variable
      :math:`i` in the pivot order. The sign will be positive if :math:`i` is a
      :math:`1\times1` pivot, and negative if :math:`i` is
      part of a :math:`2 \times 2` pivot.
   :o real d (2,n) [out]: returns the :math:`2\times2` block diagonal of
      :math:`D`. d(1,i) stores :math:`D_{ii}` and d(2,i) stores
      :math:`D_{(i+1)i}`.

.. f:subroutine:: ssids_alter(d,akeep,fkeep,options,inform)

   Alter the entries of the diagonal factor :math:`D` for a symmetric indefinite
   factorization. The pivot order remains the same.

   :p real d (2,n) [in]: New entries of :math:`D`.
      d(1,i) stores :math:`D_{ii}` and d(2,i) stores :math:`D_{(i+1)i}`.
   :p ssids_akeep akeep [in]: symbolic factorization returned by preceding
      call to :f:subr:`ssids_analyse()` or :f:subr:`ssids_analyse_coord()`.
   :p ssids_fkeep fkeep [in]: numeric factorization returned by preceding
      call to :f:subr:`ssids_factor()`.
   :p ssids_options options [in]: specifies algorithm options to be used
      (see [ssids:type:option]).
   :p ssids_inform inform [out]: returns information about the execution of the
      routine (see [ssids:type:inform]).

=============
Derived types
=============

.. f:type:: ssids_options

   The derived data type ssids\_options is used to specify the options
   used within ``SSIDS``. The components, that are automatically given
   default values in the definition of the type, are:

   :f integer print_level [default=0]: the level of printing. The different
      levels are:

      +-----+-------------------------------------------------+
      | < 0 | No printing.                                    |
      +-----+-------------------------------------------------+
      | = 0 | Error and warning messages only.                |
      +-----+-------------------------------------------------+
      | = 1 | As 0, plus basic diagnostic printing.           |
      +-----+-------------------------------------------------+
      | > 1 | As 1, plus some additional diagnostic printing. |
      +-----+-------------------------------------------------+

   :f integer unit_diagnostics [default=6]: Fortran unit number for
      diagnostics printing. Printing is suppressed if <0.
   :f integer unit_error [default=6]: Fortran unit number for printing of
      error messages. Printing is suppressed if <0.
   :f integer unit_warning [default=6]: Fortran unit number for printing of
      warning messages. Printing is suppressed if <0.
   :f integer ordering [default=1]: Ordering method to use in analyse phase:

      +-------------+---------------------------------------------------------+
      | 0           | User-supplied ordering is used (`order` argument to     |
      |             | :f:subr:`ssids_analyse()` or                            |
      |             | :f:subr:`ssids_analyse_coord()`).                       |
      +-------------+---------------------------------------------------------+
      | 1 (default) | METIS ordering with default settings.                   |
      +-------------+---------------------------------------------------------+
      | 2           | Matching-based elimination ordering is computed (the    |
      |             | Hungarian algorithm is used to identify large           |
      |             | off-diagonal entries. A restricted METIS ordering is    |
      |             | then used that forces these on to the subdiagonal).     |
      |             |                                                         |
      |             | **Note:** This option should only be chosen for         |
      |             | indefinite systems. A scaling is also computed that may |
      |             | be used in :f:subr:`ssids_factor()` (see %scaling       |
      |             | below).                                                 |
      +-------------+---------------------------------------------------------+

   :f integer nemin [default=32]: supernode amalgamation threshold. Two
      neighbours in the elimination tree are merged if they both involve fewer
      than nemin eliminations. The default is used if nemin<1.
   :f ignore_numa [default=true]: If true, all CPUs and GPUs are treated as
      belonging to a single NUMA region.
   :f use_gpu [default=true]: Use an NVIDIA GPU if present.
   :f integer(long) min_gpu_work [default=1e10]: Minimum number of flops
      in subtree before scheduling on GPU.
   :f real max_load_inbalance [default=1.5]: Maximum permissiable load
      inbalance for leaf subtree allocations. Values less than 1.0 are treated
      as 1.0.
   :f real gpu_perf_coeff [default=1.5]: GPU perfromance coefficient. How many
      times faster a GPU is than CPU at factoring a subtree.
   :f integer scaling [default=0]: scaling algorithm to use:

      +---------------+-------------------------------------------------------+
      | <=0 (default) | No scaling (if ``scale(:)`` is not present on call to |
      |               | :f:subr:`ssids_factor()`, or user-supplied scaling (if|
      |               | ``scale(:)`` is present).                             |
      +---------------+-------------------------------------------------------+
      | =1            | Compute using weighted bipartite matching via the     |
      |               | Hungarian Algorithm (``MC64`` algorithm).             |
      +---------------+-------------------------------------------------------+
      | =2            | Compute using a weighted bipartite matching via the   |
      |               | Auction Algorithm (may be lower quality than that     |
      |               | computed using the Hungarian Algorithm, but can be    |
      |               | considerably faster).                                 |
      +---------------+-------------------------------------------------------+
      | =3            | Use matching-based ordering generated during the      |
      |               | analyse phase using options%ordering=2. The scaling   |
      |               | will be the same as that generated with %scaling= 1   |
      |               | if the matrix values have not changed. This option    |
      |               | will generate an error if a matching-based ordering   |
      |               | was not used during analysis.                         |
      +---------------+-------------------------------------------------------+
      | >=4           | Compute using the norm-equilibration algorithm of     |
      |               | Ruiz.                                                 |
      +---------------+-------------------------------------------------------+

   :f integer(long) small_subtree_threshold [default=4e6]: Maximum number of
      flops in a subtree treated as a single task. See
      :ref:`method section <ssids_small_leaf>`.
   :f integer cpu_block_size [default=256]: Block size to use for
      parallelization of large nodes on CPU resources.
   :f logical action [default=.true.]: continue factorization of singular matrix
      on discovery of zero pivot if true (a warning is issued), or abort if
      false.
   :f integer pivot_method [default=1]: Pivot method to be used on CPU, one of:

      +-------------+----------------------------------------------------------+
      | 0           | Aggressive a posteori pivoting. Cholesky-like            |
      |             | communication pattern is used, but a single failed pivot |
      |             | requires restart of node factorization and potential     |
      |             | recalculation of all uneliminated entries.               |
      +-------------+----------------------------------------------------------+
      | 1 (default) | Block a posteori pivoting. A failed pivot only requires  |
      |             | recalculation of entries within its own block column.    |
      +-------------+----------------------------------------------------------+
      | 2           | Threshold Partial Pivoting. Not parallel.                |
      +-------------+----------------------------------------------------------+

   :f real small [default=1d-20]: threshold below which an entry is treated as
      equivalent to `0.0`.
   :f real u [default=0.01]: relative pivot threshold used in symmetric
      indefinite case. Values outside of the range :math:`[0,0.5]` are treated
      as the closest value in that range.

.. f:type:: ssids_inform

   Used to return information about the progress and needs of the algorithm.

   :f integer flag: exit status of the algorithm (see table below).
   :f integer matrix_dup: number of duplicate entries encountered (if
      :f:subr:`ssids_analyse()` called with check=true, or any call to
      :f:subr:`ssids_analyse_coord()`).
   :f integer matrix_missing_diag: number of diagonal entries without an
      explicit value (if :f:subr:`ssids_analyse()` called with check=true, or
      any call to :f:subr:`ssids_analyse_coord()`).
   :f integer matrix_outrange: number of out-of-range entries encountered (if
      :f:subr:`ssids_analyse()` called with check=true, or any call to
      :f:subr:`ssids_analyse_coord()`).
   :f integer matrix_rank: (estimated) rank (structural after analyse phase,
      numerical after factorize phase).
   :f integer maxdepth: maximum depth of the assembly tree.
   :f integer maxfront: maximum front size (without pivoting after analyse
      phase, with pivoting after factorize phase).
   :f integer num_delay: number of delayed pivots. That is, the total
      number of fully-summed variables that were passed to the father node
      because of stability considerations. If a variable is passed further
      up the tree, it will be counted again.
   :f integer(long) num_factor: number of entries in :math:`L` (without pivoting
      after analyse phase, with pivoting after factorize phase).
   :f integer(long) num_flops: number of floating-point operations for Cholesky
      factorization (indefinte needs slightly more). Without pivoting after
      analyse phase, with pivoting after factorize phase.
   :f integer num_neg: number of negative eigenvalues of the matrix :math:`D`
      after factorize phase.
   :f integer num_sup: number of supernodes in assembly tree.
   :f integer num_two: number of :math:`2 \times 2` pivots used by the
      factorization (i.e. in the matrix :math:`D`).
   :f integer stat: Fortran allocation status parameter in event of allocation
      error (0 otherwise).
   :f integer cublas_error: CUBLAS error code in the event of a CUBLAS error
      (0 otherwise).
   :f integer cuda_error: CUDA error code in the event of a CUDA error
      (0 otherwise). Note that due to asynchronous execution, CUDA errors may 
      not be reported by the call that caused them.

   +-------------+-------------------------------------------------------------+
   | inform%flag | Return status                                               |
   +=============+=============================================================+
   | 0           | Success.                                                    |
   +-------------+-------------------------------------------------------------+
   | -1          | Error in sequence of calls (may be caused by failure of a   |
   |             | preceding call).                                            |
   +-------------+-------------------------------------------------------------+
   | -2          | n<0 or ne<1.                                                |
   +-------------+-------------------------------------------------------------+
   | -3          | Error in ptr(:).                                            |
   +-------------+-------------------------------------------------------------+
   | -4          | CSC format: All variable indices in one or more columns are |
   |             | out-of-range.                                               |
   |             |                                                             |
   |             | Coordinate format: All entries are out-of-range.            |
   +-------------+-------------------------------------------------------------+
   | -5          | Matrix is singular and options%action=.false.               |
   +-------------+-------------------------------------------------------------+
   | -6          | Matrix found not to be positive definite.                   |
   +-------------+-------------------------------------------------------------+
   | -7          | ptr(:) and/or row(:) not present, but required as           |
   |             | :f:subr:`ssids_analyse()` was called with check=.false,.    |
   +-------------+-------------------------------------------------------------+
   | -8          | options%ordering out of range, or options%ordering=0 and    |
   |             | order parameter not provided or not a valid permutation.    |
   +-------------+-------------------------------------------------------------+
   | -9          | options%ordering=-2 but val(:) was not supplied.            |
   +-------------+-------------------------------------------------------------+
   | -10         | ldx<n or nrhs<1.                                            |
   +-------------+-------------------------------------------------------------+
   | -11         | job is out-of-range.                                        |
   +-------------+-------------------------------------------------------------+
   | -13         | Called :f:subr:`ssids_enquire_posdef()` on indefinite       |
   |             | factorization.                                              |
   +-------------+-------------------------------------------------------------+
   | -14         | Called :f:subr:`ssids_enquire_indef()` on positive-definite |
   |             | factorization.                                              |
   +-------------+-------------------------------------------------------------+
   | -15         | options%scaling=3 but a matching-based ordering was not     |
   |             | performed during analyse phase.                             |
   +-------------+-------------------------------------------------------------+
   | -50         | Allocation error. If available, the stat parameter is       |
   |             | returned in inform%stat.                                    |
   +-------------+-------------------------------------------------------------+
   | -51         | CUDA error. The CUDA error return value is returned in      |
   |             | inform%cuda_error.                                          |
   +-------------+-------------------------------------------------------------+
   | -52         | CUBLAS error. The CUBLAS error return value is returned in  |
   |             | inform%cublas_error.                                        |
   +-------------+-------------------------------------------------------------+
   | +1          | Out-of-range variable indices found and ignored in input    |
   |             | data. inform%matrix_outrange is set to the number of such   |
   |             | entries.                                                    |
   +-------------+-------------------------------------------------------------+
   | +2          | Duplicate entries found and summed in input data.           |
   |             | inform%matrix_dup is set to the number of such entries.     |
   +-------------+-------------------------------------------------------------+
   | +3          | Combination of +1 and +2.                                   |
   +-------------+-------------------------------------------------------------+
   | +4          | One or more diagonal entries of :math:`A` are missing.      |
   +-------------+-------------------------------------------------------------+
   | +5          | Combination of +4 and +1 or +2.                             |
   +-------------+-------------------------------------------------------------+
   | +6          | Matrix is found be (structurally) singular during analyse   |
   |             | phase. This will overwrite any of the above warning flags.  |
   +-------------+-------------------------------------------------------------+
   | +7          | Matrix is found to be singular during factorize phase.      |
   +-------------+-------------------------------------------------------------+
   | +8          | Matching-based scaling found as side-effect of              |
   |             | matching-based ordering ignored                             |
   |             | (consider setting options%scaling=3).                       |
   +-------------+-------------------------------------------------------------+
   | +50         | OpenMP processor binding is disabled. Consider setting      |
   |             | the environment variable OMP_PROC_BIND=true (this may       |
   |             | affect performance on NUMA systems).                        |
   +-------------+-------------------------------------------------------------+

.. _ssids_example:

=======
Example
=======

Suppose we wish to factorize the matrix

.. math::

   A = \left(\begin{array}{ccccc}
      2. & 1.                \\
      1. & 4. & 1. &    & 1. \\
         & 1. & 3. & 2.      \\
         &    & 2. & 0. &    \\
         & 1. &    &    & 2.
   \end{array}\right)

and then solve for the right-hand side

.. math::

   B = \left(\begin{array}{c}
      4.    \\
      17.   \\
      19.   \\
      6.    \\
      12.
   \end{array}\right).

The following code may be used.

.. literalinclude:: ../../examples/Fortran/ssids.f90
   :language: Fortran


This produces the following output:

::

     Warning from ssids_analyse. Warning flag =   4
     one or more diagonal entries is missing

     The computed solution is:
      1.0000000000E+00  2.0000000000E+00  3.0000000000E+00
      4.0000000000E+00  5.0000000000E+00
     Pivot order:   4    5   -2   -1    3

.. _ssids_method:

======
Method
======

Partition of work across available resources
--------------------------------------------

Once the ordering has been determined and the assembly tree determined in the
analyse phase, the tree is broken into a number of leaf subtrees rooted at a
single node, leaving a root subtree/forest consisting of all remaining nodes
above those. Each leaf subtree is pre-assigned to a particular NUMA region or
GPU for the factorization phase. Exact details of the algorithm used for
finding these subtrees and their assignment can be found in the paper [2]_.

The factorization phase has two steps. In the first, leaf subtrees are
factorized in parallel on their assigned resources. Once all leaf subtrees are
factored, the second phase begins where all CPU resources cooperate to factorize
the root subtree.

At present the solve phase is performed in serial.

Data checking
-------------

If check is set to .true. on the call to :f:subr:`ssids_analyse()` or if
:f:subr:`ssids_analyse_coord()` is called, the user-supplied matrix data is
checked for errors. The cleaned integer matrix data (duplicates are
summed and out-of-range indices discarded) is stored in akeep. The use
of checking is optional on a call to :f:subr:`ssids_analyse()` as it incurs both
time and memory overheads. However, it is recommended since the
behaviour of the other routines in the package is unpredictable if
duplicates and/or out-of-range variable indices are entered.

If the user has supplied an elimination order it is checked for errors.
Otherwise, an elimination order is generated by the package. The
elimination order is used to construct an assembly tree. On exit from
:f:subr:`ssids_analyse()` (and :f:subr:`ssids_analyse_coord()`), `order(:)` is
set so that `order(i)` holds the position of variable :math:`i` in the
elimination order. If an ordering was supplied by the user, this order may
differ, but will be equivalent in terms of fill-in.

Factorization performed
-----------------------

The factorization performed depends on the value of `posdef`.

If **posdef=true**, a Cholesky factorization is performed:

.. math:: SAS = PL(PL)^T.

Pivoting is not performed, so :math:`P` is the permutation determined in the
analysis phase. :math:`S` is a diagonal scaling matrix.

If **posdef=false**, a symmetric indefinite factorization is performed:

.. math:: SAS = PLD(PL)^T

Pivoting is performed, so :math:`P` may differ from the permutation determined
in the analysis phase, though it is kept as close as possible to minimize fill.
The exact pivoting algorithm varies depending on the particular kernel the
algorithm chooses to employ.

Full details of the algorithms used are provided in [1]_ for GPUs and [2]_ for
CPUs.

.. _ssids_small_leaf:

Small Leaf Subtrees
-------------------

For subtrees allocated to run on the CPU, the factorization of small nodes near
the leaves of the tree can be amalgamated into a single parallel task (normally
each would be treated as its own OpenMP task to be scheduled). This can reduce
scheduling overheads, especially on small problems. If the total number of
operations for a subtree root at a given node is less than
`options.small_subtree_threshold`, that subtree is treated as a single task.

References
----------

.. [1] J.D. Hogg, E. Ovtchinnikov and J.A. Scott. (2014).
   *A sparse symmetric indefinite direct solver for GPU architectures*.
   ACM Transactions on Mathematical Software 42(1), Article 1, 25 pages.
   [`DOI: 10.1145/2756548 <https://doi.org/10.1145/2756548>`_]
   [`Preprint RAL-P-2014-006 <https://epubs.stfc.ac.uk/work/12189719>`_]

.. [2] J.D. Hogg. (2016).
   *A new sparse LDL^T solver using a posteriori threshold pivoting*.
   RAL Technical Report. RAL-TR-2016-0xx, to appear.
