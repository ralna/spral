*************************************************
SSIDS - Sparse Symmetric Indefinite Direct Solver
*************************************************

.. code-block:: C

   #include <spral_ssids.h> /* or <spral.h> for all packages */

=======
Purpose
=======

This package solves one or more sets of :math:`n\times n`
sparse **symmetric** equations  :math:`AX = B` using a multifrontal method on an
**NVIDIA GPU**.
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

SSIDS returns bit-compatible results.

An option exists to scale the matrix. In this case, the factorization of
the scaled matrix  :math:`\overline{A} = S A S` is computed,
where :math:`S` is a diagonal scaling matrix.

Version history
---------------

2014-03-17 Version 1.0.0
   Initial release

[For detail, see ChangeLog]

==============
Usage overview
==============

Solving :math:`AX=B` using SSIDS is a four stage process.

1. Call :c:func:`spral_ssids_analyse()` to perform a symbolic factorization,
   stored in `akeep`.
2. Call :c:func:`spral_ssids_factor()` to perform a numeric factorization,
   stored in `fkeep`. More than one numeric factorization can refer to the same
   `akeep`.
3. Call :c:func:`spral_ssids_solve1()` or :c:func:`spral_ssids_solve()` to
   perform a solve with the factors. More than one solve can be performed with
   the same `fkeep`.
4. Once all desired solutions have been performed, free memory with
   :c:func:`spral_ssids_free()`.

In addition, advanced users may use the following routines:

* :c:func:`spral_ssids_enquire_posdef()` and
  :c:func:`spral_ssids_enquire_indef()` return
  the diagonal entries of the factors and the pivot sequence.
* :c:func:`spral_ssids_alter()` allows altering the diagonal entries of the
  factors.


.. note::

   **Bit-compatibility:**
   If used with bit-compatible BLAS and compiler options, this routine will
   return bit compatible results. That is, consecutive runs with the same data
   on the same machine produces exactly the same solution.

=================
Basic Subroutines
=================

Note: For the most efficient use of the package, CSC format should be used
without checking.

.. c:function:: void spral_ssids_default_options(struct spral_ssids_options *options)

   Intialises members of options structure to default values.

   :param options: Structure to be initialised.

.. c:function:: void spral_ssids_analyse(bool check, int n, int *order, const int *ptr, const int *row, const double *val, void **akeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform)

   Perform the analyse (symbolic) phase of the factorization for a matrix
   supplied in :doc:`CSC format<csc_format>`. The resulting symbolic factors
   stored in :c:type:`akeep` should be passed unaltered in the following call to
   :c:func:`spral_ssids_factor()`.

   :param check: if true, matrix data is checked. Out-of-range entries
      are dropped and duplicate entries are summed.
   :param n: number of columns in :math:`A`.
   :param order: may be `NULL`; otherwise must be an array of size `n`
      used on entry a user-supplied ordering
      (:c:type:`options.ordering=0 <spral_ssids_options.ordering>`).
      On return, the actual ordering used.
   :param ptr[n+1]: column pointers for :math:`A`
      (see :doc:`CSC format<csc_format>`).
   :param row[ptr[n]]: row indices for :math:`A`
      (see :doc:`CSC format<csc_format>`).
   :param val: may be `NULL`; otherwise must be an array of size `ptr[n]`
      containing non-zero values for :math:`A`
      (see :doc:`CSC format<csc_format>`).
      Only used if a matching-based ordering is requested.
   :param akeep: returns symbolic factorization, to be passed unchanged to
      subsequent routines.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).

   **Note:** If a user-supplied ordering is used, it may be altered by this
   routine, with the altered version returned in `order[]`. This version will be
   equivalent to the original ordering, except that some supernodes may have
   been amalgamated, a topographic ordering may have been applied to the tree
   and the order of columns within a supernode may have been adjusted to improve
   cache locality.

.. c:function:: void spral_ssids_analyse_coord(int n, int *order, int ne, const int *row, const int *col, const double *val, void **akeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform)

   :param n: number of columns in :math:`A`.
   :param order: may be `NULL`; otherwise must be an array of size `n`
      used on entry a user-supplied ordering
      (:c:type:`options.ordering=0 <spral_ssids_options.ordering>`).
      On return, the actual ordering used.
   :param ne: number of non-zero entries in :math:`A`.
   :param row[ne]: row indices for :math:`A`
      (see :doc:`Coordinate format<coord_format>`).
   :param col[ne]: column indices for :math:`A`
      (see :doc:`Coordinate format<coord_format>`).
   :param val: may be `NULL`; otherwise must be an array of size `ne`
      containing non-zero values for :math:`A`
      (see :doc:`Coordinate format<coord_format>`).
      Only used if a matching-based ordering is requested.
   :param akeep: returns symbolic factorization, to be passed unchanged to
      subsequent routines.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).

   **Note:** If a user-supplied ordering is used, it may be altered by this
   routine, with the altered version returned in `order[]`. This version will be
   equivalent to the original ordering, except that some supernodes may have
   been amalgamated, a topographic ordering may have been applied to the tree
   and the order of columns within a supernode may have been adjusted to improve
   cache locality.

.. c:function:: void spral_ssids_factor(bool posdef, const int *ptr, const int *row, const double *val, double *scale, void *akeep, void **fkeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform)

   :param posdef: true if matrix is matrix is positive-definite
   :param ptr: may be `NULL`; otherwise a length `n+1` array of column pointers
      for :math:`A`, only required if :f:type:`akeep` was obtained by running
      :c:func:`spral_ssids_analyse()` with `check=true`, in which case it must
      be unchanged since that call.
   :param row: may be `NULL`; otherwise a length `ptr[n]` array of row indices
      for :math:`A`, only required if :f:type:`akeep` was obtained by running
      :c:func:`spral_ssids_analyse()` with `check=true`, in which case it must
      be unchanged since that call.
   :param val[]: non-zero values for :math:`A` in same format as for
      the call to :c:func:`spral_ssids_analyse()` or
      :c:func:`spral_ssids_analyse_coord()`.
   :param scale: may be `NULL`; otherwise a length `n` array for diagonal
      scaling. `scale[i-1]` contains entry :math:`S_ii` of :math:`S`. Must be
      supplied by user on entry if
      c:member:`options%scaling=0 <spral_ssids_options.scaling>` (user-supplied
      scaling). On exit, returns scaling used.
   :param akeep: symbolic factorization returned by preceding call to
      :c:func:`ssids_analyse()` or :c:func:`ssids_analyse_coord()`.
   :param fkeep: returns numeric factorization, to be passed unchanged
      to subsequent routines.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).

.. c:function:: void spral_ssids_solve1(int job, double *x1, void *akeep, void *fkeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform)
   
   Solve (for a single right-hand side) one of the following equations:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0             | :math:`Ax=b`             |
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

   :param job: specifies equation to solve, as per above table.
   :param x1[n]: right-hand side :math:`b` on entry, solution :math:`x`
      on exit.
   :param akeep: symbolic factorization returned by preceding
      call to :c:func:`spral_ssids_analyse()` or
      :c:func:`spral_ssids_analyse_coord()`.
   :param fkeep: numeric factorization returned by preceding
      call to :c:func:`spral_ssids_factor()`.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).

.. c:function:: void spral_ssids_solve(int job, int nrhs, double *x, int ldx, void *akeep, void *fkeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform)
   
   Solve (for multiple right-hand sides) one of the following equations:

   +---------------+--------------------------+
   | `job`         | Equation solved          |
   +===============+==========================+
   | 0             | :math:`AX=B`             |
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

   :param job: specifies equation to solve, as per above table.
   :param nrhs: number of right-hand sides.
   :param x[ldx*nrhs]: right-hand sides :math:`B` on entry,
      solutions :math:`X` on exit. The `i`-th entry of right-hand side `j`
      is in position `x[j*ldx+i]`.
   :param ldx: leading dimension of `x`.
   :param akeep: symbolic factorization returned by preceding
      call to :c:func:`spral_ssids_analyse()` or
      :c:func:`spral_ssids_analyse_coord()`.
   :param fkeep: numeric factorization returned by preceding
      call to :c:func:`spral_ssids_factor()`.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).

.. c:function:: int spral_ssids_free_akeep(void **akeep)

   Frees memory and resources associated with :c:type:`akeep`.

   :param akeep: symbolic factors to be freed.
   :returns: 0 on success, or a CUDA error code on failure.

.. c:function:: int spral_ssids_free_fkeep(void **fkeep)

   Frees memory and resources associated with :c:type:`fkeep`.

   :param fkeep: numeric factors to be freed.
   :returns: 0 on success, or a CUDA error code on failure.

.. c:function:: int spral_ssids_free(void **akeep, void **fkeep)

   Frees memory and resources associated with :c:type:`akeep` and
   :c:type:`fkeep`.

   :param akeep: symbolic factors to be freed.
   :param fkeep: numeric factors to be freed.
   :returns: 0 on success, or a CUDA error code on failure.

.. warning::

   The free routine(s) must be called by the user. Merely deallocating
   :c:type:`akeep` or :c:type:`fkeep`, or allowing them to go out of scope will
   result in memory leaks. :c:type:`akeep` should only be deallocated after all
   associated numeric factorizations :c:type:`fkeep` have been freed.

====================
Advanced subroutines
====================

.. c:function:: void spral_ssids_enquire_posdef(const void *akeep, const void *fkeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform, double *d)

   Return the diagonal entries of the Cholesky factor.

   :param akeep: symbolic factorization returned by preceding
      call to :c:func:`spral_ssids_analyse()` or
      :c:func:`spral_ssids_analyse_coord()`.
   :param fkeep: numeric factorization returned by preceding
      call to :c:func:`spral_ssids_factor()`.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).
   :param d[n]: returns the diagonal of :math:`L`. d[i-1] stores the
      entry :math:`L_{ii}`.

.. c:function:: void spral_ssids_enquire_indef(const void *akeep, const void *fkeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform, int *piv_order, double *d)

   Return the pivot order and/or values of :math:`D` of the Symmetric Indefinite
   Factorization.

   :param akeep: symbolic factorization returned by preceding
      call to :c:func:`spral_ssids_analyse()` or
      :c:func:`spral_ssids_analyse_coord()`.
   :param fkeep: numeric factorization returned by preceding
      call to :c:func:`spral_ssids_factor()`.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).
   :param piv_order: may be `NULL`; otherwise a length `n` array for the pivot
      order. On return, :math:`|\,\texttt{piv_order[i]}|` gives the position
      of variable :math:`i` in the pivot order.
   :param d: may be `NULL`; otherwise a length `2*n` array for the
      :math:`2\times2` block diagonal of :math:`D`. `d[2*(i-1)+0]` stores
      :math:`D_{ii}` and `d[2*(i-1)+1]` stores :math:`D_{(i+1)i}`.

.. c:function:: void spral_ssids_alter(const double *d, const void *akeep, void *fkeep, const struct spral_ssids_options *options, struct spral_ssids_inform *inform)

   Alter the entries of the diagonal factor :math:`D` for a symmetric indefinite
   factorization. The pivot order remains the same.

   :param d[2*n]: New entries of :math:`D`.
      `d[2*(i-1)+0]` stores :math:`D_{ii}` and `d[2*(i-1)+1]` stores
      :math:`D_{(i+1)i}`.
   :param akeep: symbolic factorization returned by preceding
      call to :c:func:`spral_ssids_analyse()` or
      :c:func:`spral_ssids_analyse_coord()`.
   :param fkeep: numeric factorization returned by preceding
      call to :c:func:`spral_ssids_factor()`.
   :param options: specifies algorithm options to be used
      (see :c:type:`spral_ssids_options`).
   :param inform: returns information about the execution of the routine
      (see :c:type:`spral_ssids_inform`).

   **Note:** This routine is not compatabile with the option
   ``options.presolve=1``.

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

   :f integer nemin [default=8]: supernode amalgamation threshold. Two
      neighbours in the elimination tree are merged if they both involve fewer
      than nemin eliminations. The default is used if nemin<1.
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

   :f logical action [default=.true.]: continue factorization of singular matrix
      on discovery of zero pivot if true (a warning is issued), or abort if
      false.
   :f real u [default=0.01]: relative pivot threshold used in symmetric
      indefinite case. Values outside of the range :math:`[0,0.5]` are treated
      as the closest value in that range.
   :f logical use_gpu_solve [default=.true.]: use GPU for solve phase if true,
      or CPU if false. 
   :f integer presolve [default=0]: perform presolve operations during factorize
      to accelerate solves. It may take the following values:

      +-------------+----------------------------------------------------------+
      | 0 (default) | Minimal work is performed during :f:subr:`ssids_factor()`|
      |             | to prepare for the solve.                                |
      +-------------+----------------------------------------------------------+
      | 1           | The explicit inverse of the                              |
      |             | :math:`\mathrm{nelim}\times\mathrm{nelim}` block in each |
      |             | supernode is precalculated during                        |
      |             | :f:subr:`ssids_factor()` (where ``nelim`` is the number  |
      |             | of variables eliminated at that supernode). As the matrix|
      |             | :math:`L` is overwritten, the routine                    |
      |             | :f:subr:`ssids_alter() cannot be used.                   |
      |             |                                                          |
      |             | This option requires %use_gpu_solve=.true.               |
      +-------------+----------------------------------------------------------+

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
   :f integer cublas_error: CUDA error code in the event of a CUDA error
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
   | -12         | The combination of options%use_gpu_solve and                |
   |             | options%presolve are not compatible with the requested      |
   |             | operation.                                                  |
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

======
Method
======

:f:subr:`ssids_analyse()` and :f:subr:`ssids_analyse_coord()`
-------------------------------------------------------------

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
:f:subr:`ssids_analyse()` (and :f:subr:`ssids_analyse_coord()`), order(:) is
set so that order(i) holds the position of variable :math:`i` in the elimination
order. If an ordering was supplied by the user, this order may differ,
but will be equivalent in terms of fill-in.

:f:subr:`ssids_factor()`
------------------------

:f:subr:`ssids_factor()` optionally computes a scaling and then performs the
numerical factorization. The user must specify whether or not the matrix
is positive definite. If posdef is set to .true., no pivoting is
performed and the computation will terminate with an error if a
non-positive pivot is encountered.

The factorization uses the assembly tree that was set up by the analyse
phase. At each node, entries from :math:`A` and, if it is not a leaf
node, the generated elements and any delayed pivots from its child nodes
must be assembled. Separate kernels handle each of these.

The kernel that performs the assembly from the child nodes considers one
parent-child assembly at a time. Each generated element from a child is
divided into a number of tiles, and a thread block launched to assemble
each tile into a dense submatrix using a simple mapping array to
determine the destination row and column of each entry.
Bit-compatibility is achieved by ensuring the child entries are always
assembled in the same order.

A dense partial factorization of the fully summed columns is then
performed. The fully summed columns are split into a number of tiles
that are each handled by an associated block. Factorization proceeds one
column of tiles at a time. The pivoting condition is chosen to ensure
that all entries of :math:`L` have absolute value less than
:math:`u^{-1}`. This limits the growth of the entries of the
:math:`D` factor and ensures that any solves will be backwards stable.
The details are described in [1]_.

If a pivot candidate does not pass the pivot tests at a node, it is
delayed to its parent node, where further elimination operations may
make it acceptable. Delaying pivots leads to additional fill-in and
floating-point operations beyond that predicted by :f:subr:`ssids_analyse()`
(or :f:subr:`ssids_analyse_coord()`), and may result in additional memory
allocations being required. The number of delayed pivots can often be
reduced by using appropriate scaling.

At each non-root node, the majority of the floating-point operations
involve the formation of the generated element. This is handled by a
single dedicated kernel; again, see [1]_ for details.

At the end of the factorization, data structures for use in future calls
to :f:subr:`ssids_solve()` are prepared. If ``options%presolve=1``, the block
of :math:`L` corresponding to the eliminated variables is explicitly
inverted to accelerate future calls to :f:subr:`ssids_solve()` at the cost of
making :f:subr:`ssids_factor()` slower.

:f:subr:`ssids_solve()`
-----------------------

If ``options%use_gpu_solve=.false.``, data is moved to
the CPU if required and the BLAS calls are used to perform a solve using
the assembly tree and factors generated on previous calls.

Otherwise, the solve is conducted on the GPU in a similar fashion. If
``options%presolve=0``, custom GPU implementations of ``_trsv()`` and
``_gemv()`` are used to handle multiple independent operations. If
multiple right-hand sides are to be solved for, the single right-hand
side solve is looped over. If ``options%presolve=1``, ``_trsv()`` can be
replaced by the much more parallel (and hence faster) ``_gemv()``. In
this case multiple right-hand sides are handled at the same time.

References
----------

.. [1] J.D. Hogg, E. Ovtchinnikov and J.A. Scott. (2014).
   A sparse symmetric indefinite direct solver for GPU architectures.
   RAL Technical Report. RAL-P-2014-0xx, to appear.
