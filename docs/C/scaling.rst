********************************
SCALING - Sparse matrix scalings
********************************

.. code-block:: C

   #include <spral_scaling.h> /* or <spral.h> for all packages */

=======
Purpose
=======

This package generates various scalings (and matchings) of real sparse matrices.

Given a **symmetric** matrix :math:`A`, it finds a diagonal matrix :math:`D`
such that the scaled matrix

.. math::

   \hat{A} = DAD

has specific numerical properties.

Given an **unsymmetric** or **rectangular** matrix :math:`A`, it finds
diagonal matrices :math:`D_r` and :math:`D_c` such that the scaled matrix

.. math::
   \hat{A} = D_r A D_c

has specific numerical properties.

The specific numerical properties delivered depends on the algorithm used:

Matching-based
   algorithms scale :math:`A` such that the maximum (absolute) value in each row
   and column of :math:`\hat{A}` is exactly :math:`1.0`, where the entries of
   maximum value form a maximum cardinality matching. The
   :ref:`Hungarian algorithm <hungarian_algorithm>` delivers an optimal matching
   slowly, whereas the :ref:`auction algorithm <auction_algorithm>` delivers an
   approximate matching quickly.
Norm-equilibration
   algorithms scale :math:`A` such that the infinity norm of each row and
   column of :math:`\hat{A}` is :math:`1.0\pm \tau` (for some user specified
   tolerance :math:`\tau`).

.. _auction_algorithm:

=================
Auction Algorithm
=================

Routines
""""""""

.. c:function:: void spral_scaling_auction_default_options(struct spral_scaling_auction_options *options)

   Intialises members of options structure to default values.

   :param options: Structure to be initialised.

.. c:function:: void spral_scaling_auction_sym(int n, const int *ptr, const int *row, const double *val, double *scaling, int *match, const struct spral_scaling_auction_options *options, struct spral_scaling_auction_inform *inform)

   Find a matching-based symmetric scaling using the auction algorithm.

   The scaled matrix is such that the entry of maximum absolute value in each
   row and column is (approximately) :math:`1.0`.

   :param n: number of columns in :math:`A`.
   :param ptr[n+1]: columns pointers for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param row[ptr[n]]: row indices for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param val[ptr[n]]: non-zero values for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param scaling[n]: returns scaling found by routine.
   :param match: may be `NULL`; otherwise, an array of size `n` to output the
      matching found by routine. Row `i` is matched to column `match[i]`, or is
      unmatched if `match[i]=0`.
   :param options: controls behaviour of routine.
   :param inform: returns information on execution of routine.

.. c:function:: void spral_scaling_auction_sym_long(int n, const int64_t *ptr, const int *row, const double *val, double *scaling, int *match, const struct spral_scaling_auction_options *options, struct spral_scaling_auction_inform *inform)

   As :c:func:`spral_scaling_auction_sym`, except `ptr` has type ``int64_t``.

.. c:function:: void spral_scaling_auction_unsym(int m, int n, const int *ptr, const int *row, const double *val, double *rscaling, double *cscaling, int *match, const struct spral_scaling_auction_options *options, struct spral_scaling_auction_inform *inform)

   Find a matching-based unsymmetric scaling using the auction algorithm.

   The scaled matrix is such that the entry of maximum absolute value in each
   row and column is (approximately) :math:`1.0`.

   :param m: number of rows in :math:`A`
   :param n: number of columns in :math:`A`
   :param ptr[n+1]: columns pointers for :math:`A` (see :doc:`CSC format<csc_format>`)
   :param row[ptr[n]]: row indices for :math:`A` (see :doc:`CSC format<csc_format>`)
   :param val[ptr[n]]: non-zero values for :math:`A` (see :doc:`CSC format<csc_format>`)
   :param rscaling[m]: returns row scaling found by routine
   :param cscaling[n]: returns column scaling found by routine
   :param match: may be `NULL`; otherwise, an array of size `m` to output the
      matching found by routine. Row `i` is matched to column `match[i]`, or is
      unmatched if `match[i]=0`.
   :param options: controls behaviour of routine
   :param inform: returns information on execution of routine

.. c:function:: void spral_scaling_auction_unsym_long(int m, int n, const int64_t *ptr, const int *row, const double *val, double *rscaling, double *cscaling, int *match, const struct spral_scaling_auction_options *options, struct spral_scaling_auction_inform *inform)

   As :c:func:`spral_scaling_auction_unsym`, except `ptr` has type ``int64_t``.

Data-types
""""""""""

.. c:type:: struct spral_scaling_auction_options

   Used to specify options to the routines :c:func:`spral_scaling_auction_sym`
   and :c:func:`spral_scaling_auction_unsym`. The routine
   :c:func:`spral_scaling_auction_default_options` may be used to intialise
   with default values.

   Please refer to the :ref:`method section<auction_algorithm_method>` for
   details on how these parameters are used.

   .. c:member:: int array_base

      Indexing base for arrays. Either 0 (C indexing) or 1 (Fortran indexing).
      Default is 0.

   .. c:member:: float eps_initial

      Initial value of improvement parameter :math:`\epsilon`.
      Default is 0.01.

   .. c:member:: int max_iterations

      Maximum number of iterations.
      Default is 30000.

   .. c:member:: int max_unchanged[3]

      Together with `min_proportion[]`, specifies termination conditions.
      Default is `{10, 100, 100}`.

   .. c:member:: float min_proportion[3]

      Together with `max_unchanged[]`, specifies termination conditions.
      Default is `{0.9, 0.0, 0.0}`.

.. c:type:: struct spral_scaling_auction_inform

   Used to return information about the execution of the algorithm.

   .. c:member:: int flag

      Gives the exit status of the algorithm (see table below)

   .. c:member:: int iterations

      Number of iterations performed.

   .. c:member:: int matched

      Number of rows and columns that have been matched.

   .. c:member:: int stat

      Fortran stat parameter in the event of an allocation failure (set to 0
      otherwise).

   .. c:member:: int unmatchable

      Number of columns designated as unmatchable
      (there is no way to match it that improves the quality of the matching).

   Note: As the algorithm may terminate before a full matching is obtained,
   :c:member:`spral_scaling_auction_inform.matched` provides only a lower bound
   on the structural rank. However,
   :c:member:`spral_scaling_auction_inform.unmatchable` provides an approximate
   lower bound on the structural rank deficiency.

   +-------------+------------------------------------------------------------+
   | inform.flag | Return status                                              |
   +=============+============================================================+
   | 0           | Success.                                                   |
   +-------------+------------------------------------------------------------+
   | -1          | Allocation error. Fortran stat value is returned in        |
   |             | :c:member:`spral_scaling_auction_inform.stat`.             |
   +-------------+------------------------------------------------------------+

Example
"""""""

The following code shows an example usage of :f:subr:`auction_scale_sym`.

.. literalinclude:: ../../examples/C/scaling/auction_sym.c
   :language: C

The above code produces the following output::

   Initial matrix:
   Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
   0:   2.0000E+00   1.0000E+00
   1:   1.0000E+00   4.0000E+00   1.0000E+00                8.0000E+00
   2:                1.0000E+00   3.0000E+00   2.0000E+00
   3:                             2.0000E+00
   4:                8.0000E+00                             2.0000E+00
   Matching:          0         4         3         2         1
   Scaling:    7.07E-01  1.62E-01  2.78E-01  1.80E+00  7.72E-01
   Scaled matrix:
   Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
   0:   1.0000E+00   1.1443E-01
   1:   1.1443E-01   1.0476E-01   4.5008E-02                1.0000E+00
   2:                4.5008E-02   2.3204E-01   1.0000E+00
   3:                             1.0000E+00
   4:                1.0000E+00                             1.1932E+00

.. _auction_algorithm_method:

Method
""""""

This algorithm finds a fast approximation to the matching and scaling produced
by the HSL package MC64. If an optimal matching is required, use the
Hungarian algorithm instead. The algorithm works by solving the following
maximum product optimization problem using an auction algorithm. The scaling
is derived from the dual variables associated with the solution.

.. math::
   \max_{\sigma} & \prod_{i=1}^m\prod_{j=1}^n |a_{ij}|\sigma_{ij} & \\
   \mathrm{s.t.} & \sum_{i=1}^m\sigma_{ij} = 1, & \forall j=1,n \\
                 & \sum_{j=1}^n\sigma_{ij} = 1, & \forall i=1,m \\
                 & \sigma_{ij} \in \{0,1\}.

The array :math:`\sigma` gives a matching of rows to columns.

By using the transformation

.. math::

   w_{ij} = \log c_j - \log |a_{ij}|,

where :math:`c_j = \max_i |a_{ij}|`, the maximum product problem in
:math:`a_{ij}` is replaced by a minimum sum problem in :math:`w_{ij}` where all
entries are positive. By standard optimization theory, there exist dual
variables :math:`u` and :math:`v` corresponding to the constraints that satisfy
the first order optimality conditions

.. math::

   w_{ij} - u_i - v_j = 0, & \mbox{ if } \sigma_{ij }=1, \\
   w_{ij} - u_i - v_j \ge 0, & \mbox{ if } \sigma_{ij }=0.

To obtain a scaling we define scaling matrices :math:`D_r` and :math:`D_c` as

.. math::

   d^r_i = e^{u_i},

   d^c_i = e^{v_i}.

If a symmetric scaling is required, we average these as

.. math::

   d_i = \sqrt{d^r_id^c_i}.

By the first order optimality conditions, these scaling matrices guarantee that

.. math::
   d^r_i|a_{ij}|d^c_j = 1, && \mbox{if } \sigma_{ij}=1, \\
   d^r_i|a_{ij}|d^c_j \le 1, && \mbox{if } \sigma_{ij}=0.

To solve the minimum sum problem an auction algorithm is used. The
algorithm is **not** guaranteed to find an optimal matching. However it
can find an approximate matching very quickly. A matching is maintained along
with the row pricing vector :math:`u`. In each major iteration, we loop over
each column in turn. If the column :math:`j` is unmatched, we calculate the
value :math:`p_i = w_{ij} - u_i` for each entry and find the maximum across the
column. If this maximum is positive, the current matching can be improved by
matching column :math:`j` with row :math:`i`. This may mean that the previous
match of row :math:`i` now becomes unmatched. We update the price of row
:math:`i`, that is :math:`u_i`, to reflect this new benefit and continue to the
next column.

To prevent incremental shuffling, we insist that the value of adding a new
column is at least a threshold value :math:`\epsilon` above zero, where
:math:`\epsilon` is based on the last iteration in which row :math:`i` changed
its match. This is done by adding :math:`\epsilon` to the price :math:`u_i`,
where :math:`\epsilon = \mathrm{options.eps_initial} + \mathrm{itr} / (n+1)`,
where itr is the current iteration number.

The algorithm terminates if any of the following are satsified:

* All entries are matched.
* The number of major iterations exceeds
  :c:member:`options.max_iterations <spral_scaling_auction_options.max_iterations>`.
* At least :c:member:`options.max_unchanged[0] <spral_scaling_auction_options.max_unchanged>`
  iterations have passed without the cardinality of the matching increasing,
  and the proportion of matched columns is
  :c:member:`options.min_proportion[0] <spral_scaling_auction_options.min_proportion>`.
* At least :c:member:`options.max_unchanged[1] <spral_scaling_auction_options.max_unchanged>`
  iterations have passed without the cardinality of the matching increasing,
  and the proportion of matched columns is
  :c:member:`options.min_proportion[1] <spral_scaling_auction_options.min_proportion>`.
* At least :c:member:`options.max_unchanged[2] <spral_scaling_auction_options.max_unchanged>`
  iterations have passed without the cardinality of the matching increasing,
  and the proportion of matched columns is
  :c:member:`options.min_proportion[2] <spral_scaling_auction_options.min_proportion>`.

The different combinations given by options.max_unchanged[] and
options.min_proportion[] allow a wide range of termination
heuristics to be specified by the user depending on their particular needs.
Note that the matching and scaling produced will always be approximate as
:math:`\epsilon` is non-zero.

Further details are given in the following paper:

.. [1] J.D. Hogg and J.A. Scott. (2014). On the efficient scaling of sparse symmetric matrices using an auction algorithm. RAL Technical Report RAL-P-2014-002. [`STFC TR <https://epubs.stfc.ac.uk/work/11539865>`_]

============================
Norm-equilibration Algorithm
============================

Routines
""""""""

.. c:function:: void spral_scaling_equilib_default_options(struct spral_scaling_equilib_options *options)

   Intialises members of options structure to default values.

   :param options: Structure to be initialised.

.. c:function:: void spral_scaling_equilib_sym(int n, const int *ptr, const int *row, const double *val, double *scaling, const struct spral_scaling_equilib_options *options, struct spral_scaling_equilib_inform *inform)

   Find a matching-based symmetric scaling using the norm-equilibration
   algorithm.

   The scaled matrix is such that the infinity norm of each row and column are
   equal to :math:`1.0`.

   :param n: number of columns in :math:`A`.
   :param ptr[n]: columns pointers for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param row[ptr[n]]: row indices for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param val[ptr[n]]: non-zero values for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param scaling[n]: returns scaling found by routine.
   :param options: controls behaviour of routine.
   :param inform: returns information on execution of routine.

.. c:function:: void spral_scaling_equilib_sym_long(int n, const int64_t *ptr, const int *row, const double *val, double *scaling, const struct spral_scaling_equilib_options *options, struct spral_scaling_equilib_inform *inform)

   As :c:func:`spral_scaling_equilib_sym`, except `ptr` has type ``int64_t``.

.. c:function:: void spral_scaling_equilib_unsym(int m, int n, const int *ptr, const int *row, const double *val, double *rscaling, double *cscaling, const struct spral_scaling_equilib_options *options, struct spral_scaling_equilib_inform *inform)

   Find a matching-based unsymmetric scaling using the norm-equilibration
   algorithm.

   The scaled matrix is such that the infinity norm of each row and column are
   equal to :math:`1.0`.

   :param m: number of rows in :math:`A`.
   :param n: number of columns in :math:`A`.
   :param ptr[n+1]: columns pointers for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param row[ptr[n]]: row indices for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param val[ptr[n]]: non-zero values for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param rscaling[m]: returns row scaling found by routine.
   :param cscaling[n]: returns column scaling found by routine.
   :param options: controls behaviour of routine.
   :param inform: returns information on execution of routine.

.. c:function:: void spral_scaling_equilib_unsym_long(int m, int n, const int64_t *ptr, const int *row, const double *val, double *rscaling, double *cscaling, const struct spral_scaling_equilib_options *options, struct spral_scaling_equilib_inform *inform)

   As :c:func:`spral_scaling_equilib_unsym`, except `ptr` has type ``int64_t``.

Data-types
""""""""""

.. c:type:: struct spral_scaling_equilib_options

   Used to specify options to the routines :c:func:`spral_scaling_equilib_sym()`
   and :c:func:`spral_scaling_equilib_unsym()`. The routine
   :c:func:`spral_scaling_equilib_default_options()` may be used to intialise
   with default values.

   Please refer to the :ref:`method section<equilib_algorithm_method>` for
   details on how these parameters are used.

   .. c:member:: int array_base

      Indexing base for arrays. Either 0 (C indexing) or 1 (Fortran indexing).
      Default is 0.

   .. c:member:: int max_iterations

      Maximum number of iterations.
      Default is 10.

   .. c:member:: float tol

      Convergence tolerance :math:`\tau`.
      Default is 1e-8.

.. c:type:: struct spral_scaling_equilib_inform

   Used to return information about the execution of the algorithm.

   .. c:member:: int flag

      Gives the exit status of the algorithm (see table below).

   .. c:member:: int iterations

      Number of iteration performed.

   .. c:member:: int stat

      Holds the Fortran stat parameter in the event of an
      allocation failure (set to 0 otherwise).

   +-------------+-------------------------------------------------------------+
   | inform.flag | Return status                                               |
   +=============+=============================================================+
   | 0           | Success.                                                    |
   +-------------+-------------------------------------------------------------+
   | -1          | Allocation error. Fortran stat value is returned in         |
   |             | :c:member:`inform.stat <spral_scaling_equilib_inform.stat>`.|
   +-------------+-------------------------------------------------------------+

Example
"""""""

The following code shows an example usage of
:c:func:`spral_scaling_equilib_sym()`.

.. literalinclude:: ../../examples/C/scaling/equilib_sym.c
   :language: C

The above code produces the following output::

   Initial matrix:
   Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
   0:   2.0000E+00   1.0000E+00
   1:   1.0000E+00   4.0000E+00   1.0000E+00                8.0000E+00
   2:                1.0000E+00   3.0000E+00   2.0000E+00
   3:                             2.0000E+00
   4:                8.0000E+00                             2.0000E+00
   Scaling:    7.07e-01   3.54e-01   5.77e-01   8.66e-01   3.54e-01
   Scaled matrix:
   Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
   0:   1.0000E+00   2.5000E-01
   1:   2.5000E-01   5.0000E-01   2.0412E-01                1.0000E+00
   2:                2.0412E-01   1.0000E+00   9.9960E-01
   3:                             9.9960E-01
   4:                1.0000E+00                             2.5000E-01

.. _equilib_algorithm_method:

Method
""""""

This algorithm is very similar to that used by the HSL routine MC77.
An iterative method is used to scale the infinity norm of both rows and columns
to :math:`1.0` with an asymptotic linear rate of convergence of
:math:`\frac{1}{2}`, preserving symmetry if the matrix is symmetric.

For unsymmetric matrices, the algorithm outline is as follows:

* :math:`D_r^{(0)} = I, D_c^{(0)}=I`
* **for** (:math:`k=1,` options.max_iterations)

   * :math:`A^{(k-1)} = D_r^{(k-1)} A D_c^{(k-1)}`
   * :math:`(D_r^{(k)})_{ii} = (D_r^{(k-1)})_{ii}\; /\; \sqrt{\max_j(A^{(k-1)})_{ij}}`
   * :math:`(D_c^{(k)})_{jj} = (D_c^{(k-1)})_{jj}\; /\; \sqrt{\max_i(A^{(k-1)})_{ij}}`
   * **if** ( :math:`|1-\|A^{(k-1)}\|_{\max}|\le` options.tol ) **exit**

For symmetric matrices, :math:`A^{(k-1)}` is symmetric, so :math:`D_r^{(k)} = D_c^{(k)}`, and some operations can be skipped.

Further details are given in the following paper:

.. [2] P. Knight, D. Ruiz and B. Ucar. (2012). A symmetry preserving algorithm for matrix scaling. INRIA Research Report 7552. [`INRIA TR <https://hal.inria.fr/inria-00569250v3/document>`_]

.. _hungarian_algorithm:

===================
Hungarian Algorithm
===================

Routines
""""""""

.. c:function:: void spral_scaling_hungarian_default_options(struct spral_scaling_hungarian_options *options)

   Intialises members of options structure to default values.

   :param options: Structure to be initialised.

.. c:function:: void spral_scaling_hungarian_sym(int n, const int *ptr, const int *row, const double *val, double *scaling, int *match, const struct spral_scaling_hungarian_options *options, struct spral_scaling_hungarian_inform *inform)

   Find a matching-based symmetric scaling using the Hungarian algorithm.

   The scaled matrix is such that the entry of maximum absolute value in each
   row and column is :math:`1.0`.

   :param n: number of columns in :math:`A`.
   :param ptr[n+1]: columns pointers for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param row[ptr[n]]: row indices for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param val[ptr[n]]: non-zero values for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param scaling[n]: returns scaling found by routine.
   :param match: may be `NULL`; otherwise, an array of size `n` to output the
      matching found by routine. Row `i` is matched to column `match[i]`, or is
      unmatched if `match[i]=0`.
   :param options: controls behaviour of routine.
   :param inform: returns information on execution of routine.

.. c:function:: void spral_scaling_hungarian_sym_long(int n, const int64_t *ptr, const int *row, const double *val, double *scaling, int *match, const struct spral_scaling_hungarian_options *options, struct spral_scaling_hungarian_inform *inform)

   As :c:func:`spral_scaling_hungarian_sym`, except `ptr` has type ``int64_t``.

.. c:function:: void spral_scaling_hungarian_unsym(int m, int n, const int *ptr, const int *row, const double *val, double *rscaling, double *cscaling, int *match, const struct spral_scaling_hungarian_options *options, struct spral_scaling_hungarian_inform *inform)

   Find a matching-based symmetric scaling using the Hungarian algorithm.

   The scaled matrix is such that the entry of maximum absolute value in each
   row and column is :math:`1.0`.

   :param m: number of rows in :math:`A`.
   :param n: number of columns in :math:`A`.
   :param ptr[n+1]: columns pointers for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param row[ptr[n]]: row indices for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param val[ptr[n]]: non-zero values for :math:`A` (see :doc:`CSC format<csc_format>`).
   :param rscaling[m]: returns row scaling found by routine.
   :param cscaling[n]: returns column scaling found by routine.
   :param match: may be `NULL`; otherwise, an array of size `n` to output the
      matching found by routine. Row `i` is matched to column `match[i]`, or is
      unmatched if `match[i]=0`.
   :param options: controls behaviour of routine.
   :param inform: returns information on execution of routine.

.. c:function:: void spral_scaling_hungarian_unsym_long(int m, int n, const int64_t *ptr, const int *row, const double *val, double *rscaling, double *cscaling, int *match, const struct spral_scaling_hungarian_options *options, struct spral_scaling_hungarian_inform *inform)

   As :c:func:`spral_scaling_hungarian_unsym()`, except `ptr` has type ``int64_t``.

Data-types
""""""""""

.. c:type:: struct spral_scaling_hungarian_options

   Used to specify options to the routines :c:func:`spral_scaling_hungarian_sym`
   and :c:func:`spral_scaling_hungarian_unsym()`. The routine
   :c:func:`spral_scaling_hungarian_default_options()` may be used to intialise
   with default values.

   Please refer to the :ref:`method section<hungarian_algorithm_method>` for
   details on how these parameters are used.

   .. c:member:: int array_base

      Indexing base for arrays. Either 0 (C indexing) or 1 (Fortran indexing).
      Default is 0.

   .. c:member:: bool scale_if_singular

      Behaviour for structurally singular matrices. If `true`, a partial scaling
      corresponding to a maximum cardinality matching will be returned.
      If `false`, an identity scaling is returned with an error code.
      Default is `false`.

   .. warning::

      If `options.scale_if_singular=true`, the resulting scaling will
      only be maximal for the matched rows/columns, and extreme care should be
      taken to ensure its use is meaningful!

.. f:type:: hungarian_inform

   Used to return information about the execution of the algorithm.

   .. c:member:: int flag

      Exit status of the algorithm (see table below)

   .. c:member:: int matched

      Number of rows and columns that have been matched.

   .. c:member:: int stat

      Holds the Fortran stat parameter in the event of an allocation failure
      (set to 0 otherwise).

   **Note:** The number matched gives the structural rank of the matrix.

   +-------------+------------------------------------------------------------+
   | inform.flag | Return status                                              |
   +=============+============================================================+
   | 0           | Success.                                                   |
   +-------------+------------------------------------------------------------+
   | +1          | Warning: Matrix is structurally rank-deficient.            |
   |             | Only returned if `options.scale_if_singular=true`.         |
   +-------------+------------------------------------------------------------+
   | -1          | Error: Allocation failed.                                  |
   |             | Fortran stat value is returned in `inform.stat`.           |
   +-------------+------------------------------------------------------------+
   | -2          | Error: Matrix is structurally rank-deficient.              |
   |             | Only returned if `options.scale_if_singular=false`. Scaling|
   |             | vector(s) will be set to the identity, and a maximum       |
   |             | cardinality matching will be returned in `match[]` (if     |
   |             | present).                                                  |
   +-------------+------------------------------------------------------------+

Example
"""""""

The following code shows an example usage of
:c:func:`spral_scaling_hungarian_sym()`.

.. literalinclude:: ../../examples/C/scaling/hungarian_sym.c
   :language: C

The above code produces the following output::

   Initial matrix:
   Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
   0:   2.0000E+00   1.0000E+00
   1:   1.0000E+00   4.0000E+00   1.0000E+00                8.0000E+00
   2:                1.0000E+00   3.0000E+00   2.0000E+00
   3:                             2.0000E+00
   4:                8.0000E+00                             2.0000E+00
   Matching:          0          4          3          2          1
   Scaling:    7.07e-01   3.54e-01   5.77e-01   8.66e-01   3.54e-01
   Scaled matrix:
   Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
   0:   1.0000E+00   2.5000E-01
   1:   2.5000E-01   5.0000E-01   2.0412E-01                1.0000E+00
   2:                2.0412E-01   1.0000E+00   1.0000E+00
   3:                             1.0000E+00
   4:                1.0000E+00                             2.5000E-01


.. _hungarian_algorithm_method:

Method
""""""

This algorithm is the same as used by the HSL package MC64. A scaling
is derived from dual variables found during the solution of the below
maximum product optimization problem using the Hungarian algorithm.

.. math::
   \max_{\sigma} & \prod_{i=1}^m\prod_{j=1}^n |a_{ij}|\sigma_{ij} & \\
   \mathrm{s.t.} & \sum_{i=1}^m\sigma_{ij} = 1, & \forall j=1,n \\
                 & \sum_{j=1}^n\sigma_{ij} = 1, & \forall i=1,m \\
                 & \sigma_{ij} \in \{0,1\}.

The array :math:`\sigma` gives a matching of rows to columns.

By using the transformation

.. math::
   w_{ij} = \log c_j - \log |a_{ij}|,

where :math:`c_j = \max_i |a_{ij}|`, the maximum product problem in
:math:`a_{ij}` is replaced by a minimum sum problem in :math:`w_{ij}` where all
entries are positive. By standard optimization theory, there exist dual
variables :math:`u` and :math:`v` corresponding to the constraints that satisfy
the first order optimality conditions

.. math::
   w_{ij} - u_i - v_j = 0, & \mbox{if } \sigma_{ij }=1, \\
   w_{ij} - u_i - v_j \ge 0, & \mbox{if } \sigma_{ij }=0.

To obtain a scaling we define scaling matrices :math:`D_r` and :math:`D_c` as

.. math::
   d^r_i = e^{u_i},

   d^c_i = e^{v_i}.

If a symmetric scaling is required, we average these as

.. math::
   d_i = \sqrt{d^r_id^c_i}.

By the first order optimality conditions, these scaling matrices guarantee that

.. math::
   d^r_i|a_{ij}|d^c_j = 1, && \mbox{if } \sigma_{ij}=1, \\
   d^r_i|a_{ij}|d^c_j \le 1, && \mbox{if } \sigma_{ij}=0.

To solve the minimum sum problem, the Hungarian algorithm maintains an optimal
matching on a subset of the rows and columns. It proceeds to grow this set by
finding augmenting paths from an unmatched row to an unmatched column. The
algorithm is guaranteed to find the optimal solution in a fixed number of steps,
but can be very slow as it may need to explore the full matrix a number of
times equal to the dimension of the matrix. To minimize the solution time, a
warmstarting heuristic is used to construct an initial optimal subset matching.

Further details are given in the following paper:

.. [3] I.S. Duff and J. Koster. (1997). The design and use of algorithms for permuting large entries to the diagonal of sparse matrices. SIAM J. Matrix Anal. Applics. 20(4), pp 889--901. [`Journal <http://dx.doi.org/10.1137/S0895479897317661>`_] [`Preprint <https://epubs.stfc.ac.uk/work/33194>`_]
