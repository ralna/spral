************************************************
SSMFE - Sparse Symmetric Matrix-Free Eigensolver
************************************************

.. code-block:: C

   #include <spral_ssmfe.h> /* or <spral.h> for all packages */

=======
Purpose
=======

This package computes extreme (leftmost and/or rightmost)
eigenpairs :math:`\{\lambda_i, x_i\}` of the following eigenvalue problems:

- the standard eigenvalue problem

   .. math:: A x = \lambda x,

- the generalized eigenvalue problem

   .. math:: A x = \lambda B x,

- the buckling problem

   .. math:: B x = \lambda A x,

where :math:`A` and :math:`B` are **real symmetric** (or **Hermitian**) matrices
and :math:`B` is **positive definite**.

This package provides a user-friendly wrapper around
:f:mod:`spral_ssmfe_expert`, which in turn provides a wrapper around
:f:mod:`spral_ssmfe_core`. If more fine-tuned control of the eigensolver is
required, use those modules instead.

Version history
---------------

2015-04-20 Version 1.0.0
    Initial release

[for detail please see ChangeLog]

==============
Usage overview
==============

The eigensolver subroutines behind this package implement a block
iterative algorithm. The block nature of this algorithm allows the user
to benefit from highly optimized linear algebra subroutines and from the
ubiquitous multicore architecture of modern computers. It also makes
this algorithm more reliable than Krylov-based algorithms employed e.g.
by ARPACK in the presence of clustered eigenvalues. However, convergence
of the iterations may be slow if the density of the spectrum is high.

Thus, good performance (in terms of speed) is contingent on the
following two factors:

i. the number of desired eigenpairs must be substantial (e.g. not fewer
   than the number of CPU cores), and
ii. the employment of a convergence acceleration technique.

The acceleration techniques that can be used are shift-and-invert and
preconditioning.

The former requires the direct solution of linear systems
with the matrix :math:`A` or its linear combination with :math:`B`, for which a
sparse symmetric indefinite solver (such as HSL_MA97 or SPRAL_SSIDS)
can be employed.

The latter applies to the case of positive definite
:math:`A` and requires a matrix or an operator :math:`T`, called *a
preconditioner*, such that the vector :math:`v = T f` is an
approximation to the solution :math:`u` of the system :math:`A u = f`
(see the simple :ref:`example below <example>`). Note: This
technique is only recommended for experienced users.

===========
Subroutines
===========

.. c:function:: void spral_ssmfe_default_options(struct spral_ssmfe_options *options)

   Intialises members of options structure to default values.

   :param options: Structure to be initialised.

.. c:function:: void spral_ssmfe_standard_double(struct spral_ssmfe_rcid *rci, int left, int mep, double *lambda, int n, double *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   Computes the left-most eigenpairs of the standard eigenvalue problem

   .. math:: Ax = \lambda x

   Optionally uses preconditioning.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci.job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see :c:member:`inform.flag                 |
   |          | <spral_ssmfe_inform.flag>`.                                   |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see                     |
   |          | :c:member:`inform.flag <spral_ssmfe_inform.flag>`.            |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`Y = AX`.                                     |
   +----------+---------------------------------------------------------------+
   |  2       | Apply preconditioner :math:`Y = TX`.                          |
   +----------+---------------------------------------------------------------+

   The matrices :math:`X` and :math:`Y` are pointed to by components of `rci`.

   :param rci: Reverse communication type.
      :c:member:`rci.job <spral_ssmfe_rcid.job>` must be
      set to `0` before the first call.
   :param left: Number of left eigenpairs to find.
   :param mep: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left`.
   :param lambda[mep]: Current eigenvalue estimates in ascending
      order.
   :param n: Size of matrix :math:`A`.
   :param x[n][ldx]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      :c:member:`options.user_x>0 <spral_ssmfe_options.user_x>`.
   :param ldx: Leading dimension of `x`.
   :param keep: Internal workspace used by routine.
   :param options: specifies algorithm options to be used.
   :param inform: returns information about the exection of the
      routine.

.. c:function:: void spral_ssmfe_standard_double_complex(struct spral_ssmfe_rciz *rci, int left, int mep, double *lambda, int n, double complex *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   As :c:func:`spral_ssmfe_standard_double()`, but types of ``rci`` and ``x``
   changed to support type ``double complex``.

.. c:function:: void spral_ssmfe_standard_shift_double(struct spral_ssmfe_rcid *rci, double sigma, int left, int right, int mep, double *lambda, int n, double *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   Computes eigenpairs of the standard eigenvalue problem

   .. math:: Ax = \lambda x

   in the vicinity of a given value :math:`\sigma`.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci.job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see :c:member:`inform.flag                 |
   |          | <spral_ssmfe_inform.flag>`.                                   |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see                     |
   |          | :c:member:`inform.flag <spral_ssmfe_inform.flag>`.            |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`Y = AX`.                                     |
   +----------+---------------------------------------------------------------+
   |  9       | Solve :math:`(A-\sigma I)Y = X` for Y.                        |
   +----------+---------------------------------------------------------------+

   The matrices :math:`X` and :math:`Y` are components of `rci`.

   :param rci: Reverse communication type.
      :c:member:`rci.job <spral_ssmfe_rcid.job>` must be
      set to `0` before the first call.
   :param sigma: Shift value :math:`sigma`.
   :param left: Number of left eigenpairs to find.
   :param right: Number of right eigenpairs to find.
   :param mep: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left+right`.
   :param lambda[mep]: Current eigenvalue estimates in ascending
      order.
   :param n: Size of matrix :math:`A`.
   :param x[n][ldx]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      :c:member:`options.user_x>0 <spral_ssmfe_options.user_x>`.
   :param ldx: Leading dimension of `x`.
   :param keep: Internal workspace used by routine.
   :param options: specifies algorithm options to be used.
   :param inform: returns information about the exection of the
      routine.

.. c:function:: void spral_ssmfe_standard_shift_double_complex(struct spral_ssmfe_rciz *rci, double sigma, int left, int right, int mep, double *lambda, int n, double complex *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   As :c:func:`spral_ssmfe_standard_shift_double()`, but types of ``rci``
   and ``x`` changed to support type ``double complex``.

.. c:function:: void spral_ssmfe_generalized_double(struct spral_ssmfe_rcid *rci, int left, int mep, double *lambda, int n, double *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   Computes the left-most eigenpairs of the generalized eigenvalue problem

   .. math:: Ax = \lambda B x

   Optionally uses preconditioning.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci.job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see :c:member:`inform.flag                 |
   |          | <spral_ssmfe_inform.flag>`.                                   |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see                     |
   |          | :c:member:`inform.flag <spral_ssmfe_inform.flag>`.            |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`Y = AX`.                                     |
   +----------+---------------------------------------------------------------+
   |  2       | Apply preconditioner :math:`Y = TX`.                          |
   +----------+---------------------------------------------------------------+
   |  3       | Calculate :math:`Y = BX`.                                     |
   +----------+---------------------------------------------------------------+

   The matrices :math:`X` and :math:`Y` are components of `rci`.

   :param rci: Reverse communication type.
      :c:member:`rci.job <spral_ssmfe_rcid.job>` must be
      set to `0` before the first call.
   :param left: Number of left eigenpairs to find.
   :param mep: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left`.
   :param lambda[mep]: Current eigenvalue estimates in ascending
      order.
   :param n: Size of matrix :math:`A`.
   :param x[n][ldx]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      :c:member:`options.user_x>0 <spral_ssmfe_options.user_x>`.
   :param ldx: Leading dimension of `x`.
   :param keep: Internal workspace used by routine.
   :param options: specifies algorithm options to be used.
   :param inform: returns information about the exection of the
      routine.

.. c:function:: void spral_ssmfe_generalized_double_complex(struct spral_ssmfe_rciz *rci, int left, int mep, double *lambda, int n, double complex *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   As :c:func:`spral_ssmfe_generalized_double()`, but types of ``rci`` and
   ``x`` changed to support type ``double complex``.

.. c:function:: void spral_ssmfe_generalized_shift_double(struct spral_ssmfe_rcid *rci, double sigma, int left, int right, int mep, double *lambda, int n, double *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   Computes eigenpairs of the generalized eigenvalue problem

   .. math:: Ax = \lambda B x

   in the vicinity of a given value :math:`\sigma`.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci.job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see :c:member:`inform.flag                 |
   |          | <spral_ssmfe_inform.flag>`.                                   |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see                     |
   |          | :c:member:`inform.flag <spral_ssmfe_inform.flag>`.            |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`Y = AX`.                                     |
   +----------+---------------------------------------------------------------+
   |  3       | Calculate :math:`Y = BX`.                                     |
   +----------+---------------------------------------------------------------+
   |  9       | Solve :math:`(A-\sigma B)Y = X` for Y.                        |
   +----------+---------------------------------------------------------------+

   The matrices :math:`X` and :math:`Y` are components of `rci`.

   :param rci: Reverse communication type.
      :c:member:`rci.job <spral_ssmfe_rcid.job>` must be
      set to `0` before the first call.
   :param sigma: Shift value :math:`sigma`.
   :param left: Number of left eigenpairs to find.
   :param right: Number of right eigenpairs to find.
   :param mep: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left+right`.
   :param lambda[mep]: Current eigenvalue estimates in ascending
      order.
   :param n: Size of matrix :math:`A`.
   :param x[n][ldx]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      :c:member:`options.user_x>0 <spral_ssmfe_options.user_x>`.
   :param ldx: Leading dimension of `x`.
   :param keep: Internal workspace used by routine.
   :param options: specifies algorithm options to be used.
   :param inform: returns information about the exection of the
      routine.

.. c:function:: void spral_ssmfe_generalized_shift_double_complex(struct spral_ssmfe_rciz *rci, double sigma, int left, int right, int mep, double *lambda, int n, double complex *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   As :c:func:`spral_ssmfe_generalized_shift_double()`, but types of ``rci`` and
   ``x`` changed to support type ``double complex``.

.. c:function:: void spral_ssmfe_buckling_double(struct spral_ssmfe_rcid *rci, double sigma, int left, int right, int mep, double *lambda, int n, double *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   Computes eigenpairs of the buckling problem

   .. math:: Bx = \lambda A x

   in the vicinity of a given value :math:`\sigma`.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci.job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see :c:member:`inform.flag                 |
   |          | <spral_ssmfe_inform.flag>`.                                   |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see                     |
   |          | :c:member:`inform.flag <spral_ssmfe_inform.flag>`.            |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`Y = AX`.                                     |
   +----------+---------------------------------------------------------------+
   |  3       | Calculate :math:`Y = BX`.                                     |
   +----------+---------------------------------------------------------------+
   |  9       | Solve :math:`(B-\sigma A)Y = X` for Y.                        |
   +----------+---------------------------------------------------------------+

   The matrices :math:`X` and :math:`Y` are components of `rci`.

   :param rci: Reverse communication type.
      :c:member:`rci.job <spral_ssmfe_rcid.job>` must be
      set to `0` before the first call.
   :param sigma: Shift value :math:`sigma`.
   :param left: Number of left eigenpairs to find.
   :param right: Number of right eigenpairs to find.
   :param mep: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left+right`.
   :param lambda[mep]: Current eigenvalue estimates in ascending
      order.
   :param n: Size of matrix :math:`A`.
   :param x[n][ldx]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      :c:member:`options.user_x>0 <spral_ssmfe_options.user_x>`.
   :param ldx: Leading dimension of `x`.
   :param keep: Internal workspace used by routine.
   :param options: specifies algorithm options to be used.
   :param inform: returns information about the exection of the
      routine.

.. c:function:: void spral_ssmfe_buckling_double_complex(struct spral_ssmfe_rciz *rci, double sigma, int left, int right, int mep, double *lambda, int n, double complex *x, int ldx, void **keep, const struct spral_ssmfe_options *options, struct spral_ssmfe_inform *inform)

   As :c:func:`spral_ssmfe_buckling_double()`, but types of ``rci`` and
   ``x`` changed to support type ``double complex``.

.. c:function:: void spral_ssmfe_free_double(void **keep, struct spral_ssmfe_inform *inform)

   Free memory allocated in `keep` and `inform`.

   :param keep: Workspace to be freed.
   :param inform: Information type to be freed.

   .. warning::

      As memory in ``keep`` and ``inform`` has been allocated using Fortran
      functions, this routine **must** be called to avoid a memory leak.

.. c:function:: void spral_ssmfe_free_double_complex(void **keep, struct spral_ssmfe_inform *inform)

   As :c:func:`spral_ssmfe_free_double()`, but for `double complex` versions
   of types.

=============
Derived types
=============

.. c:type:: struct spral_ssmfe_rcid

   Real-valued reverse communication interface (RCI) type.

   .. c:member:: int job

      Reverse-communication task to perform.

   .. c:member:: int nx

      Number of columns in `x` and `y`.

   .. c:member:: double x[nx][n]

      Vector to be transformed by RCI task. Allocated by routine.

   .. c:member:: double y[nx][n]

      Vector to store result of RCI task. Allocated by routine.

.. c:type:: struct spral_ssmfe_rciz

   Complex-valued reverse communication interface (RCI) type.

   .. c:member:: int job

      Reverse-communication task to perform.

   .. c:member:: int nx

      Number of columns in `x` and `y`.

   .. c:member:: double complex x[nx][n]

      Vector to be transformed by RCI task. Allocated by routine.

   .. c:member:: double complex y[nx][n]

      Vector to store result of RCI task. Allocated by routine.

.. c:type:: struct spral_ssmfe_options

   Options that control the algorithm.

   .. c:member:: double abs_tol_lambda.

      Absolute tolerance for estimated
      eigenvalue convergence test, see :ref:`method section <method>`.
      Negative values are treated as the default.
      Default is `0.0`.

   .. c:member:: double abs_tol_residual

      Absolute tolerance for residual
      convergence test, see :ref:`method section <method>`.
      Negative values are treated as the default.
      Default is `0.0`.

   .. c:member:: int max_iterations

      Maximum number of iterations.
      Default is `100`.

   .. c:member:: double rel_tol_lambda

      Relative tolerance for estimated eigenvalue error convergence test, see
      :ref:`method section <method>`. Negative
      values are treated as the default.
      Default is `0.0`.

   .. c:member:: double rel_tol_residual

      Relative tolerance for residual
      convergence test, see :ref:`method section <method>`. If both
      `abs_tol_residual` and `rel_tol_residual` are 0.0, then the
      residual norms are not taken into consideration by the convergence
      test. Negative values are treated as the default.
      Default is `0.0`.

   .. c:member:: double tol_x

      Tolerance for estimated eigenvector error
      convergence test, see :ref:`method section <method>`.
      If tol_x is set to `0.0`, the eigenvector error is not estimated. If
      a negative value is assigned, the tolerance is set to
      `sqrt(DBL_EPSILON)`.
      Default is `-1.0`.

   .. c:member:: int print_level

      Amount of printing. Possible values are:

      +----+------------------------------------------------------------------+
      | <0 | no printing                                                      |
      +----+------------------------------------------------------------------+
      |  0 | error and warning messages only                                  |
      +----+------------------------------------------------------------------+
      |  1 | the type (standard or generalized) and the size of the problem,  |
      |    | the number of eigenpairs requested, the error tolerances and the |
      |    | size of the subspace are printed before the iterations start     |
      +----+------------------------------------------------------------------+
      |  2 | as above but, for each eigenpair tested for convergence, the     |
      |    | iteration number, the index of the eigenpair, the eigenvalue,    |
      |    | whether it has converged, the residual norm, and the error       |
      |    | estimates are also printed                                       |
      +----+------------------------------------------------------------------+
      | >2 | as 1 but with all eigenvalues, whether converged, residual norms |
      |    | and eigenvalue/eigenvector error estimates printed on each       |
      |    | iteration.                                                       |
      +----+------------------------------------------------------------------+

      Note that for eigenpairs that are far from convergence, ‘rough’ error
      estimates are printed (the estimates that are actually used by the
      stopping criteria, see Section [ssmfe:method], only become available on
      the last few iterations).

      Default is `0`.

   .. c:member:: int unit_error

      Fortran unit number for error messages. Printing
      suppressed if negative.
      Default is `6` (stdout).

   .. c:member:: int unit_diagnostic

      Fortran unit number for diagnostic messages.
      Printing suppressed if negative.
      Default is `6` (stdout).


   .. c:member:: int unit_warning

      Fortran unit number for warning messages.
      Printing suppressed if negative.
      Default is `6` (stdout).

   .. c:member:: double left_gap

      Minimal acceptable distance between last
      computed left eigenvalue and rest of spectrum.
      For :c:func:`spral_ssmfe_standard_double` and
      :c:func:`spral_ssmfe_generalized_double` the
      last computed left eigenvalue is the rightmost of those computed.
      For other routines it is the leftmost.
      If set to a negative value :math:`\delta`, the minimal distance is taken
      as :math:`|\delta|` times the average distance between the computed
      eigenvalues. Note that for this option to have any effect, the value of
      `mep` must be larger than `left+right`.
      See :ref:`method section <method>` for further explanation.
      Default is `0.0`.

   .. c:member:: int max_left

      Number of eigenvalues to left of
      :math:`\sigma`, or a negative value if not known.
      Default is `-1`.

   .. c:member:: int max_right

      Number of eigenvalues to right of
      :math:`\sigma`, or a negative value if not known.
      Default is `-1`.

   .. c:member:: double right_gap

      As :c:member:`left_gap <spral_ssmfe_options.left_gap>`, but for right
      eigenvalues.
      Default is `0.0`.

   .. c:member:: int user_x

      Number of eigenvectors for which an initial
      guess is supplied in `x(:,:)` on the first call. Such eigenvectors must
      be lineraly independent.
      Default is `0`.

.. c:type:: spral_ssmfe_inform

   Information on progress of the algorithm.

   .. c:member:: int flag

      Return status of algorithm. See table below.

   .. c:member:: int iteration

      Number of iterations.

   .. c:member:: int left

      Number of converged left eigenvalues.

   .. c:member:: double next_left

      Upon completion, next left eigenvalue in spectrum
      (see :c:member:`options.left_gap <spral_ssmfe_options.left_gap>`).

   .. c:member:: double next_right

      Upon completion, next right eigenvalue in spectrum
      (see :c:member:`options.right_gap <spral_ssmfe_options.right_gap>`).

   .. c:member:: int non_converged

      Number of non-converged eigenpairs.

   .. c:member:: int right

      Number of converged right eigenvalues.

   .. c:member:: int stat

      Fortran allocation status in event of failure

   +--------------+-----------------------------------------------------------+
   | `inform.flag`|                                                           |
   +==============+===========================================================+
   |   -1         | rci.job is out-of-range.                                  |
   +--------------+-----------------------------------------------------------+
   |   -9         | n is out-of-range.                                        |
   +--------------+-----------------------------------------------------------+
   |  -10         | ldx is out-of-range.                                      |
   +--------------+-----------------------------------------------------------+
   |  -11         | left is out-of-range.                                     |
   +--------------+-----------------------------------------------------------+
   |  -12         | right is out-of-range.                                    |
   +--------------+-----------------------------------------------------------+
   |  -13         | mep is less than the number of desired eigenpairs.        |
   +--------------+-----------------------------------------------------------+
   | -100         | Not enough memory; `inform.stat` contains the value of the|
   |              | Fortran stat parameter.                                   |
   +--------------+-----------------------------------------------------------+
   | -200         | :math:`B` is not positive definite or `user_x>0` and      |
   |              | linearly dependent initial guesses were supplied.         |
   +--------------+-----------------------------------------------------------+
   |   +1         | The iterations have been terminated because no further    |
   |              | improvement in accuracy is possible (this may happen if   |
   |              | :math:`B` or the preconditioner is not positive definite, |
   |              | or if the components of the residual vectors are so small |
   |              | that the round-off errors make them essentially random).  |
   |              | The value of `inform.non_converged` is set to the number  |
   |              | of non-converged eigenpairs.                              |
   +--------------+-----------------------------------------------------------+
   |   +2         | The maximum number of iterations `max_iterations` has been|
   |              | exceeded. The value of `inform.non_converged` is set to   |
   |              | the number of non-converged eigenpairs.                   |
   +--------------+-----------------------------------------------------------+
   |   +3         | The solver had run out of storage space for the converged |
   |              | eigenpairs before the gap in the spectrum required by     |
   |              | `options.left_gap` and/or `options.right_gap` was reached.|
   |              | The value of `inform.non_converged` is set to the number  |
   |              | of non-converged eigenpairs.                              |
   +--------------+-----------------------------------------------------------+

   If the computation is terminated with the error code 2 or 3, the computation
   is not complete, but may be restarted with larger values of `max_iterations`
   and/or `mep`. In this case the user should set `options.user_x` to
   `info.left + info.right` and restart the reverse communication loop. An
   alternative option is to use one of the advanced solver procedures from
   :doc:`ssmfe_expert` or :doc:`ssmfe_core` that delegate the
   storage of computed eigenpairs and the termination of the computation to the
   user.

.. _example:

========
Examples
========

Preconditioning example
-----------------------

The following code computes the 5 leftmost eigenpairs of the matrix
:math:`A` of order 100 that approximates the two-dimensional Laplacian
operator on a 20-by-20 grid. One forward and one backward Gauss-Seidel
update are used for preconditioning, which halves the number of
iterations compared with solving the same problem without
preconditioning. The header `laplace2d.h`
(examples/C/ssmfe/laplace2d.h) supplies a subroutine
`apply_laplacian()` that multiplies a block of vectors by :math:`A`, and
a subroutine `apply_gauss_seidel_step()` that computes :math:`y = T x`
for a block of vectors :math:`x` by applying one forward and one
backward update of the Gauss-Seidel method to the system
:math:`A y = x`.

.. literalinclude:: ../../examples/C/ssmfe/precond_ssmfe.c
   :language: C

This code produces the following output::

   6 eigenpairs converged in 19 iterations
    lambda[0] = 4.4676695e-02
    lambda[1] = 1.1119274e-01
    lambda[2] = 1.1119274e-01
    lambda[3] = 1.7770878e-01
    lambda[4] = 2.2040061e-01
    lambda[5] = 2.2040061e-01

Note that the code computed one extra eigenpair because of the
insufficient gap between the 5th and 6th eigenvalues.

Shift-and-invert example
------------------------

The following code computes the eigenpairs of the matrix of order 64
that approximates the two-dimensional Laplacian operator on 8-by-8 grid
with eigenvalues near the shift `sigma=1.0`. For the shifted
solve, LAPACK subroutines DSYTRS and DSYTRF are used, which perform the
LDLT-factorization and the solution of the factorized system
respectively. The matrix of the discretized Laplacian is computed by the
subroutine `set_2d_laplacian_matrix()` from the `laplace2d.h` header
(examples/C/ssmfe/laplace2d.h). The header `ldltf.h`
(examples/C/ssmfe/ldltf.h) supplies the function
`num_neg_D()` that counts the number of negative eigenvalues of the
D-factor.

.. literalinclude:: ../../examples/C/ssmfe/shift_invert.c
   :language: C

This code produces the following output::

   Eigenvalues near 1.000000e+00 (took 5 iterations)
    lambda[0] = 2.4122952e-01
    lambda[1] = 5.8852587e-01
    lambda[2] = 5.8852587e-01
    lambda[3] = 9.3582223e-01
    lambda[4] = 1.1206148e+00
    lambda[5] = 1.1206148e+00
    lambda[6] = 1.4679111e+00
    lambda[7] = 1.4679111e+00
    lambda[8] = 1.7733184e+00

Hermitian example
-----------------

The following code computes the 5 leftmost eigenpairs of the
differential operator :math:`i \frac{d}{dx}` acting in the space of
periodic functions discretized by central differences on a uniform mesh
of 80 steps.

.. literalinclude:: ../../examples/C/ssmfe/hermitian.c
   :language: C

This code produces the following output::

   5 eigenpairs converged in 25 iterations
    lambda[0] = -2.0000000e+00
    lambda[1] = -1.9938347e+00
    lambda[2] = -1.9938347e+00
    lambda[3] = -1.9753767e+00
    lambda[4] = -1.9753767e+00

.. _method:

======
Method
======

:f:mod:`spral_ssmfe_core``, upon which :f:mod:`spral_ssmfe` is built,
implements a block iterative algorithm based on the Jacobi-conjugate
preconditioned gradients (JCPG) method [2]_, [3]_. This algorithm
simultaneously computes :math:`m < n` approximate eigenpairs, where the block
size :math:`m` exceeds the number :math:`n_e` of desired eigenpairs for the
sake of better convergence, namely, :math:`m = n_e + \min(10, 0.1 n_e)`.

An approximate eigenpair :math:`\{x,\lambda\}` is considered to have
converged if the following three conditions are all satisfied:

#. if `options.abs_tol_lambda` and `options.rel_tol_lambda` are not both
   equal to zero, then the estimated error in the approximate eigenvalue
   must be less than or equal to
   :math:`\max(\mathrm{options.abs\_tol\_lambda}, \delta*\mathrm{options.rel\_tol\_lambda})`,
   where :math:`\delta` is the estimated average distance between
   eigenvalues.

#. if `options.tol_x` is not zero, then the estimated sine of the angle
   between the approximate eigenvector and the invariant subspace
   corresponding to the eigenvalue approximated by :math:`\lambda` must
   be less than or equal to `options.tol_x`.

#. if `options.abs_tol_residual` and `options.rel_tol_residual` are not
   both equal to zero, then the Euclidean norm of the residual,
   :math:`\|A x - \lambda B x\|_2`, must be less than or equal to
   :math:`\max(\mathrm{options.abs\_tol\_residual}, \mathrm{options.rel\_tol\_residual}*\|\lambda B x\|_2)`.

The extra eigenpairs are not checked for convergence, as their role is
purely auxiliary.

If the gap between the last computed eigenvalue and the rest of the
spectrum is small, then the accuracy of the corresponding eigenvector
may be very low. To prevent this from happening, the user should set the
eigenpairs storage size mep to a value that is larger than the number of
desired eigenpairs, and set the options `options.left_gap` and
`options.right_gap` to non-zero values :math:`\delta_l` and
:math:`\delta_r`. These values determine the size of the minimal
acceptable gaps between the computed eigenvalues and the rest of the
spectrum, :math:`\delta_l` referring to either leftmost eigenvalues (for
:f:subr:`ssmfe_standard()` and :f:subr:`ssmfe_generalized()` only) or those
to the left of the shift `sigma`, and :math:`\delta_r` to those to the right of
the shift `sigma`. Positive values of :math:`\delta_l` and :math:`\delta_r` set
the gap explicitely, and negative values require the gap to be not less than
their absolute value times the average distance between the computed
eigenvalues. A recommended value of :math:`\delta_l` and
:math:`\delta_r` is `-0.1`. The value of `mep` has little effect on
the speed of computation, hence it might be set to any reasonably large
value. The larger the value of `mep`, the larger the size of an eigenvalue
cluster for which accurate eigenvectors can be computed, notably: to
safeguard against clusters of size up to :math:`k`, it is sufficient to
set `mep` to the number of desired eigenpairs plus :math:`k - 1`.

When using the solver procedures that employ the shift-and-invert
technique, it is very important to ensure that the numbers of desired
eigenvalues each side of the shift do not exceed the actual numbers of
these eigenvalues, as the eigenpairs ‘approximating’ non-existing
eigenpairs of the problem will not converge. It is therefore strongly
recommended that the user employs a linear system solver that performs
the :math:`LDL^T` factorization of the shifted system, e.g. `HSL_MA97` or
`SPRAL_SSIDS`. The :math:`LDL^T` factorization of the matrix
:math:`A - \sigma B` consists in finding a lower triangular matrix :math:`L`, a
block-diagonal matrix :math:`D` with :math:`1\times 1` and
:math:`2\times 2` blocks on the diagonal and a permutation matrix
:math:`P` such that :math:`P^T(A - \sigma B)P = L D L^T`. By the inertia
theorem, the number of eigenvalues to the left and right from the shift
:math:`\sigma` is equal to the number of negative and positive
eigenvalues of :math:`D`, which allows quick computation of the
eigenvalue numbers each side of the shift.

References
----------

.. [1] E. E. Ovtchinnikov and J. Reid (2010).
   *A preconditioned block conjugate gradient algorithm for computing extreme
   eigenpairs of symmetric and Hermitian problems*.
   Technical Report RAL-TR-2010-19.

.. [2] E. E. Ovtchinnikov (2008).
   *Jacobi correction equation, line search and conjugate gradients in
   Hermitian eigenvalue computation I: Computing an extreme eigenvalue*.
   SIAM J. Numer. Anal., 46:2567–2592.

.. [3] E. E. Ovtchinnikov (2008).
   *Jacobi correction equation, line search and conjugate gradients in
   Hermitian eigenvalue computation II: Computing several extreme eigenvalues*.
   SIAM J. Numer. Anal., 46:2593–2619.
