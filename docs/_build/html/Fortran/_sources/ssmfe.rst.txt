***************************************************************
:f:mod:`spral_ssmfe` - Sparse Symmetric Matrix-Free Eigensolver
***************************************************************
.. f:module:: spral_ssmfe
   :synopsis: Sparse Symmetric Matrix-Free Eigensolver

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

The eigensolver subroutines behind :f:mod:`spral_ssmfe` implement a block
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

.. f:subroutine:: ssmfe_standard(rci,left,mep,lambda,n,x,ldx,keep,options,inform)

   Computes the left-most eigenpairs of the standard eigenvalue problem

   .. math:: Ax = \lambda x

   Optionally uses preconditioning.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci%job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see `inform%flag`.                         |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see `inform%flag`.      |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`Y = AX`.                                     |
   +----------+---------------------------------------------------------------+
   |  2       | Apply preconditioner :math:`Y = TX`.                          |
   +----------+---------------------------------------------------------------+

   The matrices :math:`X` and :math:`Y` are components of `rci`.

   :p ssmfe_rcid rci [inout]: Reverse communication type. `rci%job` must be
      set to `0` before the first call. (Type :f:type:`ssmfe_rciz` in complex
      version).
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer n [in]: Size of matrix :math:`A`.
   :p real x (ldx, n) [inout]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      `options%user_x>0`. (Type complex in complex version).
   :p integer ldx [in]: Leading dimension of `x`.
   :p ssmfe_keepd keep [inout]: Internal workspace used by routine. (Type
      :f:type:`ssmfe_keepz` in complex version).
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of the
      routine.

.. f:subroutine:: ssmfe_standard_shift(rci,sigma,left,right,mep,lambda,n,x,ldx,keep,options,inform)

   Computes eigenpairs of the standard eigenvalue problem

   .. math:: Ax = \lambda x

   in the vicinity of a given value :math:`\sigma`.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci%job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see `inform%flag`.                         |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see `inform%flag`.      |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`Y = AX`.                                     |
   +----------+---------------------------------------------------------------+
   |  9       | Solve :math:`(A-\sigma I)Y = X` for Y.                        |
   +----------+---------------------------------------------------------------+

   The matrices :math:`X` and :math:`Y` are components of `rci`.

   :p ssmfe_rcid rci [inout]: Reverse communication type. `rci%job` must be
      set to `0` before the first call. (Type :f:type:`ssmfe_rciz` in complex
      version).
   :p real sigma [in]: Shift value :math:`sigma`.
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer right [in]: Number of right eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left+right`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer n [in]: Size of matrix :math:`A`.
   :p real x (ldx, n) [inout]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      `options%user_x>0`. (Type complex in complex version).
   :p integer ldx [in]: Leading dimension of `x`.
   :p ssmfe_keepd keep [inout]: Internal workspace used by routine. (Type
      :f:type:`ssmfe_keepz` in complex version).
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of the
      routine.

.. f:subroutine:: ssmfe_generalized(rci,left,mep,lambda,n,x,ldx,keep,options,inform)

   Computes the left-most eigenpairs of the generalized eigenvalue problem

   .. math:: Ax = \lambda B x

   Optionally uses preconditioning.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci%job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see `inform%flag`.                         |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see `inform%flag`.      |
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

   :p ssmfe_rcid rci [inout]: Reverse communication type. `rci%job` must be
      set to `0` before the first call. (Type :f:type:`ssmfe_rciz` in complex
      version).
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer n [in]: Size of matrix :math:`A`.
   :p real x (ldx, n) [inout]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      `options%user_x>0`. (Type complex in complex version).
   :p integer ldx [in]: Leading dimension of `x`.
   :p ssmfe_keepd keep [inout]: Internal workspace used by routine. (Type
      :f:type:`ssmfe_keepz` in complex version).
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of the
      routine.

.. f:subroutine:: ssmfe_generalized_shift(rci,sigma,left,right,mep,lambda,n,x,ldx,keep,options,inform)

   Computes eigenpairs of the generalized eigenvalue problem

   .. math:: Ax = \lambda B x

   in the vicinity of a given value :math:`\sigma`.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci%job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see `inform%flag`.                         |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see `inform%flag`.      |
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

   :p ssmfe_rcid rci [inout]: Reverse communication type. `rci%job` must be
      set to `0` before the first call. (Type :f:type:`ssmfe_rciz` in complex
      version).
   :p real sigma [in]: Shift value :math:`sigma`.
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer right [in]: Number of right eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left+right`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer n [in]: Size of matrix :math:`A`.
   :p real x (ldx, n) [inout]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      `options%user_x>0`. (Type complex in complex version).
   :p integer ldx [in]: Leading dimension of `x`.
   :p ssmfe_keepd keep [inout]: Internal workspace used by routine. (Type
      :f:type:`ssmfe_keepz` in complex version).
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of the
      routine.

.. f:subroutine:: ssmfe_buckling(rci,sigma,left,right,mep,lambda,n,x,ldx,keep,options,inform)

   Computes eigenpairs of the buckling problem

   .. math:: Bx = \lambda A x

   in the vicinity of a given value :math:`\sigma`.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci%job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see `inform%flag`.                         |
   +----------+---------------------------------------------------------------+
   | -2       | Restart computation. Non-fatal error, see `inform%flag`.      |
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

   :p ssmfe_rcid rci [inout]: Reverse communication type. `rci%job` must be
      set to `0` before the first call. (Type :f:type:`ssmfe_rciz` in complex
      version).
   :p real sigma [in]: Shift value :math:`sigma`.
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer right [in]: Number of right eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs.
      See :ref:`method section <method>` for guidance on selecting a good
      value. Must be at least `left+right`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer n [in]: Size of matrix :math:`A`.
   :p real x (ldx, n) [inout]: Current eigenvector estimates corresponding to
      eigenvalues in `lambda`. Used to supply initial estimates if
      `options%user_x>0`. (Type complex in complex version).
   :p integer ldx [in]: Leading dimension of `x`.
   :p ssmfe_keepd keep [inout]: Internal workspace used by routine. (Type
      :f:type:`ssmfe_keepz` in complex version).
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of the
      routine.

.. f:subroutine:: ssmfe_free(keep,inform)

   Free memory allocated in `keep` and `inform`. Unnecessary if both are going
   out of scope.

   :p ssmfe_keepd keep [inout]: Workspace to be freed.
      (Type :f:type:`ssmfe_keepz` in complex version).
   :p ssmfe_inform inform [inout]: Information type to be freed.

=============
Derived types
=============

.. f:type:: ssmfe_rcid

   Real-valued reverse communication interface (RCI) type.

   :f integer job: Reverse-communication task to perform.
   :f integer nx: Number of columns in `x` and `y`.
   :f real x (n, nx): Vector to be transformed by RCI task.
   :f real y (n, nx): Vector to store result of RCI task.

.. f:type:: ssmfe_rciz

   Complex-valued reverse communication interface (RCI) type.

   :f integer job: Reverse-communication task to perform.
   :f integer nx: Number of columns in `x` and `y`.
   :f complex x (n, nx): Vector to be transformed by RCI task.
   :f complex y (n, nx): Vector to store result of RCI task.

.. f:type:: ssmfe_options

   Options that control the algorithm.

   :f real abs_tol_lambda [default=0.0]: absolute tolerance for estimated
      eigenvalue convergence test, see Section [ssmfe:method].
      Negative values are treated as the default.
   :f real abs_tol_residual [default=0.0]: absolute tolerance for residual
      convergence test, see Section [ssmfe:method]. Negative values are
      treated as the default.
   :f integer max_iterations [default=100]: maximum number of iterations.
   :f real rel_tol_lambda [default=0.0]: relative tolerance for estimated
      eigenvalue error convergence test, see Section [ssmfe:method]. Negative
      values are treated as the default.
   :f real rel_tol_residual [default=0.0]: relative tolerance for residual
      convergence test, see Section [ssmfe:method]. If both
      `abs_tol_residual` and `rel_tol_residual` are 0.0, then the
      residual norms are not taken into consideration by the convergence
      test, see Section [ssmfe:method]. Negative values are treated as the
      default.
   :f real tol_x [default=-1.0]: tolerance for estimated eigenvector error
      convergence test, see Section [ssmfe:method].
      If tol_x is set to `0.0`, the eigenvector error is not estimated. If
      a negative value is assigned, the tolerance is set to
      `sqrt(epsilon(lambda))`.
   :f integer print_level [default=0]: amount of printing. Possible values
      are:

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

   :f integer unit_error [default=6]: unit number for error messages. Printing
      suppressed if negative.
   :f integer unit_diagnostic [default=6]: unit number for diagnostic messages.
      Printing suppressed if negative.
   :f integer unit_warning [default=6]: unit number for warning messages.
      Printing suppressed if negative.
   :f real left_gap [default=0.0]: minimal acceptable distance between last
      computed left eigenvalue and rest of spectrum.
      For :f:subr:`ssmfe_standard()` and :f:subr:`ssmfe_generalized()` the
      last computed left eigenvalue is the rightmost of those computed.
      For other routines it is the leftmost.
      If set to a negative value :math:`\delta`, the minimal distance is taken
      as :math:`|\delta|` times the average distance between the computed
      eigenvalues. Note that for this option to have any effect, the value of
      `mep` must be larger than `left+right`. See Section [ssmfe:method] for
      further explanation.
   :f integer max_left [default=-1]: number of eigenvalues to left of
      :math:`\sigma`, or a negative value if not known.
   :f integer max_right [default=-1]: number of eigenvalues to right of
      :math:`\sigma`, or a negative value if not known.
   :f real right_gap [default=0.0]: as `left_gap`, but for right eigenvalues.
   :f integer user_x [default=0]: number of eigenvectors for which an initial
      guess is supplied in `x(:,:)` on the first call. Such eigenvectors must
      be lineraly independent.

.. f:type:: ssmfe_inform

   Information on progress of the algorithm.

   :f integer flag: return status of algorithm. See table below.
   :f integer iteration: number of iterations.
   :f integer left: number of converged left eigenvalues.
   :f real next_left: upon completion, next left eigenvalue in spectrum
      (see `options%left_gap`).
   :f real next_right: upon completion, next right eigenvalue in spectrum
      (see `options%right_gap`).
   :f integer non_converged: number of non-converged eigenpairs.
   :f integer right: number of converged right eigenvalues.
   :f integer stat: allocation status in event of failure

   +--------------+-----------------------------------------------------------+
   | `inform%flag`|                                                           |
   +==============+===========================================================+
   |   -1         | rci%job is out-of-range.                                  |
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
   | -100         | Not enough memory; `inform%stat` contains the value of the|
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
   |              | The value of `inform%non_converged` is set to the number  |
   |              | of non-converged eigenpairs.                              |
   +--------------+-----------------------------------------------------------+
   |   +2         | The maximum number of iterations `max_iterations` has been|
   |              | exceeded. The value of `inform%non_converged` is set to   |
   |              | the number of non-converged eigenpairs.                   |
   +--------------+-----------------------------------------------------------+
   |   +3         | The solver had run out of storage space for the converged |
   |              | eigenpairs before the gap in the spectrum required by     |
   |              | `options%left_gap` and/or `options%right_gap` was reached.|
   |              | The value of `inform%non_converged` is set to the number  |
   |              | of non-converged eigenpairs.                              |
   +--------------+-----------------------------------------------------------+

   If the computation is terminated with the error code 2 or 3, the computation
   is not complete, but may be restarted with larger values of `max_iterations`
   and/or `mep`. In this case the user should set `options%user_x` to
   `info%left + info%right` and restart the reverse communication loop. An
   alternative option is to use one of the advanced solver procedures from
   :f:mod:`spral_ssmfe_expert` or :f:mod:`spral_ssmfe_core` that delegate the
   storage of computed eigenpairs and the termination of the computation to the
   user.

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
preconditioning. The module `laplace2d`
(examples/Fortran/ssmfe/laplace2d.f90) supplies a subroutine
`apply_laplacian()` that multiplies a block of vectors by :math:`A`, and
a subroutine `apply_gauss_seidel_step()` that computes :math:`y = T x`
for a block of vectors :math:`x` by applying one forward and one
backward update of the Gauss-Seidel method to the system
:math:`A y = x`.

.. literalinclude:: ../../examples/Fortran/ssmfe/precond_ssmfe.f90
   :language: Fortran

This code produces the following output::

      6 eigenpairs converged in 19 iterations
     lambda( 1) = 4.4676695E-02
     lambda( 2) = 1.1119274E-01
     lambda( 3) = 1.1119274E-01
     lambda( 4) = 1.7770878E-01
     lambda( 5) = 2.2040061E-01
     lambda( 6) = 2.2040061E-01

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
subroutine `set_2d_laplacian_matrix()` from the `laplace2d` module
(examples/Fortran/ssmfe/laplace2d.f90). The module `ldltf`
(examples/Fortran/ssmfe/ldltf.f90) supplies the function
`num_neg_D()` that counts the number of negative eigenvalues of the
D-factor.

.. literalinclude:: ../../examples/Fortran/ssmfe/shift_invert.f90
   :language: Fortran

This code produces the following output:

::

     Eigenvalues near  1.00E+00 (took  5 iterations)
     lambda( 1) = 2.4122952E-01
     lambda( 2) = 5.8852587E-01
     lambda( 3) = 5.8852587E-01
     lambda( 4) = 9.3582223E-01
     lambda( 5) = 1.1206148E+00
     lambda( 6) = 1.1206148E+00
     lambda( 7) = 1.4679111E+00
     lambda( 8) = 1.4679111E+00
     lambda( 9) = 1.7733184E+00

Hermitian example
-----------------

The following code computes the 5 leftmost eigenpairs of the
differential operator :math:`i \frac{d}{dx}` acting in the space of
periodic functions discretized by central differences on a uniform mesh
of 80 steps.

.. literalinclude:: ../../examples/Fortran/ssmfe/hermitian.f90
   :language: Fortran

This code produces the following output::

      5 eigenpairs converged in 25 iterations
     lambda( 1) = -2.0000000E+00
     lambda( 2) = -1.9938347E+00
     lambda( 3) = -1.9938347E+00
     lambda( 4) = -1.9753767E+00
     lambda( 5) = -1.9753767E+00

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

#. if `options%abs_tol_lambda` and `options%rel_tol_lambda` are not both
   equal to zero, then the estimated error in the approximate eigenvalue
   must be less than or equal to
   :math:`\max(\mathrm{options\%abs\_tol\_lambda}, \delta*\mathrm{options\%rel\_tol\_lambda})`,
   where :math:`\delta` is the estimated average distance between
   eigenvalues.

#. if `options%tol_x` is not zero, then the estimated sine of the angle
   between the approximate eigenvector and the invariant subspace
   corresponding to the eigenvalue approximated by :math:`\lambda` must
   be less than or equal to `options%tol_x`.

#. if `options%abs_tol_residual` and `options%rel_tol_residual` are not
   both equal to zero, then the Euclidean norm of the residual,
   :math:`\|A x - \lambda B x\|_2`, must be less than or equal to
   :math:`\max(\mathrm{options\%abs\_tol\_residual}, \mathrm{options\%rel\_tol\_residual}*\|\lambda B x\|_2)`.

The extra eigenpairs are not checked for convergence, as their role is
purely auxiliary.

If the gap between the last computed eigenvalue and the rest of the
spectrum is small, then the accuracy of the corresponding eigenvector
may be very low. To prevent this from happening, the user should set the
eigenpairs storage size mep to a value that is larger than the number of
desired eigenpairs, and set the options `options%left_gap` and
`options%right_gap` to non-zero values :math:`\delta_l` and
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
