*************************************************************************************
:f:mod:`spral_ssmfe_core` - Sparse Symmetric Matrix-Free Eigensolver (Core Algorithm)
*************************************************************************************
.. f:module:: spral_ssmfe_core
   :synopsis: Sparse Symmetric Matrix-Free Eigensolver (Core Algorithm)

=======
Purpose
=======

This package computes extreme (leftmost and/or rightmost)
eigenpairs :math:`\{\lambda_i, x_i\}` of the following eigenvalue problems:

- the standard eigenvalue problem

   .. math:: A x = \lambda x,

- the generalized eigenvalue problem

   .. math:: A x = \lambda B x,

- the eigenvalue problem

   .. math:: AB x = \lambda x

where :math:`A` and :math:`B` are **real symmetric** (or **Hermitian**) matrices
and :math:`B` is **positive definite**.

The modules :f:mod:`spral_ssmfe` and :f:mod:`spral_ssmfe_expert` provide a
simplified interface to this routine, and should be used if the user does not
require access to low level features provided in this package.


Major version history
---------------------

2014-11-20 Version 1.0.0
    Initial release


==============
Usage overview
==============

:f:mod:`spral_ssmfe_core` implements a block iterative algorithm for
simultaneous computation of several eigenpairs for the problems above.
The block nature of this
algorithm allows the user to benefit from highly optimized linear
algebra subroutines and from the ubiquitous multicore architecture of
modern computers. It also makes this algorithm more reliable than
Krylov-based algorithms employed by e.g. ARPACK in the presence of
clustered eigenvalues. However, convergence of the iterations may be
slow if the density of the spectrum is high.

Thus, good performance (in terms of speed) is contingent on the
following two factors:

i. the number of desired eigenpairs must be substantial (e.g. not fewer than
   the number of CPU cores), and
ii. the employment of a convergence acceleration technique.

The acceleration techniques that can be used are shift-and-invert and
preconditioning.

The former rewrites the eigenvalue problem for a matrix :math:`M` as
:math:`Ax = \lambda x` with :math:`A = (M - \sigma I)^{-1}`, where
:math:`I` is the identity matrix and :math:`\sigma` is a real value near
eigenvalues of interest, and the generalized problem :math:`M x = \mu B x` as
the :math:`Ax = \lambda B x` with :math:`A = (M - \sigma B)^{-1}`.

The latter applies to the case of positive definite
:math:`A` and requires a matrix or an operator :math:`T`, called *a
preconditioner*, such that the vector :math:`v = T f` is an
approximation to the solution :math:`u` of the system :math:`A u = f`
(see the simple example :ref:`below <ssmfe_core_example>`). Note: This
technique is only recommended for experienced users.

For futher detail on the algorithm, see the outline in the
:ref:`method section <ssmfe_core_method>` below, and associated references.

The routine :f:mod:`spral_ssmfe` provides a user-friendly interface to this
algorithm, whilst :f:mod:`spral_ssmfe_expert` provides an interface that allows
users to manage their own memory. If this routine is used instead of
:f:mod:`spral_ssmfe_expert`, the user is additionally responsible for deciding
when a ufficient number of eigenpairs have been computed to sufficient
accuracy. The amount of computation performed by the solver subroutines
in this package and the memory they use are negligible. These features
facilitate the use of these subroutines for shared-memory, out-of-core and
hybrid computation.

===========
Subroutines
===========

To use the solver procedures, the user must maintain a workspace of `(kw+1)`
blocks each containing `m` vectors of size `n`. For notational convienience
we refer to this workspace as a Fortran array ``W(n,m,0:kw)``, but the user
is free to store it as they wish. Note the block dimension is indexed from
zero, not from one. The following table provides minimum values of `kw` for
each setup:

 +-------------------+------------+------------+------------+------------+
 |                   |        minAprod=T       |       minAprod=F        |
 +-------------------+------------+------------+------------+------------+
 | Problem           | minBprod=T | minBprod=F | minBprod=T | minBprod=F |
 +===================+============+============+============+============+
 | standard          |     7      |     5      |     3      |     3      |
 +-------------------+------------+------------+------------+------------+
 | standard_shift    |     7      |     5      |     N/A    |     N/A    |
 +-------------------+------------+------------+------------+------------+
 | generalized       |     7      |     5      |     5      |     3      |
 +-------------------+------------+------------+------------+------------+
 | generalized_shift |     7      |     7      |     N/A    |     N/A    |
 +-------------------+------------+------------+------------+------------+
 | buckling          |     7      |     7      |     N/A    |     N/A    |
 +-------------------+------------+------------+------------+------------+

Further, the user must also store the converged eigenvectors :math:`X`, and
(for generalised problems) their :math:`B`-images :math:`BX` using
separate storage, e.g. ``X(n,mep), BX(n,mep)``.
In addition to being output, the routine may need to
reorthagonalise against these from time to time.

The first (zero-indexed) block holds the eigenvector approximations: the user
must fill this block with :math:`m` linearly independent vectors before the
first call to a solver procedure.

The number of desired eigenpairs may exceed :math:`m`: whenever converged
eigenpairs have been detected, a solver procedure reports the indices of
these eigenpairs and they must be moved by the user to separate storage
(``X(:)``).

.. f:subroutine:: ssmfe(rci,problem,left,right,m,lambda,rr,ind,keep,options,inform)

   Computes specified number of leftmost and rightmost eigenvalues and
   corresponding eigenvectors.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are:

   +----------+---------------------------------------------------------------+
   | `rci%job`| Task to be performed                                          |
   +==========+===============================================================+
   | -3       | None. Fatal error, see `inform%flag`.                         |
   +----------+---------------------------------------------------------------+
   | -2       | Failed to converge, see `inform%flag`.                        |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`\bar{V} = AU`.                               |
   +----------+---------------------------------------------------------------+
   |  2       | Apply preconditioner :math:`\bar{V} = TU`. (Copy if T=I).     |
   +----------+---------------------------------------------------------------+
   |  3       | Compute :math:`\bar{V} = BU`                                  |
   +----------+---------------------------------------------------------------+
   |  4       | Test convergence for each of :math:`m` eigenvalues:           |
   |          |                                                               |
   |          | * If eigenpair `i` has converged, set inform%converged(i) to  |
   |          |   a positive number.                                          |
   |          | * Otherwise, leave at current value.                          |
   |          |                                                               |
   |          | Tests may use the estimated eigenvalue errors                 |
   |          | ``inform%err_lambda(i)`` and eigenvector errors               |
   |          | ``inform%err_x(i)``.                                          |
   +----------+---------------------------------------------------------------+
   |  5       | Copy converged eigenvectors :math:`X` to user storage:        |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%kx)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%kx)``.       |
   |          |                                                               |
   |          | Optionally save their :math:`B`-images (if :math:`B\ne I`)    |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%ky)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%ky)``.       |
   +----------+---------------------------------------------------------------+
   | 11       | If `rci%i.eq.0`, copy :math:`\bar{V} = U`.                    |
   |          |                                                               |
   |          | Otherwise, reorder columns of block `rci%kx` such that column |
   |          | `ind(j)` becomes the new column `j` for `j=1, ..., rci%nx`    |
   |          |                                                               |
   |          | Note: if ``rci%kx.eq.rci%ky``, only reorder once.             |
   +----------+---------------------------------------------------------------+
   | 12       | Compute the dot products                                      |
   |          |                                                               |
   |          | .. math:: r_{ii} = U_i \cdot \bar{V}_i                        |
   +----------+---------------------------------------------------------------+
   | 13       | Perform the scalings                                          |
   |          |                                                               |
   |          | .. math:: U_i = U_i/\sqrt{(U_i\cdot \bar{V}_i)}               |
   |          |                                                               |
   |          | and                                                           |
   |          |                                                               |
   |          | .. math:: \bar{V}_i = \bar{V}_i/\sqrt{(U_i\cdot \bar{V}_i)}   |
   |          |                                                               |
   |          | for each column :math:`U_i` and :math:`\bar{V}_i` of :math:`U`|
   |          | and :math:`\bar{V}`.                                          |
   |          |                                                               |
   |          | Note: if ``rci%kx.eq.rci%ky``, only scale once.               |
   +----------+---------------------------------------------------------------+
   | 14       | Perform the updates                                           |
   |          |                                                               |
   |          | .. math:: \bar{V}_i = \bar{V}_i + r_{ii} U_i                  |
   |          |                                                               |
   |          | for each column :math:`\bar{V}_i` of :math:`\bar{V}`          |
   +----------+---------------------------------------------------------------+
   | 15       | Perform the update                                            |
   |          |                                                               |
   |          | .. math:: R = \alpha U^T V + \beta R                          |
   +----------+---------------------------------------------------------------+
   | 16       | Perform the update                                            |
   |          |                                                               |
   |          | .. math:: V = \alpha U R + \beta V                            |
   +----------+---------------------------------------------------------------+
   | 17       | Perform the update                                            |
   |          |                                                               |
   |          | .. math:: U = \alpha U R                                      |
   |          |                                                               |
   |          | Note: :math:`V` may be used as a workspace                    |
   +----------+---------------------------------------------------------------+
   | 21       | :math:`B`-orthogonalize columns of :math:`V` to all vectors   |
   |          | :math:`X` by solving                                          |
   |          |                                                               |
   |          | .. math:: (X^TBX) Q = X^T \bar{V}                             |
   |          |                                                               |
   |          | for :math:`Q` and updating                                    |
   |          |                                                               |
   |          | .. math::                                                     |
   |          |                                                               |
   |          |    U       & = & U - XQ \\                                    |
   |          |                                                               |
   |          | If :math:`B\ne I`, the :\math:`\bar{V}` must also be updated  |
   |          | as                                                            |
   |          |                                                               |
   |          | .. math::                                                     |
   |          |                                                               |
   |          |    \bar{V} & = & \bar{V} - BXQ,                               |
   |          |                                                               |
   |          | or                                                            |
   |          |                                                               |
   |          | .. math:: \bar{V} = BU                                        |
   +----------+---------------------------------------------------------------+
   | 22       | :math:`B`-orthogonalize columns of :math:`U` to all vectors   |
   |          | :math:`X` by solving                                          |
   |          |                                                               |
   |          | .. math:: (X^TBX) Q = X^T U                                   |
   |          |                                                               |
   |          | for :math:`Q` and updating                                    |
   |          |                                                               |
   |          | .. math:: U = U - BXQ                                         |
   +----------+---------------------------------------------------------------+
   | 999      | Restart:                                                      |
   |          |                                                               |
   |          | If `rci%k>0`: Restart suggested with block size               |
   |          | `m >= rci%nx + rci%i + rci%j`, adjusting workspace size       |
   |          | to match. Set `rci%i=0` and `rci%j=0` and recall the routine. |
   |          | If a restart is not desirable, routine may be recalled with   |
   |          | no change to parameters.                                      |
   |          |                                                               |
   |          | If `rci%k=0`: Restart required with the same block size.      |
   |          |                                                               |
   |          | In both cases, the first block ``W(:,:,0)`` should retain     |
   |          | vectors ``rci%jx:rci%jx+rci%nx-1``, filling remaining vectors |
   |          | randomly such that the entire set of columns is linearly      |
   |          | independent from each other and also from the converged       |
   |          | eigenvectors.                                                 |
   +----------+---------------------------------------------------------------+

   The matrices are defined as follows:

   * :math:`U` = ``W(:, rci%jx:rci%jx+rci%nx-1, rci%kx)``
   * :math:`V` = ``W(:, rci%jy:rci%jy+rci%ny-1, rci%ky)``
   * :math:`\bar{V}` = ``W(:, rci%jy:rci%jy+rci%nx-1, rci%ky)``
   * :math:`R` = ``rr(rci%i:rci%i+rci%nx-1, rci%j:rci%j+rci%ny-1, rci%k)``

   and :math:`\alpha` and :math:`\beta` are given by ``rci%alpha`` and
   ``rci%beta`` respectively. We use the notation :math:`r_{ii}` to refer
   to the :math:`i`-th diagonal element of :math:`R`, being
   ``rr(rci%i+i-1,rci%j+i-1,rci%k)``.

   :p ssmfe_rcid rci [inout]: Reverse communication type. `rci%job` must be
      set to `0` before the first call. (Type :f:type:`ssmfe_rciz` in complex
      version).
   :p integer problem [in]: Problem to be solved, one of:

      +----+-----------------------+
      | 0  | :math:`Ax=\lambda x`  |
      +----+-----------------------+
      | >0 | :math:`Ax=\lambda Bx` |
      +----+-----------------------+
      | <0 | :math:`ABx=\lambda x` |
      +----+-----------------------+

   :p integer left [in]: Number of left eigenpairs to find. On return with
      ``rci%job=5``, can be set to zero if sufficient have been found.
   :p integer right [in]: Number of right eigenpairs to find. On return with
      ``rci%job=5``, can be set to zero if sufficient have been found.
   :p real lambda (m) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer m [in]: Block size of workspace `W`. Must be at least `2`.
   :p real rr (2*m,2*m,3) [inout]: reverse communication workspace.
      (Type `complex` in complex version).
   :p integer ind (m) [inout]: reverse communication workspace.
   :p ssmfe_core_keep keep [inout]: Internal workspace used by routine.
   :p ssmfe_core_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of
      the routine.

.. f:subroutine:: ssmfe_largest(rci,problem,nep,m,lambda,rr,ind,keep,options,inform)

   Computes specified number of eigenvalues of largest magnitude and
   corresponding eigenvectors.

   Uses reverse-communication. Upon return the user must perform a task
   specified by the `rci` parameter and recall the routine. Possible values of
   `rci` and associated tasks are as for :f:subr:`ssmfe()` above.

   :p ssmfe_rcid rci [inout]: Reverse communication type. `rci%job` must be
      set to `0` before the first call. (Type :f:type:`ssmfe_rciz` in complex
      version).
   :p integer problem [in]: Problem to be solved, one of:

      +----+-----------------------+
      | 0  | :math:`Ax=\lambda x`  |
      +----+-----------------------+
      | >0 | :math:`Ax=\lambda Bx` |
      +----+-----------------------+
      | <0 | :math:`ABx=\lambda x` |
      +----+-----------------------+

   :p integer nep [in]: Number of eigenpairs to find.
   :p real lambda (m) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer m [in]: Block size of workspace `W`. Must be at least `2`.
   :p real rr (2*m,2*m,3) [inout]: reverse communication workspace.
      (Type `complex` in complex version).
   :p integer ind (m) [inout]: reverse communication workspace.
   :p ssmfe_core_keep keep [inout]: Internal workspace used by routine.
   :p ssmfe_core_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of
      the routine.

.. f:subroutine:: ssmfe_free([keep,inform])

   Free memory allocated in `keep` and `inform`. Unnecessary if both are going
   out of scope.

   :o ssmfe_expert_keep keep [inout]: Workspace to be freed.
   :o ssmfe_inform inform [inout]: Information type to be freed.

=============
Derived types
=============

.. f:type:: ssmfe_core_options

   Options that control the algorithm.

   :f integer err_est [default=2]: error estimation scheme, one of:

      +-------------+---------------------------------------------------------+
      | 1           | Residual error bounds: modified Davis-Kahan estimate for|
      |             | eigenvector error and Lehmann bounds for eigenvale error|
      |             | (see method section).                                   |
      +-------------+---------------------------------------------------------+
      | 2 (default) | Convergence curve-based estimate.                       |
      +-------------+---------------------------------------------------------+

   :f integer extra_left [default=0]: number of extra approximate eigenvectors
      corresponding to leftmost eigenvalues used to enhance convergence.
   :f integer extra_right [default=0]: number of extra approximate eigenvectors
      corresponding to rightmost eigenvalues used to enhance convergence.
   :f logical minAprod [default=.true.]: If true, minimize number of
      multiplications with :math:`A` by requiring 2 additional blocks of memory
      for the workspace ``W(:,:,:)``. If false, three returns with `rci%job=1`
      occur per iteration instead of one. Must be true if ``problem < 0``.
   :f logical minBprod [default=.true.]: If true, minimize number of
      multiplications with :math:`B` by requiring 2 additional blocks of memory
      for the workspace ``W(:,:,:)``. If false, at least three returns with
      `rci%job=3` occur per iteration instead of one.
   :f real min_gap [default=0.0]: Restart sensitivity: if the relative
      distance between the last eigenvalue of interest on either margin of
      the spectrum and the rest of the spectrum is smaller than `min_gap`,
      the solver procedure suggests restart (`rci%job=999`).
      The values `rci%i` and `rci%j` are set to the numbers of eigenvalues on
      the left and right margin of the spectrum that are too close to the
      eigenvalues of interest, causing slow convergence. The default value
      of 0.0 means no restart is ever suggested. Must be in the range
      :math:`[0.0,1.0]`.
   :f real cf_max[default=1.0]: Stagnation sensitivity: if the value
      :math:`q_{ij}` (see method section) is greater than `cf_max` for
      :math:`i > 5`, the eigenpair is marked as stagnated by setting
      `inform%converged(j)` to a negative value. The default value of
      1.0 indicates that the estimated asymptotic convergence
      factor is not used for stagnation detection.
      Must be in the range :math:`[0.5, 1.0]`.

.. f:type:: ssmfe_infrom

   Information on progress of the algorithm.

   :f integer converged (mep) [allocatable]: Convergence status.

      * If ``converged(j)>0``, the eigenpair `(lambda(j), X(j))` converged
        on iteration `converged(j)`.
      * If ``converged(j)=0``, the eigenpair `(lambda(j), X(j))` is still
        converging.
      * If ``converged(j)<0``, the eigenpair `(lambda(j), X(j))` stagnated
        at iteration `converged(j)`.

   :f real err_lambda (mep) [allocatable]: estimated eigenvalue errors for
      converged and stagnated eigenvalues.
   :f real err_x (mep) [allocatable]: estimated eigenvector errors for
      converged and stagnated eigenvectors.
   :f integer flag: return status of algorithm. See table below.
   :f integer iteration: number of iterations.
   :f real residual_norms (mep) [allocatable]: Euclidean norms of residuals
      for `(lambda(:), X(:))` on return with ``rci%job=4, 5``.
   :f integer stat: allocation status in event of failure

   +--------------+-----------------------------------------------------------+
   | `inform%flag`|                                                           |
   +==============+===========================================================+
   |   -1         | m is out-of-range.                                        |
   +--------------+-----------------------------------------------------------+
   |   -2         | rci%job is out-of-range.                                  |
   +--------------+-----------------------------------------------------------+
   |   -3         | options%err_est is out-of-range.                          |
   +--------------+-----------------------------------------------------------+
   |   -4         | options%minAprod is incompatible with selected routine.   |
   +--------------+-----------------------------------------------------------+
   |   -5         | options%extra_left or options%extra_right is out-of-range.|
   +--------------+-----------------------------------------------------------+
   |   -6         | options%min_gap is out-of-range.                          |
   +--------------+-----------------------------------------------------------+
   |   -7         | options%cf_max is out-of-range.                           |
   +--------------+-----------------------------------------------------------+
   |  -11         | left is out-of-range (:f:subr:`ssmfe`) or                 |
   |              | nep is out-of-range (:f:subr:`ssmfe_largest`).            |
   +--------------+-----------------------------------------------------------+
   |  -12         | right is out-of-range.                                    |
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

.. _ssmfe_core_example:

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
(examples/Fortran/ssmfe/laplace2d.f90) supplies the subroutine
`apply_laplacian()` that multiplies a block of vectors by :math:`A`, and
the subroutine `apply_gauss_seidel_step()` that computes :math:`y = T x`
for a block of vectors :math:`x` by applying one forward and one
backward update of the Gauss-Seidel method to the system
:math:`A y = x`.

.. literalinclude:: ../../examples/Fortran/ssmfe/precond_core.f90
   :language: Fortran

This code produces the following output:

::

      5 eigenpairs converged in  72 iterations
     lambda(1) = 4.4676695E-02
     lambda(2) = 1.1119274E-01
     lambda(3) = 1.1119274E-01
     lambda(4) = 1.7770878E-01
     lambda(5) = 2.2040061E-01

.. _ssmfe_core_method:

======
Method
======

The algorithm
-------------

The algorithm implemented by :f:subr:`ssmfe()` and :f:subr`ssmfe_largest()`
is based on the Jacobi-conjugate preconditioned gradients (JCPG) method. The
outline of the algorithm, assuming for simplicity that :math:`0 < left \le m`
and :math:`right = 0` and using Matlab notation, is as follows.

-  **Initialization.** Perform the Rayleigh-Ritz procedure in the trial
   subspace spanned by the columns of an :math:`{\tt n} \times {\tt m}`
   matrix :math:`X` i.e. compute

   -  :math:`L = X’*A*X`, if ``problem >= 0``, and :math:`L = X’*B*A*B*X`
      otherwise,

   -  :math:`M = X’*B*X` (:math:`B=I` if ``problem = 0``),

   and solve the generalized eigenvalue problem for the matrix pencil
   :math:`L - t*M`, i.e. compute an :math:`m\times m`
   matrix :math:`Q` such that :math:`Q’*M*Q` is the identity matrix and
   :math:`D = Q’*L*Q` is a diagonal matrix. Compute :math:`X = X*Q` and set
   :math:`Z = F = []`.

-  **Main loop.**

   DO

   #. If ``problem>=0``, compute the residual matrix :math:`R = A*X - B*X*D`,
      where :math:`D` is a diagonal matrix with entries

      .. math::

         D(j,j) = \frac{X(:,j)'*A*X(:,j)}{X(:,j)'*B*X(:,j)}

      else compute :math:`R = A*B*X - X*D`, where :math:`D` is a diagonal
      matrix with entries

      .. math::

         D(j,j) = \frac{X(:,j)'*B*A*B*X(:,j)}{X(:,j)'*B*X(:,j)}.

   #. Perform the orthogonalization of :math:`R` to constraints :math:`C` by
      updating :math:`R = R - B*C*(C’*R)`.

   #. If ``options%err_est= 1``, compute :math:`R’*B*R` and use it to compute
      error bounds; otherwise only compute the diagonal of this matrix and use
      to compute error bounds. Test for converged eigenpairs and move
      converged eigenvectors from :math:`X` to :math:`C` and reduce :math:`m`
      accordingly. Exit the main loop if :math:`X = []`.

   #. If problem is non-negative, compute the preconditioned gradient
      matrix :math:`Y = T*R`.

   #. If :math:`Z` is not empty, conjugate :math:`Y` to :math:`Z`, i.e.

      a. if problem is non-negative, then compute :math:`P = Z’*A*Y`,
         otherwise compute :math:`P = Z’*B*A*B*Y`;

      b. compute :math:`S = Z’*B*Y`;

      c. update :math:`Y = Y + Z*H`, where

         .. math::

            H(i,j) = \frac{P(i,j) - S(i,j)*D(j,j)}{F(i,i) - D(j,j)}

         where :math:`F` is described in the final step below..

   #. Perform orthogonalization of the search direction matrix :math:`Y` to
      constraints :math:`C` by updating :math:`Y = Y - C*(C’*B*Y)`.

   #. :math:`B`-normalize the columns of :math:`Y` and reorder them so that
      the :math:`B`-norm of the projection of :math:`Y(:,j)` onto the linear
      span of :math:`X` and :math:`Y(:,1:j-1)` is an increasing function of
      :math:`j`. Compute :math:`M = [X Y]’*B*[X Y]`. If the condition number
      of :math:`M` is greater than the allowed limit of :math:`10^4` then
      start removing the last columns of :math:`Y` and :math:`M` and
      respective rows of :math:`M` until the condition number falls below the
      limit.

   #. Perform the Rayleigh-Ritz procedure in the linear span of columns of
      :math:`[X Y]`. Update :math:`X` by selecting Ritz vectors corresponding
      to the leftmost :math:`m` Ritz values, and place the remaining Ritz
      vectors into :math:`Z` and corresponding Ritz values onto the diagonal
      of :math:`F`.

   END DO

The orthogonalization to constraints on step 2 ensures that the
algorithm deals with the residuals for the constrained problem rather
than with those for the unconstrained one, which may not go to zeros if
the constraints are not close enough to exact eigenvectors.

The definition of :math:`H(i,j)` in step 5c ensures optimality of the new
search direction vector :math:`Y(:,j)`, notably, the search for the minimum of
the Rayleigh quotient :math:`D(j,j)` in the direction of this vector produces
the asymptotically smallest possible value.

The search directions cleanup procedure employed on step 7 ensures that the
convergence of iterations is not damaged by the computational errors of
the LAPACK eigensolver _SYGV, in the Rayleigh-Ritz procedure of step 8.
The accuracy of _SYGV is affected by the condition number of :math:`M`.
If the latter is very large, poor accuracy in the computed Ritz vectors may
lead to convergence failure. The ordering of :math:`Y(:,j)` ensures that the
‘least important’ search direction is dropped if the condition number of
:math:`M` is unacceptably large. In practice, the loss of search direction is
rare and does not lead to convergence problems.

If the number of sought eigenpairs exceeds :math:`m`, then :math:`m` is not
reduced on step 3. Instead, the approximate eigenvectors moved to :math:`C`
are replaced with vectors from :math:`Z`.

Error estimation
----------------

Standard problem
~~~~~~~~~~~~~~~~

If ``options%err_est = 1``, the error estimates for the eigenvalues are
based on the eigenvalues of a matrix of the form

.. math::

   \hat A = \tilde\Lambda_k - S_k^T S_k,

where :math:`\tilde\Lambda_k` is a diagonal matrix with the :math:`k-1`
leftmost Ritz values :math:`\tilde\lambda_j` on the diagonal, and the
columns of :math:`S_k` are the respective residual vectors
:math:`r_j = A \tilde x_j - \tilde\lambda_j \tilde x_j` divided by
:math:`\sqrt{\lambda_k - \tilde\lambda_j}`. If :math:`k` is such that
:math:`\tilde\lambda_{k-1} < \lambda_k`, then the eigenvalues of
:math:`\hat A` are the left-hand side bounds for eigenvalues
:math:`\lambda_i`, and thus the difference
:math:`\tilde\lambda_j - \hat\lambda_j` estimates the eigenvalue error
:math:`\tilde\lambda_j - \lambda_j`. The unknown :math:`\lambda_k` is
replaced by :math:`\tilde\lambda_k`, and select the maximal
:math:`k \le m` for which the distance between
:math:`\tilde\lambda_{k-1}` and :math:`\tilde\lambda_k` exceeds the sum
of the absolute error tolerance for eigenvalues and the Frobenius norm
of the matrix formed by the residuals :math:`r_j, j = 1,\ldots,k-1`. If
:math:`\tilde\lambda_j - \hat\lambda_j` is close to the machine
accuracy, it may be too polluted by round-off errors to rely upon. In
such case, we use instead

.. math::

   \tilde\lambda_j - \lambda_j \le \delta_j \approx
   \frac{\|r_j\|^2}{\tilde\lambda_k - \lambda_j}.

The eigenvector errors are estimated based on the Davis-Kahan
inequality:

.. math::

   \min_{x \in \mathcal{X}_{k-1}}
   \sin\{\tilde x_j; x\} \le
   \frac{\|r_j\|}{\lambda_k - \tilde\lambda_j} \approx
   \frac{\|r_j\|}{\tilde\lambda_k - \tilde\lambda_j},

where :math:`\mathcal{X}_{k-1}` is the invariant subspace corresponding
to :math:`k-1` leftmost eigenvalues.

If ``options%err_est = 2`` the errors are estimated based on the
eigenvalue decrements history, which produces an estimate for the
asymptotic convergence facotr, the geometrical average of the eigenvalue
error reduction per iteration:

.. math::

   q_{ij} = \left|
         \frac{\lambda_j^i - \lambda_j^{i-1}} {\lambda_j^i - \lambda_j^0}
      \right|^{\frac{1}{i}}
      \approx\left|
         \frac{\lambda_j - \lambda_j^{i-1}} {\lambda_j - \lambda_j^0}
      \right|^{\frac{1}{i}} = \left|
         \frac{\lambda_j - \lambda_j^{i-1}} {\lambda_j - \lambda_j^{i-2}}
         \cdots
         \frac{\lambda_j - \lambda_j^1} {\lambda_j - \lambda_j^0}
      \right|^{\frac{1}{i}}

where :math:`\lambda_j^i` is the approximation to :math:`\lambda_j` on
:math:`i`-th iteration (see Technical Report [1]_ for further
details). Unlike the residual estimates mentioned in this section, such
‘kinematic’ error estimates are not guaranteed to be upper bounds for
the actual errors. However, the numerical tests have demonstrated that
kinematic error estimates are significantly more accurate, i.e. closer
to the actual error, than the residual-based estimates. Furthermore,
they straightforwardly apply to the generalized case as well.

Generalized problems
~~~~~~~~~~~~~~~~~~~~

In the case of the generalized eigenvalue problem , all of the residual
norms in the previous section must be replaced with
:math:`\|\cdot\|_{B^{-1}}`-norm of the residual
:math:`r_j = A \tilde x_j - \tilde\lambda_j B \tilde x_j`
(:math:`\|r_j\|_{B^{-1}}^2 = r_j^* B^{-1} r_j`) or its upper estimate,
e.g. :math:`\beta_1^{-1/2}\|\cdot\|`, where :math:`\beta_1` is the
smallest eigenvalue of :math:`B`. Hence, if :math:`\beta_1` is known,
then the error tolerances for eigenvalues and eigenvectors must be
multiplied by :math:`\beta_1` and :math:`\sqrt{\beta_1}` respectively.
If no estimate for :math:`\|\cdot\|_{B^{-1}}`-norm is available, then
the use of non-zero residual tolerances and
``options%err_est = 1`` is not recommended. In the case of
problem , the residuals are computed as
:math:`r_j = A B \tilde x_j - \tilde \lambda_j \tilde x_j`,
:math:`B`-norms of :math:`r_j` are used in and , and Lehmann matrix
becomes :math:`\hat A = \tilde\Lambda_k - S_k^T B\ S_k`.

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
