*****************************************************************************************
:f:mod:`spral_ssmfe_expert` - Sparse Symmetric Matrix-Free Eigensolver (Expert interface)
*****************************************************************************************
.. f:module:: spral_ssmfe_expert
   :synopsis: Sparse Symmetric Matrix-Free Eigensolver (Expert interface)

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

The module :f:mod:`spral_ssmfe` provides a more user-friendly wrapper around
this code. Conversely, :f:mod:`spral_ssmfe_core` provides a lower level
implementation of the core solver, which this package provides a wrapper for.


Major version history
---------------------

2014-11-20 Version 1.0.0
    Initial release

[for detail please see ChangeLog]

==============
Usage overview
==============

The eigensolver subroutines behind :f:mod:`spral_ssmfe_expert` implement a
block iterative algorithm. The block nature of this algorithm allows the
user to benefit from highly optimized linear algebra subroutines and
from the ubiquitous multicore architecture of modern computers. It also
makes this algorithm more reliable than Krylov-based algorithms employed
e.g. by ARPACK in the presence of clustered eigenvalues. However,
convergence of the iterations may be slow if the density of the spectrum
is high.

Thus, good performance (in terms of speed) is contingent on the
following two factors:

i. the number of desired eigenpairs must be substantial (e.g. not fewer than
   the number of CPU cores), and
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
(see a simple example in Section [ssmfe:example:precond]). Note: This
technique is only recommended for experienced users.

In this expert interface, the user must handle storage of all vectors,
facilitating advanced memory handling techniques required for parallel,
hybrid and/or out-of-core execution. If there is no requirement to store
these vectors, consider using the simplified interface of :f:mod:`spral_ssmfe`
instead.

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

.. f:subroutine:: ssmfe_standard(rci,left,mep,lambda,m,rr,ind,keep,options,inform)

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
   | -2       | Failed to converge, see `inform%flag`.                        |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  1       | Calculate :math:`\bar{V} = AU`.                               |
   +----------+---------------------------------------------------------------+
   |  2       | Apply preconditioner :math:`\bar{V} = TU`. (Copy if T=I).     |
   +----------+---------------------------------------------------------------+
   |  5       | Copy converged eigenvectors :math:`X` to user storage:        |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%kx)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%kx)``.       |
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
   | 21       | Orthogonalize columns of :math:`V` to all vectors :math:`X`   |
   |          | by solving                                                    |
   |          |                                                               |
   |          | .. math:: (X^TX) Q = X^T \bar{V}                              |
   |          |                                                               |
   |          | for :math:`Q` and updating                                    |
   |          |                                                               |
   |          | .. math:: U = U - XQ                                          |
   +----------+---------------------------------------------------------------+
   | 22       | Orthogonalize columns of :math:`U` to all vectors :math:`X`   |
   |          | by solving                                                    |
   |          |                                                               |
   |          | .. math:: (X^TX) Q = X^T U                                    |
   |          |                                                               |
   |          | for :math:`Q` and updating                                    |
   |          |                                                               |
   |          | .. math:: U = U - XQ                                          |
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
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs. See [method] section for
      guidance on selecting a good value. Must be at least `left`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer m [in]: Block size of workspace `W`. Must be at least `2`.
   :p real rr (2*m,2*m,3) [inout]: reverse communication workspace.
      (Type `complex` in complex version).
   :p integer ind (m) [inout]: reverse communication workspace.
   :p ssmfe_expert_keep keep [inout]: Internal workspace used by routine.
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of
      the routine.

.. f:subroutine:: ssmfe_standard_shift(rci,sigma,left,right,mep,lambda,m,rr,ind,keep,options,inform)

   Computes eigenpairs of the standard eigenvalue problem

   .. math:: Ax = \lambda x

   in the vicinity of a given value :math:`sigma`.

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
   |  2       | Apply preconditioner :math:`\bar{V} = TU`. (Copy if T=I).     |
   +----------+---------------------------------------------------------------+
   |  5       | Copy converged eigenvectors :math:`X` to user storage:        |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%kx)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%kx)``.       |
   +----------+---------------------------------------------------------------+
   |  9       | Compute :math:`V = (A-\sigma I)^{-1} U`                       |
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
   | 21       | Orthogonalize columns of :math:`V` to all vectors :math:`X`   |
   |          | by solving                                                    |
   |          |                                                               |
   |          | .. math:: (X^TX) Q = X^T \bar{V}                              |
   |          |                                                               |
   |          | for :math:`Q` and updating                                    |
   |          |                                                               |
   |          | .. math:: U = U - XQ                                          |
   +----------+---------------------------------------------------------------+
   | 22       | Orthogonalize columns of :math:`U` to all vectors :math:`X`   |
   |          | by solving                                                    |
   |          |                                                               |
   |          | .. math:: (X^TX) Q = X^T U                                    |
   |          |                                                               |
   |          | for :math:`Q` and updating                                    |
   |          |                                                               |
   |          | .. math:: U = U - XQ                                          |
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
   :p real sigma [in]: Shift value :math:`sigma`.
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer right [in]: Number of right eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs. See [method] section for
      guidance on selecting a good value. Must be at least `left+right`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer m [in]: Block size of workspace `W`. Must be at least `2`.
   :p real rr (2*m,2*m,3) [inout]: reverse communication workspace.
      (Type `complex` in complex version).
   :p integer ind (m) [inout]: reverse communication workspace.
   :p ssmfe_expert_keep keep [inout]: Internal workspace used by routine.
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of
      the routine.

.. f:subroutine:: ssmfe_generalized(rci,left,mep,lambda,m,rr,ind,keep,options,inform)

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
   |  5       | Copy converged eigenvectors :math:`X` to user storage:        |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%kx)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%kx)``.       |
   |          |                                                               |
   |          | Optionally save their :math:`B`-images:                       |
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
   |          |    \bar{V} & = & \bar{V} - BXQ                                |
   |          |                                                               |
   |          | The update of :math:`\bar{V}` may be replaced by              |
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
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs. See [method] section for
      guidance on selecting a good value. Must be at least `left`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer m [in]: Block size of workspace `W`. Must be at least `2`.
   :p real rr (2*m,2*m,3) [inout]: reverse communication workspace.
      (Type `complex` in complex version).
   :p integer ind (m) [inout]: reverse communication workspace.
   :p ssmfe_expert_keep keep [inout]: Internal workspace used by routine.
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of
      the routine.

.. f:subroutine:: ssmfe_generalized_shift(rci,sigma,left,right,mep,lambda,m,rr,ind,keep,options,inform)

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
   | -2       | Failed to converge, see `inform%flag`.                        |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  3       | Compute :math:`\bar{V} = BU`                                  |
   +----------+---------------------------------------------------------------+
   |  5       | Copy converged eigenvectors :math:`X` to user storage:        |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%kx)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%kx)``.       |
   |          |                                                               |
   |          | Optionally save their :math:`B`-images:                       |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%ky)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%ky)``.       |
   +----------+---------------------------------------------------------------+
   |  9       | Compute :math:`V = (A-\sigma B)^{-1} U`                       |
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
   |          |    \bar{V} & = & \bar{V} - BXQ                                |
   |          |                                                               |
   |          | The update of :math:`\bar{V}` may be replaced by              |
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
   :p real sigma [in]: Shift value :math:`sigma`.
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer right [in]: Number of right eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs. See [method] section for
      guidance on selecting a good value. Must be at least `left+right`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer m [in]: Block size of workspace `W`. Must be at least `2`.
   :p real rr (2*m,2*m,3) [inout]: reverse communication workspace.
      (Type `complex` in complex version).
   :p integer ind (m) [inout]: reverse communication workspace.
   :p ssmfe_expert_keep keep [inout]: Internal workspace used by routine.
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of
      the routine.

.. f:subroutine:: ssmfe_buckling(rci,sigma,left,right,mep,lambda,m,rr,ind,keep,options,inform)

   Computes the eigenpairs of the buckling problem

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
   | -2       | Failed to converge, see `inform%flag`.                        |
   +----------+---------------------------------------------------------------+
   | -1       | None. Computation complete.                                   |
   +----------+---------------------------------------------------------------+
   |  3       | Compute :math:`\bar{V} = BU`                                  |
   +----------+---------------------------------------------------------------+
   |  5       | Copy converged eigenvectors :math:`X` to user storage:        |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%kx)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%kx)``.       |
   |          |                                                               |
   |          | Optionally save their :math:`B`-images:                       |
   |          |                                                               |
   |          | * If `rci%i>0`: ``W(:,rci%jx:rci%jx+rci%nx-1,rci%ky)``.       |
   |          | * Else:         ``W(:,rci%jx-rci%nx+1:rci%jx,rci%ky)``.       |
   +----------+---------------------------------------------------------------+
   |  9       | Compute :math:`V = (B-\sigma A)^{-1} U`                       |
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
   |          |    \bar{V} & = & \bar{V} - BXQ                                |
   |          |                                                               |
   |          | The update of :math:`\bar{V}` may be replaced by              |
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
   :p real sigma [in]: Shift value :math:`sigma`.
   :p integer left [in]: Number of left eigenpairs to find.
   :p integer right [in]: Number of right eigenpairs to find.
   :p integer mep [in]: Number of working eigenpairs. See [method] section for
      guidance on selecting a good value. Must be at least `left+right`.
   :p real lambda (mep) [inout]: Current eigenvalue estimates in ascending
      order.
   :p integer m [in]: Block size of workspace `W`. Must be at least `2`.
   :p real rr (2*m,2*m,3) [inout]: reverse communication workspace.
      (Type `complex` in complex version).
   :p integer ind (m) [inout]: reverse communication workspace.
   :p ssmfe_expert_keep keep [inout]: Internal workspace used by routine.
   :p ssmfe_options options [in]: specifies algorithm options to be used.
   :p ssmfe_inform inform [inout]: returns information about the exection of
      the routine.

.. f:subroutine:: ssmfe_free(keep,inform)

   Free memory allocated in `keep` and `inform`. Unnecessary if both are going
   out of scope.

   :p ssmfe_expert_keep keep [inout]: Workspace to be freed.
   :p ssmfe_inform inform [inout]: Information type to be freed.

=============
Derived types
=============

.. f:type:: ssmfe_rcid

   Real-valued reverse communication interface (RCI) type.

   :f integer job: Reverse-communication task to perform.
   :f integer jx: First column of :math:`U` in block.
   :f integer kx: Block to which :math:`U` belongs.
   :f integer nx: Number of columns in :math:`U` and :math:`\bar{V}`, and
      number of rows in :math:`R`.
   :f integer jy: First column of :math:`V` in block.
   :f integer ky: Block to which :math:`V` belongs.
   :f integer ny: Number of columns in :math:`V` and :math:`R`.
   :f integer i: First row of :math:`R` in ``rr(:,:,:)``.
   :f integer j: First column of :math:`R` in ``rr(:,:,:)``.
   :f integer k: Block of :math:`R` in ``rr(:,:,:)``.
   :f real alpha: Coefficient for matrix multiplication.
   :f real beta: Coefficient for matrix multiplication.

.. f:type:: ssmfe_rciz

   Real-valued reverse communication interface (RCI) type.

   :f integer job: Reverse-communication task to perform.
   :f integer jx: First column of :math:`U` in block.
   :f integer kx: Block to which :math:`U` belongs.
   :f integer nx: Number of columns in :math:`U` and :math:`\bar{V}`, and
      number of rows in :math:`R`.
   :f integer jy: First column of :math:`V` in block.
   :f integer ky: Block to which :math:`V` belongs.
   :f integer ny: Number of columns in :math:`V` and :math:`R`.
   :f integer i: First row of :math:`R` in ``rr(:,:,:)``.
   :f integer j: First column of :math:`R` in ``rr(:,:,:)``.
   :f integer k: Block of :math:`R` in ``rr(:,:,:)``.
   :f complex alpha: Coefficient for matrix multiplication.
   :f complex beta: Coefficient for matrix multiplication.

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
   :f logical minAprod [default=.true.]: If true, minimize number of
      multiplications with :math:`A` by requiring 2 additional blocks of memory
      for the workspace ``W(:,:,:)``. Must be true for calls to
      :f:subr:`ssmfe_standard_shift()`, :f:subr:`ssmfe_generalized_shift()`,
      and :f:subr:`ssmfe_buckling()`.
   :f logical minBprod [default=.true.]: If true, minimize number of
      multiplications with :math:`B` by requiring 2 additional blocks of memory
      for the workspace ``W(:,:,:)``.
   :f real right_gap [default=0.0]: as `left_gap`, but for right eigenvalues.
   :f integer user_x [default=0]: number of eigenvectors for which an initial
      guess is supplied in `x(:,:)` on the first call. Such eigenvectors must
      be lineraly independent.


.. f:type:: ssmfe_inform

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
   :f integer left: number of converged left eigenvalues.
   :f real next_left: upon completion, next left eigenvalue in spectrum
      (see `options%left_gap`).
   :f real next_right: upon completion, next right eigenvalue in spectrum
      (see `options%right_gap`).
   :f real residual_norms (mep) [allocatable]: Euclidean norms of residuals
      for `(lambda(:), X(:))` on return with ``rci%job=5``.
   :f integer non_converged: number of non-converged eigenpairs.
   :f integer right: number of converged right eigenvalues.
   :f integer stat: allocation status in event of failure

   +--------------+-----------------------------------------------------------+
   | `inform%flag`|                                                           |
   +==============+===========================================================+
   |   -1         | rci%job is out-of-range.                                  |
   +--------------+-----------------------------------------------------------+
   |   -2         | m is out-of-range.                                        |
   +--------------+-----------------------------------------------------------+
   |   -3         | options%err_est is out-of-range.                          |
   +--------------+-----------------------------------------------------------+
   |   -4         | options%minAprod is incompatible with selected routine.   |
   +--------------+-----------------------------------------------------------+
   |   -5         | options%extra_left or options%extra_right is out-of-range.|
   +--------------+-----------------------------------------------------------+
   |   -6         | options%min_gap is out-of-range.                          |
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
for a block of vectors :math:`x` by applying one forward
and one backward update of the Gauss-Seidel method to the system
:math:`A y = x`.

.. literalinclude:: ../../examples/Fortran/ssmfe/precond_expert.f90
   :language: Fortran

This code produces the following output:

::

      6 eigenpairs converged in 129 iterations
     lambda( 1) = 4.4676695E-02
     lambda( 2) = 1.1119274E-01
     lambda( 3) = 1.1119274E-01
     lambda( 4) = 1.7770878E-01
     lambda( 5) = 2.2040061E-01
     lambda( 6) = 2.2040061E-01

Note that the code computed one extra eigenpair because of the
insufficient gap between the 5th and 6th eigenvalues.

======
Method
======

The algorithm
-------------

The solver procedures of :f:mod:`spral_ssmfe_expert` are interfaces to solver
procedures of :f:mod:`spral_ssmfe_core`, which implement a block iterative
algorithm based on the Jacobi-conjugate preconditioned gradients method
[2]_, [3]_. Further information on the algorithm used by
:f:mod:`spral_ssmfe_expert` can be found in the specification document for
:f:mod:`spral_ssmfe_core` and in [1]_.

Stopping criteria
-----------------

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


Improving eigenvector accuracy
------------------------------

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
the gap explicitly, and negative values require the gap to be not less than
their absolute value times the average distance between the computed
eigenvalues. A recommended value of :math:`\delta_l` and
:math:`\delta_r` is :math:`-0.1`. The value of `mep` has little effect on
the speed of computation, hence it might be set to any reasonably large
value. The larger the value of `mep`, the larger the size of an eigenvalue
cluster for which accurate eigenvectors can be computed, notably: to
safeguard against clusters of size up to :math:`k`, it is sufficient to
set mep to the number of desired eigenpairs plus :math:`k - 1`.

The use of shifted matrix factorization
---------------------------------------

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

Error estimation
----------------

Standard problem
~~~~~~~~~~~~~~~~

If ``options%err_est=1``, the error estimates for the eigenvalues are
based on the eigenvalues of a matrix of the form

.. math::

   \begin{aligned}
   \label{L.mx}
   \hat A = \tilde\Lambda_k - S_k^T S_k,\end{aligned}

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
of the matrix formed by the residuals :math:`r_j, j = 1, \ldots, k-1`.
If :math:`\tilde\lambda_j - \hat\lambda_j` is close to the machine
accuracy, it may be too polluted by round-off errors to rely upon. In
such case, we use instead

.. math::

   \begin{aligned}
   \tilde\lambda_j - \lambda_j \le \delta_j \approx
   \frac{\|r_j\|^2}{\tilde\lambda_k - \lambda_j}.\end{aligned}

The eigenvector errors are estimated based on the Davis-Kahan
inequality:

.. math::

   \begin{aligned}
   \min_{x \in \mathcal{X}_{k-1}}
   \sin\{\tilde x_j; x\} \le
   \frac{\|r_j\|}{\lambda_k - \tilde\lambda_j} \approx
   \frac{\|r_j\|}{\tilde\lambda_k - \tilde\lambda_j},\end{aligned}

where :math:`\mathcal{X}_{k-1}` is the invariant subspace corresponding
to :math:`k-1` leftmost eigenvalues.

If ``options%err_est=2`` the errors are estimated based on the
eigenvalue decrements history, which produces an estimate for the
average eigenvalue error reduction per iteration, which in turn yields
error estimates for both eigenvalues and eigenvectors. Unlike the
residual estimates mentioned in this section, such ‘kinematic’ error
estimates are not guaranteed to be upper bounds for the actual errors.
However, the numerical tests have demonstrated that kinematic error
estimates are significantly more accurate, i.e. closer to the actual
error, than the residual-based estimates. Furthermore, they
straightforwardly apply to the generalized case as well.

Generalized problem
~~~~~~~~~~~~~~~~~~~

In the case of the generalized eigenvalue problem solved by iterations
with preconditioning, all of the residual norms in the previous section
must be replaced with :math:`\|\cdot\|_{B^{-1}}`-norm of the residual
:math:`r_j = A \tilde x_j - \tilde\lambda_j B \tilde x_j`
(:math:`\|r_j\|_{B^{-1}}^2 = r_j^* B^{-1} r_j`) or its upper estimate,
e.g. :math:`\beta_1^{-1/2}\|\cdot\|`, where :math:`\beta_1` is the
smallest eigenvalue of :math:`B`. Hence, if :math:`\beta_1` is known,
then the error tolerances for eigenvalues and eigenvectors must be
multiplied by :math:`\beta_1` and :math:`\sqrt{\beta_1}` respectively.
If no estimate for :math:`\|\cdot\|_{B^{-1}}`-norm is available, then
the use of non-zero residual tolerances and
``options%err_est=1`` is not recommended. In the case of
problems solved by iterations with shift-and-invert and the problem ,
the residuals are computed as
:math:`r_j = T B \tilde x_j - \tilde \lambda_j \tilde x_j` where
:math:`T = (A - \sigma B)^{-1}` for and :math:`T = (B - \sigma A)^{-1}`
for , and :math:`B`-norms of :math:`r_j` are used, so that Lehmann
matrix becomes :math:`\hat A = \tilde\Lambda_k - S_k^T B\ S_k`. 0 Note
that the residual estimates may considerably overestimate the actual
error of direct iterations because of the use of the Euclidean norm of
the residual, which is too strong a norm for it when :math:`A` is the
discretization of a differential operator.

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
