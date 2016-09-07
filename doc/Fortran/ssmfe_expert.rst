Major version history
=====================

2014-11-20 Version 1.0.0
    Initial release

Installation
============

Please see the SPRAL install documentation. In particular note that:

-  A BLAS library is required.

-  A LAPACK library is required.

Usage overview
==============

The eigensolver subroutines behind ``SPRAL_SSMFE_EXPERT`` implement a
block iterative algorithm. The block nature of this algorithm allows the
user to benefit from highly optimized linear algebra subroutines and
from the ubiquitous multicore architecture of modern computers. It also
makes this algorithm more reliable than Krylov-based algorithms employed
e.g. by ARPACK in the presence of clustered eigenvalues. However,
convergence of the iterations may be slow if the density of the spectrum
is high.

Thus, good performance (in terms of speed) is contingent on the
following two factors: (i) the number of desired eigenpairs must be
substantial (e.g. not less than the number of CPU cores), and (ii) the
employment of a convergence acceleration technique. The acceleration
techniques that can be used are shift-and-invert and preconditioning.
The former requires the direct solution of linear systems with the
matrix :math:`A` or its linear combination with :math:`B`, for which a
sparse symmetric indefinite solver (such as HSL\_MA97 or SPRAL\_SSIDS)
can be employed. The latter applies to the case of positive definite
:math:`A` and requires a matrix or an operator (that is, an algorithm
producing a vector :math:`v = T u` for a given vector :math:`u`)
:math:`T`, called *a preconditioner*, such that the vector
:math:`v = T f` is an approximation to the solution :math:`u` of the
system :math:`A u = f`. This technique is more sophisticated and is
likely to be of interest only to experienced users.

Further information on the algorithm used by ``SPRAL_SSMFE_EXPERT`` can
be found in Section [ssmfe\ :sub:`e`\ xpert:method].

``SPRAL_SSMFE_EXPERT`` delegates a considerable part of the computation
to the user. The user’s code stores all vectors of size equal to the
problem size :math:`n`. ``SPRAL_SSMFE_EXPERT`` is not “aware” of
:math:`n` or how these vectors are stored; all operations on these
vectors are performed by the user. The amount of computation performed
by the solver subroutines of ``SPRAL_SSMFE_EXPERT`` and the memory they
use are negligible. These features facilitate the use of these
subroutines for shared-memory, out-of-core and hybrid computation. A
simpler but less flexible interface to ``SPRAL_SSMFE_EXPERT`` is offered
by ``SPRAL_SSMFE``.

Calling sequences
-----------------

| Access to the package requires a USE statement

| 

The following procedures are available to the user:

-  ssmfe\_standard() computes leftmost eigenpairs of , optionally using
   preconditioning if :math:`A` is positive definite

-  ssmfe\_standard\_shift() computes eigenpairs of near a given shift
   using the shift-and-invert technique

-  ssmfe\_generalized() computes leftmost eigenpairs of , optionally
   using preconditioning if :math:`A` is positive definite

-  ssmfe\_generalized\_shift() computes eigenpairs of near a given shift
   using the shift-and-invert technique

-  ssmfe\_buckling() computes eigenpairs of near a given shift using the
   shift-and-invert technique

-  ssmfe\_free() should be called after all other calls are complete. It
   frees memory referenced by ``keep`` and ``inform``.

The main solver procedures must be called repeatedly using a reverse
communication interface. The procedure ``ssmfe_free()`` should be called
once after the final call to a solver procedure to deallocate all arrays
that have been allocated by the solver procedure.

0 Several problems can be solved simultaneously, i.e. the package does
not require the solution of one problem to be finished before the
solution of the next starts, as long as for each problem a separate set
of arguments for the above subroutines is used. However, if two or more
problems of the same type need to be solved, it is reasonable to solve
them one after another to reduce memory requirements.

Package types
-------------

``INTEGER`` denotes default ``INTEGER``, and ``REAL`` denotes
``REAL(kind=kind(0d0))``. The term **call type** means
``REAL(kind=kind(0d0))`` for calls to the double precision real
interface, and ``COMPLEX(kind=kind(0d0))`` for calls to the double
precision complex interface.

Derived types
-------------

For each problem, the user must employ the derived types defined by the
module to declare scalars of the types ssmfe\_rcid (real version) or
ssmfe\_rciz (complex version), ssmfe\_expert\_keep, ssmfe\_options and
ssmfe\_inform. The following pseudocode illustrates this.

::

          use SPRAL_SSMFE_EXPERT
          ...
          type (ssmfe_rcid       ) :: rcid
          type (ssmfe_expert_keep) :: keep
          type (ssmfe_options    ) :: options
          type (ssmfe_inform     ) :: inform
          ...

Argument lists
==============

 ``ssmfe_standard()``, ``ssmfe_standard_shift()``, ``ssmfe_generalized()``,
``ssmfe_generalized_shift()``, and ``ssmfe_buckling()`` 
---------------------------------------------------------------------------

**To compute the leftmost eigenpairs of , optionally using
preconditioning, the following call must be made repeatedly: **

To use the solver procedures, the user needs to maintain a workspace W
containing kw + 1 blocks of m vectors of size :math:`n`. A value kw = 7
is always sufficient. However, if options%minAprod :math:`=` .false. and
either options%minBprod :math:`=` .false. or the standard eigenvalue
problem is solved, then kw = 3 is sufficient; if options%minAprod
:math:`=` .true. and either options%minBprod :math:`=` .true. or
ssmfe\_generalized\_shift or ssmfe\_buckling are used, then kw must be
at least 7; otherwise kw = 5 is sufficient. Solver procedures use
indices 1 to m to refer to vectors inside each block and indices 0 to kw
to refer to particular blocks. The first (zero-indexed) block holds the
eigenvector approximations: the user must fill this block with m
linearly independent vectors before the first call to a solver
procedure.

The number of desired eigenpairs may exceed m: whenever converged
eigenpairs have been detected, a solver procedure reports the indices of
these eigenpairs and they must be moved by the user to a separate
eigenvectors’ storage X.

When :math:`B \ne I`, it is expedient to have storage BX for the
:math:`B`-images of the converged eigenvectors, i.e. BX = B\*X.

To simplify the description of the reverse communication interface,
below we assume that an array W(n,m,0:kw) of package type is used as a
workspace, and that arrays X(n, mep) and BX(n, mep) of package type are
used for storing the computed eigenvectors and their :math:`B`-images.
The transpose (real or complex, depending on the package type) of a
matrix H is denoted by H\ :math:`^\prime`.

The meaning of the arguments of the solver procedures is as follows.

is an  scalar of type ssmfe\_rcid in the real version and ssmfe\_rciz in
the complex version. Before the first call, rci%job must be set to 0. No
other values may be assigned to rci%job by the user. After each call,
the value of rci%job must be inspected by the user’s code and the
appropriate action taken:

: fatal error return, the computation must be terminated;

: not all desired eigenpairs converged to required accuracy, see
Section [ssmfe\ :sub:`e`\ xpert:errors];

: the computation is complete and successful.

: (ssmfe\_standard() and ssmfe\_generalized() only) the user must
compute :math:`V = A U`, where

:math:`U=` W(:, ix:jx, rci%kx),  with  ix :math:`=` rci%jx  and  jx
:math:`=` ix + rci%nx - 1,

:math:`V=` W(:, iy:jy, rci%ky),  with  iy :math:`=` rci%jy  and  jy
:math:`=` iy + rci%nx - 1.

: (ssmfe\_standard() and ssmfe\_generalized() only) the user must
compute :math:`V = T U` if preconditioning is used or copy :math:`U` to
:math:`V` otherwise, where :math:`U` and :math:`V` are as for rci%job =
1.

: (ssmfe\_generalized(), ssmfe\_generalized\_shift() and
ssmfe\_buckling() only) the user must compute :math:`V = B U` where
:math:`U` and :math:`V` are as for rci%job = 1.

: the user must save the converged eigenvectors to the eigenvector
storage X and, optionally, for problems and , save their
:math:`B`-images. The converged eigenvectors are columns of the
:math:`{\tt n}\times {\tt m}` matrix W(:,:,rci%kx) and their
:math:`B`-images are respective columns of W(:,:,rci%ky) that are
identified by rci%i, rci%jx and rci%nx as follows: if rci%i > 0, then
the column numbers run from rci%jx to rci%jx + rci%nx - 1, and if rci%i
< 0, then they run from rci%jx - rci%nx + 1 to rci%jx.

: (ssmfe\_standard\_shift(), ssmfe\_generalized\_shift() and
ssmfe\_buckling() only) the user must compute
:math:`V = A_\sigma^{-1} U`, where :math:`A_\sigma = A - \sigma I` and
:math:`I` is :math:`n\times n` identity, for problem ,
:math:`A_\sigma = A - \sigma B` for problem , and
:math:`A_\sigma = B - \sigma A` for problem .

: if rci%i = 0, then the user must perform a copy
:math:`V \leftarrow U`, where :math:`U` and :math:`V` are as for rci%job
= 1, otherwise the columns of W(:,:,rci%kx) and W(:,:,rci%ky) (if rci%kx
:math:`\not=` rci%ky) must be reordered using the index array ind so
that the column ind(j) becomes column j for j = 1, …, rci%nx.

: for each i = 0, 1,..., rci%nx - 1, the user must compute the dot
product of the columns

and

and place it in

.

: if rci%kx :math:`=` rci%ky, then for each i = 0, 1,..., rci%nx - 1,
the user must perform the scaling

,

where :math:`s_i` is the 2-norm of the column W(:, rci%jx + i, rci%kx),
otherwise the user must perform the scalings

,

where :math:`s_i` is the square root of the dot product of the columns
W(:, rci%jx + i, rci%kx) and W(:, rci%jy + i, rci%ky). No scaling is to
be applied to zero columns.

14: for each i = 0, 1,..., rci%nx - 1, the user must perform
axpy-updates:

.

15: the user must perform the matrix multiplication:

.

16: the user must perform the matrix multiplication:

.

17: the user must perform the multiplication:

.

W(:, rci%jy :math:`:` rci%jy + rci%ny - 1, rci%ky) can be used as a
workspace.

21: the user must :math:`B`-orthogonalize the columns of W specified by
rci%nx, rci%jx and rci%kx to all vectors stored in X by solving the
system

for Q and updating

.

For problems and , the respective columns of W(:,:,rci%ky), which store
:math:`B`-images of the respective columns of W(:,:,rci%kx), must be
updated accordingly, either by applying B to these vectors or using the
columns of BX, i.e.

;

22: the user must solve the system

for Q and perform the update

,

where X and BX are same as in the case rci%job = 21 (in the case of
problem , rci%job = 21 and 22 require exactly the same computation).

999: If rci%k > 0, then a restart, normally with a larger block size m,
is suggested with the aim of achieving better convergence. If the
suggestion is accepted, the user must compute the new block size as m =
rci%nx + k + l, where k :math:`\ge` rci%i and l :math:`\ge` rci%j,
reallocate the workspace array W if the new block size is different from
the old one, and set rci%i = 0 and rci%j = 0. If the restart is not
acceptable (e.g. the new block size exceeds a certain limit set by the
user), then nothing needs to be done. If rci%k == 0, then the restart
with the same block size m is required. In both restart cases, the first
block W(:,:,0) of the new workspace should retain the vectors
W(:,i:j,0), where i = rci%jx and j = i + rci%nx - 1, from the old
workspace. The remaining m - rci%nx columns of W(:,:,0) must be filled
with arbitrary vectors that are linearly independent from the converged
eigenvectors and such that the entire set of the columns of W(:,:,0) is
linearly independent.

**Restriction:** rci%job = 0, rci%i = 0 and rci%j = 0 are the only
assignments to the components of rci that can be done by the user. The
first one can only be done before the first call. The other two can only
be done if rci%job = 999 and rci%k > 0.

is an  scalar of type ``REAL`` that holds the shift, a value around
which the desired eigenvalues are situated.

is an  scalar of type default ``INTEGER`` that holds the number of
desired eigenvalues to the left of sigma. **Restriction:** :math:`0 < `
left + right :math:`\le` min(mep, n/2), where right is zero for
ssmfe\_standard() and ssmfe\_generalized().

is an  scalar of type default ``INTEGER`` that holds the number of
desired eigenvalues to the right of sigma. **Restriction:** :math:`0 < `
left + right :math:`\le` min(mep, n/2).

is an  scalar of type default ``INTEGER`` that holds the size of the
array lambda. See Section [ssmfe\ :sub:`e`\ xpert:method] for guidance
on setting mep. **Restriction:** mep is not less than the number of
desired eigenpairs.

is an array of type ``REAL`` and size mep that is used to store the
computed eigenvalues. After a successful completion of the computation
it contains eigenvalues in ascending order. This array must not be
changed by the user.

| is an  scalar of type INTEGER that holds the block size of the user’s
workspace W. **Restriction:**
| 2 :math:`\le` m :math:`<` n.

is an work array of package type, and dimensions 2\*m, 2\*m and 3. It
can only be changed by the user when instructed to do so by rci%job.

is an   array of default integer type, and size at least m. It must not
be changed by the user. It is used for reordering the columns of some
blocks of W.

is an  scalar of type ssmfe\_expert\_keep that holds private data.

is an   scalar of type ssmfe\_options. Its components offer the user a
range of options, see Section [ssmfe\ :sub:`e`\ xpert:type:options]. It
must not be changed by the user between calls.

is an  scalar of type ssmfe\_inform. Its components provide information
about the execution of the subroutine, see
Section [ssmfe\ :sub:`e`\ xpert:type:inform]. It must not be changed by
the user.

``ssmfe_free()``
----------------

**At the end of the computation, the memory allocated by the solver
procedures should be released by making the following subroutine call:
**

``keep``
    is an  scalar of type ssmfe\_expert\_keep, optional. On exit, its
    components that are allocatable arrays will have been deallocated.

``inform``
    is an  scalar of type ssmfe\_inform, optional. On exit, its
    components that are allocatable arrays will have been deallocated.

Derived types
=============

``type(ssmfe_options)``
-----------------------

The derived data type ssmfe\_options has the following components.

``abs_tol_lambda``
    is a scalar of type ``REAL`` that holds an absolute tolerance used
    when testing the estimated eigenvalue error, see
    Section [ssmfe\ :sub:`e`\ xpert:method]. The default value is 0.
    Negative values are treated as the default.

``abs_tol_residual``
    is a scalar of type ``REAL`` that holds an absolute tolerance used
    when testing the residual, see
    Section [ssmfe\ :sub:`e`\ xpert:method]. The default value is 0.
    Negative values are treated as the default.

``max_iterations``
    is a scalar of type default ``INTEGER`` that contains the maximum
    number of iterations to be performed. The default value is 100.
    **Restriction:** max\_it :math:`\ge` 0.

``rel_tol_lambda``
    is a scalar of type ``REAL`` that holds a relative tolerance used
    when testing the estimated eigenvalue error, see
    Section [ssmfe\ :sub:`e`\ xpert:method]. The default value is 0.
    Negative values are treated as the default.

``rel_tol_residual``
    is a scalar of type ``REAL`` that holds a relative tolerance used
    when testing the residual, see
    Section [ssmfe\ :sub:`e`\ xpert:method]. If both abs\_tol\_residual
    and rel\_tol\_residual are set to 0, then the residual norms are not
    taken into consideration by the convergence test, see
    Section [ssmfe\ :sub:`e`\ xpert:method]. The default value is 0.
    Negative values are treated as the default.

``tol_x``
    is a scalar of type ``REAL`` that holds a tolerance used when
    testing the estimated eigenvector error, see
    Section [ssmfe\ :sub:`e`\ xpert:method]. If tol\_x is set to zero,
    the eigenvector error is not estimated. If a negative value is
    assigned, the tolerance is set to 10\*epsilon(lambda). The default
    value is -1.0.

``print_level``
    | is a scalar of type default ``INTEGER`` that determines the amount
    of printing. Possible values are:

    | r@ : p0.85 :math:`<0` & no printing;
    | :math:`0` & error and warning messages only;
    | :math:`1` & the type (standard or generalized) and the size of the
    problem, the number of eigenpairs requested, the error tolerances
    and the size of the subspace are printed before the iterations
    start;
    | :math:`2` & as :math:`1` but, for each eigenpair tested for
    convergence (see Section [ssmfe\ :sub:`e`\ xpert:method]), the
    iteration number, the index of the eigenpair, the eigenvalue,
    whether it has converged, the residual norm, and the error estimates
    are printed;
    | :math:`>2` & as :math:`1` but with all eigenvalues, whether
    converged, residual norms and eigenvalue/eigenvector error estimates
    printed on each iteration.

    The default value is 0. Note that for eigenpairs that are far from
    convergence, ‘rough’ error estimates are printed (the estimates that
    are actually used by the stopping criteria, see
    Section [ssmfe\ :sub:`e`\ xpert:method], only become available on
    the last few iterations).

``unit_error``
    is a scalar of type default ``INTEGER`` that holds the unit number
    for error messages. Printing is suppressed if unit\_error < 0. The
    default value is 6.

``unit_diagnostic``
    is a scalar of type default ``INTEGER`` that holds the unit number
    for messages monitoring the convergence. Printing is suppressed if
    unit\_diagnostics < 0. The default value is 6.

``unit_warning``
    is a scalar of type default ``INTEGER`` that holds the unit number
    for warning messages. Printing is suppressed if unit\_warning < 0.
    The default value is 6.

``err_est``
    is a scalar of type default ``INTEGER`` that defines which error
    estimation scheme for eigenvalues and eigenvectors is to be used by
    the stopping criterion. Two schemes are implemented. If err\_est =
    1, residual error bounds are used, namely, a modified Davis-Kahan
    estimate for the eigenvector error and the Lehmann bounds for the
    eigenvalue error. (see Section [ssmfe\ :sub:`e`\ xpert:errors:est]).
    If err\_est = 2, then the eigenvector and eigenvalue errors are
    estimated by analyzing the convergence curve for the eigenvalues
    (see Section [ssmfe\ :sub:`e`\ xpert:errors:est]). The default is
    err\_est = 2. **Restriction:** err\_est = 1 or 2.

``extra_left``
    is a scalar of type default ``INTEGER`` that holds the number of
    extra approximate eigenvectors corresponding to leftmost eigenvalues
    that are of no interest to the user and are iterated solely to
    enhance convergence. The default is extra\_left = 0.
    **Restriction:** extra\_left :math:`\ge` 0.

``extra_right``
    is a scalar of type default ``INTEGER`` that holds the number of
    extra approximate eigenvectors corresponding to rightmost
    eigenvalues that are of no interest to the user and are iterated
    solely to enhance convergence. The default is extra\_right = 0.
    **Restriction:** extra\_right :math:`\ge` 0.

``left_gap``
    is a scalar of type ``REAL`` that is only used when left is
    non-zero, and specifies the minimal acceptable distance between the
    last computed left eigenvalue and the rest of the spectrum. For
    ssmfe\_standard() and ssmfe\_generalized(), the last computed left
    eigenvalue is the rightmost of the computed ones, and for the other
    procedures it is the leftmost. If set to a negative value
    :math:`\delta`, the minimal distance is taken as :math:`|\delta|`
    times the average distance between the computed eigenvalues. Note
    that for this option to have any effect, the value of mep must be
    larger than left + right: see
    Section [ssmfe\ :sub:`e`\ xpert:method] for further explanation. The
    default value is 0.

``max_left``
    is a scalar of type default ``INTEGER`` that holds the number of
    eigenvalues to the left from :math:`\sigma`, or a negative value, if
    this number is not known (cf.
    Section [ssmfe\ :sub:`e`\ xpert:sec:si]). The default is max\_left =
    -1.

``max_right``
    is a scalar of type default ``INTEGER`` that holds the number of
    eigenvalues to the right from :math:`\sigma`, or a negative value,
    if this number is not known. (cf.
    Section [ssmfe\ :sub:`e`\ xpert:sec:si]). The default is max\_right
    = -1.

``minAprod``
    is a scalar of type default ``LOGICAL`` that determines whether the
    number of multiplications by :math:`A` is to be reduced at the
    expense of memory. If :math:`{\tt minAprod = .false.}`, on each
    iteration three returns to the user with rci%job = 1 are made for
    multiplications of rci%nx vectors by :math:`A`. Otherwise, only one
    such return is made at each iteration but the number kw of blocks in
    the user’s work array W must be increased by 2. The default is
    minAprod = .true.. **Restriction:** minAprod = .true. for
    ssmfe\_standard\_shift(), ssmfe\_generalized\_shift() and
    ssmfe\_buckling().

``minBprod``
    is a scalar of type default ``LOGICAL`` that determines whether the
    number of multiplications by :math:`B` is to be reduced at the
    expense of memory. If :math:`{\tt minBprod = .false.}`, on each
    iteration at least three returns to the user with rci%job = 3 are
    made for multiplications of rci%nx vectors by :math:`B`. Otherwise,
    only one such return is made at each iteration but the number kw of
    blocks in the user’s work array W must be increased by 2. The
    default is minBprod = .true..

``right_gap``
    is a scalar of type ``REAL`` that is only used by ssmfe\_shift,
    ssmfe\_gen\_shift and ssmfe\_buckling with a non-zero right, and has
    the same meaning as options%left\_gap but for the rightmost computed
    eigenvalue. The default value is 0.

``user_x``
    is a scalar of type default ``INTEGER``. If user\_x > 0 then the
    first user\_x columns of x(:,:) will be used as initial guesses for
    eigenvectors. Hence, if the user has good approximations to some of
    the required eigenvectors, the computation time may be reduced by
    putting these approximations into the first user\_x columns of
    x(:,:). The default value is 0, i.e. the columns of x(:,:) are
    overwritten by the solver. **Restriction:** 0 :math:`\le` user\_x
    :math:`\le` m, the first user\_x columns in x(:,:) must be linearly
    independent.

``type(ssmfe_inform)``
----------------------

The derived data type ssmfe\_inform is used to hold information from the
execution of the solver procedures. The components are:

``converged``
    is a rank-1 allocatable array of type default ``INTEGER`` that is
    allocated to have size mep on a call with rci%job = 0 or 999. If, on
    some iteration i, an eigenpair (lambda(j), X(j)) has been accepted
    as converged, then converged(j) = i; if the convergence stagnated
    then converged(j) = -i; otherwise converged(j) = 0.

``err_lambda``
    is a rank-1 allocatable array of type REAL that is allocated to have
    size mep on a call with rci%job = 0 or 999. err\_lmd(i) contains the
    estimated eigenvalue error for the approximate eigenvalue lambda(i)
    if info%converged(i) is non-zero, and -1.0 otherwise.

``err_x``
    is a rank-1 allocatable array of type REAL. This array is allocated
    to have size mep on a call with rci%job = 0 or 999, and is used for
    storing the eigenvector errors in the same way as err\_lmd is used
    for storing the eigenvalue errors.

``flag``
    is a scalar of type default ``INTEGER`` that is used as an error
    flag. If a call is successful, flag has value 0. A nonzero value of
    flag indicates an error or a warning (see
    Section [ssmfe\ :sub:`e`\ xpert:errors]).

``iteration``
    is a scalar of type default ``INTEGER`` that holds the number of
    iterations since the previous rci%job = 0 or rci%job = 999 call.

``left``
    is a scalar of type default ``INTEGER`` that holds the number of
    converged eigenvalues on the left, i.e. the total number of
    converged eigenpairs of or the number of the converged eigenvalues
    of or to the left of sigma.

``next_left``
    is a scalar of type default REAL that holds the non-converged
    eigenvalue next to the last converged on the left (cf.
    options%left\_gap).

``next_right``
    is a scalar of type default REAL that holds the non-converged
    eigenvalue next to the last converged on the right (cf.
    options%right\_gap).

``non_converged``
    is a scalar of type default ``INTEGER`` that holds the number of
    non-converged eigenpairs (see
    Section [ssmfe\ :sub:`e`\ xpert:errors]).

``residual_norms``
    is a rank-1 allocatable array of type default REAL that is allocated
    to have size mep on a call with rci%job = 0 or 999. On returns with
    rci%job = 5, residual\_norms(i) contains the Euclidean norm of the
    residual for lambda(i), X(i).

``right``
    is a scalar of type default ``INTEGER`` that holds the number of
    converged eigenvalues of or to the right of sigma.

``stat``
    is a scalar of type default ``INTEGER`` that holds the allocation
    status (see Section [ssmfe\ :sub:`e`\ xpert:errors]).

Error flags
===========

A successful return from a solver procedure is indicated by
inform%flag\ :math:`=`\ 0. A negative value indicates an error, a
positive value indicates a warning.

Possible negative values of inform%flag are as follows:

  -1 Incorrect value of rci%job.

  -2 Block size m is out-of-range.

  -3 Incorrect value of options%err\_est.

  -4 Incorrect value of options%minAprod.

  -5 Incorrect value of options%extra\_left or options%extra\_right.

  -6 Incorrect value of options%min\_gap.

 -11 Incorrect value of left.

 -12 Incorrect value of right.

 -13 mep is less than the number of desired eigenpairs.

-100 Not enough memory; inform%stat contains the value of the Fortran
stat parameter.

-200 :math:`B` is not positive definite or initial eigenvectors are
linearly dependent.

Possible positive values are:

1 The iterations have been terminated because no further improvement in
accuracy is possible (this may happen if the preconditioner is not
positive definite, or if the components of the residual vectors are so
small that the round-off errors make them essentially random). The value
of inform%non\_converged is set to the number of non-converged
eigenpairs.

2 The maximal number of iterations has been exceeded. The value of
inform%non\_converged is set to the number of non-converged eigenpairs.

3 Out of storage space for the converged eigenpairs. The value of
inform%non\_converged is set to the number of non-converged eigenpairs.

Method
======

The algorithm
-------------

The solver procedures of ``SPRAL_SSMFE_EXPERT`` are interfaces to solver
procedures of ``SPRAL_SSMFE_CORE``, which implement a block iterative
algorithm based on the Jacobi-conjugate preconditioned gradients method
[2,3]. Further information on the algorithm used by
``SPRAL_SSMFE_EXPERT`` can be found in the specification document for
``SPRAL_SSMFE_CORE`` and in [1].

Stopping criteria
-----------------

An approximate eigenpair :math:`\{x,\lambda\}` is considered to have
converged if the following three conditions are all satisfied:

#. if options%abs\_tol\_lambda and options%rel\_tol\_lambda are not both
   equal to zero, then the estimated error in the approximate eigenvalue
   must be less than or equal to

   max(options%abs\_tol\_lambda,
   :math:`\delta`\ \*options%rel\_tol\_lambda),

   where :math:`\delta` is the estimated average distance between
   eigenvalues.

#. if options%tol\_x is not zero, then the estimated sine of the angle
   between the approximate eigenvector and the invariant subspace
   corresponding to the eigenvalue approximated by :math:`\lambda` must
   be less than or equal to options%tol\_x.

#. if options%abs\_tol\_residual and options%rel\_tol\_residual are not
   both equal to zero, then the Euclidean norm of the residual,
   :math:`\|A x - \lambda B x\|_2`, must be less than or equal to

   max(options%abs\_tol\_residual,
   options%rel\_tol\_residual\*\ :math:`\|\lambda B x\|_2`).

The extra eigenpairs are not checked for convergence, as their role is
purely auxiliary.

Improving eigenvector accuracy
------------------------------

If the gap between the last computed eigenvalue and the rest of the
spectrum is small, then the accuracy of the corresponding eigenvector
may be very low. To prevent this from happening, the user should set the
eigenpairs storage size mep to a value that is larger than the number of
desired eigenpairs, and set the options options%left\_gap and
options%right\_gap to non-zero values :math:`\delta_l` and
:math:`\delta_r`. These values determine the size of the minimal
acceptable gaps between the computed eigenvalues and the rest of the
spectrum, :math:`\delta_l` referring to either leftmost eigenvalues (for
ssmfe\_standard() and ssmfe\_generalized() only) or those to the left of
the shift sigma, and :math:`\delta_r` to those to the right of the shift
sigma. Positive values of :math:`\delta_l` and :math:`\delta_r` set the
gap explicitly, and negative values require the gap to be not less than
their absolute value times the average distance between the computed
eigenvalues. A recommended value of :math:`\delta_l` and
:math:`\delta_r` is :math:`-0.1`. The value of mep has little effect on
the speed of computation, hence it might be set to any reasonably large
value. The larger the value of mep, the larger the size of an eigenvalue
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
the LDLT factorization of the shifted system, e.g. HSL\_MA97 or
SPRAL\_SSIDS. The LDLT factorization of the matrix :math:`A - \sigma B`
consists in finding a unit lower triangular matrix :math:`L`, a
block-diagonal matrix :math:`D` with :math:`1\times 1` and
:math:`2\times 2` blocks on the main diagonal and a permutation matrix
:math:`P` such that :math:`P^T(A - \sigma B)P = L D L^T`. By inertia
theorem, the number of eigenvalues to the left and right from the shift
:math:`\sigma` is equal to the number of negative and positive
eigenvalues of :math:`D`, which allows quick computation of the
eigenvalue numbers each side of the shift.

Error estimation
----------------

Standard problem
~~~~~~~~~~~~~~~~

If options%err\_est = 1, the error estimates for the eigenvalues are
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

If options%err\_est\ :math:`=`\ 2 the errors are estimated based on the
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
options%err\_est\ :math:`=`\ 1 is not recommended. In the case of
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

[1] E. E. Ovtchinnikov and J. Reid. A preconditioned block conjugate
gradient algorithm for computing extreme eigenpairs of symmetric and
Hermitian problems. Technical Report RAL-TR-2010-19, 2010.

E. E. Ovtchinnikov, *Jacobi correction equation, line search and
conjugate gradients in Hermitian eigenvalue computation I: Computing an
extreme eigenvalue*, SIAM J. Numer. Anal., **46**:2567–2592, 2008.

E. E. Ovtchinnikov, *Jacobi correction equation, line search and
conjugate gradients in Hermitian eigenvalue computation II: Computing
several extreme eigenvalues*, SIAM J. Numer. Anal., **46**:2593–2619,
2008.

Examples
========

Preconditioning example
-----------------------

The following code computes the 5 leftmost eigenpairs of the matrix
:math:`A` of order 100 that approximates the two-dimensional Laplacian
operator on a 20-by-20 grid. One forward and one backward Gauss-Seidel
update are used for preconditioning, which halves the number of
iterations compared with solving the same problem without
preconditioning. The module laplace2d
(``examples/Fortran/ssmfe/laplace2d.f90``) supplies the subroutine
apply\_laplacian() that multiplies a block of vectors by :math:`A`, and
the subroutine apply\_gauss\_seidel\_step() that computes
:math:`y = T x` for a block of vectors :math:`x` by applying one forward
and one backward update of the Gauss-Seidel method to the system
:math:`A y = x`. This code produces the following output:

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
