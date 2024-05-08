******************************************************
:f:mod:`spral_lsmr` - Sparse Least Squares LSMR Solver
******************************************************
.. f:module:: spral_lsmr
   :synopsis: Sparse Least Squares LSMR Solver

=======
Purpose
=======

This package uses the LSMR iterative method to solve sparse linear
equations and sparse least-squares problems of the form:

.. math::

   \begin{array}{ll}
      \mbox{1. Nonsymmetric equations:} &
         \mbox{minimize } \|x\|_2 \mbox{ subject to }Ax = b, \\
      \mbox{2. Linear least squares:} & \mbox{minimize  } \|Ax - b\|_2^2,\\
      \mbox{3. Regularized least squares:} &
         \mbox{minimize  } \|Ax - b\|_2^2 + \lambda^2\|x\|_2^2,
   \end{array}

where the :math:`m \times n` matrix :math:`A` may be square or rectangular, and
may have any rank. The scalar :math:`\lambda` is a damping parameter. If
:math:`\lambda > 0`, the solution is regularized in the sense that a unique
soluton always exists, and :math:`\|x\|_2` is always bounded.

Preconditioning may be used to try to reduce the number of iterations.
A suitable choice for the preconditioner depends on the user's knowledge of
:math:`A`. For a user-chosen :math:`n \times n` nonsingular matrix :math:`P`,
LSMR solves

.. math::

   \begin{array}{ll}
      \mbox{1. Nonsymmetric equations:} &
         \mbox{minimize  } \|Py\|_2 \mbox{  subject to }APy = b, \\
      \mbox{2. Linear least squares:} & \mbox{minimize  } \|APy - b\|_2^2 ,\\
      \mbox{3. Regularized least squares:} &
         \mbox{minimize  } \|APy - b\|_2^2 + \lambda^2\|Py\|_2^2 ,
   \end{array}

The user must then recover the final solution :math:`x` by computing
:math:`x=Py`. :math:`P` will be a good preconditioner if :math:`AP` is
significantly better conditioned than :math:`A`.

Reverse communication is used for preconditioning operations :math:`Pz` and
:math:`P^Tz` and matrix-vector products of the form :math:`Av` and :math:`A^Tu`.

The method used is based on the Golub-Kahan bidiagonalization process.
It is algebraically equivalent to applying MINRES to the normal
equation :math:`(A^TA+\lambda^2I)x=A^Tb` (or
:math:`((AP)^T(AP)+\lambda^2I)y=(AP)^Tb`), but has better numerical properties,
especially if :math:`A` is ill-conditioned.

Notation
--------

In the rest of this documentation, we use the following notation:

.. math::

   \bar{A} = \left( \begin{array}{c}
            A \\
            \lambda I
         \end{array} \right)P, \hspace{1cm} \bar{b} = \left( \begin{array}{c}
            b \\
            0
         \end{array} \right),

and

.. math::

   r = b - APy , \hspace{1cm}  \bar{r} = \bar{b} - \bar{A}y,
         \hspace{1cm}  x = Py.


Version history
---------------

2016-05-10 Version 1.0.0.
    Initial release.

[For detailed history, see ChangeLog]

===========
Subroutines
===========

.. f:subroutine:: lsmr_solve(action, m, n, u, v, y, keep, options, inform[,damp])

   Solve the least squares problem. Uses reverse communication. Upon return the
   user must perform a task specified by the `action` parameter and recall the
   routine. Possible values of action and associated tasks are:

   +----------+----------------------------------------------------------------+
   | `action` | Task to be performed                                           |
   +==========+================================================================+
   | 0        | None. Computation has terminated. Check inform%flag for reason.|
   +----------+----------------------------------------------------------------+
   | 1        | The user must compute                                          |
   |          |                                                                |
   |          | .. math:: v = v + P^TA^Tu                                      |
   |          |                                                                |
   |          | without altering :math:`u`.                                    |
   +----------+----------------------------------------------------------------+
   | 2        | The user must compute                                          |
   |          |                                                                |
   |          | .. math:: u = u + APv                                          |
   |          |                                                                |
   |          | without altering :math:`v`.                                    |
   +----------+----------------------------------------------------------------+
   | 3        | The user may test for convergence.                             |
   |          |                                                                |
   |          | To continue, recall without changing any of the arguments.     |
   +----------+----------------------------------------------------------------+


   :p integer action [inout]: reverse communication task (see table above).
      Must be 0 on first call.
   :p integer m [in]: number of rows in :math:`A`.
   :p integer n [in]: number of columns in :math:`A`.
   :p real u (m) [inout]: the vector :math:`u`. Must contain :math:`b` on first
      call.
   :p real v (n) [inout]: the vector :math:`v`.
   :p real y (n) [inout]: the current solution :math:`y` of the preconditioned
      problem (note :math:`x=Py`).
   :p lsmr_keep keep [inout]: private internal data for LSMR.
   :p lsmr_options options [in]: controls execution of algorithm.
   :p lsmr_inform inform [inout]: information about execution of algorithm.
      Must not be changed by the user.
   :o real damp [in]: Damping parameter :math:`\lambda`.

.. f:subroutine:: lsmr_free(keep,stat)

   Free memory allocated in `keep` (unnecessary if `keep` is going out of
   scope).

   :p lsmr_keep [inout]: private data to be freed.
   :p integer stat [out]: return 0 on success, or Fortran stat parameter on
      failed deallocation.

   .. note::

      If a series of problems is being solved sequentially, the same keep may be
      used without calling :f:subr:`lsmr_free()` between each solution.

=============
Derived types
=============

.. f:type:: lsmr_options

   Specify options used by LSMR algorithm.

   :f integer print_freq_head [default=20]: frequency of printing heading
      information (that is, how many lines are printed before the heading
      information is reprinted).
   :f integer print_freq_itn [default=10]: frequency of printing status.
      There is printing on each of the first `print_freq_itn` iterations and
      then printing every `print_freq_itn` iterations.
   :f integer unit_diagnostics [default=6]: Fortran unit for diagnostic
      printing. Printing is suppressed if negative.
   :f integer unit_error [default=6]: Fortran unit for printing error messages.
      Printing is suppressed if negative.
   :f real atol [default=sqrt(epsilon(1.0d0))]: Relative error in :math:`A`.
      i.e. if :math:`A` is accurate to about 6 digits, set atol to 1.0d-6.
      Only used if `options%ctest=3`.
   :f real btol [default=sqrt(epsilon(1.0d0))]: Relative error in :math:`b`.
      i.e. if :math:`b` is accurate to about 6 digits, set btol to 1.0d-6.
      Only used if `options%ctest=3`.
   :f real conlim [default=1/(10*sqrt(epsilon(1.0d0))]: Upper limit on
      :math:`cond(\bar{A})`, apparent condition number of :math:`\bar{A}`.
      Only used if `options%ctest=3`.
   :f integer ctest [default=3]: Convergence test to use, one of:

      +-------------+----------------------------------------------------------+
      | 1           | User to test convergence (`action=3`).                   |
      |             | *Without*  computation of norms in inform.               |
      +-------------+----------------------------------------------------------+
      | 2           | User to test convergence (`action=3`).                   |
      |             | *With*  computation of norms in inform.                  |
      +-------------+----------------------------------------------------------+
      | 3 (default) | Convergence test of Fond and Saunders (in preconditioner |
      |             | norm).                                                   |
      +-------------+----------------------------------------------------------+

   :f integer itnlim [default=-1]: Maximum number of iterations. If negative,
      the limit used is :math:`4n`.
   :f integer itn_test [default=-1]: Number of iterations between user
      convergence tests. If negative, use :math:`min(n,10)`.
   :f integer localSize [default=0]: Number of historical vectors to use for
      reorthogonalization.

.. f:type:: lsmr_inform

   Information about progress of algorithm.

   :f real condAP: Estimate of :math:`cond(\bar{A})`. A very high value of
      condAP may again indicate an error in the products with :math:`A`,
      :math:`A^T`, :math:`P`, or :math:`P^T`. A negative value indicates that
      no estimate is currently available.
      This component is not used if `options%ctest=1`.
   :f integer flag: exit status of algorithm. See table below.
   :f integer itn: number of iterations performed
   :f real normAP: Estimate of Frobenius norm of :math:`\bar{A}`. If
      :math:`\lambda` is small and the columns of :math:`AP` have all been
      scaled to have length 1.0, normAP should increase to roughly
      :math:`\sqrt{n}`. A radically different value for normAP may indicate an
      error in the user-supplied products with :math:`A`, :math:`A^T`,
      :math:`P`, or :math:`P^T`. A negative value indicates that no estimate is
      currently available.
      This component is not used if `options%ctest=1`.
   :f real normAP: Estimate of :math:`\|\bar{A}^T\bar{r}\|_2`, (normal
      equations residual). This should be small in all cases. Note that
      normAPr will often be smaller than the true value.
      A negative value indicates that no estimate is currently available.
      This component is not used if `options%ctest=1`.
   :f real normr: Estimate of :math:`\|\bar{r}\|_2`. This will be small
      if :math:`Ax = b` has a solution. A negative value indicates that no
      estimate is currently available.
      This component is not used if `options%ctest=1`.
   :f real normy: Estimate of :math:`\|y\|_2`. A negative value indicates that
      no estimate is currently available.
      This component is not used if `options%ctest=1`.
   :f integer stat: The Fortran stat parameter in the event of a failed
      allocation.

   +-------------+-------------------------------------------------------------+
   | inform%flag | Interpretation                                              |
   +=============+=============================================================+
   |  0          | :math:`x = 0.0` is the exact solution.                      |
   |             | No iterations were performed.                               |
   +-------------+-------------------------------------------------------------+
   |  1          | The equations :math:`Ax = b` are probably compatible.       |
   |             | :math:`\|Ax - b\|_2` is sufficiently small, given the values|
   |             | of `options%atol` and `options%btol`.                       |
   |             | (`options%ctest=3` only).                                   |
   +-------------+-------------------------------------------------------------+
   |  2          | If damp is not present or is zero then the system           |
   |             | :math:`Ax = b` is probably not compatible. A least-squares  |
   |             | solution has been obtained that is sufficiently accurate,   |
   |             | given the value of `options%atol`. Otherwise, damped        |
   |             | least-squares solution has been obtained that is            |
   |             | sufficiently accurate, given the value of `options%atol`.   |
   |             | (`options%ctest=3` only).                                   |
   +-------------+-------------------------------------------------------------+
   |  3          | An estimate of cond(:math:`\bar{A}`) has exceeded           |
   |             | `options%conlim`. The system :math:`Ax = b` appears to be   |
   |             | ill-conditioned, or there could be an error in the products |
   |             | with :math:`A`, :math:`A^T`, :math:`P`, or :math:`P^T`.     |
   |             | (`options%ctest=3` only).                                   |
   +-------------+-------------------------------------------------------------+
   |  4          | :math:`\|APy - b \|_2` is small enough for this machine.    |
   |             | (`options%ctest=3` only).                                   |
   +-------------+-------------------------------------------------------------+
   |  5          | The least-squares solution is good enough for this machine. |
   |             | (options%ctest=3 only).                                     |
   +-------------+-------------------------------------------------------------+
   |  6          | The estimate `inform%condAP` appears to be too large for    |
   |             | this machine. (`options%ctest=3` only).                     |
   +-------------+-------------------------------------------------------------+
   |  7          | The iteration limit `options%itnlim` has been reached.      |
   +-------------+-------------------------------------------------------------+
   |  8          | An array allocation failed.                                 |
   +-------------+-------------------------------------------------------------+
   |  9          | An array deallocation failed.                               |
   +-------------+-------------------------------------------------------------+
   | 10          | Either `m<0` or `n<0`.                                      |
   +-------------+-------------------------------------------------------------+

=======
Example
=======

The following code illustrates the use of LMSR

.. literalinclude:: ../../examples/Fortran/lsmr.f90
   :language: Fortran

This returns the following output:

::

     Exit LSMR with inform%flag =   2 and inform%itn =   3
     LS solution is:
          0.18      0.26      0.17

======
Method
======

Algorithm
---------

The method used is based on the Golub-Kahan bidiagonalization process.
It is algebraically equivalent to applying MINRES to the normal equation
:math:`(A^TA+\lambda^2I)x=A^Tb` (or :math:`((AP)^T(AP)+\lambda^2I)y=(AP)^Tb`,
:math:`Py = x`, if preconditioning is used), but has better numerical
properties, especially if :math:`A` is ill-conditioned. Full details may be
found in [1]_.

Scaling
-------

LSMR uses an iterative method to approximate the solution. The number of
iterations required to reach a certain accuracy depends strongly on the
scaling of the problem. Poor scaling of the rows or columns of
:math:`A` should therefore be avoided where possible. For example, in
problem 1 the solution is unaltered by row-scaling. If a row of
:math:`A` is very small or large compared to the other rows of
:math:`A`, the corresponding row of :math:`( A\;  b )` should be scaled
up or down.

In problems 1 and 2, the solution :math:`x` is easily recovered
following column-scaling. Unless better information is known, the
nonzero columns of :math:`A` should be scaled so that they all have the
same Euclidean norm (e.g., 1.0). In problem 3, there is no freedom to
re-scale if `damp` is nonzero. However, the value of `damp` should be
assigned only after attention has been paid to the scaling of :math:`A`.

The parameter `damp` is intended to help regularize ill-conditioned
systems, by preventing the true solution from being very large. Another
aid to regularization is provided by the `inform%condAP`, which may be
used to terminate iterations before the computed solution becomes very
large.

Initial estimate
----------------

Note that :math:`x` (or :math:`y` for the preconditioned problem) is not
an input parameter. If some initial estimate :math:`x_0` of :math:`x` is
known and if :math:`\lambda = 0`, one could proceed as follows:

1. Compute a residual vector :math:`r_0 = b - Ax_0`.

2. Use LSMR to solve the system :math:`A \delta x = r_0`.

3. Add the correction :math:`\delta x` to obtain a final solution
   :math:`x = x_0 + \delta x`.

This can be generalized for the preconditioned case. The guess :math:`x_0` has
to be available before and after the calls to :f:subr:`lsmr_solve()`. To judge
the benefits, suppose :f:subr:`lsmr_solve()` takes :math:`k_1` iterations to
solve :math:`Ax = b` and :math:`k_2` iterations to solve
:math:`A \delta x = r_0`. If :math:`x_0` is "good", :math:`\|r_0\|_2` will be
smaller than :math:`\|b\|_2`. If the same stopping tolerances `options%atol`
and `options%btol` are used for each system, :math:`k_1` and :math:`k_2` will
be similar, but the final solution :math:`x = x_0 + \delta x` should be more
accurate. The only way to reduce the total work is to use a larger stopping
tolerance for the second system. If some value `options%btol` is suitable for
:math:`Ax=b`, the larger value
:math:`\mathrm{options\%btol}*\|b\|_2 / \|r_0\|_2` should be suitable for
:math:`A \delta x = r_0`.

References
----------

.. [1] D.C.-L. Fong and M.A. Saunders (2011).
   *LSMR: An iterative algorithm for sparse least-squares problems*.
   SIAM J. Sci. Comput. 33:5, 2950-2971.
   [`DOI: 10.1137/10079687X <https://doi.org/10.1137/10079687X>`_]
