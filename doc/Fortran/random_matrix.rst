2014-03-06 Version 1.0.0
    Initial release

Installation
============

Please see the SPRAL install documentation.

Usage overview
==============

Calling sequences
-----------------

Access to the package requires a USE statement:

::

       use spral_random_matrix

The following procedure is available to the user:

-  random\_matrix\_generate() generates a random matrix to the supplied
   specification.

Derived types
-------------

The user must supply a random number generator state using the type
defined in the ``spral_random`` module. The following pseudo-code
illustrates how such a generator may be defined.

::

          use spral_random, only : random_state
          ...
          type (random_state) :: state
          ...

Further details are given in the documentation for the ``spral_random``
module. The user may wish to use the routines ``get_random_seed()`` and
``set_random_seed()`` to control the matrix generated.

Optional arguments
------------------

We use square brackets to indicate *optional* arguments. In each call,
optional arguments follow the argument info. Since we reserve the right
to add additional optional arguments in future releases of the code,
**we strongly recommend that all optional arguments be called by
keyword, not by position**.

Subroutines
===========

``random_matrix_generate()``
----------------------------

If ``matrix_type`` specifies a symmetric or skew symmetric matrix, only
the lower half matrix will be returned to the user.

``state``
    is an  scalar of type random\_state. It is used as the state of the
    pseudo-random number generator used.

``matrix_type``
    is an  scalar of type default INTEGER. It specifies the matrix type
    to be generated, and must be one of:

    -  Undefined (matrix generated will be unsymmetric/rectangular)

    -  Rectangular matrix (:math:`\texttt{m}\ne\texttt{n}`)

    -  Unsymmetric (:math:`\texttt{m}=\texttt{n}`)

    -  Symmetric positive-definite (:math:`\texttt{m}=\texttt{n}`,
       ``nnz`` entries in lower triangle returned, matrix made
       diagonally dominant if ``val`` is present, non-singularity
       assured regardless of the value of the optional argument
       nonsingular)

    -  Symmetric indefinite (:math:`\texttt{m}=\texttt{n}`, ``nnz``
       entries in lower triangle returned)

    -  Skew symmetric (:math:`\texttt{m}=\texttt{n}`, ``nnz`` entries in
       lower triangle returned)

``m``
    is an  scalar of type default INTEGER that specifies the number of
    rows in the matrix. **Restriction:** m\ :math:`\geq`\ 1.

``n``
    is an  scalar of type default INTEGER that specifies the number of
    columns in the matrix. **Restriction:** n\ :math:`\geq`\ 1, and
    consistent with ``matrix_type``.

``nnz``
    is an  scalar of type default INTEGER that specifies the number of
    non-zeroes in the matrix. **Restriction:** nnz\ :math:`\geq`\ 1 (and
    ``nnz``\ :math:`\geq\min(\texttt{m},\texttt{n})` if non-singularity
    requested, or positive-definite).

``ptr(:)``
    is an  array of type default INTEGER and size n+1. On exit, ptr(j)
    specifies the position in row(:) of the first entry in column j and
    ptr(n+1)\ :math:`=`\ nnz+1.

``row(:)``
    is an  array of type default INTEGER and size nnz. On exit, row(j)
    specifies the row to which the j-th entry belongs.

``flag``
    is an  scalar of type default INTEGER. On exit, it specifies a
    return code that indicates success with the value 0 or failure with
    a negative value detailed in
    Section [random\ :sub:`m`\ atrix:errors].

``stat``
    is an optional  scalar of type default INTEGER. If present, on exit
    it returns the value of the Fortran ``stat`` parameter on the last
    ``allocate()`` call. In particular, if ``flag`` indicates an
    allocation failure, it returns further information on the failure.

``val(:)``
    is an optional  array of type REAL(wp) and size ``nnz``. On exit,
    ``val(j)`` gives the value of the ``j``-th entry. Entries are
    generated from a uniform distribution on the interval
    :math:`[-1,1]`. In the positive-definite case only, diagonal entries
    are given a value equal to the sum of the off-diagonal entries in
    the row plus a value chosen uniformally at random from the interval
    :math:`(0,1]`.

``nonsingular``
    is an optional  scalar of type LOGICAL. If present with the value
    .true., the generated matrix is guaranteed to have a transversal of
    size :math:`\min({\tt m}, {\tt n})`. In the symmetric or skew
    symmetric case this will be the natural diagonal. In the unsymmetric
    and rectangular cases a random matching is used. In the symmetric
    positive-definite case, this value is ignored (it is treated as
    .true.). In all other cases, if nonsingular is not present, or is
    present with the value .false., a maximum transversal is not
    guaranteed and the generated matrix may be structurally rank
    deficient.

``sort``
    is an optional  scalar of type LOGICAL. If present with the value
    .true., the row entries of the generated matrix will be sorted into
    ascending order within each column. Otherwise, if sort is not
    present, or is present with the value .false., entries may be
    returned in a random order within each column.

Return codes
============

A successful return is indicated by flag having the value zero. A
negative value is associated with an error message.

Possible negative values are:

:math:`-`\ 1 An allocation error has occurred. If present, stat will
contain the Fortran stat value returned by the failed allocate() call.

:math:`-`\ 2 An invalid value of matrix\_type was supplied.

:math:`-`\ 3 At least one of m, n, or nnz was less than :math:`1`.

:math:`-`\ 4 The (in)equality of m and n was inconsistent with
matrix\_type.

:math:`-`\ 5 A non-singular matrix was requested, but
:math:`\texttt{nnz}<\min(\texttt{m},\texttt{n})`.

Method
======

If structural non-singularity is requested, first
:math:`\min({\tt m}, {\tt n})` entries are generated as follows:

Unsymmetric or Rectangular
    Random permutations of the rows and columns are generated. The first
    :math:`\min({\tt m}, {\tt n})` entries of these permutations are
    used to specify the entries of a maximum transversal.

Symmetric
    The diagonal is added to the matrix explicitly.

The remaining non-zero entries are then assigned to columns uniformally
at random. In the symmetric case, a weighting is used in proportion to
the number of entries below the diagonal. If the selected column for a
given non-zero is already full, a new random sample is drawn.

Once the number of entries in each column has been determined, and any
required maximum transversal inserted, row indices are determined
uniformally at random. Should a non-zero in that row already be present
in the column, a new random sample is drawn.

In all cases, values are drawn uniformally at random from the range
:math:`(-1,1)`. In the positive-definite case, a post-processing step
sums the absolute values of all the entries in each column and replaces
the diagonal with this value.

Example
=======

The following code generates a random :math:`4 \times 5` matrix with
:math:`8` non-zeroes that is non-singular. This produces the following
output:

::

    Generating a   4 x  5 non-singular matrix with   8 non-zeroes
    Generated matrix:
    Matrix of undefined type, dimension 4x5 with 8 entries.
    1:                                         -1.0744E-01   9.1000E-01
    2:                             9.5364E-01                1.0912E-01
    3:                                          1.1631E-01  -5.8957E-01
    4:  -9.0631E-01                                          7.7313E-01

