*************************************************************
:f:mod:`spral_random_matrix` - Pseudo-random Matrix Generator
*************************************************************
.. f:module:: spral_random_matrix
   :synopsis: Pseudo-random Matrix Generator

=======
Purpose
=======

This package generates a random sparse matrix of specified size and density in
compressed sparse column format. Either the pattern or both the pattern and
values can be generated. Both symmetric and unsymmetric matrices can be
generated, and structural non-degeneracy can optionally be ensured, and the
row indices can be sorted within columns.

Version history
---------------

2016-09-08 Version 1.1.0
   Add long support

2014-03-06 Version 1.0.0
    Initial release

========
Routines
========

.. f:function:: random_matrix_generate(state,matrix_type,m,n,nnz,ptr,row,flag[,stat,val,nonsingular,sort])

   Generate an :math:`m\times n` random matrix with :math:`nnz` non-zero
   entries.

   If `matrix_type` specifies a symmetric or skew symmetric matrix, only
   the lower half matrix will be returned to the user.

   :p random_state state [inout]: State of the pseudo-random number generator
      to use.
   :p integer matrix_type [in]: Type of matrix to generate. One of:

      +---+-------------------------------------------------------------------+
      | 0 | Undefined (matrix generated will be unsymmetric/retangular).      |
      +---+-------------------------------------------------------------------+
      | 1 | Rectangular matrix.                                               |
      +---+-------------------------------------------------------------------+
      | 2 | Unsymmetric.                                                      |
      +---+-------------------------------------------------------------------+
      | 3 | Symmetric positive-definite.                                      |
      |   |                                                                   |
      |   | `nnz` entries in lower triangle returned, matrix made diagonally  |
      |   | dominant if ``val`` is present, non-singularity assured regardless|
      |   | of argument `nonsingular`.                                        |
      +---+-------------------------------------------------------------------+
      | 4 | Symmetric indefinite.                                             |
      |   |                                                                   |
      |   | `nnz` entries in lower triangle returned.                         |
      +---+-------------------------------------------------------------------+
      | 5 | Skew symmetric.                                                   |
      |   |                                                                   |
      |   | `nnz` entries in lower triangle returned.                         |
      +---+-------------------------------------------------------------------+

   :p integer m [in]: Number of rows in the matrix.
   :p integer n [in]: Number of columns in the matrix.
   :p integer(long) nnz [in]: Number of non-zeroes in the matrix.
   :p integer(long) ptr (n+1) [out]: Column pointers of the matrix
      (see :doc:`CSC format<csc format>`).
   :p integer row (nnz) [out]: Row indices of the matrix
      (see :doc:`CSC format<csc_format>`).
   :p integer flag [out]: Exit status of the algorithm, 0 on success, otherwise
      an error occured, see table below.
   :o integer stat [out]: Stat parameter of last ``allocate()`` call.
   :o integer val (nnz) [out]: non-zero values of the matrix
      (see :doc:`CSC format<csc_format>`).
   :o logical nonsingular [in]: Ensure matrix is non-singular if present with
      value ``.true.``. Such a matrix is guaranteed to have a transversal of
      size :math:`\min({\tt m}, {\tt n})`. In the symmetric or skew
      symmetric case this will be the natural diagonal. In the unsymmetric
      and rectangular cases a random matching is used. Otherwise the matrix
      may or may not be structurally singular.
      Note that symmetric positive-definite matrices are always non-singular.
   :o logical sort [in]: Sort entries in each column into ascending order if
      present with value ``.true.``.
      Otherwise entries may be in any order within a column.

   Possible exit status values are:

   +--------+-----------------------------------------------------------------+
   | `flag` | Status                                                          |
   +========+=================================================================+
   | -1     | An allocation error has occurred. If present, `stat` will       |
   |        | contain the Fortran stat value returned by the failed           |
   |        | allocate() call.                                                |
   +--------+-----------------------------------------------------------------+
   | -2     | An invalid value of matrix_type was supplied.                   |
   +--------+-----------------------------------------------------------------+
   | -3     | At least one of m, n, or nnz was less than 1.                   |
   +--------+-----------------------------------------------------------------+
   | -4     | The (in)equality of m and n was inconsistent with matrix_type.  |
   +--------+-----------------------------------------------------------------+
   | -5     |  A non-singular matrix was requested, but :math:`nnz<\min(m,n)`.|
   +--------+-----------------------------------------------------------------+

   .. note::

      A version is also provided in which ``nnz`` and ``ptr`` have type default
      integer, however users are encouraged to use 64-bit integers to ensure
      code can handle large matrices.

=======
Example
=======

The following code generates a random :math:`4 \times 5` matrix with
:math:`8` non-zeroes that is non-singular.

.. literalinclude:: ../../examples/Fortran/random_matrix.f90
   :language: Fortran

This produces the following output:

::

    Generating a   4 x  5 non-singular matrix with   8 non-zeroes
    Generated matrix:
    Matrix of undefined type, dimension 4x5 with 8 entries.
    1:                                         -1.0744E-01   9.1000E-01
    2:                             9.5364E-01                1.0912E-01
    3:                                          1.1631E-01  -5.8957E-01
    4:  -9.0631E-01                                          7.7313E-01

======
Method
======

If structural non-singularity is requested, first
:math:`\min(m, n)` entries are generated as follows:

Unsymmetric or Rectangular
    Random permutations of the rows and columns are generated. The first
    :math:`\min(m, n)` entries of these permutations are
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
