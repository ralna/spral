**********************************************
RANDOM_MATRIX - Pseudo-random Matrix Generator
**********************************************

.. code-block:: C

   #include <spral_random_matrix.h> /* or <spral.h> for all packages */

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

===================
Seed Initialization
===================

Prior to first use, the random number generator state must be initialised:

.. code-block:: C

   int state = SPRAL_RANDOM_INITIAL_SEED;

See :doc:<random> documentation for further information.

========
Routines
========

.. c:function:: int spral_random_matrix_generate(int *state, enum spral_matrix_type matrix_type, int m, int n, int nnz, int ptr[n+1], int row[nnz], double *val, int flags)

   Generate an :math:`m\times n` random matrix with `math:`nnz` non-zero
   entries.

   If `matrix_type` specifies a symmetric or skew symmetric matrix, only
   the lower half matrix will be returned to the user.

   :param state: State of the pseudo-random number generator to use.
   :param matrix_type: Type of matrix to generate, see
      :c:type:`spral_matrix_type`.
   :param m: Number of rows in the matrix.
   :param n: Number of columns in the matrix.
   :param nnz: Number of non-zeroes in the matrix.
   :param ptr: Column pointers of the matrix
      (see :doc:`CSC format<csc_format>`).
   :param row: Row indices of the matrix
      (see :doc:`CSC format<csc_format>`).
   :param val: If not `NULL`, array of size `nnz` for non-zero values of the
      matrix (see :doc:`CSC format<csc_format>`).
   :param flags: Logical combination (i.e. bitwise-or) of the following
      possible values:

         SPRAL_RANDOM_MATRIX_FINDEX
            Entries in the arrays `ptr[]` and `row[]` should be numbered from
            1 (i.e. Fortran indexing). If this flag is not set, these arrays
            will be numbered from 0 (i.e. C indexing).

         SPRAL_RANDOM_MATRIX_NONSINGULAR
            Matrix will have a transversal of size
            :math:`\min({\tt m}, {\tt n})`. In the symmetric or skew
            symmetric case this will be the natural diagonal. In the unsymmetric
            and rectangular cases a random matching is used. If this flag is not
            set, the matrix may or may not be structurally singular.
            Note that symmetric positive-definite matrices are always
            non-singular.

         SPRAL_RANDOM_MATRIX_SORT
            Matrix will have entries sorted into ascending order within columns.
            If this flag is not set, entries may occur in any order.

   :returns: 0 on success, otherwise refer to table below for error code.

   Possible exit status values are:

   +--------+-----------------------------------------------------------------+
   |  0     | Success                                                         |
   +--------+-----------------------------------------------------------------+
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

.. c:function:: int spral_random_matrix_generate_long(int *state, enum spral_matrix_type matrix_type, int m, int n, int64_t nnz, int64_t ptr[n+1], int row[nnz], double *val, int flags)

   As :c:func:`spral_random_matrix_generate`, except ``nnz`` and ``ptr`` are
   ``int64_t``.

======
Macros
======

The following preprocessor macros are defined:

.. c:macro:: SPRAL_RANDOM_MATRIX_FINDEX 1

   Flag to use Fortran indexing on call to
   :c:func:`spral_random_matrix_generate()`.

.. c:macro:: SPRAL_RANDOM_MATRIX_NONSINGULAR 2

   Flag to generate non-singular matrix on call to
   :c:func:`spral_random_matrix_generate()`.

.. c:macro:: SPRAL_RANDOM_MATRIX_SORT 4

   Flag to sort row indices on call to
   :c:func:`spral_random_matrix_generate()`.

=======
Example
=======

The following code generates a random :math:`4 \times 5` matrix with
:math:`8` non-zeroes that is non-singular.

.. literalinclude:: ../../examples/C/random_matrix.c
   :language: C

This produces the following output:

::

    Generating a 4 x  5 non-singular matrix with 8 non-zeroes
    Generated matrix:
    Matrix of undefined type, dimension 4x5 with 8 entries.
    0:                                         -1.0744E-01   9.1000E-01
    1:                             9.5364E-01                1.0912E-01
    2:                                          1.1631E-01  -5.8957E-01
    3:  -9.0631E-01                                          7.7313E-01

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
