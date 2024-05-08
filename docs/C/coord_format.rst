=========================
Coordinate (Coord) Format
=========================

This standard data format consists of the following data:

.. code-block:: C

   int    m;         /* number of rows (unsymmetric only) */
   int    n;         /* number of columns */
   int    ne;        /* number of entries in matrix (may have type int64_t) */
   int    row[ne];   /* row indices */
   int    col[ne];   /* column indices */
   double val[ne];   /* numerical values */

The arrays should be set such that the ``k``-th entry is in row
``row[k]`` and column ``col[k]`` with value ``val[k]``. Entries that are
zero, including those on the diagonal, need not be specified.

For **symmetric matrices**, only the lower *or* upper triangular entries of
:math:`A` should be supplied. For **unsymmetric matrices**, all entries in the
matrix should be supplied. Duplicate entries will be summed and out-of-range
entries will be ignored.

Some SPRAL routines offer only input in :doc:`CSC format<csc_format>`, you
may need to use routines from the :doc:`matrix_util` package to convert
data from Coordinate to CSC format.

To illustrate the Coordinate format, the matrix

.. math::

   \left( \begin{array}{ccccc}
      1.1 & 2.2 &     & 3.3 &     \\
      2.2 &     & 4.4 &     &     \\
          & 4.4 & 5.5 &     & 6.6 \\
      3.3 &     &     & 7.7 & 8.8 \\
          &     & 6.6 & 8.8 & 9.9
   \end{array} \right)

is described by the following data:

.. code-block:: C

   int    n     = 5;
   int    ne    = 9;
   int    row[] = { 1,   2,   3,   4,   3,   5,   4,   5,   5 };
   int    col[] = { 1,   1,   2,   1,   3,   3,   4,   4,   5 };
   double val[] = { 1.1, 2.2, 4.4, 3.3, 5.5, 6.6, 7.7, 8.8, 9.9 };
