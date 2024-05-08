=====================================
Compressed Sparse Column (CSC) Format
=====================================

This standard data format consists of the following data:

.. code-block:: C

   int    m;                  /* number of rows (unsymmetric only) */
   int    n;                  /* number of columns */
   int    ptr[n+1];           /* column pointers (may have type int64_t) */
   int    row[ ptr[n]-1 ];    /* row indices */
   double val[ ptr[n]-1 ];    /* numerical values */

Non-zero matrix entries are ordered by increasing column index and stored in
the arrays `row[]` and `val[]` such that `row[k]` holds
the row number and `val[k]` holds the value of the k-th entry.
The `ptr[]` array stores column pointers such that `ptr[i]` is
the position in `row[]` and `val[]` of
the first entry in the i-th column, and `ptr[n]` is
the total number of entries. `ptr[]` may be either `int` or `int64_t`.
There must be no duplicate or out of range entries.
Entries that are zero, including those on the diagonal, need not be specified.

For **symmetric matrices**, only the lower triangular entries of :math:`A`
should be supplied. For **unsymmetric matrices**, all entries in the matrix
should be supplied.

Note that most SPRAL routines offer **no checking** of user data, and the
behaviour of these routines with misformatted data is undefined. You may use
routines from the :doc:`matrix_util` package to convert data to and
check data stored in this format.

To illustrate the CSC format, the matrix

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
   int    ptr[] = { 1,             4,   5,        7,        9,    10 };
   int    row[] = { 1,   2,   4,   3,   3,   5,   4,   5,   5 };
   double val[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9 };
