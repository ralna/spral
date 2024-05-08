****************************************************
RUTHERFORD_BOEING - RB File Utilities
****************************************************

.. code-block:: C

   #include <spral_rutherford_boeing.h> /* or <spral.h> for all packages */

=======
Purpose
=======

This package provides routines to read and write matrices stored in
files using the :ref:`Rutherford-Boeing format <rb_format>`.

At present, reading and writing of supplementary data (e.g. right-hand sides)
is not supported. If it is required to read and write these files, we
recommend the use of the HSL routines
`MC55 <http://www.hsl.rl.ac.uk/catalogue/mc55.html>`_ and
`MC56 <http://www.hsl.rl.ac.uk/catalogue/mc56.html>`_ that are available
without charge (although redistribution is not permitted).

Version history
---------------

2016-09-08 Version 1.0.0
   Initial release.

[For detailed history, see ChangeLog]

===========
Subroutines
===========

.. c:function:: void spral_rb_default_read_options(struct spral_rb_read_options *options)

   Intialises members of :c:type:`spral_rb_read_options` structure to default
   values.

   :param options: Structure to be initialised.

.. c:function:: void spral_rb_default_write_options(struct spral_rb_write_options *options)

   Intialises members of :c:type:`spral_rb_write_options` structure to default
   values.

   :param options: Structure to be initialised.

.. c:function:: int rb_peek(const char *filename, int *m, int *n, int64_t *nelt, int64_t *nvar, int64_t *nval, enum spral_matrix_type *matrix_type, char *type_code, char *title, char *identifier)

   Returns information about a matrix :math:`A` stored in a specified file using
   Rutherford Boring format (only information from the file header is accessed).

   :param filename: File to read.
   :param m: If not `NULL`, set to number of rows in :math:`A`.
   :param n: If not `NULL`, set to number of columns in :math:`A`.
   :param nelt: If not `NULL`, set to number of elements in file. Set to 0 if
      the matrix is not in elemental format.
   :param nvar: If not `NULL`, set to number of row indices in file.
   :param nval: If not `NULL`, set to number of reals in file.
   :param matrix_type: If not `NULL`, set to type of matrix
      (see :c:type:`spral_matrix_type`). Note that RB files do not
      distguish between symmetric positive-definite and symmetric indefinite
      matrices, so the latter matrix type is used for all symmetric matrices.
   :param type_code: If not `NULL`, must point to a length 4
      character buffer. Set to the three letter type code from the
      file (see :ref:`table <type_code>`).
   :param title: If not `NULL`, must point to a length 73
      character buffer. Set to the title field of file.
   :param identifier: If not `NULL`, must point to a length 9
      character buffer. Identifier field of file.
   :returns: Exit status, see :ref:`table below <exit_status>`.

.. c:function:: int spral_rb_read(const char *filename, void **handle, enum spral_matrix_type *matrix_type, int *m, int *n, int64_t **ptr, int **row, double **val, const struct spral_rb_read_options *options, char *title, char *identifier, int *state)

   Reads a CSC matrix from a file stored in RB format.

   .. warning::
      Memory is allocated using Fortran routines and stored in `handle`.
      The arrays describing the matrix are just pointers into this data
      structure. It can only be freed using the
      :c:func:`spral_rb_free_handle()` routine, and not by `free()`.

   :param filename: File to read.
   :param handle: Handle for underlying memory holding matrix. Must be freed
      using :c:func:`spral_rb_free_handle()` after last access to `ptr`, `row`
      and `val`.
   :param matrix_type: Type of matrix to read
      (see :c:type:`spral_matrix_type`). Note that RB files do not
      distguish between symmetric positive-definite and symmetric indefinite
      matrices, so the latter matrix type is used for all symmetric matrices.
   :param m: Number of rows in :math:`A`.
   :param n: Number of columns in :math:`A`.
   :param ptr: Column pointers (see :doc:`CSC format <csc_format>`).
      Set to point to an array allocated by the routine. Invalidated when
      `handle` is freed.
   :param row: Row indices
      (see :doc:`CSC format <csc_format>`).
      Set to point to an array allocated by the routine. Invalidated when
      `handle` is freed.
   :param val: If not `NULL`, values of non-zero entries
      (see :doc:`CSC format <csc_format>`).
      Set to point to an array allocated by the routine. Invalidated when
      `handle` is freed.
   :param options: Options for reading matrix (see
      :c:type:`spral_rb_read_options`).
   :param title: If not `NULL`, must point to a length 73
      character buffer. Set to the title field of file.
   :param identifier: If not `NULL`, must point to a length 9
      character buffer. Identifier field of file.
   :param state: If not `NULL`, the random state to use for random number
      generation (see :doc:`random`).
   :returns: Exit status, see :ref:`table below <exit_status>`.

.. c:function:: int spral_rb_read_ptr32(const char *filename, void **handle, enum spral_matrix_type *matrix_type, int *m, int *n, int **ptr, int **row, double **val, const struct spral_rb_read_options *options, char *title, char *identifier, int *state)

   As :c:func:`spral_rb_read()` except ``ptr`` has type ``int``.

   .. note::

      This is just a wrapper around the 64-bit call.

      Users are encouraged to prefer the 64-bit version.

.. c:function:: int spral_rb_write(const char *filename, enum spral_matrix_type matrix_type, int m, int n, const int64_t *ptr, const int *row, const double * val, const struct spral_rb_write_options *options, const char *title, const char *identifier)

   Writes a CSC format matrix to the specified file in RB format.

   :param filename: File to write. Existing files will be overwritten.
   :param matrix_type: Type of matrix to write, see
      :c:type:`enum spral_matrix_type`. (will be converted into the second
      character of the :ref:`type code <type_code>`).
   :param m: Number of rows in :math:`A`.
   :param n: Number of columns in :math:`A`.
   :param ptr[n]: Column pointers (see :doc:`CSC format <csc_format>`).
   :param row[ptr[n]]: Row indices (see :doc:`CSC format <csc_format>`).
   :param val[ptr[n]]: Values of non-zero entries
      (see :doc:`CSC format <csc_format>`). If `NULL` a pattern only matrix
      is written.
   :param options: Options for writing matrix (see
      :c:type:`spral_rb_write_options`).
   :param title: The title field for the file. Maximum length is 72
      characters. Defaults to ``"Matrix"`` if `NULL`.
   :param identifier: Identifier field of file. Maximum
      length is 8 characters. Defaults to ``"0"`` if `NULL`.
   :returns: Exit status, see :ref:`table below <exit_status>`.

.. c:function:: int spral_rb_write_ptr32(const char *filename, enum spral_matrix_type matrix_type, int m, int n, const int *ptr, const int *row, const double * val, const struct spral_rb_write_options *options, const char *title, const char *identifier)

   As :c:func:`spral_rb_write()` except ``ptr`` has type ``int``.

   .. note::

      This is just a wrapper around the 64-bit call.

      Users are encouraged to prefer the 64-bit version to ensure support
      for large matrices.


Return codes
------------

.. table:: Return codes
   :name: exit_status

   +------------+-------------------------------------------------------------+
   | `inform`   | Status                                                      |
   +============+=============================================================+
   | 0          | Success                                                     |
   +------------+-------------------------------------------------------------+
   | -1         | Failed to open file.                                        |
   +------------+-------------------------------------------------------------+
   | -2         | Not a valid Rutherford-Boeing file.                         |
   +------------+-------------------------------------------------------------+
   | -3         | Error on i/o operation.                                     |
   +------------+-------------------------------------------------------------+
   | -4         | Attempted to read data type not supported by routine.       |
   +------------+-------------------------------------------------------------+
   | -5         | Attempted to read matrix in elemental format.               |
   +------------+-------------------------------------------------------------+
   | -6         | Invalid matrix type.                                        |
   +------------+-------------------------------------------------------------+
   | -10        | `options%extra_space<1.0`.                                  |
   +------------+-------------------------------------------------------------+
   | -11        | `options%lwr_upr_full` has invalid value.                   |
   +------------+-------------------------------------------------------------+
   | -12        | `options%values` has invalid value.                         |
   +------------+-------------------------------------------------------------+
   | -20        | Memory allocation failed.                                   |
   +------------+-------------------------------------------------------------+
   | +1         | Values are stored in an auxiliary file (not read)           |
   +------------+-------------------------------------------------------------+


=============
Derived types
=============

.. c:type:: struct spral rb_read_options

   Specify options for reading matrices.

   .. c:member:: bool add_diagonal

      Add any diagonal entries that are missing if true.
      Default is false

   .. c:member:: float extra_space

      Proportion of extra space to allow in `row[]` and `val[]` arrays.
      They are allocated to have size `options.extra_space * ptr[n]`.
      Default is 1.0.

   .. c:member:: int lwr_upr_full

      Return lower triangle, upper triangle
      or both for symmetric and skew-symmetric matrices. One of:

      +-------------+---------------------------------------------------------+
      | 1 (default) | Lower triangular entries only.                          |
      +-------------+---------------------------------------------------------+
      | 2           | Upper triangular entries only.                          |
      +-------------+---------------------------------------------------------+
      | 3           | Both lower and upper triangular entries.                |
      +-------------+---------------------------------------------------------+

      Default is 1.

   .. c:member:: int values

      Whether to read and/or generate values. One of:

      +-------------+---------------------------------------------------------+
      | 0 (default) | Read values from file, only if present.                 |
      +-------------+---------------------------------------------------------+
      | 1           | Do not read or generate value, return pattern only.     |
      +-------------+---------------------------------------------------------+
      | 2           | Read values from file. If no values are present,        |
      |             | randomly generate symmetric values.                     |
      +-------------+---------------------------------------------------------+
      | 3           | Read values from file. If no values are present,        |
      |             | randomly generate symmetric, diagonally dominant        |
      |             | values.                                                 |
      +-------------+---------------------------------------------------------+
      | 4           | Read values from file. If no values are present,        |
      |             | randomly generate (unsymmetric) values.                 |
      +-------------+---------------------------------------------------------+
      | -2          | Randomly generate symmetric values. Any values in file  |
      |             | are ignored.                                            |
      +-------------+---------------------------------------------------------+
      | -3          | Randomly generate symmetric, diagonally dominant        |
      |             | values. Any values in file are ignored.                 |
      +-------------+---------------------------------------------------------+
      | -4          | Randomly generate (unsymmetric) values. Any values in   |
      |             | file are ignored.                                       |
      +-------------+---------------------------------------------------------+

      Default is 0.

.. c:type:: struct spral_rb_write_options

   Specify options for writing matrices.

   .. c:member:: char val_format[20]

      Fortran format string to use when writing values. Should not exceed 80
      characters per line.
      Default is "(3e24.16)".

   .. note::

      Formats for integer data will be automatically determined based on the
      maximum values to be represented.

=======
Example
=======

Reading a matrix
----------------

The following code reads a matrix from the file "matrix.rb":

.. literalinclude:: ../../examples/C/rutherford_boeing/rb_read.c
   :language: C

This produces the following output, when run on the file generated by the
`rb_write.c` example in the next section::

Matrix 'SPRAL_RUTHERFORD_BOEING test matrix'
Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
0:   2.0000E+00   1.0000E+00
1:   1.0000E+00   4.0000E+00   1.0000E+00                8.0000E+00
2:                1.0000E+00   3.0000E+00   2.0000E+00
3:                             2.0000E+00
4:                8.0000E+00                             2.0000E+00

Writing a matrix
----------------

The following code writes a matrix to the file "matrix.rb":

.. literalinclude:: ../../examples/C/rutherford_boeing/rb_write.c
   :language: C

This produces the following file::

   SPRAL_RUTHERFORD_BOEING test matrix                                     0
                5             1             1             3
   rsa                        5             5             8             0
   (40i2)          (40i2)          (3e24.16)
    1 3 6 8 8 9
    1 2 2 3 5 3 4 5
     0.2000000000000000E+01  0.1000000000000000E+01  0.4000000000000000E+01
     0.1000000000000000E+01  0.8000000000000000E+01  0.3000000000000000E+01
     0.2000000000000000E+01  0.2000000000000000E+01

======
Method
======

Generation of random values
---------------------------

Values are generated uniformly at random from :math:`[-1,1]`.
If a diagonally dominant matrix is requested, the diagonal entry in each
row is set to :math:`\max(100, 10k)`, where `k` is the number of entries
in the column.

If a random `state` is not provided by the user, the default initial state
from the :f:mod:`spral_random` module is used.

.. _rb_format:

Rutherford Boeing File format
-----------------------------

The Rutherford Boeing file format is described in the following report [1]_.
A file may either contain a sparse matrix or supplementary data. However, this
package only supports the former.

Sparse matrices are stored in either :doc:`CSC format <csc_format>` or
`Elemental format`. A three letter type code is used to encodes this and
additional information, as per the following table:

.. table:: Type code
   :name: type_code

   +-----------------+-----------------------+-------------------------+
   | First character | Second character      | Third character         |
   +=================+=======================+=========================+
   | **r**: real     | **s**: symmetric      | **a**: CSC format       |
   +-----------------+-----------------------+-------------------------+
   | **c**: complex  | **u**: unsymmetric    | **e**: Elemental format |
   +-----------------+-----------------------+-------------------------+
   | **i**: integer  | **h**: Hermitian      |                         |
   +-----------------+-----------------------+-------------------------+
   | **p**: pattern  | **z**: skew symmetric |                         |
   +-----------------+-----------------------+-------------------------+
   | **q**: pattern  | **r**: rectangular    |                         |
   +-----------------+-----------------------+-------------------------+

The difference between the **p** and **q** pattern types is that the latter
indicates that values are supplied in an auxiliary file (this package does not
support reading such files).

Further information may be found in:

.. [1] I.S. Duff, R.G. Grimes and J.G. Lewis (1997).
   *The Rutherford-Boeing Sparse Matrix Collection*.
   RAL Technical Report `RAL-TR-97-031 <http://purl.org/net/epubs/work/26879>`_.
