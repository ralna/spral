****************************************************
RUTHERFORD_BOEING - RB File Utilities
****************************************************

.. code-block:: C

   #include <spral_rutherford_boeing.h> /* or <spral.h> for all packages */

=======
Purpose
=======

This package provides routines to read and write matrices stored in
Rutherford-Boeing format.

At present, reading and writing of supplementary data (e.g. right-hand sides)
is not supported. If it is required to read and write these files, we
recommend the use of the HSL routines
`MC55 <http://www.hsl.rl.ac.uk/catalogue/mc55.html>`_ and
`MC56 <http://www.hsl.rl.ac.uk/catalogue/mc56.html>`_ that are available
for free (though redistribution is not permitted).

Version history
---------------

2016-09-08 Version 1.0.0
   Initial release.

[For detailed history, see ChangeLog]

===========================
Data formats and type codes
===========================

A Rutherford Boeing (RB) file may either contain a sparse matrix or
supplementary data. This package only supports the former.

The RB format support storage of sparse matrices in either
:doc:`CSC format <csc_format>` or in :doc:`Element format <element_format>`.
Additionally a three letter type code is supplied that encodes additional
information, as per the following table:

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

The difference between the `p` and `q` pattern types is that the latter
indicates that values are supplied in an auxiliary file (this package does not
support reading such files).

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

.. c:function:: int rb_peek(const char *filename, int *m, int *n, long *nelt, long *nvar, long *nval, enum spral_matrix_type *matrix_type, char *type_code, char *title, char *identifier)

   Returns information about a matrix :math:`A` stored in specified file (only
   information from the file header is accessed).

   :param filename: File to read.
   :param m: If not `NULL`, set to number of rows in :math:`A`.
   :param n: If not `NULL`, set to number of columns in :math:`A`.
   :param nelt: If not `NULL`, set to number of elements in file. If the matrix
      is assembled, 0 is returned.
   :param nvar: If not `NULL`, set to number of row indices in file.
   :param nval: If not `NULL`, set to number of reals in file.
   :param matrix_type: If not `NULL`, set to type of matrix
      (see :c:type:`spral_matrix_type`). Note that RB files do not
      distguish between symmetric positive-definite and symmetric indefinite
      matrices, so the latter matrix type is used for type code `'s'`.
   :param type_code: If not `NULL`, must point to a length 4
      character buffer. Set to the three letter type code from the
      file (see :ref:`table <type_code>`).
   :param title: If not `NULL`, must point to a length 73
      character buffer. Set to the title field of file.
   :param identifier: If not `NULL`, must point to a length 9
      character buffer. Identifier field of file.
   :returns: Exit status, see :ref:`table below <exit_status>`.

.. c:function:: int spral_rb_read(const char *filename, void **handle, enum spral_matrix_type *matrix_type, int *m, int *n, long **ptr, int **row, double **val, struct spral_rb_read_options *options, char *title, char *identifier, int *state)

   Reads a CSC format matrix from the specified file.
   
   .. warning::
      Memory is allocated using Fortran routines and stored in `handle`.
      The arrays describing the matrix are just pointers into this data
      structure. It can only be freed using the
      :c:func:`spral_rb_free_handle()` routine, and not by `free()`.

   :param filename: File to read.
   :param handle: Handle for underlying memory holding matrix. Must be freed
      using :c:func:`spral_rb_free_handle()` after last access to `ptr`, `row`
      and `val`.
   :param matrix_type: Type of matrix to write
      (see :c:type:`spral_matrix_type`). Note that RB files do not
      distguish between symmetric positive-definite and symmetric indefinite
      matrices, so the latter matrix type is used for type code `'s'`.
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

.. c:function:: int spral_rb_read_ptr32(const char *filename, void **handle, enum spral_matrix_type *matrix_type, int *m, int *n, int **ptr, int **row, double **val, struct spral_rb_read_options *options, char *title, char *identifier, int *state)

   As :c:func:`spral_rb_read()` except ``ptr`` has type ``int``.
   Users are encouraged to prefer the 64-bit version.

.. f:subroutine:: rb_write(filename,matrix_type,m,n,ptr,row,val,options,inform[,title,identifier])

   Writes a CSC format matrix to the specified file.

   :p character(len=*) filename [in]: File to write. Existing files will be
      overwritten.
   :p integer matrix_type [in]: Type of matrix to write, one of the values
      defined in :f:mod:`spral_matrix_util` (will be converted into the second
      character of the :ref:`type code <type_code>`).
   :p integer m [in]: Number of rows in matrix.
   :p integer n [in]: Number of columns in matrix.
   :p integer ptr(n+1) [in]: Column pointers
      (see :doc:`CSC format <csc_format>`).
   :p integer row(ptr(n+1)-1) [in]: Row indices
      (see :doc:`CSC format <csc_format>`).
   :p real val(ptr(n+1)-1) [in]: Values of non-zero entries
      (see :doc:`CSC format <csc_format>`).
   :p rb_write_options options [in]: Options for writing matrix (see
      :f:type:`rb_write_options`).
   :p integer inform [out]: Exit status, see :ref:`table below <exit_status>`.
   :o character(len=*) title [in]: Title field of file. Maximum length is 72
      character. Defaults to ``"Matrix"`` if not present.
   :o character(len=*) identifier [in]: Identifier field of file. Maximum
      length is 8 characters. Defaults to ``"0"`` if not present.

Exit status codes
-----------------

.. table:: Exit status codes
   :name: exit_status

   +------------+-------------------------------------------------------------+
   | ``inform`` | Status                                                      |
   +============+=============================================================+
   | 0          | Success                                                     |
   +------------+-------------------------------------------------------------+
   | -1         | Failed to open file.                                        |
   +------------+-------------------------------------------------------------+
   | -2         | Not a valid for Rutherford-Boeing file.                     |
   +------------+-------------------------------------------------------------+
   | -3         | Error on i/o operation.                                     |
   +------------+-------------------------------------------------------------+
   | -4         | Attempted to read data type not supported by routine.       |
   +------------+-------------------------------------------------------------+
   | -5         | Attempted to read element as assembled or vice versa.       |
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

.. f:type:: rb_read_options

   Specify options for reading matrices.

   :f logical add_diagonal [default=.false.]: Add any diagonal entries that are
      missing.
   :f real extra_space [default=1.0]: Proportion of extra space to allow in
      `row(:)` and `val(:)` arrays. They are allocated to have size
      `options%extra_space * (ptr(n+1)-1)`.
   :f integer lwr_upr_full [default=1]: Return lower triangle, upper triangle
      or both for symmetric and skew-symmetric matrices. One of:

      +-------------+---------------------------------------------------------+
      | 1 (default) | Lower triangular entries only.                          |
      +-------------+---------------------------------------------------------+
      | 2           | Upper triangular entries only.                          |
      +-------------+---------------------------------------------------------+
      | 3           | Both lower and upper triangular entries.                |
      +-------------+---------------------------------------------------------+

   :f integer values [default=0]: Whether to read and/or generate values. One
      of:

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

.. f:type:: rb_write_options

   Specify options for writing matrices.
   
   :f character(len=20) val_format [default="(3e24.16)"]: Fortran format
      string to use when writing values. Should not exceed 80 characters per
      line.

=======
Example
=======

Reading a matrix
----------------

The following code reads a matrix from the file "matrix.rb":

.. literalinclude:: ../../examples/Fortran/rutherford_boeing/rb_read.f90
   :language: Fortran

This produces the following output, when run on the file generated by the
`rb_read.f90` example in the next section::

   Matrix 'SPRAL_RUTHERFORD_BOEING test matrix'
   Real symmetric indefinite matrix, dimension 5x5 with 8 entries.
   1:   2.0000E+00   1.0000E+00                                       
   2:   1.0000E+00   4.0000E+00   1.0000E+00                8.0000E+00
   3:                1.0000E+00   3.0000E+00   2.0000E+00             
   4:                             2.0000E+00                          
   5:                8.0000E+00                             2.0000E+00


Writing a matrix
----------------

The following code writes a matrix to the file "matrix.rb":

.. literalinclude:: ../../examples/Fortran/rutherford_boeing/rb_write.f90
   :language: Fortran

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

If a random `state` isn't provided by the user, the default initial state
from the :f:mod:`spral_random` module is used.

File format
-----------

The Rutherford Boeing file format is described in the following report:

.. [1] I.S. Duff, R.G. Grimes and J.G. Lewis (1997).
   *The Rutherford-Boeing Sparse Matrix Collection*.
   RAL Technical Report `RAL-TR-97-031 <http://purl.org/net/epubs/work/26879>`_.
