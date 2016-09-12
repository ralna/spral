****************************************************
:f:mod:`spral_rutherford_boeing` - RB File Utilities
****************************************************
.. f:module spral_rutherford_boeing
   :synopsis: Rutherford Boeing Rile Utilities

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

.. f:subroutine:: rb_peek(filename,inform[,m,n,nelt,nvar,nval,type_code,title,identifier])

   Returns information about a matrix :math:`A` stored in specified file (only
   information from the file header is accessed).

   :p character(len=*) filename [in]: File to read.
   :p integer inform [out]: Exit status, see :ref:`table below <exit_status>`.
   :o integer m [out]: Number of rows in :math:`A`.
   :o integer n [out]: Number of columns in :math:`A`.
   :o integer(long) nelt [out]: Number of elements in file. If the matrix
      is assembled, 0 is returned.
   :o integer(long) nvar [out]: Number of row indices in file.
   :o integer(long) nval [out]: Number of reals in file.
   :o character(len=3) type_code [out]: The three letter type code from the
      file (see :ref:`table <type_code>`).
   :o character(len=72) title [out]: Title field of file.
   :o character(len=8) identifier [out]: Identifier field of file.

.. f:subroutine:: rb_peek_unit(unit,info[,m,n,nelt,nvar,nval,type_code,title,identifier,no_rewind])

   Returns information about a matrix :math:`A` accessible on the specified
   Fortran unit. The unit will be rewound to its original location unless
   ``no_rewind=.true.`` is specified.

   Arguments are as for :f:subr:`rb_peek()`, except:

   :p integer unit [in]: Fortran unit from which to read file.
   :o logical no_rewind [in]: If present and `.true.`, the file is not rewound
      after reading the header.

.. f:subroutine:: rb_read(filename,m,n,ptr,row,val,options,inform[,type_code,title,identifier,state])

   Reads a CSC format matrix from the specified file.

   :p character(len=*) filename [in]: File to read.
   :p integer m [out]: Number of rows in :math:`A`.
   :p integer n [out]: Number of columns in :math:`A`.
   :p integer ptr(\:) [allocatable, out]: Column pointers
      (see :doc:`CSC format <csc_format>`). Will be allocated by the routine
      to have sufficient size.
   :p integer row(\:) [allocatable, out]: Row indices
      (see :doc:`CSC format <csc_format>`). Will be allocated by the routine
      to have sufficient size.
   :p real val(\:) [allocatable, out]: Values of non-zero entries
      (see :doc:`CSC format <csc_format>`) if present and requested.
      Will be allocated by the routine to have sufficient size.
   :p rb_read_options options [in]: Options for reading matrix (see
      :f:type:`rb_read_options`).
   :p integer inform [out]: Exit status, see :ref:`table below <exit_status>`.
   :o character(len=3) type_code [out]: The three letter type code from the
      file (see :ref:`table <type_code>`).
   :o character(len=72) title [out]: Title field of file.
   :o character(len=8) identifier [out]: Identifier field of file.
   :o random_state state [inout]: Random state to use for random number
      generation (if required).

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
