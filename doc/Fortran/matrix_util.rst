*********************************************
:f:mod:`spral_matrix_util` - Matrix utilities
*********************************************
.. f:module:: spral_matrix_util
   :synopsis: Matrix utilities

=======
Purpose
=======

This packages contains assorted utility routines and datatypes for:
   * matrix data storage and format conversion
   * printing of matrices

Version history
---------------

2016-09-07 Version 0.1.0
   API still open to redefinition: only documented entries are considered at
   all "fixed".


=================
Matrix type codes
=================

.. f:variable:: integer SPRAL_MATRIX_UNSPECIFIED = 0

   User doesn't wish to specify matrix type, use default behaviour for
   routine.

.. f:variable:: integer SPRAL_MATRIX_REAL_RECT = 1

   Rectangular real-valued matrix, :math:`m\ne n`.

.. f:variable:: integer SPRAL_REAL_UNSYM = 2

   Square real-valued unsymmetric matrix, :math:`m=n`.

.. f:variable:: integer SPRAL_REAL_SYM_PSDEF = 3

   Symmetric real-valued positive-definite matrix.

.. f:variable:: integer SPRAL_REAL_SYM_INDEF = 4

   Symmetric real-valued indefinite matrix.

.. f:variable:: integer SPRAL_REAL_SKEW = 6

   Skew-symmetric real-valued matrix.

.. f:variable:: integer SPRAL_MATRIX_CPLX_RECT = -1

   Rectangular complex-valued matrix, :math:`m\ne n`.

.. f:variable:: integer SPRAL_CPLX_UNSYM = -2

   Square complex-valued unsymmetric matrix, :math:`m=n`.

.. f:variable:: integer SPRAL_CPLX_HERM_PSDEF = -3

   Hermitian complex-valued positive-definite matrix.

.. f:variable:: integer SPRAL_CPLX_SYM_INDEF = -4

   Hermitian complex-valued indefinite matrix.

.. f:variable:: integer SPRAL_CPLX_SYM = -5

   Symmetric complex-valued matrix.

.. f:variable:: integer SPRAL_CPLX_SKEW = -6

   Skew-symmetric complex-valued matrix.
