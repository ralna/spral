******************************
MATRIX_UTIL - Matrix utilities
******************************

.. code-block:: C

   #include <spral_matrix_util.h> /* or <spral.h> for all packages */

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

==========
Data types
==========

.. c:type:: enum spral_matrix_type

   .. c:member:: SPRAL_MATRIX_UNSPECIFIED

      User doesn't wish to specify matrix type, use default behaviour for
      routine.

   .. c:member:: SPRAL_MATRIX_REAL_RECT

      Rectangular real-valued matrix, :math:`m\ne n`.

   .. c:member:: SPRAL_REAL_UNSYM

      Square real-valued unsymmetric matrix, :math:`m\eq n`.

   .. c:member:: SPRAL_REAL_SYM_PSDEF

      Symmetric real-valued positive-definite matrix.

   .. c:member:: SPRAL_REAL_SYM_INDEF

      Symmetric real-valued indefinite matrix.

   .. c:member:: SPRAL_REAL_SKEW

      Skew-symmetric real-valued matrix.

   .. c:member:: SPRAL_MATRIX_CPLX_RECT

      Rectangular complex-valued matrix, :math:`m\ne n`.

   .. c:member:: SPRAL_CPLX_UNSYM

      Square complex-valued unsymmetric matrix, :math:`m\eq n`.

   .. c:member:: SPRAL_CPLX_HERM_PSDEF

      Hermitian complex-valued positive-definite matrix.

   .. c:member:: SPRAL_CPLX_SYM_INDEF

      Hermitian complex-valued indefinite matrix.

   .. c:member:: SPRAL_CPLX_SYM

      Symmetric complex-valued matrix.

   .. c:member:: SPRAL_CPLX_SKEW

      Skew-symmetric complex-valued matrix.
