**************************
Data types and conventions
**************************

Data types
----------

In this documentation, we use:

   * **integer(long)** to denote 64 bit integers, where
     **long = selected_int_kind(18)**.
   * **real** to mean either 32-bit or 64-bit floating point. Most packages
     only implement the latter (i.e. double precision).


Optional arguments
------------------

Optional arguments are indicated by square brackets [] in argument lists.
We strongly recommend the use of keywords in argument lists when using
optional arguments, rather than relying on position to indicate which you are
using. For example,

.. code-block:: Fortran

   call rb_peek("matrix.rb", inform, matrix_type=mat_type, title=mat_title)
