******************************************************
:f:mod:`spral_random` - Pseudo-random number generator
******************************************************
.. f:module:: spral_random
   :synopsis: Random number generator

=======
Purpose
=======

This package generates pseudo-random numbers using a linear congruential
generator. It should generate the same random numbers using any standards
compliant Fortran compiler on any architecture so long as the default
integer and real kinds are the same.

The seed can optionally be observed or specified by the user. Otherwise
a default seed of 486502 is used.

Version history
---------------

2016-09-08 Version 1.1.0
   Add support for long integers

2014-04-07 Version 1.0.0
   Initial release

========
Routines
========

Random Number Generation
------------------------

.. f:function:: random_real(state[, positive])

   Return a real uniformly at random from the interval :math:`(-1,1)`
   (positive =.true.) or :math:`(0,1)` (positive =.false.).

   :p random_state state [inout]: current state of RNG.
   :o logical positive [in,default=.false.]:
      if .true., sample from :math:`(0,1)`;
      otherwise, sample from :math:`(-1,1)`.
   :r random_real: Sampled value.
   :rtype random_real: real

.. f:function:: random_integer(state, n)

   Return an integer uniformly at random from the interval :math:`[1,n]`.

   :p random_state state [inout]: current state of the RNG.
   :p integer(kind) n [in]: largest value in range to be sampled. `kind` may be
      either default or long integer (return type will match).
   :r random_integer: Sampled value.
   :rtype random_integer: integer(kind)

.. f:function:: random_logical(state)

   Return a logical with equal probability of being .true. or .false..

   :p random_state state [inout]: current state of the RNG.
   :r random_logical: Sampled value.
   :rtype random_logical: logical

Get/Set Random Seed
-------------------

.. f:function:: random_get_seed(state)

   Return the current random seed stored in state.

   The stream of random numbers generated after this call can be reproduced
   through the same sequence of calls after seed has been passed to
   :f:subr:`random_set_seed`.

   :p random_state state [in]: state variable to extract state from
   :r random_get_seed: current random seed
   :rtype: integer

.. f:subroutine:: random_set_seed(state, seed)

   Set the random seed stored in the state variable

   :p random_state state [inout]: state variable to set seed for.
   :p integer seed [in]: new seed.

==========
Data Types
==========

.. f:type:: random_state

   State of the random number generator.
   Compontents are not available to the user, but may be examined and
   altered through calls to :f:func:`random_get_seed` and
   :f:func:`random_set_seed` respectively.

=======
Example
=======
The following code:

.. literalinclude:: ../../examples/Fortran/random.f90
   :language: Fortran

Produces the following output::

   Some random values
   Sample Unif(-1,1)               =   0.951878630556
   Sample Unif(0,1)                =   0.395779648796
   Sample Unif(1, ..., 20)         =                3
   Sample Unif(1, ..., 20*huge(0)) =      33572664025
   Sample B(1,0.5)                 =                F

   The same random values again
   Sample Unif(-1,1)               =   0.951878630556
   Sample Unif(0,1)                =   0.395779648796
   Sample Unif(1, ..., 20)         =                3
   Sample Unif(1, ..., 20*huge(0)) =      33572664025
   Sample B(1,0.5)                 =                F

======
Method
======
We use a linear congruential generator of the following form:

.. math::

   X_{n+1} = (aX_n + c)\quad \mathrm{mod}\; m

with the following constants

.. math::

   a = 1103515245, \\
   c = 12345, \\
   m = 2^{31}.

According to Wikipedia, this is the same as used in glibc.

The LCG is evolved before each sample is taken, and the sample is based on the
new value.

The routines :f:func:`random_get_seed` and :f:func:`random_set_seed` allow the
user to get and set the current value of :math:`X_n`. The default seed is
:math:`X_0 = 486502`.

In :f:func:`random_real`
------------------------

Samples from :math:`\mathrm{Unif}(0,1)` are generated as

.. math::

   \frac{\mathrm{real}(X_n)}{\mathrm{real}(m)},

and samples from :math:`\mathrm{Unif}(-1,1)` are generated as

.. math::

   1−\frac{\mathrm{real}(2X_n)}{\mathrm{real}(m)}.

In :f:func:`random_integer`
---------------------------

Samples from :math:`\mathrm{Unif}(1,\ldots,n)` are generated as

.. math::

   \mathrm{int}\left(X_n\frac{\mathrm{real}(n)}{\mathrm{real}(m)}\right) + 1

In :f:func:`random_logical`
----------------------------

Returns the value of the Fortran expression

.. code-block:: Fortran

   (1 .eq. random_integer(state,2))
