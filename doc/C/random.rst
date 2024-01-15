***************************************
RANDOM - Pseudo-random number generator
***************************************

.. code-block:: C

   #include <spral_random.h> /* or <spral.h> for all packages */

=======
Purpose
=======

This package generates pseudo-random numbers using a linear congruential
generator. It should generate the same random numbers using any standards
compliant Fortran compiler on any architecture so long as the default
integer and real kinds are the same.

The seed can optionally be observed or specified by the user.

===========
Random Seed
===========
The random number generator's state is stored in the variable `state` that
is common to all calls. Before its first use, `state` should be initialized
by the user to an initial seed value. For convienience the following
preprocessor macro is defined.

.. c:macro:: SPRAL_RANDOM_INITIAL_SEED 486502

   Convience macro for initial random seed.

   Example of use:

   .. code-block:: C

      int state = SPRAL_RANDOM_INITIAL_SEED;

The user may change the seed at any time, for example to restore a previous
value. The same sequence of calls following the (re-)initialization of `seed`
to the same value will produce the same sequence of pseudo-random values.

========
Routines
========

.. c:function:: double spral_random_real(int *state, bool positive)

   Return a real uniformly at random from the interval :math:`(-1,1)`
   (`positive=true`) or :math:`(0,1)` (`positive=false`).

   :param state: current state of RNG.
   :param positive: if true, sample from :math:`(0,1)`;
      otherwise, sample from :math:`(-1,1)`.
   :returns: Sampled value.

.. c:function:: int spral_random_integer(int *state, int n)

   Return an int uniformly at random from the interval :math:`[1,n]`.

   :param state: current state of the RNG.
   :param n: largest value in range to be sampled.
   :returns: Sampled value.

.. c:function:: int64_t spral_random_long(int *state, int64_t n)

   Return a int64_t uniformly at random from the interval :math:`[1,n]`.

   :param state: current state of the RNG.
   :param n: largest value in range to be sampled.
   :returns: Sampled value.

.. c:function:: bool spral_random_logical(int *state)

   Return a logical with equal probability of being `true` or `false`.

   :param state: current state of the RNG.
   :returns: Sampled value.

=======
Example
=======
The following code:

.. literalinclude:: ../../examples/C/random.c
   :language: C

Produces the following output::

   Some random values
   Sample Unif(-1,1)               =   0.951878630556
   Sample Unif(0,1)                =   0.395779648796
   Sample Unif(1, ..., 20)         =                3
   Sample Unif(1, ..., 20*INT_MAX) =      33572664025
   Sample B(1,0.5)                 =            false

   The same random values again
   Sample Unif(-1,1)               =   0.951878630556
   Sample Unif(0,1)                =   0.395779648796
   Sample Unif(1, ..., 20)         =                3
   Sample Unif(1, ..., 20*INT_MAX) =      33572664025
   Sample B(1,0.5)                 =            false

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

The variable `state` stores the current value of :math:`X_n`.

In :c:func:`spral_random_real`
------------------------------

Samples from :math:`\mathrm{Unif}(0,1)` are generated as

.. math::

   \frac{\mathrm{real}(X_n)}{\mathrm{real}(m)},

and samples from :math:`\mathrm{Unif}(-1,1)` are generated as

.. math::

   1−\frac{\mathrm{real}(2X_n)}{\mathrm{real}(m)}.

In :c:func:`spral_random_integer`
---------------------------------

Samples from :math:`\mathrm{Unif}(1,\ldots,n)` are generated as

.. math::

   \mathrm{int}\left(X_n\frac{\mathrm{real}(n)}{\mathrm{real}(m)}\right) + 1

In :c:func:`spral_random_logical`
---------------------------------

Returns the value of the Fortran expression

.. code-block:: Fortran

   (1 .eq. spral_random_integer(state,2))
