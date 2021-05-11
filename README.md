![STFC logo](http://www.stfc.ac.uk/stfc/includes/themes/MuraSTFC/assets/legacy/2473_web_2.png)

# SPRAL: The Sparse Parallel Robust Algorithm Library

[![travis status](https://travis-ci.org/ralna/spral.svg?branch=master)](https://travis-ci.org/ralna/spral)

An open-source (BSD) library for sparse linear algebra and associated
algorithms. It is primarily developed by the Numerical Analysis group at
STFC Rutherford Appleton Laboratory ([hsl@stfc.ac.uk](mailto:hsl@stfc.ac.uk)).

## Documentation

For detailed information about the SPRAL packages and API see [Fortran
documentation](http://www.numerical.rl.ac.uk/spral/doc/latest/Fortran/)
or [C
documentation](http://www.numerical.rl.ac.uk/spral/doc/latest/C/).

## Packages

- **LSMR** - Solves sparse least squares problems using LSMR
  algorithm.
- **RANDOM** - Pseudo-random number generator.
- **RANDOM_MATRIX** - Generates random matrices for testing purposes.
- **RUTHERFORD_BOEING** - Read and write matrices in Rutherford-Boeing
  format.
- **SCALING** - Calculates matrix scalings through a variety of
  algorithms
- **SSIDS** - Sparse Symmetric Indefinite Direct Solver.
- **SSMFE** - Sparse Symmetric Matrix-Free Eigensolver. Uses
                      Jacobi-conjugate preconditioned gradients
                      method.

If the functionality you are looking for is not support, it may be offered by
our proprietary licenced [HSL Library](http://www.hsl.rl.ac.uk/)
(free to academics).

## Installation
We use a standard autotools-based build:
```bash
./autogen.sh # If compiling from fresh git checkout
mkdir build
cd build
../configure
make
make install
```

## Usage at a Glance
When using SSIDS, ensure the following environment variables are set:
```bash
export OMP_CANCELLATION=TRUE
export OMP_NESTED=TRUE
export OMP_PROC_BIND=TRUE
```
