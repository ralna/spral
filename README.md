# SPRAL: The Sparse Parallel Robust Algorithm Library

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

If the functionality you are looking for is not supported, it may be offered by
our proprietary licenced [HSL Library](http://www.hsl.rl.ac.uk/)
(free to academics).

## Installation
Please note that we require [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) 
and [hwloc](https://www.open-mpi.org/projects/hwloc/) to be installed 
(hwloc should be compiled with CUDA support if building for GPU).

We use a standard autotools-based build:
```bash
./autogen.sh # If compiling from fresh git checkout
mkdir build
cd build
../configure --with-metis="-L/path/to/metis -lmetis"
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
