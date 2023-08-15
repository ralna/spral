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
cp nvcc_arch_sm.c build/ # If building for GPU
cd build
../configure --with-metis="-L/path/to/metis -lmetis"
make
make install
```

## Usage at a Glance
When using SSIDS, ensure the following environment variables are set:
```bash
export OMP_CANCELLATION=TRUE
export OMP_PROC_BIND=TRUE
```

## Generate shared library

SPRAL must be compiled with `-fPIC` to be able to generate a shared library.
The static library `libspral.a` can be converted to a shared library using one of the following commands:
```bash
# Linux
gfortran -fPIC -shared -Wl,--whole-archive libspral.a -Wl,--no-whole-archive -lgomp -lblas -llapack -lhwloc -lmetis -lstdc++ -o libspral.so

# Windows
gfortran -fPIC -shared -Wl,--whole-archive libspral.a -Wl,--no-whole-archive -lgomp -lopenblas -lhwloc -lmetis -lstdc++ -o libspral.dll

# Mac
gfortran -fPIC -shared -Wl,-all_load libspral.a -Wl,-noall_load -lgomp -lopenblas -lhwloc -lmetis -lstdc++ -o libspral.dylib
```

## Citing
To cite SSIDS, please use the following reference:

J. Hogg, E. Ovtchinnikov, and J. Scott, (2016). A sparse symmetric indefinite direct solver for GPU architectures. ACM Transactions on Mathematical Software (TOMS), 42(1), 1-25, https://dx.doi.org/10.1145/275654

In BibTeX, the citation is:

@article{hogg2016sparse,
  title={A sparse symmetric indefinite direct solver for GPU architectures},
  author={Hogg, Jonathan D and Ovtchinnikov, Evgueni and Scott, Jennifer A},
  journal={ACM Transactions on Mathematical Software (TOMS)},
  volume={42},
  number={1},
  pages={1--25},
  year={2016},
  publisher={ACM New York, NY, USA}
}