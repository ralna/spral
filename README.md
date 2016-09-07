![STFC logo](http://www.stfc.ac.uk/stfc/includes/themes/MuraSTFC/assets/legacy/2473_web_2.png)
# SPRAL: The Sparse Parallel Robust Algorithm Library
An open-source (BSD) library for sparse linear algebra and associated
algorithms. It is primarily developed by the Numerical Analysis group at
STFC Rutherford Appleton Laboratory [hsl@stfc.ac.uk](mailto:hsl@stfc.ac.uk).

## Packages

- **LSMR**          - Solves sparse least squares problems using LSMR algorithm.
- **RANDOM**        - Pseudo-random number generator.
- **RANDOM_MATRIX** - Generates random matrices for testing purposes.
- **SCALING**       - Calculates matrix scalings through a variety of algorithms
- **SSIDS**         - Sparse Symmetric Indefinite Direct Solver. Requires an
                      NVIDIA GPU!
- **SSMFE**         - Sparse Symmetric Matrix-Free Eigensolver. Uses
                      Jacobi-conjugate preconditioned gradients method.

If the functionality you are looking for is not support, it may be offered by
our proprietary licenced [HSL Library](http://www.hsl.rl.ac.uk/)
(free to academics).

## Installation
We use a standard autotools-based build:
```bash
./configure
make
make install
```
