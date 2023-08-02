## Dependencies
### Software Directory
First, create a directory where the dependencies will be compiled, e.g.,
```bash
mkdir -p ${HOME}/Software
```
The remainder of this guide assumes such a directory has been created.

### METIS
Next, compile the [COIN-OR version of METIS](https://github.com/coin-or-tools/ThirdParty-Metis), e.g.,
```bash
cd ${HOME}/Software
git clone https://github.com/coin-or-tools/ThirdParty-Metis.git
cd ThirdParty-Metis && ./get.Metis
mkdir build
cd build
../configure --prefix=${PWD}
make && make install
export METISDIR=${PWD}
```

### hwloc
Next, ensure the [hardware locality library](https://www.open-mpi.org/projects/hwloc/) is installed on your system.
For example, on Ubuntu 20.04 LTS, this can be accomplished via
```bash
sudo apt-get install hwloc libhwloc-dev
```

### BLAS
A number of BLAS libraries can be used for compilation, including the proprietary [Intel MKL](https://software.intel.com/en-us/mkl).
One open-source option is [OpenBLAS](https://www.openblas.net), which can be installed on Ubuntu 20.04 LTS via
```bash
sudo apt install libopenblas-dev
```
This is the library used throughout the remainder of this guide.

### CUDA (optional)
If you're installing with [NVIDIA CUDA](https://developer.nvidia.com/cuda-downloads) GPU support, ensure it is installed and that the following environment variables are set:
```
export CUDA_HOME="/usr/local/cuda" # Change this to your system-specific path.
export PATH="${PATH}:${CUDA_HOME}/bin"
export LIBRARY_PATH="${LIBRARY_PATH}:${CUDA_HOME}/lib64"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${CUDA_HOME}/lib64"
export C_INCLUDE_PATH="${CPLUS_INCLUDE_PATH}:${CUDA_HOME}/include"
export CPLUS_INCLUDE_PATH="${CPLUS_INCLUDE_PATH}:${CUDA_HOME}/include"
export NVCC_INCLUDE_FLAGS="${NVCC_INCLUDE_FLAGS}:-I${CUDA_HOME}/include"
```

## Compilation
### Multicore CPUs Only
To compile with only multicore CPU support, execute
```bash
cd ${HOME}/Software
git clone https://github.com/ralna/spral.git
cd spral
./autogen.sh # If compiling from scratch.
mkdir build
cd build
CFLAGS=-fPIC CPPFLAGS=-fPIC CXXFLAGS=-fPIC FFLAGS=-fPIC \
   FCFLAGS=-fPIC ../configure --prefix=${PWD}/build \
   --with-blas="-lopenblas" --with-lapack="-llapack" \
   --with-metis="-L${METISDIR}/lib -lcoinmetis" \
   --with-metis-inc-dir="${METISDIR}/include/coin-or/metis"
make && make install
```

### Multicore CPUs and NVIDIA GPUs (optional)
To compile with multicore CPU and NVIDIA GPU support, execute
```bash
cd ${HOME}/Software
git clone https://github.com/ralna/spral.git
cd spral
./autogen.sh # If compiling from scratch.
mkdir build
cd build
CFLAGS=-fPIC CPPFLAGS=-fPIC CXXFLAGS=-fPIC FFLAGS=-fPIC \
   FCFLAGS=-fPIC NVCCFLAGS="-shared -Xcompiler -fPIC" \
   ../configure --prefix=${PWD}/build \
   --with-blas="-lopenblas" --with-lapack="-llapack" \
   --with-metis="-L${METISDIR}/lib -lcoinmetis" \
   --with-metis-inc-dir="${METISDIR}/include/coin-or/metis"
make && make install
```
If your GPUs are not being recognized, consider uncommenting line 91 and commenting line 92 of `src/hw_topology/hwloc_wrapper.hxx`, then recompiling.
The methods for recognizing GPUs do not seem to function independently of the type of system being used.

## Usage
For future use, set the SPRAL directory environment variable via
```bash
export SPRALDIR=${PWD}/build
```
from the `${HOME}/Software/spral` directory.
Also, ensure the following environment variables are set when using the library:
```bash
export OMP_CANCELLATION=TRUE
export OMP_PROC_BIND=TRUE
```
