subdir('kernels')

if build_gpu and nvcc_available and libcuda.found()
  libspral_src += files('alloc.f90',
                        'cpu_solve.f90',
                        'datatypes.f90',
                        'dense_factor.f90',
                        'factor.f90',
                        'interfaces.f90',
                        'smalloc.f90',
                        'solve.f90',
                        'subtree.f90')
else
  libspral_src += files('subtree_no_cuda.f90')
endif
