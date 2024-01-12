subdir('cuda')
subdir('hw_topology')
subdir('ssids')
subdir('ssmfe')

libspral_src += files('blas_iface.f90',
                      'core_analyse.f90',
                      'lapack_iface.f90',
                      'lsmr.f90',
                      'match_order.f90',
                      'matrix_util.f90',
                      'pgm.f90',
                      'random.f90',
                      'random_matrix.f90',
                      'rutherford_boeing.f90',
                      'scaling.f90',
                      'timer.f90')

libspral_cpp_src += files('compat.cxx',
                          'omp.cxx')

if libmetis_version == '4'
  libspral_src += files('metis4_wrapper.F90')
else
  libspral_src += files('metis5_wrapper.F90')
endif
