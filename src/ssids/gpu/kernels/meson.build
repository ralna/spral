if build_gpu and libcuda.found()
  libspral_nvcc_src += files('assemble.cu',
                             'dense_factor.cu',
                             'reorder.cu',
                             'solve.cu',
                             'syrk.cu')
endif
