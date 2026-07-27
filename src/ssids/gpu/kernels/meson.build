if build_gpu_cuda
  libspral_nvcc_src += files('assemble.cu',
                             'dense_factor.cu',
                             'reorder.cu',
                             'solve.cu',
                             'syrk.cu')
elif build_gpu_amd
  libspral_hip_src += files('assemble.cu',
                            'dense_factor.cu',
                            'reorder.cu',
                            'solve.cu',
                            'syrk.cu')
endif
