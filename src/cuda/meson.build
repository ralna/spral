if build_gpu and nvcc_available and libcuda.found()
  libspral_src += files('cuda.f90')
  libspral_nvcc_src += files('api_wrappers.cu')
else
  libspral_src += files('cuda_nocuda.f90')
endif
