if build_gpu_dev
  libspral_src += files('cuda.f90')
  if build_gpu_cuda
    libspral_nvcc_src += files('api_wrappers.cu')
  else
    libspral_hip_src += files('api_wrappers.cu')
  endif
else
  libspral_src += files('cuda_nocuda.f90')
endif
