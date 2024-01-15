if build_gpu and nvcc_available and libcuda.found()
  binspral_src += files('cuda_helper_gpu.f90')
else
  binspral_src += files('cuda_helper_nogpu.f90')
endif

binspral_src += files('spral_ssids.F90')
