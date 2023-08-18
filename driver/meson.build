if build_gpu and libcuda.found()
  libspral_src += files('cuda_helper_gpu.f90')
else
  libspral_src += files('cuda_helper_nogpu.f90')
endif

libspral_src += files('spral_ssids.F90')
