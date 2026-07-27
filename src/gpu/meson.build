# AMD ROCm backend: the shim exports the CUDA-named runtime entry points
# (cudaMalloc, ...) as hip* wrappers so the Fortran iso_c_binding layer links
# unchanged. Not needed on the CUDA backend (those symbols come from libcudart).
if build_gpu_amd
  libspral_hip_src += files('cuda_hip_shim.cxx')
endif
