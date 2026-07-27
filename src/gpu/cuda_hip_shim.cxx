/** \file
 *  \copyright 2026 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *
 *  \brief
 *  CUDA-named entry points for the Fortran iso_c_binding layer, on HIP builds.
 *
 *  The Fortran module spral_cuda (src/cuda/cuda.f90) binds directly, by C name,
 *  to a subset of the CUDA runtime API (bind(C, name="cudaMalloc"), ...). Those
 *  symbols do not exist in the HIP runtime, so on an AMD build we provide them
 *  here as thin extern "C" wrappers around the equivalent hip* calls. The
 *  Fortran layer is then usable verbatim, with no per-symbol #ifdef.
 *
 *  This translation unit deliberately does NOT include spral_gpu_runtime.hxx:
 *  that header #defines cudaMalloc -> hipMalloc etc. for single-source device
 *  code, whereas here we need to *define* symbols literally named cuda*.
 *
 *  Only compiled when building SSIDS for HIP. On a CUDA build these symbols
 *  come from libcudart and this file is not used.
 *
 *  NOTE (WIP): compile-validated under hipcc (ROCm 7.2); not yet run on AMD
 *  hardware. See src/gpu/README.hip.md.
 */

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

extern "C" {

/* --- device management --- */
int cudaSetDevice(int device) {
   return (int) hipSetDevice(device);
}
int cudaGetDeviceCount(int *count) {
   return (int) hipGetDeviceCount(count);
}
int cudaDeviceSynchronize(void) {
   return (int) hipDeviceSynchronize();
}
int cudaDeviceEnablePeerAccess(int peerDevice, int flags) {
   return (int) hipDeviceEnablePeerAccess(peerDevice, (unsigned int) flags);
}
int cudaMemGetInfo(size_t *free, size_t *total) {
   return (int) hipMemGetInfo(free, total);
}

/* --- error handling --- */
int cudaGetLastError(void) {
   return (int) hipGetLastError();
}
const char *cudaGetErrorString(int error) {
   return hipGetErrorString((hipError_t) error);
}

/* --- memory --- */
int cudaMalloc(void **devPtr, size_t size) {
   return (int) hipMalloc(devPtr, size);
}
int cudaFree(void *devPtr) {
   return (int) hipFree(devPtr);
}
int cudaMemset(void *devPtr, int value, size_t count) {
   return (int) hipMemset(devPtr, value, count);
}
int cudaMemcpy(void *dst, const void *src, size_t count, int kind) {
   return (int) hipMemcpy(dst, src, count, (hipMemcpyKind) kind);
}
int cudaMemcpy2D(void *dst, size_t dpitch, const void *src, size_t spitch,
      size_t width, size_t height, int kind) {
   return (int) hipMemcpy2D(dst, dpitch, src, spitch, width, height,
         (hipMemcpyKind) kind);
}

/* --- shared-memory bank config: no AMD analogue, accept and ignore --- */
int cudaDeviceSetSharedMemConfig(int config) {
   (void) config;
   return (int) hipSuccess;
}
int cudaDeviceGetSharedMemConfig(int *pConfig) {
   if (pConfig) *pConfig = 0; /* cudaSharedMemBankSizeDefault */
   return (int) hipSuccess;
}

} /* extern "C" */
