/** \file
 *  \copyright 2026 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *
 *  \brief
 *  Single-source GPU runtime compatibility layer for SPRAL/SSIDS.
 *
 *  The SSIDS device code was written against the CUDA runtime and cuBLAS.
 *  This header lets the very same .cu / .cxx sources compile unchanged under
 *  either nvcc (NVIDIA) or hipcc (AMD ROCm): when built for HIP it pulls in
 *  the HIP runtime and hipBLAS and aliases the handful of CUDA symbols the
 *  device code actually uses onto their HIP equivalents; otherwise it just
 *  includes the genuine CUDA headers and is a no-op.
 *
 *  Scope: this only covers the CUDA surface exercised by the SSIDS kernels
 *  and their C wrappers (see src/ssids/gpu/kernels and src/cuda). It is
 *  deliberately NOT a general CUDA->HIP shim.
 *
 *  NOTE (WIP): compile-validated under hipcc (ROCm 7.2) and nvcc (CUDA 12),
 *  but not yet run on AMD hardware. See src/gpu/README.hip.md.
 */
#pragma once

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIP__) || defined(SPRAL_USE_HIP)

/* ------------------------------------------------------------------ *
 *  AMD / HIP build
 * ------------------------------------------------------------------ */
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hip/hip_complex.h>
#include <hipblas/hipblas.h>

/* --- runtime: error handling --- */
using cudaError = hipError_t;
using cudaError_t = hipError_t;
#define cudaSuccess                    hipSuccess
#define cudaErrorInsufficientDriver    hipErrorInsufficientDriver
#define cudaErrorNoDevice              hipErrorNoDevice
#define cudaGetLastError               hipGetLastError
#define cudaGetErrorString             hipGetErrorString
#define cudaThreadSynchronize          hipDeviceSynchronize
#define cudaDeviceSynchronize          hipDeviceSynchronize

/* --- runtime: device management --- */
#define cudaSetDevice                  hipSetDevice
#define cudaGetDeviceCount             hipGetDeviceCount
#define cudaDeviceEnablePeerAccess     hipDeviceEnablePeerAccess
#define cudaMemGetInfo                 hipMemGetInfo

/* --- runtime: memory --- */
#define cudaMalloc                     hipMalloc
#define cudaFree                       hipFree
#define cudaMemset                     hipMemset
#define cudaMemsetAsync                hipMemsetAsync
#define cudaMemcpy                     hipMemcpy
#define cudaMemcpy2D                   hipMemcpy2D
#define cudaMemcpyAsync                hipMemcpyAsync
#define cudaMemcpy2DAsync              hipMemcpy2DAsync

#define cudaMemcpyKind                 hipMemcpyKind /* #define: used as `enum cudaMemcpyKind` */
#define cudaMemcpyHostToHost           hipMemcpyHostToHost
#define cudaMemcpyHostToDevice         hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost         hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice       hipMemcpyDeviceToDevice
#define cudaMemcpyDefault              hipMemcpyDefault

/* --- runtime: streams --- */
using cudaStream_t = hipStream_t;
#define cudaStreamCreate               hipStreamCreate
#define cudaStreamDestroy              hipStreamDestroy
#define cudaStreamSynchronize          hipStreamSynchronize

/* --- runtime: events --- */
using cudaEvent_t = hipEvent_t;
#define cudaEventCreateWithFlags       hipEventCreateWithFlags
#define cudaEventDestroy               hipEventDestroy
#define cudaEventRecord                hipEventRecord
#define cudaEventSynchronize           hipEventSynchronize
#define cudaEventDefault               hipEventDefault
#define cudaEventBlockingSync          hipEventBlockingSync
#define cudaEventDisableTiming         hipEventDisableTiming

/* --- shared-memory bank config: no AMD analogue (LDS has no configurable
 *     bank width), so these become harmless no-ops. --- */
using cudaSharedMemConfig = unsigned int;
#define cudaSharedMemBankSizeDefault   0u
#define cudaSharedMemBankSizeFourByte  1u
#define cudaSharedMemBankSizeEightByte 2u
static inline hipError_t cudaDeviceSetSharedMemConfig(unsigned int) {
   return hipSuccess;
}
static inline hipError_t cudaDeviceGetSharedMemConfig(unsigned int *cfg) {
   if (cfg) *cfg = cudaSharedMemBankSizeDefault;
   return hipSuccess;
}

/* --- complex --- */
using cuDoubleComplex = hipDoubleComplex;
using cuComplex       = hipComplex;

/* --- cuBLAS -> hipBLAS (source-compatible layer) --- */
using cublasHandle_t    = hipblasHandle_t;
using cublasStatus_t    = hipblasStatus_t;
using cublasOperation_t = hipblasOperation_t;
#define CUBLAS_STATUS_SUCCESS          HIPBLAS_STATUS_SUCCESS
#define CUBLAS_OP_N                    HIPBLAS_OP_N
#define CUBLAS_OP_T                    HIPBLAS_OP_T
#define CUBLAS_OP_C                    HIPBLAS_OP_C
#define cublasCreate                   hipblasCreate
#define cublasDestroy                  hipblasDestroy
#define cublasSetStream                hipblasSetStream
#define cublasDgemm                    hipblasDgemm
#define cublasDsyrk                    hipblasDsyrk

#else /* ---------------- NVIDIA / CUDA build ---------------- */

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuComplex.h>
#include <cublas_v2.h>

#endif

/* Per-backend kernel launch bounds. NVIDIA honours the
 * (maxThreads, minBlocksPerSM) occupancy hint, tuned for its SMs. On AMD that
 * second term maps to waves-per-EU and the NVIDIA-tuned values misfire (the
 * ROCm compiler warns it cannot meet them), so drop the min-blocks term and
 * let the compiler choose occupancy. Re-tuning per AMD arch is future work. */
#if defined(__HIP_PLATFORM_AMD__) || defined(__HIP__) || defined(SPRAL_USE_HIP)
#define SPRAL_LAUNCH_BOUNDS(maxThreads, minBlocks) __launch_bounds__(maxThreads)
#else
#define SPRAL_LAUNCH_BOUNDS(maxThreads, minBlocks) \
   __launch_bounds__(maxThreads, minBlocks)
#endif

/* Intra-warp barrier over a lane mask. HIP requires a 64-bit mask (wavefronts
 * are up to 64 lanes); CUDA uses a 32-bit mask. Callers build the mask as a
 * 64-bit value and this casts to the width the backend expects. */
#if defined(__HIP_PLATFORM_AMD__) || defined(__HIP__) || defined(SPRAL_USE_HIP)
#define SPRAL_SYNCWARP(mask) __syncwarp((unsigned long long)(mask))
#else
#define SPRAL_SYNCWARP(mask) __syncwarp((unsigned)(mask))
#endif
