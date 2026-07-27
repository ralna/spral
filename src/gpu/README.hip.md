# AMD GPU (HIP/ROCm) port — work in progress

This directory holds the compatibility layer that lets the SSIDS device code
(originally CUDA-only) build for AMD GPUs via HIP/ROCm, **single-source**: the
exact same `.cu` sources compile under either `nvcc` or `hipcc`.

## Files

| File | Role |
|------|------|
| `spral_gpu_runtime.hxx` | Single-source header. Under HIP, includes the HIP runtime + hipBLAS and aliases the CUDA symbols used by the kernels onto their `hip*` / `hipblas*` equivalents. Under CUDA, includes the genuine CUDA headers (no-op). |
| `cuda_hip_shim.cxx` | Defines the CUDA-named runtime entry points (`cudaMalloc`, `cudaMemcpy`, ...) as thin `extern "C"` wrappers around `hip*`, so the Fortran `bind(C, name="cuda...")` layer in `src/cuda/cuda.f90` links unchanged. HIP build only. |

The device sources were adapted to include `gpu/spral_gpu_runtime.hxx` instead of
`<cuda_runtime.h>`/`<cublas_v2.h>` directly:
`src/cuda/api_wrappers.cu`, and all of `src/ssids/gpu/kernels/*.cu` plus
`dtrsv.h` (PTX `%smid` → AMD `s_getreg`, guarded; vestigial `<cuComplex.h>`
include guarded).

## Status — compile-validated

Every device translation unit compiles cleanly under **both** toolchains
(no regression to the CUDA build):

```
# ROCm 7.2 / gfx906           # CUDA 12.0
hipcc -x hip -c <file> -Isrc  nvcc -c <file> -Isrc
```
- `src/cuda/api_wrappers.cu`   ✅ hipcc  ✅ nvcc
- `src/gpu/cuda_hip_shim.cxx`  ✅ hipcc  (exports cudaMalloc, cudaMemcpy, ...)
- `src/ssids/gpu/kernels/syrk.cu`         ✅ ✅
- `src/ssids/gpu/kernels/solve.cu`        ✅ ✅
- `src/ssids/gpu/kernels/assemble.cu`     ✅ ✅
- `src/ssids/gpu/kernels/reorder.cu`      ✅ ✅
- `src/ssids/gpu/kernels/dense_factor.cu` ✅ ✅

**Not yet done / not yet validated** — this is compile-level only, no GPU was
available to run on:

1. **Meson wiring.** Add `hipcc`/ROCm detection and a `gpu=amd` option; compile
   the `.cu` with `hipcc`; build+link `cuda_hip_shim.cxx`; link
   `hipblas` + `amdhip64` instead of `cudart`+`cublas`. Nothing here is wired
   into the build yet.
2. **Warp size 32 vs 64 (correctness risk).** `dtrsv.h` fixes
   `TRSV_NB_TASK = 32` "= warpSize" and does warp-synchronous work. On AMD a
   wavefront is 64 lanes; a 32-thread block is half a wavefront. This *may* be
   safe (32 lanes are still lockstep within one wavefront) but must be verified
   on hardware; likewise the reductions in `solve.cu`.
3. **Inter-block synchronisation.** The batched trsv uses spin-locks with
   `__threadfence_system()` + `atomicAdd`. The HIP/AMD memory model differs from
   CUDA's; the busy-wait sync must be validated on-device.
4. **hipBLAS semantics.** `spral_cublasDgemm` passes `alpha`/`beta` as host
   pointers — verify hipBLAS pointer-mode default matches.
5. **`__launch_bounds__` tuning.** Occupancy hints (e.g. `(64,14)`, `(256,8)`)
   were tuned for NVIDIA; re-tune for the target AMD arch.
6. **`cudaDeviceSetSharedMemConfig`** is a no-op on AMD (LDS has no configurable
   bank width) — harmless, but confirm no perf assumption depends on it.

## Manual compile check

```sh
ROCM=/opt/rocm
hipcc -x hip -c src/ssids/gpu/kernels/solve.cu -o /tmp/solve.o \
      -Isrc -I$ROCM/include
```
