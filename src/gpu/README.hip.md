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

## Building for AMD

```sh
meson setup build -Dgpu=true -Dgpu_backend=amd \
      -Drocm_path=/opt/rocm -Dgpu_arch=gfx90a   # gfx942 MI300, gfx906 ...
ninja -C build
```

The meson wiring (`meson.build`, `meson_options.txt`, and the `src/**/meson.build`
gates) detects `hipcc`, finds `amdhip64`/`hipblas` under `rocm_path`, compiles the
`.cu` + shim with `hipcc` (each via a `custom_target`, `--offload-arch=<gpu_arch>`)
and links the objects into `libspral`. The NVIDIA path is unchanged and selected
with `-Dgpu_backend=cuda` (the default). GPU discovery for the hardware topology
uses `hipGetDeviceCount` (see `src/hw_topology`); precise GPU<->NUMA affinity is
used when hwloc is built with the ROCm SMI backend (`HAVE_HWLOC_RSMI`), otherwise
all visible GPUs are attached to the first NUMA region.

Verified here (ROCm 7.2, gfx906): `meson setup -Dgpu_backend=amd` + `ninja`
compile all six kernels + the shim to device objects **and link the full
`libspral.so`** end-to-end (host compilers gcc/gfortran).

> **Build tip.** Use gcc/gfortran as the host compilers (`CC=gcc CXX=g++
> FC=gfortran`). With ROCm's clang as the host C++ compiler, OpenMP linking
> currently fails (`__kmpc_*` undefined) because meson links `-lgomp`; that is a
> pre-existing clang-on-Linux OpenMP quirk, not specific to the AMD backend.

**Not yet done / not yet validated** — everything below needs actual AMD
hardware, which was not available here (work so far is configure + compile only):

1. **Run-time correctness.** No SSIDS factorise/solve has been executed on a GPU.
2. **Warp-synchronous trsv (fixed, needs on-device check).** `dblkSolve` /
   `dblkSolve_trans` in `dtrsv.h` shared a `volatile` scalar across lanes
   relying on 32-lane lockstep; they now use explicit `SPRAL_SYNCWARP` barriers
   (correct on AMD wavefronts and NVIDIA Volta+). The fix compiles on both
   backends but its numerical result must still be verified on hardware.
   A full audit of `dtrsv.h` found the rest properly synchronised (`slvinv`,
   `slvinv_trans`, `tocache`, `transpose`, `nextRow` all use `__syncthreads`),
   **except `slv21`** (the inner forward-substitution loop, ~L287-294) which has
   the *same* warp-synchronous `xs` sharing. It is deliberately left unpatched:
   its threads diverge (per-`y` `continue`, early `return`), so a wrong
   `__syncwarp` mask would deadlock — it needs a careful rework validated on
   hardware, not a blind barrier.
3. **Inter-block synchronisation (unresolved).** The batched trsv and the
   assembly kernels use spin-locks (`while(sync[...] < ...)`) with
   `__threadfence*` + `atomicAdd`. This assumes all blocks are co-resident
   (forward progress) and the HIP/AMD memory model differs from CUDA's; it must
   be validated on-device and may need a cooperative-launch rewrite. Note also
   `dtrsv.h` mixes `__threadfence_system()` and `__threadfence()` for the same
   pattern (device scope suffices for single-GPU).
4. **hipBLAS semantics.** `spral_cublasDgemm` passes `alpha`/`beta` as host
   pointers — verify hipBLAS pointer-mode default matches.
5. **`__launch_bounds__` re-tuning.** The per-backend `SPRAL_LAUNCH_BOUNDS` macro
   drops the NVIDIA-tuned min-blocks hint on AMD (so the compiler picks
   occupancy); tuning proper AMD values needs profiling on the target arch.
6. **`cudaDeviceSetSharedMemConfig`** is a no-op on AMD (LDS has no configurable
   bank width) — harmless, but confirm no perf assumption depends on it.

## Bugs fixed during code review

- **int64 index truncation** (would corrupt large factorizations, pre-existing):
  `nnz` in `assemble.cu::cu_load_nodes_sc`, and `offc`/`offa` in
  `syrk.cu` (`multisyrk_type::offc` was read into an `int`).
- **Warp-synchronous races** in `dtrsv.h` `dblkSolve`/`dblkSolve_trans` — see (2).
- **`SM_3X` relied on `__CUDA_ARCH__`** (undefined under hipcc, silently 0): now
  selected explicitly per backend in `syrk.cu` and `reorder.cu`.
- **`min`/`max` macros** guarded with `#ifndef` to avoid clashing with the
  toolchain's `std::min`/`std::max`.

## Manual compile check

```sh
ROCM=/opt/rocm
hipcc -x hip -c src/ssids/gpu/kernels/solve.cu -o /tmp/solve.o \
      -Isrc -I$ROCM/include
```
