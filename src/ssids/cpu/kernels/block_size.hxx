/** \file
 *  \copyright 2026 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *
 *  Front-size-adaptive outer block size for the dense node kernels.
 */
#pragma once

#include <algorithm>
#include <cstdint>

/* Tuning constants for the front-size-adaptive block size. These are fixed at
 * build time via the Meson options ssids_block_div / ssids_block_tiles_per_thread
 * / ssids_block_min / ssids_block_max (which the build system turns into the
 * macros below). SPRAL_SSIDS_BLOCK_MIN in particular is, by default, determined
 * by a configure-time BLAS efficiency probe (see meson/blas_block_probe.cxx).
 * The fallback defaults here keep the header usable in a standalone compile. */
#ifndef SPRAL_SSIDS_BLOCK_DIV
#define SPRAL_SSIDS_BLOCK_DIV 0            /* 0 => thread-scaled (see below) */
#endif
#ifndef SPRAL_SSIDS_BLOCK_TILES_PER_THREAD
#define SPRAL_SSIDS_BLOCK_TILES_PER_THREAD 2
#endif
#ifndef SPRAL_SSIDS_BLOCK_MIN
#define SPRAL_SSIDS_BLOCK_MIN 128
#endif
#ifndef SPRAL_SSIDS_BLOCK_MAX
#define SPRAL_SSIDS_BLOCK_MAX 512
#endif

namespace spral { namespace ssids { namespace cpu {

/** Inner block size of the dense LDLT kernel: the granularity of its recursive
 *  diagonal-block factorization, and the unit to which adaptive_block_size()
 *  rounds. Defined here so the LDLT (ldlt_app.cxx) and Cholesky (cholesky.cxx)
 *  kernels share a single definition rather than each hard-coding 32. */
static const int INNER_BLOCK_SIZE = 32;

/** Front-size-adaptive outer block size.
 *
 *  The optimal block size for the dense factorization of a supernode is not a
 *  global constant: it tracks both the front's size and how many threads are
 *  working it.
 *
 *  - Absolute clamps [MIN,MAX] bound single-block efficiency (BLAS shape,
 *    cache footprint, per-task overhead) and are independent of thread count.
 *    MIN is, by default, chosen at build time by a BLAS efficiency probe: the
 *    smallest square tile at which this machine's BLAS runs near peak.
 *  - The ramp rate governs how finely a front is tiled, i.e. how much parallel
 *    work it exposes -- so it scales with the number of threads. We aim for
 *    roughly TILES_PER_THREAD block-rows per thread, giving a divisor of
 *    TILES_PER_THREAD * nthreads. Small fronts and/or many threads therefore
 *    pull the block size down (more tiles to keep everyone busy) until the MIN
 *    floor; very large fronts and/or few threads push it up to the MAX ceiling.
 *
 *  \param m front row count
 *  \param granularity block size is rounded to a multiple of this; both kernels
 *         pass INNER_BLOCK_SIZE (the LDLT inner block size)
 *  \param nthreads number of threads collaborating on this front, i.e. the size
 *         of the OpenMP team assigned to this subtree by SSIDS' topology model
 *         (captured once at subtree entry). Serial callers should pass 1, for
 *         which the rule selects large tiles. Values \f$\le 0\f$ are treated as
 *         1 defensively.
 *
 *  A build may instead pin a fixed divisor (thread-independent) by setting the
 *  Meson option ssids_block_div to a positive value; 0 (the default) selects
 *  the thread-scaled divisor above.
 */
inline int adaptive_block_size(int m, int granularity, int nthreads) {
   int const min_blk = std::max(granularity, SPRAL_SSIDS_BLOCK_MIN);
   int const max_blk = std::max(min_blk, SPRAL_SSIDS_BLOCK_MAX);

   int const nt = std::max(1, nthreads);

   // Ramp divisor: a positive build-time override, else thread-scaled.
   int const fixed_div = SPRAL_SSIDS_BLOCK_DIV;
   int64_t const div = (fixed_div > 0)
      ? static_cast<int64_t>(fixed_div)
      : static_cast<int64_t>(SPRAL_SSIDS_BLOCK_TILES_PER_THREAD) * nt;

   int b = static_cast<int>(static_cast<int64_t>(m) / div);
   b = ((b + granularity/2) / granularity) * granularity; // round to nearest
   if(b < min_blk) b = min_blk;
   if(b > max_blk) b = max_blk;
   return b;
}

}}} /* namespaces spral::ssids::cpu */
