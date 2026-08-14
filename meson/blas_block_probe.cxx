/* blas_block_probe.cxx
 *
 * Configure-time probe that chooses the SSIDS front-size-adaptive block-size
 * floor (SPRAL_SSIDS_BLOCK_MIN) from this machine's BLAS. It measures the
 * single-threaded double-precision GEMM rate for a range of small square tile
 * sizes and reports the smallest tile at which the BLAS runs near peak -- i.e.
 * the smallest tile for which the tiled dense factorization is BLAS-efficient.
 * Below that size the kernels run below peak and per-task/per-call overhead
 * dominates, so it is a sensible lower clamp for the adaptive ramp.
 *
 * The program prints a single integer (a multiple of 32, clamped to a sane
 * range) on stdout; the Meson build captures it via compiler.run(). It is
 * deliberately self-contained and quick (well under a second). If anything is
 * unexpected it prints a conservative fallback and still exits 0, so the build
 * is never blocked by the probe.
 *
 * Note: because it times the machine at configure time, the chosen value can
 * vary slightly between configurations run under different load. Builds that
 * need a fixed, reproducible floor should set the Meson option
 * -Dssids_block_min=<value> to a positive number, which bypasses this probe.
 */
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <algorithm>

// Reference-BLAS Fortran symbol (OpenBLAS, Netlib, Accelerate, MKL all export
// this). Column-major, no transpose.
extern "C" void dgemm_(const char* transa, const char* transb,
                       const int* m, const int* n, const int* k,
                       const double* alpha, const double* a, const int* lda,
                       const double* b, const int* ldb,
                       const double* beta, double* c, const int* ldc);

namespace {

// GEMM rate (flop/s) for a square b-by-b-by-b product, best of a few trials.
double gemm_rate(int b) {
   std::vector<double> A(static_cast<size_t>(b) * b, 1.0);
   std::vector<double> B(static_cast<size_t>(b) * b, 1.0);
   std::vector<double> C(static_cast<size_t>(b) * b, 0.0);
   const char N = 'N';
   const double one = 1.0, zero = 0.0;
   const double work = 2.0 * b * b * b; // flop per gemm

   // Enough repetitions that each timed burst lasts a few milliseconds.
   int reps = static_cast<int>(3.0e8 / work);
   if (reps < 5) reps = 5;

   // Warm up (first call may allocate BLAS buffers / spin up threads).
   dgemm_(&N, &N, &b, &b, &b, &one, A.data(), &b, B.data(), &b, &zero, C.data(), &b);

   double best = 0.0;
   for (int trial = 0; trial < 3; ++trial) {
      auto t0 = std::chrono::steady_clock::now();
      for (int r = 0; r < reps; ++r)
         dgemm_(&N, &N, &b, &b, &b, &one, A.data(), &b, B.data(), &b,
                &zero, C.data(), &b);
      auto t1 = std::chrono::steady_clock::now();
      double secs = std::chrono::duration<double>(t1 - t0).count();
      if (secs > 0.0) best = std::max(best, work * reps / secs);
   }
   return best;
}

} // anonymous namespace

int main() {
   // Measure single-threaded BLAS: the tiled kernels each run on one thread
   // inside the task scheduler, so single-thread efficiency is what sets the
   // useful floor. Set before the first BLAS call so the runtime picks it up.
   setenv("OPENBLAS_NUM_THREADS", "1", 1);
   setenv("OMP_NUM_THREADS", "1", 1);
   setenv("MKL_NUM_THREADS", "1", 1);
   setenv("VECLIB_MAXIMUM_THREADS", "1", 1); // Apple Accelerate

   const int gran = 32;   // report a multiple of the LDLT inner block size
   const int lo   = 64;   // never floor below this
   const int hi   = 256;  // never floor above this
   const double frac = 0.90; // "near peak" threshold

   // Sample square tile sizes from small to a size comfortably past the knee.
   std::vector<int> sizes;
   for (int b = 32; b <= 384; b += 32) sizes.push_back(b);

   double peak = 0.0;
   std::vector<double> rate(sizes.size(), 0.0);
   for (size_t i = 0; i < sizes.size(); ++i) {
      rate[i] = gemm_rate(sizes[i]);
      peak = std::max(peak, rate[i]);
   }

   // Smallest tile reaching frac of peak.
   int chosen = hi;
   if (peak > 0.0) {
      for (size_t i = 0; i < sizes.size(); ++i) {
         if (rate[i] >= frac * peak) { chosen = sizes[i]; break; }
      }
   } else {
      chosen = 128; // BLAS produced no timing; fall back conservatively
   }

   // Round to a multiple of the granularity and clamp to [lo, hi].
   chosen = ((chosen + gran / 2) / gran) * gran;
   if (chosen < lo) chosen = lo;
   if (chosen > hi) chosen = hi;

   std::printf("%d\n", chosen);
   return 0;
}
