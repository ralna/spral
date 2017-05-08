/* nvcc -DNDEBUG nvcc_arch_sm.c -o nvcc_arch_sm -lcuda */
/* If nvcc is not functioning properly, use the native compiler; e.g.,
 * on macOS Sierra with CUDA 8 and the latest Clang releases, call:
 * clang -DNDEBUG -I$CUDADIR/include nvcc_arch_sm.c -o nvcc_arch_sm -L$CUDADIR/lib -lcuda */

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

int main(int argc, char *argv[])
{
  CUresult result = CUDA_SUCCESS;
  int count = -1;
  int *major = (int*)NULL;
  int *minor = (int*)NULL;
  CUdevice dev = (CUdevice)-1;
  int major_min = 0;
  int major_max = 0;
  int minor_min = 0;
  int minor_max = 0;
  int chosen = 0;

  switch (result = cuInit(0u)) {
  case CUDA_ERROR_INVALID_VALUE:
    (void)fprintf(stderr, "cuInit: CUDA_ERROR_INVALID_VALUE\n");
    return EXIT_FAILURE;
  case CUDA_ERROR_INVALID_DEVICE:
    (void)fprintf(stderr, "cuInit: CUDA_ERROR_INVALID_DEVICE\n");
    return EXIT_FAILURE;
  case CUDA_SUCCESS:
    break;
  default:
    (void)fprintf(stderr, "cuInit: %d\n", result);
    return EXIT_FAILURE;
  }

  switch (result = cuDeviceGetCount(&count)) {
  case CUDA_ERROR_DEINITIALIZED:
    (void)fprintf(stderr, "cuDeviceGetCount: CUDA_ERROR_DEINITIALIZED\n");
    return EXIT_FAILURE;
  case CUDA_ERROR_NOT_INITIALIZED:
    (void)fprintf(stderr, "cuDeviceGetCount: CUDA_ERROR_NOT_INITIALIZED\n");
    return EXIT_FAILURE;
  case CUDA_ERROR_INVALID_CONTEXT:
    (void)fprintf(stderr, "cuDeviceGetCount: CUDA_ERROR_INVALID_CONTEXT\n");
    return EXIT_FAILURE;
  case CUDA_ERROR_INVALID_VALUE:
    (void)fprintf(stderr, "cuDeviceGetCount: CUDA_ERROR_INVALID_VALUE\n");
    return EXIT_FAILURE;
  case CUDA_SUCCESS:
    if (count < 0) {
      (void)fprintf(stderr, "cuDeviceGetCount = %d\n", count);
      return EXIT_FAILURE;
    }
    break;
  default:
    (void)fprintf(stderr, "cuDeviceGetCount: %d\n", result);
    return EXIT_FAILURE;
  }
#ifndef NDEBUG
  (void)fprintf(stdout, "cuDeviceGetCount = %d\n", count);
  (void)fflush(stdout);
#endif /* !NDEBUG */
  if (!count) /* signal no GPUs by returning the failure code */
    return EXIT_FAILURE;

  if (!(major = (int*)calloc((size_t)(2 * count), sizeof(int))))
    return EXIT_FAILURE;
  minor = major + count;

  for (dev = 0; dev < count; ++dev) {
    switch (result = cuDeviceGetAttribute((major + dev), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, dev)) {
    case CUDA_ERROR_DEINITIALIZED:
      (void)fprintf(stderr, "cuDeviceGetAttribute(major): CUDA_ERROR_DEINITIALIZED\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_NOT_INITIALIZED:
      (void)fprintf(stderr, "cuDeviceGetAttribute(major): CUDA_ERROR_NOT_INITIALIZED\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_INVALID_CONTEXT:
      (void)fprintf(stderr, "cuDeviceGetAttribute(major): CUDA_ERROR_INVALID_CONTEXT\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_INVALID_VALUE:
      (void)fprintf(stderr, "cuDeviceGetAttribute(major): CUDA_ERROR_INVALID_VALUE\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_INVALID_DEVICE:
      (void)fprintf(stderr, "cuDeviceGetAttribute(major): CUDA_ERROR_INVALID_DEVICE\n");
      return EXIT_FAILURE;
    case CUDA_SUCCESS:
      if (major[dev] <= 0) {
        (void)fprintf(stderr, "cuDeviceGetAttribute(major) = %d\n", major[dev]);
        return EXIT_FAILURE;
      }
      break;
    default:
      (void)fprintf(stderr, "cuDeviceGetAttribute(major): %d\n", result);
      return EXIT_FAILURE;
    }
#ifndef NDEBUG
    (void)fprintf(stdout, " %d=%d.", dev, major[dev]);
    (void)fflush(stdout);
#endif /* !NDEBUG */
    switch (result = cuDeviceGetAttribute((minor + dev), CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, dev)) {
    case CUDA_ERROR_DEINITIALIZED:
      (void)fprintf(stderr, "cuDeviceGetAttribute(minor): CUDA_ERROR_DEINITIALIZED\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_NOT_INITIALIZED:
      (void)fprintf(stderr, "cuDeviceGetAttribute(minor): CUDA_ERROR_NOT_INITIALIZED\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_INVALID_CONTEXT:
      (void)fprintf(stderr, "cuDeviceGetAttribute(minor): CUDA_ERROR_INVALID_CONTEXT\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_INVALID_VALUE:
      (void)fprintf(stderr, "cuDeviceGetAttribute(minor): CUDA_ERROR_INVALID_VALUE\n");
      return EXIT_FAILURE;
    case CUDA_ERROR_INVALID_DEVICE:
      (void)fprintf(stderr, "cuDeviceGetAttribute(minor): CUDA_ERROR_INVALID_DEVICE\n");
      return EXIT_FAILURE;
    case CUDA_SUCCESS:
      if (minor[dev] < 0) {
        (void)fprintf(stderr, "cuDeviceGetAttribute(minor) = %d\n", minor[dev]);
        return EXIT_FAILURE;
      }
      break;
    default:
      (void)fprintf(stderr, "cuDeviceGetAttribute(minor): %d\n", result);
      return EXIT_FAILURE;
    }
#ifndef NDEBUG
    (void)fprintf(stdout, "%d", minor[dev]);
    (void)fflush(stdout);
#endif /* !NDEBUG */
  }
#ifndef NDEBUG
  (void)fprintf(stdout, "\n");
  (void)fflush(stdout);
#endif /* !NDEBUG */

  major_max = major_min = major[0];
  minor_max = minor_min = minor[0];
  for (dev = 1; dev < count; ++dev) {
    if (major[dev] < major_min)
      major_min = major[dev];
    if (major[dev] > major_max)
      major_max = major[dev];
    if (minor[dev] < minor_min)
      minor_min = minor[dev];
    if (minor[dev] > minor_max)
      minor_max = minor[dev];
  }

#if defined(__CUDA_API_VERSION) && __CUDA_API_VERSION >= 8000
  /* CC <= 1.x not supported in CUDA >= 8. */
  if (major_max <= 1) /* signal no GPUs by returning the failure code */
    return EXIT_FAILURE;
#endif /* __CUDA_API_VERSION */

  if (major_max == major_min) {
    /* uniform GPUs */
    if (minor_max == minor_min) {
      (void)fprintf(stdout, "-arch=sm_%d%d\n", major_max, minor_max);
      return EXIT_SUCCESS;
    }
    else if (major_max == 2) {
      (void)fprintf(stdout, "-arch=compute_20 -code=sm_20,sm_21\n");
      return EXIT_SUCCESS;
    }
  }

  for (dev = 0; dev < count; ++dev) {
    int seen = 0;
    CUdevice prev = (CUdevice)-1;
#if defined(__CUDA_API_VERSION) && __CUDA_API_VERSION >= 8000
    /* CC <= 1.x not supported in CUDA >= 8. */
    if (major[dev] <= 1)
      continue;
#endif /* __CUDA_API_VERSION */
    for (prev = dev - 1; prev >= 0; --prev)
      if (seen = ((major[prev] == major[dev]) && (minor[prev] == minor[dev])))
        break;
    if (seen)
      continue;
    if (chosen++)
      (void)fprintf(stdout, " ");
    (void)fprintf(stdout, "-gencode arch=compute_%d%d,code=sm_%d%d", major[dev], minor[dev], major[dev], minor[dev]);
    (void)fflush(stdout);
  }
  free(major);
  if (chosen) {
    (void)fprintf(stdout, "\n");
    return EXIT_SUCCESS;
  }
  /* signal no GPUs by returning the failure code */
  return EXIT_FAILURE;
}
