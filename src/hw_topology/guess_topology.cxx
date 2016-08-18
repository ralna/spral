/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 */
/** \file */
#include "guess_topology.hxx"

#include "config.h"
#include <omp.h>
#ifdef HAVE_NVCC
#include <cuda_runtime_api.h>
#endif /* HAVE_NVCC */

#include "hw_topology/hwloc_wrapper.hxx"

#include <cstdio> // debug

using namespace spral::hw_topology;

/**
 * \brief Guess hardware topology (using hwloc if available)
 * \param nregions Number of regions.
 * \param regions[nregions] Array of region descriptors, allocated by this
 *        routine. To free, call spral_hw_topology_free().
 */
extern "C"
void spral_hw_topology_guess(int* nregions, NumaRegion** regions) {
#if HAVE_HWLOC
   // Compiled with hwloc support
   HwlocTopology topology;
   auto numa_nodes = topology.get_numa_nodes();
   *nregions = numa_nodes.size();
   *regions = new NumaRegion[*nregions];
   for(int i=0; i<*nregions; ++i) {
      NumaRegion& region = (*regions)[i];
      region.nproc = topology.count_cores(numa_nodes[i]);
      region.ngpu = 0;
      region.gpus = nullptr;
   }
#else /* HAVE_HWLOC */
   // Compiled without hwloc support, just put everything in one region
   *nregions = 1;
   *regions = new NumaRegion[*nregions];
   NumaRegion& region = (*regions)[0];
   region.nproc = omp_get_max_threads();
#if HAVE_NVCC
   cudaError_t cuda_error = cudaGetDeviceCount(&region.ngpu);
   if(cuda_error != 0) {
      printf("CUDA Failed, working without GPU\n");
      region.ngpu = 0;
   }
   region.gpus = (region.ngpu > 0) ? new int[region.ngpu] : nullptr;
   for(int i=0; i<region.ngpu; ++i)
      region.gpus[i] = i;
#else /* HAVE_NVCC */
   region.ngpu = 0;
   region.gpus = nullptr;
#endif /* HAVE_NVCC */
#endif /* HAVE_HWLOC */
}

/**
 * \brief Free hardware topology allocated by spral_hw_topology_guess().
 * \param nregions Number of regions.
 * \param regions[nregions] Array of region descriptors to free.
 */
extern "C"
void spral_hw_topology_free(int nregions, NumaRegion* regions) {
   for(int i=0; i<nregions; ++i) {
      if(regions[i].gpus)
         delete[] regions[i].gpus;
   }
   delete[] regions;
}
