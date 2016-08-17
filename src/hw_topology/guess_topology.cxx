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

#include <omp.h>

using namespace spral::hw_topology;

/**
 * \brief Guess hardware topology (FIXME: using hwloc if available)
 * \param nregions Number of regions.
 * \param regions[nregions] Array of region descriptors, allocated by this
 *        routine. To free, call spral_hw_topology_free().
 */
extern "C"
void spral_hw_topology_guess(int* nregions, NumaRegion** regions) {
   *nregions = 1;
   *regions = new NumaRegion[*nregions];
   for(int i=0; i<*nregions; ++i) {
      NumaRegion& region = (*regions)[i];
      region.nproc = omp_get_max_threads();
      region.ngpu = 0;
      region.gpus = nullptr;
   }
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
