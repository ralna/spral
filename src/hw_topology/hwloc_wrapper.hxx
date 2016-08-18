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
#pragma once
#include "config.h"
#ifdef HAVE_HWLOC

#include <vector>

#include <hwloc.h>

namespace spral { namespace hw_topology {

class HwlocTopology {
public:
   // \{
   // Not copyable
   HwlocTopology(HwlocTopology const&) =delete;
   HwlocTopology operator=(HwlocTopology const&) =delete;
   // \}
   /** \brief Constructor */
   HwlocTopology() {
      hwloc_topology_init(&topology_);
      hwloc_topology_load(topology_);
   }
   /** \brief Destructor */
   ~HwlocTopology() {
      hwloc_topology_destroy(topology_);
   }

   /** \brief Return vector of Numa nodes or just machine object */
   std::vector<hwloc_obj_t> get_numa_nodes() const {
      std::vector<hwloc_obj_t> regions;
      int nregions = hwloc_get_nbobjs_by_type(topology_, HWLOC_OBJ_NODE);
      if(nregions==0) {
         // No regions, just give machine
         regions.push_back(
               hwloc_get_obj_by_type(topology_, HWLOC_OBJ_MACHINE, 0)
               );
         return regions;
      } else {
         // Iterate over regions, adding them
         regions.reserve(nregions);
         for(int i=0; i<nregions; ++i)
            regions.push_back(
                  hwloc_get_obj_by_type(topology_, HWLOC_OBJ_NODE, i)
                  );
         return regions;
      }
   }

   /** \brief Return number of cores associated to object. */
   int count_cores(hwloc_obj_t const& obj) {
      if(obj->type == HWLOC_OBJ_CORE) return 1;
      int count = 0;
      for(unsigned int i=0; i<obj->arity; ++i)
         count += count_cores(obj->children[i]);
      return count;
   }
private:
   hwloc_topology_t topology_; ///< Underlying topology object
};

}} /* namespace spral::hw_topology */

#endif /* HAVE_HWLOC */
