/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 * This file exists only for Doxygen documentation purposes and should never
 * be included for compilation.
 */
#error Documentation only

/** \brief Sparse Parallel Robust Algorithms Library
 *  \details Parent namespace for all SPRAL modules.
 */
namespace spral {
   /** \brief Sparse Symmetric Indefinite Direct Solver */
   namespace ssids {
      /** \brief CPU-specific code */
      namespace cpu {}
      /** \brief GPU-specific code */
      namespace gpu {}
   } /* namespace ssids */
} /* namespace spral */
