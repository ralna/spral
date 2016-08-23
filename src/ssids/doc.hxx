/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 *
 *  \brief Doxygen documentation
 *
 *  \warning This file exists only for Doxygen documentation purposes and should
 *  never be included for compilation.
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
