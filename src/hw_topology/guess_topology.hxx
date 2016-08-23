/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 *
 * \brief
 * Defines NumaRegion struct.
 */
#pragma once

namespace spral {
/** \brief Hardware topology module */
namespace hw_topology {

struct NumaRegion {
   int nproc;
   int ngpu;
   int *gpus;
};


}} /* namespace spral::hw_topology */
