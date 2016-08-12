/* Copyright 2016 The Science and Technology Facilities Council (STFC)
 *
 * Authors: Jonathan Hogg (STFC)
 *
 * IMPORTANT: This file is NOT licenced under the BSD licence. If you wish to
 * licence this code, please contact STFC via hsl@stfc.ac.uk
 * (We are currently deciding what licence to release this code under if it
 * proves to be useful beyond our own academic experiments)
 *
 * This header file implements compatability functions depending on the value
 * of autoconf macros.
 */
#pragma once

#include <cstdint>
#include <cstdio>

#include "config.h"

#ifndef HAVE_STD_ALIGN
// Older versions of g++ (and intel that relies on equivalent -lstdc++) don't
// define std::align, so we do it ourselves.
namespace std {
void* align(std::size_t alignment, std::size_t size, void*& ptr, std::size_t& space);
} /* namespace std */
#endif
