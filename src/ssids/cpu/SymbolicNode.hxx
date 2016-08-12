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

#include <vector>

namespace spral { namespace ssids { namespace cpu {

/** Symbolic representation of a node */
struct SymbolicNode {
   bool insmallleaf;
   int idx; //< Index of node
   int nrow; //< Number of rows
   int ncol; //< Number of columns
   SymbolicNode* first_child; //< Pointer to first child in linked list
   SymbolicNode* next_child; //< Pointer to second child in linked list
   int const* rlist; //< Pointer to row lists
   int num_a; //< Number of entries mapped from A to L
   int const* amap; //< Pointer to map from A to L locations
   int parent; //< index of parent node
   std::vector<int> contrib; //< index of expected contribution(s)
};

}}} /* end of namespace spral::ssids::cpu */
