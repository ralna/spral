/** \file
 *  \copyright 2016 The Science and Technology Facilities Council (STFC)
 *  \licence   BSD licence, see LICENCE file for details
 *  \author    Jonathan Hogg
 */
#pragma once

namespace spral { namespace ssids { namespace cpu {

class SymbolicNode;

template <typename T>
class NumericNode {
public:
   /* Symbolic node associate with this one */
   SymbolicNode const* symb;

   /* Fixed data from analyse */
   NumericNode<T>* first_child; // Pointer to our first child
   NumericNode<T>* next_child; // Pointer to parent's next child

   /* Data that changes during factorize */
   int ndelay_in; // Number of delays arising from children
   int ndelay_out; // Number of delays arising to push into parent
   int nelim; // Number of columns succesfully eliminated
   T *lcol; // Pointer to start of factor data
   int *perm; // Pointer to permutation
   T *contrib; // Pointer to contribution block
};

}}} /* namespaces spral::ssids::cpu */
