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

#include <cmath>

#include "common.hxx"

namespace spral { namespace ssids { namespace cpu {

/** Calculates LD from L and D */
template <enum operation op, typename T>
void calcLD(int m, int n, const T *l, int ldl, const T *d, T *ld, int ldld) {
   for(int col=0; col<n; ) {
      if(col+1==n || std::isfinite(d[2*col+2])) {
         // 1x1 pivot
         T d11 = d[2*col];
         if(d11 != 0.0) d11 = 1/d11; // Zero pivots just cause zeroes
         for(int row=0; row<m; row++)
            ld[col*ldld+row] = d11 * ((op==OP_N) ? l[col*ldl+row]
                                                 : l[row*ldl+col]);
         col++;
      } else {
         // 2x2 pivot
         T d11 = d[2*col];
         T d21 = d[2*col+1];
         T d22 = d[2*col+3];
         T det = d11*d22 - d21*d21;
         d11 = d11/det;
         d21 = d21/det;
         d22 = d22/det;
         for(int row=0; row<m; row++) {
            T a1 = (op==OP_N) ? l[col*ldl+row]     : l[row*ldl+col];
            T a2 = (op==OP_N) ? l[(col+1)*ldl+row] : l[row*ldl+(col+1)];
            ld[col*ldld+row]     =  d22*a1 - d21*a2;
            ld[(col+1)*ldld+row] = -d21*a1 + d11*a2;
         }
         col += 2;
      }
   }
}

}}} /* namespaces spral::ssids::cpu */
