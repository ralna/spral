#ifndef SPRAL_COMPLEX_H
#define SPRAL_COMPLEX_H

// We have to be careful with C and C++ as they have different complex types
// the following uses typedefs to work around this.
#ifdef __cplusplus
#include <complex>
typedef std::complex<double> spral_double_complex;
#else /* ifdef __cplusplus */
#include <complex.h>
typedef double complex spral_double_complex;
#endif /* ifdef __cplusplus */

#endif /* SPRAL_COMPLEX_H */
