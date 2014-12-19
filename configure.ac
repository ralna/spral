AC_INIT(spral, 2014.12.19-dev, jonathan.hogg@stfc.ac.uk)
AC_CONFIG_SRCDIR([src/ssids/ssids.f90])
AM_INIT_AUTOMAKE

# Check for compilers
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_FC
AC_PROG_RANLIB
ACX_PROG_NVCC
AM_CONDITIONAL([HAVE_NVCC], [test -n "$NVCC"])

# Establish linking flags C->Fortran
AC_FC_LIBRARY_LDFLAGS

# Check for required libraries
ACX_BLAS(,[AC_MSG_ERROR([No BLAS library found.])])
ACX_LAPACK(,[AC_MSG_ERROR([No LAPACK library found.])])
SPRAL_METIS(,
   [AC_MSG_ERROR([No MeTiS library found.])]
   )

# Output data
AC_CONFIG_FILES(Makefile)
AC_OUTPUT
