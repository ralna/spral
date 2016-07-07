AC_INIT(spral, 2016.06.24, jonathan.hogg@stfc.ac.uk)
AC_CONFIG_SRCDIR([src/ssids/ssids.f90])
AM_INIT_AUTOMAKE

# Allow disabling of GPU if desired
AC_ARG_ENABLE([gpu],
   AS_HELP_STRING([--disable-gpu], [Do not attempt to compile GPU code (otherwise presence of nvcc is autodetected)])
   )

# Check for compilers
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_FC
AC_PROG_RANLIB
AS_IF([test "x$enable_gpu" != "xno"], [
   ACX_PROG_NVCC
   ])
AM_CONDITIONAL([HAVE_NVCC], [test -n "$NVCC"])

# Add OpenMP support, set language standards
AC_LANG_PUSH(C)
AC_OPENMP
AC_PROG_CC_C99
AC_LANG(C++)
AC_OPENMP
AC_LANG_POP(C++)

# Establish linking flags C->Fortran
AC_FC_LIBRARY_LDFLAGS

# Establish -lstdc++ or equivalent
if test "x$CXXLIB" == "x"; then
   if test "x$CXX" == "xifort"; then
      CXXLIB="-cxxlib"
   else
      CXXLIB="-lstdc++"
   fi
fi
AC_SUBST(CXXLIB)

# Check for required libraries
ACX_BLAS(,[AC_MSG_ERROR([No BLAS library found.])])
ACX_LAPACK(,[AC_MSG_ERROR([No LAPACK library found.])])
SPRAL_METIS(,
   [AC_MSG_ERROR([No MeTiS library found.])]
   )

# Output data
AC_CONFIG_FILES(Makefile)
AC_OUTPUT
