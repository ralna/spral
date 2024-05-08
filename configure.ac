AC_INIT(spral, 2024.05.08, hsl@stfc.ac.uk)
AC_CONFIG_SRCDIR([src/ssids/ssids.f90])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Allow disabling of OpenMP
AC_ARG_ENABLE([openmp],
   AS_HELP_STRING([--disable-openmp], [Disable OpenMP parallelism. Not recommended.])
   )

# Allow disabling of GPU if desired
AC_ARG_ENABLE([gpu],
   AS_HELP_STRING([--disable-gpu], [Do not attempt to compile GPU code (otherwise presence of nvcc is autodetected).])
   )

# Allow user to enable profiling
AC_ARG_ENABLE([profile],
   AS_HELP_STRING([--enable-profile], [Enable support for profile generation. Require gtg library.])
   )

# Allow debugging of CUDA code
AC_ARG_ENABLE([gpudbg],
   AS_HELP_STRING([--enable-gpudbg], [Enable debugging of CUDA code.])
   )

# Allow debugging of SSIDS analyse phase
AC_ARG_ENABLE([analdbg],
   AS_HELP_STRING([--enable-analdbg], [Enable debugging of SSIDS analyse phase.])
   )

# Check for compilers
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_FC
AC_PROG_RANLIB

# Add OpenMP support
AS_IF([test "x$enable_openmp" != "xno"], [
   AC_LANG_PUSH(C)
   AC_OPENMP
   AC_LANG_POP(C)
   AC_LANG_PUSH(C++)
   AC_OPENMP
   AC_LANG_POP(C++)
   AC_LANG_PUSH(Fortran)
   AC_OPENMP
   AC_LANG_POP(Fortran)
   ])
AM_CONDITIONAL([HAVE_OPENMP], [test -n "$OPENMP_CXXFLAGS"])

# Check for NVCC
AS_IF([test "x$enable_gpu" != "xno"], [
   SPRAL_PROG_NVCC
   ])
AM_CONDITIONAL([HAVE_NVCC], [test -n "$NVCC"])

AS_IF([test "x$enable_gpudbg" = "xyes"], [
   GPUDBG=yes
   ])
AM_CONDITIONAL([HAVE_GPUDBG], [test -n "$GPUDBG"])

AS_IF([test "x$enable_analdbg" = "xyes"], [
   AC_DEFINE(ANALDBG,1,[Define to 1 to enable analyse debugging])
   ])

# Establish linking flags C->Fortran and how to link C with fortran
AC_FC_LIBRARY_LDFLAGS
SPRAL_NO_FORT_MAIN

# Establish -lstdc++ or equivalent
if test "x$CXXLIB" == "x"; then
   if test "x$CXX" == "xifort"; then
      CXXLIB="-cxxlib"
   else
      CXXLIB="-lstdc++"
   fi
fi
AC_SUBST(CXXLIB)

# Check for features we want
AC_MSG_CHECKING(for sched_getcpu())
AC_TRY_LINK([#include <sched.h>], [sched_getcpu();],
   AC_MSG_RESULT(yes); AC_DEFINE(HAVE_SCHED_GETCPU, 1, [Define to 1 if you have sched_getcpu().]),
   AC_MSG_RESULT(no)
   )


# Check for required libraries
AX_BLAS(,[AC_MSG_ERROR([No BLAS library found.])])
AX_LAPACK(,[AC_MSG_ERROR([No LAPACK library found.])])
SPRAL_METIS(,
   [AC_MSG_ERROR([No MeTiS library found.])]
   )

# Check for hwloc
PKG_PROG_PKG_CONFIG # initialise $PKG_CONFIG
PKG_CONFIG="$PKG_CONFIG --static" # we will be linking statically
PKG_CHECK_MODULES([HWLOC], [hwloc],
   AC_DEFINE(HAVE_HWLOC,1,[Define if you have hwloc library]),
   AC_MSG_WARN([hwloc not supplied: cannot detect NUMA regions])
   )

AS_IF([test "x$NVCC" != x], [
   SPRAL_NVCC_LIB
])

# Check for profiling library if desired
AS_IF([test "x$enable_profile" == "xyes"], [
   SPRAL_GTG(,[AC_MSG_ERROR([GTG library not found, cannot enable profiling])])
   AC_DEFINE(PROFILE,1,[Define to 1 to enable profiling])
   echo "Bite me $enable_profile"
   ])

# Output data
AC_CONFIG_FILES(Makefile)
AC_OUTPUT
