#
# SYNOPSIS
#
#     SPRAL_NVCC_LIB()
#
# DESCRIPTION
#
#  This macro looks for the CUDA include directories. On success, the following
#  variables are defined:
#     $NVCC_INCLUDE_FLAGS
#  The following preprocessor macros are defined HAVE_NVCC
#
AC_DEFUN([SPRAL_NVCC_LIB], [
AC_REQUIRE([SPRAL_PROG_NVCC])

# Check in default path first
NVCC_INCLUDE_FLAGS="$NVCC_INCLUDE_FLAGS"
AC_CHECK_HEADER([cuda_runtime_api.h], [spral_nvcc_inc_ok=yes])

# Check in CUDA_HOME/include
if test x"$spral_nvcc_inc_ok" != xyes; then
   NVCC_INCLUDE_FLAGS="$NVCC_INCLUDE_FLAGS -I$CUDA_HOME/include"
   save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $NVCC_INCLUDE_FLAGS"
   AC_CHECK_HEADER([cuda_runtime_api.h], [spral_nvcc_inc_ok=yes])
   CPPFLAGS=$save_CPPFLAGS
fi

# Handle success or failure
AC_SUBST(NVCC_INCLUDE_FLAGS)
if test x"$spral_nvcc_inc_ok" = xyes; then
   AC_DEFINE(HAVE_NVCC,1,[Define to 1 if you are compiling against NVCC])
else
   AC_MSG_ERROR([NVCC include path not found])
fi

])dnl SPRAL_NVCC_LIB
