#
# SYNOPSIS
#
#     AX_PROG_NVCC([compiler-search-list])
#
# DESCRIPTION
#
#  This macro looks for the CUDA C compiler to use. On success, the following
#  variables are defined:
#     $NVCC          The CUDA compiler command
#     $NVCCFLAGS     The flags pased to the CUDA compiler
#

AC_DEFUN([ACX_PROG_NVCC], [

AC_ARG_VAR(NVCC,[CUDA compiler command])
AC_ARG_VAR(NVCCFLAGS,[CUDA compiler flags])

test "x$NVCC" = x && AC_CHECK_PROGS(NVCC,nvcc)
test "x$NVCCFLAGS" = x && NVCCFLAGS="-g -arch=compute_20 -code=compute_20,sm_20,sm_35"

])dnl AX_NVCC
