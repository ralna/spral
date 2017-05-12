#
# SYNOPSIS
#
#     SPRAL_PROG_NVCC([compiler-search-list])
#
# DESCRIPTION
#
#  This macro looks for the CUDA C compiler to use. On success, the following
#  variables are defined:
#     $NVCC          The CUDA compiler command
#     $NVCCFLAGS     The flags pased to the CUDA compiler
#     $NVCC_ARCH_SM  The CUDA architectures to compile for
#

AC_DEFUN([SPRAL_PROG_NVCC], [

AC_ARG_VAR(NVCC,[CUDA compiler command])
AC_ARG_VAR(NVCCFLAGS,[CUDA compiler flags])

test "x$NVCC" = x && AC_CHECK_PROGS(NVCC,nvcc)
$NVCC -DNDEBUG nvcc_arch_sm.c -o nvcc_arch_sm -lcuda
test "x$NVCC_ARCH_SM" = x && NVCC_ARCH_SM=`./nvcc_arch_sm`
test "x$NVCCFLAGS" = x && NVCCFLAGS="-std=c++11 -g $NVCC_ARCH_SM"

])dnl SPRAL_PROG_NVCC
