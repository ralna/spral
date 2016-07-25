AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

SUFFIXES = .cu

# NVCC setup
PTX_FLAGS = -v
#NVCCFLAGS = -Iinclude -arch=sm_20 -g -Xptxas="${PTX_FLAGS}"
OPENMP_LIB = -lgomp # FIXME: autoconf this
AM_NVCC_FLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src
AM_LD_FLAGS = -lcuda
NVCCLINK = \
	$(NVCC) $(NVCCFLAGS) $(AM_NVCC_FLAGS) $(AM_LD_FLAGS) $(LDFLAGS) $(OPENMP_LIB) -o $@
.cu.o:
	$(NVCC) $(NVCCFLAGS) $(AM_NVCC_FLAGS) -dc -o $@ $<

# Include directory for standard C
AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src $(OPENMP_CFLAGS)
AM_CXXFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src -std=c++11 $(OPENMP_CXXFLAGS) $(GTG_INCLUDE)
AM_FCFLAGS = $(OPENMP_CXXFLAGS) # assume CXX and FC use same flags...

# FIXME: make below configurable
LDADD = -L. -lspral $(METIS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(GTG_LIBS) $(FCLIBS) -lrt
if HAVE_NVCC
SPRAL_LINK_LIBS = -lcublas
SPRALLINK = $(NVCCLINK)
else
SPRAL_LINK_LIBS = $(CXXLIB)
SPRALLINK = $(FCLINK)
endif

lib_LIBRARIES = libspral.a
include_HEADERS = include/spral.h

# BLAS_IFACE
# Note: first so not +=
libspral_a_SOURCES = \
	src/blas_iface.f90

# CORE_ANALYSE
libspral_a_SOURCES += \
	src/core_analyse.f90

# CUDA
if HAVE_NVCC
libspral_a_SOURCES += \
	src/cuda/api_wrappers.cu \
	src/cuda/cuda.f90 \
	src/cuda/cuda_check.h
endif

# LAPACK_IFACE
libspral_a_SOURCES += \
	src/lapack_iface.f90

# LSMR
include_HEADERS += include/spral_lsmr.h
libspral_a_SOURCES += \
	src/lsmr.f90 \
	interfaces/C/lsmr.f90
check_PROGRAMS = \
	lsmr_test \
	examples/Fortran/lsmr \
	examples/C/lsmr
lsmr_test_SOURCES = tests/lsmr.f90
examples_Fortran_lsmr_SOURCES = examples/Fortran/lsmr.f90
tests/lsmr.$(OBJEXT): libspral.a
examples/C/lsmr.$(OBJEXT): libspral.a
examples/Fortran/lsmr.$(OBJEXT): libspral.a
TESTS = lsmr_test

# MATCH_ORDER
libspral_a_SOURCES += \
	src/match_order.f90

# MATRIX_UTIL
include_HEADERS += include/spral_matrix_util.h
libspral_a_SOURCES += \
	src/matrix_util.f90 \
	interfaces/C/matrix_util.f90

# METIS_WRAPPER
libspral_a_SOURCES += \
	src/metis$(METIS_VERSION)_wrapper.f90
EXTRA_DIST = src/metis4_wrapper.f90 src/metis5_wrapper.f90

# PGM
libspral_a_SOURCES += \
	src/pgm.f90

# RUTHERFORD_BOEING
include_HEADERS += include/spral_rutherford_boeing.h
libspral_a_SOURCES += \
	src/rutherford_boeing.f90 \
	interfaces/C/rutherford_boeing.f90

# SCALING
include_HEADERS += include/spral_scaling.h
libspral_a_SOURCES += \
  	src/scaling.f90 \
	interfaces/C/scaling.f90
check_PROGRAMS += \
	scaling_test \
	examples/Fortran/scaling/auction_sym \
	examples/Fortran/scaling/auction_unsym \
	examples/Fortran/scaling/equilib_sym \
	examples/Fortran/scaling/equilib_unsym \
	examples/Fortran/scaling/hungarian_sym \
	examples/Fortran/scaling/hungarian_unsym \
	examples/C/scaling/auction_sym \
	examples/C/scaling/auction_unsym \
	examples/C/scaling/equilib_sym \
	examples/C/scaling/equilib_unsym \
	examples/C/scaling/hungarian_sym \
	examples/C/scaling/hungarian_unsym
scaling_test_SOURCES = tests/scaling.f90
examples_Fortran_scaling_auction_sym_SOURCES = examples/Fortran/scaling/auction_sym.f90
examples_Fortran_scaling_auction_unsym_SOURCES = examples/Fortran/scaling/auction_unsym.f90
examples_Fortran_scaling_equilib_sym_SOURCES = examples/Fortran/scaling/equilib_sym.f90
examples_Fortran_scaling_equilib_unsym_SOURCES = examples/Fortran/scaling/equilib_unsym.f90
examples_Fortran_scaling_hungarian_sym_SOURCES = examples/Fortran/scaling/hungarian_sym.f90
examples_Fortran_scaling_hungarian_unsym_SOURCES = examples/Fortran/scaling/hungarian_unsym.f90
examples_C_scaling_auction_sym_SOURCES = examples/C/scaling/auction_sym.c
examples_C_scaling_auction_unsym_SOURCES = examples/C/scaling/auction_unsym.c
examples_C_scaling_equilib_sym_SOURCES = examples/C/scaling/equilib_sym.c
examples_C_scaling_equilib_unsym_SOURCES = examples/C/scaling/equilib_unsym.c
examples_C_scaling_hungarian_sym_SOURCES = examples/C/scaling/hungarian_sym.c
examples_C_scaling_hungarian_unsym_SOURCES = examples/C/scaling/hungarian_unsym.c
tests/scaling.$(OBJEXT): libspral.a
examples/Fortran/scaling/auction_sym.$(OBJEXT): libspral.a
examples/Fortran/scaling/auction_unsym.$(OBJEXT): libspral.a
examples/Fortran/scaling/equilib_sym.$(OBJEXT): libspral.a
examples/Fortran/scaling/equilib_unsym.$(OBJEXT): libspral.a
examples/Fortran/scaling/hungarian_sym.$(OBJEXT): libspral.a
examples/Fortran/scaling/hungarian_unsym.$(OBJEXT): libspral.a
examples/C/scaling/auction_sym.$(OBJEXT): libspral.a
examples/C/scaling/auction_unsym.$(OBJEXT): libspral.a
examples/C/scaling/equilib_sym.$(OBJEXT): libspral.a
examples/C/scaling/equilib_unsym.$(OBJEXT): libspral.a
examples/C/scaling/hungarian_sym.$(OBJEXT): libspral.a
examples/C/scaling/hungarian_unsym.$(OBJEXT): libspral.a
TESTS += scaling_test

# TIMER
libspral_a_SOURCES += \
  	src/timer.f90

# RANDOM
include_HEADERS += include/spral_random.h
libspral_a_SOURCES += \
	src/random.f90 \
	interfaces/C/random.f90
check_PROGRAMS += \
	random_test \
	examples/Fortran/random \
	examples/C/random
random_test_SOURCES = tests/random.f90
examples_Fortran_random_SOURCES = examples/Fortran/random.f90
examples_C_random_SOURCES = examples/C/random.c
tests/random.$(OBJEXT): libspral.a
examples/Fortran/random.$(OBJEXT): libspral.a
examples/C/random.$(OBJEXT): libspral.a
TESTS += random_test

# RANDOM_MATRIX
include_HEADERS += include/spral_random_matrix.h
libspral_a_SOURCES += \
	src/random_matrix.f90 \
	interfaces/C/random_matrix.f90
check_PROGRAMS += \
	random_matrix_test \
	examples/Fortran/random_matrix \
	examples/C/random_matrix
random_matrix_test_SOURCES = tests/random_matrix.f90
examples_Fortran_random_matrix_SOURCES = examples/Fortran/random_matrix.f90
examples_C_random_matrix_SOURCES = examples/C/random_matrix.c
tests/random_matrix.$(OBJEXT): libspral.a
examples/Fortran/random_matrix.$(OBJEXT): libspral.a
examples/C/random_matrix.$(OBJEXT): libspral.a
TESTS += random_matrix_test

# SSIDS
include_HEADERS += include/spral_ssids.h
if HAVE_NVCC
libspral_a_SOURCES += \
	src/ssids/type_select_gpu.f90 \
	src/ssids/gpu/alloc.f90 \
	src/ssids/gpu/akeep.f90 \
	src/ssids/gpu/cpu_solve.f90 \
	src/ssids/gpu/datatypes.f90 \
	src/ssids/gpu/dense_factor.f90 \
	src/ssids/gpu/factor.f90 \
	src/ssids/gpu/fkeep.f90 \
	src/ssids/gpu/inform.f90 \
	src/ssids/gpu/interfaces.f90 \
	src/ssids/gpu/solve.f90 \
	src/ssids/gpu/kernels/assemble.cu \
	src/ssids/gpu/kernels/datatypes.h \
	src/ssids/gpu/kernels/dense_factor.cu \
	src/ssids/gpu/kernels/dtrsv.h \
	src/ssids/gpu/kernels/node_solve.cu \
	src/ssids/gpu/kernels/node_solve.h \
	src/ssids/gpu/kernels/reorder.cu \
	src/ssids/gpu/kernels/solve.cu \
	src/ssids/gpu/kernels/syrk.cu \
else
libspral_a_SOURCES += \
	src/ssids/gpu/inform_nogpu.f90 \
	src/ssids/type_select_nogpu.f90
endif
libspral_a_SOURCES += \
	src/ssids/alloc.f90 \
	src/ssids/akeep.f90 \
	src/ssids/analyse.f90 \
	src/ssids/contrib.f90 \
	src/ssids/datatypes.f90 \
	src/ssids/fkeep.f90 \
	src/ssids/inform.f90 \
	src/ssids/ssids.f90 \
	src/ssids/subtree.f90 \
	src/ssids/cpu/AlignedAllocator.hxx \
	src/ssids/cpu/BlockPool.hxx \
	src/ssids/cpu/cpu_iface.f90 \
	src/ssids/cpu/cpu_iface.hxx \
	src/ssids/cpu/factor.hxx \
	src/ssids/cpu/NumericNode.hxx \
	src/ssids/cpu/NumericSubtree.cxx \
	src/ssids/cpu/NumericSubtree.hxx \
	src/ssids/cpu/smalloc.cxx \
	src/ssids/cpu/smalloc.hxx \
	src/ssids/cpu/subtree.f90 \
	src/ssids/cpu/profile.cxx \
	src/ssids/cpu/profile.hxx \
	src/ssids/cpu/PoolAlloc.hxx \
	src/ssids/cpu/SmallLeafNumericSubtree.hxx \
	src/ssids/cpu/SmallLeafSymbolicSubtree.hxx \
	src/ssids/cpu/SymbolicNode.hxx \
	src/ssids/cpu/SymbolicSubtree.cxx \
	src/ssids/cpu/SymbolicSubtree.hxx \
	src/ssids/cpu/kernels/assemble.hxx \
	src/ssids/cpu/kernels/common.hxx \
	src/ssids/cpu/kernels/block_ldlt.hxx \
	src/ssids/cpu/kernels/cholesky.cxx \
	src/ssids/cpu/kernels/cholesky.hxx \
	src/ssids/cpu/kernels/ldlt_app.hxx \
	src/ssids/cpu/kernels/ldlt_nopiv.cxx \
	src/ssids/cpu/kernels/ldlt_nopiv.hxx \
	src/ssids/cpu/kernels/ldlt_tpp.cxx \
	src/ssids/cpu/kernels/ldlt_tpp.hxx \
	src/ssids/cpu/kernels/SimdVec.hxx \
	src/ssids/cpu/kernels/wrappers.cxx \
	src/ssids/cpu/kernels/wrappers.hxx \
	interfaces/C/ssids.f90
bin_PROGRAMS = spral_ssids
spral_ssids_SOURCES = \
	driver/spral_ssids.f90
if HAVE_NVCC
spral_ssids_SOURCES += \
	driver/cuda_helper_gpu.f90
else
spral_ssids_SOURCES += \
	driver/cuda_helper_nogpu.f90
endif
check_PROGRAMS += \
	ssids_test \
	ssids_kernel_test \
	examples/Fortran/ssids \
	examples/C/ssids
ssids_test_SOURCES = tests/ssids/ssids.f90
ssids_kernel_test_SOURCES = tests/ssids/kernels.cxx \
									 tests/ssids/kernels/block_ldlt.cxx \
									 tests/ssids/kernels/block_ldlt.hxx \
									 tests/ssids/kernels/cholesky.cxx \
									 tests/ssids/kernels/cholesky.hxx \
									 tests/ssids/kernels/ldlt_app.cxx \
									 tests/ssids/kernels/ldlt_app.hxx \
									 tests/ssids/kernels/framework.cxx \
									 tests/ssids/kernels/framework.hxx \
									 tests/ssids/kernels/ldlt_nopiv.cxx \
									 tests/ssids/kernels/ldlt_nopiv.hxx \
									 tests/ssids/kernels/ldlt_tpp.cxx \
									 tests/ssids/kernels/ldlt_tpp.hxx
examples_Fortran_ssids_SOURCES = examples/Fortran/ssids.f90
examples/Fortran/ssids.$(OBJEXT): libspral.a
examples_C_ssids_SOURCES = examples/C/ssids.c
tests/ssids/ssids.$(OBJEXT): libspral.a
tests/ssids/kernels.$(OBJEXT): libspral.a
examples/C/ssids.$(OBJEXT): libspral.a
TESTS += ssids_test ssids_kernel_test
spral_ssids_LDADD = $(LDADD) $(SPRAL_LINK_LIBS)
ssids_test_LDADD = $(LDADD) $(SPRAL_LINK_LIBS)
examples_Fortran_ssids_LDADD = $(LDADD) $(SPRAL_LINK_LIBS)
examples_C_ssids_LDADD = $(LDADD) $(SPRAL_LINK_LIBS)
spral_ssids_LINK = $(SPRALLINK)
ssids_test_LINK = $(SPRALLINK)
examples_Fortran_ssids_LINK = $(SPRALLINK)
examples_C_ssids_LINK = $(SPRALLINK)

# SSMFE
include_HEADERS += include/spral_ssmfe.h
libspral_a_SOURCES += \
	src/ssmfe/ssmfe.f90 \
	src/ssmfe/core.f90 \
	src/ssmfe/expert.f90 \
	interfaces/C/ssmfe.f90 \
	interfaces/C/ssmfe_core.f90 \
	interfaces/C/ssmfe_expert.f90
check_PROGRAMS += \
   ssmfe_test \
	ssmfe_ciface_test \
	examples/C/ssmfe/hermitian \
	examples/C/ssmfe/precond_core \
	examples/C/ssmfe/precond_expert \
	examples/C/ssmfe/precond_ssmfe \
	examples/C/ssmfe/shift_invert \
	examples/Fortran/ssmfe/hermitian \
	examples/Fortran/ssmfe/precond_core \
	examples/Fortran/ssmfe/precond_expert \
	examples/Fortran/ssmfe/precond_ssmfe \
	examples/Fortran/ssmfe/shift_invert
ssmfe_test_SOURCES = tests/ssmfe/ssmfe.f90
ssmfe_ciface_test_SOURCES = tests/ssmfe/ssmfe_ciface.c \
									 tests/ssmfe/laplace2d.h \
									 tests/ssmfe/ldltf.h
tests/ssmfe/ssmfe.$(OBJEXT): libspral.a
tests/ssmfe/ssmfe_ciface.$(OBJEXT): libspral.a
TESTS += ssmfe_test \
			ssmfe_ciface_test
examples_C_ssmfe_hermitian_SOURCES = \
	examples/C/ssmfe/hermitian.c
examples_C_ssmfe_precond_core_SOURCES = \
	examples/C/ssmfe/precond_core.c \
	examples/C/ssmfe/laplace2d.h
examples_C_ssmfe_precond_expert_SOURCES = \
	examples/C/ssmfe/precond_expert.c \
	examples/C/ssmfe/laplace2d.h
examples_C_ssmfe_precond_ssmfe_SOURCES = \
	examples/C/ssmfe/precond_ssmfe.c \
	examples/C/ssmfe/laplace2d.h
examples_C_ssmfe_shift_invert_SOURCES = \
	examples/C/ssmfe/shift_invert.c \
	examples/C/ssmfe/laplace2d.h \
	examples/C/ssmfe/ldltf.h
examples_Fortran_ssmfe_hermitian_SOURCES = \
	examples/Fortran/ssmfe/hermitian.f90
examples_Fortran_ssmfe_precond_ssmfe_SOURCES = \
	examples/Fortran/ssmfe/precond_ssmfe.f90 \
	examples/Fortran/ssmfe/laplace2d.f90
examples_Fortran_ssmfe_precond_expert_SOURCES = \
	examples/Fortran/ssmfe/precond_expert.f90 \
	examples/Fortran/ssmfe/laplace2d.f90
examples_Fortran_ssmfe_precond_core_SOURCES = \
	examples/Fortran/ssmfe/precond_core.f90 \
	examples/Fortran/ssmfe/laplace2d.f90
examples_Fortran_ssmfe_shift_invert_SOURCES = \
	examples/Fortran/ssmfe/shift_invert.f90 \
	examples/Fortran/ssmfe/laplace2d.f90 \
	examples/Fortran/ssmfe/ldltf.f90
examples/C/ssmfe/hermitian.$(OBJEXT): libspral.a
examples/C/ssmfe/precond_core.$(OBJEXT): libspral.a
examples/C/ssmfe/precond_expert.$(OBJEXT): libspral.a
examples/C/ssmfe/precond_ssmfe.$(OBJEXT): libspral.a
examples/C/ssmfe/shift_invert.$(OBJEXT): libspral.a
examples/Fortran/ssmfe/hermitian.$(OBJEXT): libspral.a
examples/Fortran/ssmfe/precond_core.$(OBJEXT): libspral.a \
   examples/Fortran/ssmfe/laplace2d.$(OBJEXT)
examples/Fortran/ssmfe/precond_expert.$(OBJEXT): libspral.a \
   examples/Fortran/ssmfe/laplace2d.$(OBJEXT)
examples/Fortran/ssmfe/precond_ssmfe.$(OBJEXT): libspral.a \
   examples/Fortran/ssmfe/laplace2d.$(OBJEXT)
examples/Fortran/ssmfe/shift_invert.$(OBJEXT): libspral.a \
   examples/Fortran/ssmfe/laplace2d.$(OBJEXT) \
   examples/Fortran/ssmfe/ldltf.$(OBJEXT)

# Automake doesn't know much about modules, so get rid of them ourself
clean-local:
	-rm *.mod

# Fortran 90 dependencies
driver/ssids_spral.$(OBJEXT): libspral.a
src/match_order.$(OBJEXT): src/metis$(METIS_VERSION)_wrapper.$(OBJEXT) \
                           src/scaling.$(OBJEXT)
src/random_matrix.$(OBJEXT): src/random.$(OBJEXT)
src/rutherford_boeing.$(OBJEXT): src/matrix_util.$(OBJEXT) \
                                 src/random.$(OBJEXT)
src/scaling.$(OBJEXT): src/matrix_util.$(OBJEXT)
src/ssids/alloc.$(OBJEXT): src/ssids/datatypes.$(OBJEXT)
src/ssids/akeep.$(OBJEXT): src/ssids/datatypes.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT) \
									src/ssids/subtree.$(OBJEXT)
src/ssids/analyse.$(OBJEXT): src/core_analyse.$(OBJEXT) \
                             src/pgm.$(OBJEXT) \
                             src/ssids/akeep.$(OBJEXT) \
                             src/ssids/datatypes.$(OBJEXT) \
                             src/ssids/inform.$(OBJEXT) \
									  src/ssids/cpu/subtree.$(OBJEXT)
src/ssids/contrib.$(OBJEXT): src/ssids/datatypes.$(OBJEXT)
src/ssids/datatypes.$(OBJEXT): src/scaling.$(OBJEXT)
src/ssids/fkeep.$(OBJEXT): src/ssids/akeep.$(OBJEXT) \
                           src/ssids/alloc.$(OBJEXT) \
                           src/ssids/datatypes.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT) \
                           src/ssids/cpu/subtree.$(OBJEXT)
src/ssids/inform.$(OBJEXT): src/scaling.$(OBJEXT) \
                            src/ssids/datatypes.$(OBJEXT)
src/ssids/subtree.$(OBJEXT): src/ssids/contrib.$(OBJEXT) \
                             src/ssids/datatypes.$(OBJEXT) \
									  src/ssids/inform.$(OBJEXT)
src/ssids/cpu/cpu_iface.$(OBJEXT): src/ssids/datatypes.$(OBJEXT) \
                                   src/ssids/inform.$(OBJEXT)
src/ssids/cpu/subtree.$(OBJEXT): src/ssids/contrib.$(OBJEXT) \
                                 src/ssids/datatypes.$(OBJEXT) \
                                 src/ssids/inform.$(OBJEXT) \
                                 src/ssids/subtree.$(OBJEXT) \
                                 src/ssids/cpu/cpu_iface.$(OBJEXT)
src/ssids/gpu/akeep.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                               src/ssids/akeep.$(OBJEXT) \
                               src/ssids/datatypes.$(OBJEXT) \
                               src/ssids/inform.$(OBJEXT) \
                               src/ssids/gpu/inform.$(OBJEXT)
src/ssids/gpu/cpu_solve.$(OBJEXT): src/ssids/akeep.$(OBJEXT) \
                                   src/ssids/datatypes.$(OBJEXT) \
                                   src/ssids/fkeep.$(OBJEXT) \
                                   src/ssids/inform.$(OBJEXT)
src/ssids/gpu/datatypes.$(OBJEXT): src/cuda/cuda.$(OBJEXT)
src/ssids/gpu/interfaces.$(OBJEXT): src/ssids/gpu/datatypes.$(OBJEXT)
src/ssids/gpu/dense_factor.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                                      src/ssids/gpu/alloc.$(OBJEXT) \
                                      src/ssids/gpu/datatypes.$(OBJEXT) \
                                      src/ssids/gpu/interfaces.$(OBJEXT) \
                                      src/ssids/datatypes.$(OBJEXT)
src/ssids/gpu/factor.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                                src/ssids/alloc.$(OBJEXT) \
                                src/ssids/gpu/alloc.$(OBJEXT) \
                                src/ssids/gpu/datatypes.$(OBJEXT) \
                                src/ssids/gpu/interfaces.$(OBJEXT) \
                                src/ssids/datatypes.$(OBJEXT) \
                                src/ssids/gpu/dense_factor.$(OBJEXT) \
                                src/ssids/gpu/solve.$(OBJEXT)
src/ssids/gpu/fkeep.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                               src/ssids/akeep.$(OBJEXT) \
                               src/ssids/alloc.$(OBJEXT) \
                               src/ssids/datatypes.$(OBJEXT) \
                               src/ssids/fkeep.$(OBJEXT) \
                               src/ssids/inform.$(OBJEXT) \
                               src/ssids/gpu/akeep.$(OBJEXT) \
                               src/ssids/gpu/cpu_solve.$(OBJEXT) \
                               src/ssids/gpu/datatypes.$(OBJEXT) \
                               src/ssids/gpu/interfaces.$(OBJEXT) \
                               src/ssids/gpu/factor.$(OBJEXT) \
                               src/ssids/gpu/inform.$(OBJEXT) \
                               src/ssids/gpu/solve.$(OBJEXT)
src/ssids/gpu/inform.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                                src/ssids/datatypes.$(OBJEXT) \
                                src/ssids/inform.$(OBJEXT)
src/ssids/gpu/inform_nogpu.$(OBJEXT): src/ssids/datatypes.$(OBJEXT) \
                                      src/ssids/inform.$(OBJEXT)
src/ssids/gpu/solve.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                               src/ssids/gpu/alloc.$(OBJEXT) \
                               src/ssids/gpu/datatypes.$(OBJEXT) \
                               src/ssids/gpu/interfaces.$(OBJEXT) \
                               src/ssids/datatypes.$(OBJEXT)
if HAVE_NVCC
src/ssids/ssids.$(OBJEXT): src/match_order.$(OBJEXT) \
                           src/matrix_util.$(OBJEXT) \
                           src/metis$(METIS_VERSION)_wrapper.$(OBJEXT) \
                           src/scaling.$(OBJEXT) \
                           src/ssids/akeep.$(OBJEXT) \
                           src/ssids/alloc.$(OBJEXT) \
                           src/ssids/analyse.$(OBJEXT) \
                           src/ssids/datatypes.$(OBJEXT) \
                           src/ssids/fkeep.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT) \
                           src/ssids/type_select_gpu.$(OBJEXT)
else
src/ssids/ssids.$(OBJEXT): src/match_order.$(OBJEXT) \
                           src/matrix_util.$(OBJEXT) \
                           src/metis$(METIS_VERSION)_wrapper.$(OBJEXT) \
                           src/scaling.$(OBJEXT) \
                           src/ssids/akeep.$(OBJEXT) \
                           src/ssids/alloc.$(OBJEXT) \
                           src/ssids/analyse.$(OBJEXT) \
                           src/ssids/datatypes.$(OBJEXT) \
                           src/ssids/fkeep.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT) \
                           src/ssids/type_select_nogpu.$(OBJEXT)
endif
src/ssids/type_select_gpu.$(OBJEXT): src/ssids/gpu/akeep.$(OBJEXT) \
                                     src/ssids/gpu/fkeep.$(OBJEXT) \
                                     src/ssids/gpu/inform.$(OBJEXT)
src/ssids/type_select_nogpu.$(OBJEXT): src/ssids/akeep.$(OBJEXT) \
                                       src/ssids/fkeep.$(OBJEXT) \
                                       src/ssids/inform.$(OBJEXT)
src/ssmfe/core.$(OBJEXT): src/blas_iface.$(OBJEXT) \
                          src/lapack_iface.$(OBJEXT)
src/ssmfe/expert.$(OBJEXT): src/ssmfe/core.$(OBJEXT)
src/ssmfe/ssmfe.$(OBJEXT): src/ssmfe/expert.$(OBJEXT)
if HAVE_NVCC
driver/spral_ssids.$(OBJEXT): libspral.a \
   driver/cuda_helper_gpu.$(OBJEXT)
else
driver/spral_ssids.$(OBJEXT): libspral.a \
   driver/cuda_helper_nogpu.$(OBJEXT)
endif

# CUDA header deps
src/ssids/gpu/kernels/solve.$(OBJEXT): src/ssids/gpu/kernels/dtrsv.h

# C interface deps
interfaces/C/random.$(OBJEXT): src/random.$(OBJEXT)
interfaces/C/random_matrix.$(OBJEXT): src/random_matrix.$(OBJEXT)
interfaces/C/rutherford_boeing.$(OBJEXT): src/rutherford_boeing.$(OBJEXT)
interfaces/C/scaling.$(OBJEXT): src/scaling.$(OBJEXT)
if HAVE_NVCC
interfaces/C/ssids.$(OBJEXT): src/ssids/ssids.$(OBJEXT) \
                              src/ssids/inform.$(OBJEXT) \
                              src/ssids/gpu/inform.$(OBJEXT)
else
interfaces/C/ssids.$(OBJEXT): src/ssids/ssids.$(OBJEXT) \
                              src/ssids/inform.$(OBJEXT) \
                              src/ssids/gpu/inform_nogpu.$(OBJEXT)
endif
interfaces/C/ssmfe.$(OBJEXT): src/ssmfe/ssmfe.$(OBJEXT) \
                              interfaces/C/ssmfe_expert.$(OBJEXT)
interfaces/C/ssmfe_core.$(OBJEXT): src/ssmfe/core.$(OBJEXT)
interfaces/C/ssmfe_expert.$(OBJEXT): src/ssmfe/core.$(OBJEXT) \
                                     interfaces/C/ssmfe_core.$(OBJEXT)
