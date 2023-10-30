AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

SUFFIXES = .cu

# NVCC setup
PTX_FLAGS = -v
if HAVE_OPENMP
OPENMP_LIB = -Xcompiler="$(OPENMP_CXXFLAGS)"
else
OPENMP_LIB = #-lgomp # FIXME: autoconf this
endif
if HAVE_GPUDBG
AM_NVCC_FLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src -Xptxas="${PTX_FLAGS}" -G -src-in-ptx
else
AM_NVCC_FLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src -Xptxas="${PTX_FLAGS}"
endif
AM_NVCC_FLAGS += -gencode arch=compute_52,code=sm_52
AM_NVCC_FLAGS += -gencode arch=compute_60,code=sm_60
AM_NVCC_FLAGS += -gencode arch=compute_61,code=sm_61
AM_NVCC_FLAGS += -gencode arch=compute_62,code=sm_62
AM_NVCC_FLAGS += -gencode arch=compute_70,code=sm_70
AM_NVCC_FLAGS += -gencode arch=compute_72,code=sm_72
AM_NVCC_FLAGS += -gencode arch=compute_75,code=sm_75
AM_NVCC_FLAGS += -gencode arch=compute_80,code=sm_80
AM_NVCC_FLAGS += -gencode arch=compute_86,code=sm_86
AM_NVCC_FLAGS += -gencode arch=compute_87,code=sm_87
AM_NVCC_FLAGS += -gencode arch=compute_86,code=compute_86
AM_LD_FLAGS = -lcuda
NVCCLINK = \
	$(NVCC) $(NVCCFLAGS) $(AM_NVCC_FLAGS) $(AM_LD_FLAGS) $(LDFLAGS) $(OPENMP_LIB) -o $@
.cu.o:
	$(NVCC) $(NVCCFLAGS) $(AM_NVCC_FLAGS) -rdc false -c -o $@ $<

# Include directory for standard C
AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src $(OPENMP_CFLAGS) \
				$(HWLOC_CFLAGS) -std=c99
AM_CXXFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src -std=c++11 \
				  $(OPENMP_CXXFLAGS) $(GTG_INCLUDE) $(HWLOC_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)

# FIXME: make below configurable
if HAVE_NVCC
LDADD = -L. -lspral $(METIS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(GTG_LIBS) -L$(CUDA_HOME)/lib64 $(HWLOC_LIBS) $(FCLIBS)
SPRAL_LINK_LIBS = -lcublas
SPRALLINK = $(NVCCLINK)
else
LDADD = -L. -lspral $(METIS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(GTG_LIBS) $(HWLOC_LIBS) $(FCLIBS)
SPRAL_LINK_LIBS = $(CXXLIB)
SPRALLINK = $(FCLINK)
endif

lib_LIBRARIES = libspral.a
# Note: first so not +=
include_HEADERS = \
   include/spral.h \
   include/spral_complex.h

# BLAS_IFACE
# Note: first so not +=
libspral_a_SOURCES = \
	src/blas_iface.f90

# CORE_ANALYSE
libspral_a_SOURCES += \
	src/core_analyse.f90

# COMPAT (compatability header)
libspral_a_SOURCES += \
   src/compat.cxx

# CUDA
if HAVE_NVCC
libspral_a_SOURCES += \
	src/cuda/api_wrappers.cu \
	src/cuda/cuda.f90 \
	src/cuda/cuda_check.h
else
libspral_a_SOURCES += \
	src/cuda/cuda_nocuda.f90
endif

# HW_TOPOLOGY
libspral_a_SOURCES += \
   src/hw_topology/hw_topology.f90 \
	src/hw_topology/guess_topology.cxx \
	src/hw_topology/guess_topology.hxx \
	src/hw_topology/hwloc_wrapper.hxx

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
	src/metis$(METIS_VERSION)_wrapper.F90
EXTRA_DIST = src/metis4_wrapper.F90 src/metis5_wrapper.F90

# OMP
libspral_a_SOURCES += \
	src/omp.cxx \
	src/omp.hxx

# PGM
libspral_a_SOURCES += \
	src/pgm.f90

# RUTHERFORD_BOEING
include_HEADERS += include/spral_rutherford_boeing.h
libspral_a_SOURCES += \
	src/rutherford_boeing.f90 \
	interfaces/C/rutherford_boeing.f90
check_PROGRAMS += \
	rutherford_boeing_test \
	examples/Fortran/rutherford_boeing/rb_read \
	examples/Fortran/rutherford_boeing/rb_write \
	examples/C/rutherford_boeing/rb_read \
	examples/C/rutherford_boeing/rb_write
rutherford_boeing_test_SOURCES = tests/rutherford_boeing.f90
examples_Fortran_rutherford_boeing_rb_read_SOURCES = \
	examples/Fortran/rutherford_boeing/rb_read.f90
examples_Fortran_rutherford_boeing_rb_write_SOURCES = \
	examples/Fortran/rutherford_boeing/rb_write.f90
examples_C_rutherford_boeing_rb_read_SOURCES = \
	examples/C/rutherford_boeing/rb_read.c
examples_C_rutherford_boeing_rb_write_SOURCES = \
	examples/C/rutherford_boeing/rb_write.c
tests/rutherford_boeing.$(OBJEXT): libspral.a
examples/Fortran/rutherford_boeing/rb_read.$(OBJEXT): libspral.a
examples/Fortran/rutherford_boeing/rb_write.$(OBJEXT): libspral.a
examples/C/rutherford_boeing/rb_read.$(OBJEXT): libspral.a
examples/C/rutherford_boeing/rb_write.$(OBJEXT): libspral.a
TESTS += rutherford_boeing_test

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
	src/ssids/gpu/alloc.f90 \
	src/ssids/gpu/cpu_solve.f90 \
	src/ssids/gpu/datatypes.f90 \
	src/ssids/gpu/dense_factor.f90 \
	src/ssids/gpu/factor.f90 \
	src/ssids/gpu/interfaces.f90 \
	src/ssids/gpu/solve.f90 \
	src/ssids/gpu/smalloc.f90 \
	src/ssids/gpu/subtree.f90 \
	src/ssids/gpu/kernels/assemble.cu \
	src/ssids/gpu/kernels/datatypes.h \
	src/ssids/gpu/kernels/dense_factor.cu \
	src/ssids/gpu/kernels/dtrsv.h \
	src/ssids/gpu/kernels/reorder.cu \
	src/ssids/gpu/kernels/solve.cu \
	src/ssids/gpu/kernels/syrk.cu
else
libspral_a_SOURCES += \
	src/ssids/gpu/subtree_no_cuda.f90
endif
libspral_a_SOURCES += \
	src/ssids/akeep.f90 \
	src/ssids/anal.f90 \
	src/ssids/contrib.f90 \
	src/ssids/contrib.h \
	src/ssids/contrib_free.f90 \
	src/ssids/datatypes.f90 \
	src/ssids/doc.hxx \
	src/ssids/fkeep.f90 \
	src/ssids/inform.f90 \
	src/ssids/profile.cxx \
	src/ssids/profile.hxx \
	src/ssids/profile_iface.f90 \
	src/ssids/ssids.f90 \
	src/ssids/subtree.f90 \
	src/ssids/cpu/AppendAlloc.hxx \
	src/ssids/cpu/BlockPool.hxx \
	src/ssids/cpu/BuddyAllocator.hxx \
	src/ssids/cpu/cpu_iface.f90 \
	src/ssids/cpu/cpu_iface.hxx \
	src/ssids/cpu/factor.hxx \
	src/ssids/cpu/NumericNode.hxx \
	src/ssids/cpu/NumericSubtree.cxx \
	src/ssids/cpu/NumericSubtree.hxx \
	src/ssids/cpu/subtree.f90 \
	src/ssids/cpu/SmallLeafNumericSubtree.hxx \
	src/ssids/cpu/SmallLeafSymbolicSubtree.hxx \
	src/ssids/cpu/SymbolicNode.hxx \
	src/ssids/cpu/SymbolicSubtree.cxx \
	src/ssids/cpu/SymbolicSubtree.hxx \
	src/ssids/cpu/ThreadStats.cxx \
	src/ssids/cpu/ThreadStats.hxx \
	src/ssids/cpu/Workspace.hxx \
	src/ssids/cpu/kernels/assemble.hxx \
	src/ssids/cpu/kernels/common.hxx \
	src/ssids/cpu/kernels/block_ldlt.hxx \
	src/ssids/cpu/kernels/calc_ld.hxx \
	src/ssids/cpu/kernels/cholesky.cxx \
	src/ssids/cpu/kernels/cholesky.hxx \
	src/ssids/cpu/kernels/ldlt_app.cxx \
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
	driver/spral_ssids.F90
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
									 tests/ssids/kernels/AlignedAllocator.hxx \
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
examples_C_ssids_LINK = $(SPRALLINK) $(NO_FORT_MAIN)

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
src/ssids/akeep.$(OBJEXT): src/ssids/datatypes.$(OBJEXT) \
                           src/hw_topology/hw_topology.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT) \
									src/ssids/subtree.$(OBJEXT)
if HAVE_NVCC
src/ssids/anal.$(OBJEXT): src/core_analyse.$(OBJEXT) \
                          src/cuda/cuda.$(OBJEXT) \
                          src/hw_topology/hw_topology.$(OBJEXT) \
                          src/pgm.$(OBJEXT) \
                          src/ssids/akeep.$(OBJEXT) \
                          src/ssids/datatypes.$(OBJEXT) \
                          src/ssids/inform.$(OBJEXT) \
                          src/ssids/cpu/subtree.$(OBJEXT) \
                          src/ssids/gpu/subtree.$(OBJEXT)
else
src/ssids/anal.$(OBJEXT): src/core_analyse.$(OBJEXT) \
                          src/cuda/cuda_nocuda.$(OBJEXT) \
                          src/hw_topology/hw_topology.$(OBJEXT) \
                          src/pgm.$(OBJEXT) \
                          src/ssids/akeep.$(OBJEXT) \
                          src/ssids/datatypes.$(OBJEXT) \
                          src/ssids/inform.$(OBJEXT) \
                          src/ssids/cpu/subtree.$(OBJEXT) \
                          src/ssids/gpu/subtree_no_cuda.$(OBJEXT)
endif
src/ssids/contrib.$(OBJEXT): src/ssids/datatypes.$(OBJEXT)
src/ssids/datatypes.$(OBJEXT): src/scaling.$(OBJEXT)
src/ssids/fkeep.$(OBJEXT): src/ssids/akeep.$(OBJEXT) \
                           src/ssids/datatypes.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT) \
                           src/ssids/profile_iface.$(OBJEXT) \
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
src/ssids/gpu/cpu_solve.$(OBJEXT): src/ssids/datatypes.$(OBJEXT)
src/ssids/gpu/datatypes.$(OBJEXT): src/cuda/cuda.$(OBJEXT)
src/ssids/gpu/interfaces.$(OBJEXT): src/ssids/gpu/datatypes.$(OBJEXT)
src/ssids/gpu/dense_factor.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                                      src/ssids/gpu/alloc.$(OBJEXT) \
                                      src/ssids/gpu/datatypes.$(OBJEXT) \
                                      src/ssids/gpu/interfaces.$(OBJEXT) \
                                      src/ssids/datatypes.$(OBJEXT)
src/ssids/gpu/factor.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
										  src/ssids/contrib.$(OBJEXT) \
                                src/ssids/datatypes.$(OBJEXT) \
                                src/ssids/profile_iface.$(OBJEXT) \
                                src/ssids/gpu/alloc.$(OBJEXT) \
                                src/ssids/gpu/datatypes.$(OBJEXT) \
                                src/ssids/gpu/interfaces.$(OBJEXT) \
                                src/ssids/gpu/dense_factor.$(OBJEXT) \
                                src/ssids/gpu/solve.$(OBJEXT) \
                                src/ssids/gpu/smalloc.$(OBJEXT)
src/ssids/gpu/solve.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                               src/ssids/gpu/alloc.$(OBJEXT) \
                               src/ssids/gpu/datatypes.$(OBJEXT) \
                               src/ssids/gpu/interfaces.$(OBJEXT) \
                               src/ssids/datatypes.$(OBJEXT)
src/ssids/gpu/smalloc.$(OBJEXT): src/ssids/datatypes.$(OBJEXT)
src/ssids/gpu/subtree.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
                                 src/ssids/contrib.$(OBJEXT) \
                                 src/ssids/datatypes.$(OBJEXT) \
                                 src/ssids/inform.$(OBJEXT) \
                                 src/ssids/subtree.$(OBJEXT) \
                                 src/ssids/gpu/datatypes.$(OBJEXT) \
                                 src/ssids/gpu/factor.$(OBJEXT) \
                                 src/ssids/gpu/solve.$(OBJEXT) \
                                 src/ssids/gpu/smalloc.$(OBJEXT)
src/ssids/gpu/subtree_no_cuda.$(OBJEXT): src/ssids/contrib.$(OBJEXT) \
                                         src/ssids/datatypes.$(OBJEXT) \
                                         src/ssids/subtree.$(OBJEXT) \
                                         src/ssids/inform.$(OBJEXT)
if HAVE_NVCC
src/ssids/ssids.$(OBJEXT): src/match_order.$(OBJEXT) \
                           src/matrix_util.$(OBJEXT) \
                           src/metis$(METIS_VERSION)_wrapper.$(OBJEXT) \
                           src/scaling.$(OBJEXT) \
                           src/ssids/akeep.$(OBJEXT) \
                           src/ssids/anal.$(OBJEXT) \
                           src/ssids/datatypes.$(OBJEXT) \
                           src/ssids/fkeep.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT)
else
src/ssids/ssids.$(OBJEXT): src/hw_topology/hw_topology.$(OBJEXT) \
                           src/match_order.$(OBJEXT) \
                           src/matrix_util.$(OBJEXT) \
                           src/metis$(METIS_VERSION)_wrapper.$(OBJEXT) \
                           src/scaling.$(OBJEXT) \
                           src/ssids/akeep.$(OBJEXT) \
                           src/ssids/anal.$(OBJEXT) \
                           src/ssids/datatypes.$(OBJEXT) \
                           src/ssids/fkeep.$(OBJEXT) \
                           src/ssids/inform.$(OBJEXT)
endif
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
                              src/ssids/inform.$(OBJEXT)
else
interfaces/C/ssids.$(OBJEXT): src/ssids/ssids.$(OBJEXT) \
                              src/ssids/inform.$(OBJEXT)
endif
interfaces/C/ssmfe.$(OBJEXT): src/ssmfe/ssmfe.$(OBJEXT) \
                              interfaces/C/ssmfe_expert.$(OBJEXT)
interfaces/C/ssmfe_core.$(OBJEXT): src/ssmfe/core.$(OBJEXT)
interfaces/C/ssmfe_expert.$(OBJEXT): src/ssmfe/core.$(OBJEXT) \
                                     interfaces/C/ssmfe_core.$(OBJEXT)
