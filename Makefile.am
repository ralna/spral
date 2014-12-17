AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

SUFFIXES = .cu

# NVCC setup
PTX_FLAGS = -v
#NVCCFLAGS = -Iinclude -arch=sm_20 -g -Xptxas="${PTX_FLAGS}"
AM_NVCC_FLAGS = -I$(top_srcdir)/include
AM_LD_FLAGS = -lcuda
NVCCLINK = \
	$(NVCC) $(NVCCFLAGS) $(AM_NVCC_FLAGS) $(AM_LD_FLAGS) $(LDFLAGS) -o $@
.cu.o:
	$(NVCC) $(NVCCFLAGS) $(AM_NVCC_FLAGS) -dc -o $@ $<

# Include directory for standard C
AM_CFLAGS = -I$(top_srcdir)/include -std=c99

# FIXME: make below configurable
LDADD = -L. -lspral $(METIS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) -lrt
CUDA_LIBS = -lcublas

lib_LIBRARIES = libspral.a
include_HEADERS = include/spral.h

# CORE_ANALYSE
# Note: first so not +=
libspral_a_SOURCES = \
	src/core_analyse.f90

# CUDA
if HAVE_NVCC
libspral_a_SOURCES += \
	src/cuda/api_wrappers.cu \
	src/cuda/cuda.f90 \
	src/cuda/cuda_check.h
endif

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
libspral_a_SOURCES += \
	src/rutherford_boeing.f90

# SCALING
libspral_a_SOURCES += \
  	src/scaling.f90 \
	interfaces/C/scaling.f90
check_PROGRAMS = \
	scaling_test \
	examples/Fortran/scaling/auction_sym \
	examples/Fortran/scaling/auction_unsym \
	examples/Fortran/scaling/equilib_sym \
	examples/Fortran/scaling/equilib_unsym \
	examples/Fortran/scaling/hungarian_sym \
	examples/Fortran/scaling/hungarian_unsym \
	examples/C/scaling/auction_sym \
	examples/C/scaling/equilib_sym \
	examples/C/scaling/hungarian_sym
#examples/C/scaling/auction_unsym \
#examples/C/scaling/equilib_sym \
#examples/C/scaling/hungarian_unsym
scaling_test_SOURCES = tests/scaling.f90
examples_Fortran_scaling_auction_sym_SOURCES = examples/Fortran/scaling/auction_sym.f90
examples_Fortran_scaling_auction_unsym_SOURCES = examples/Fortran/scaling/auction_unsym.f90
examples_Fortran_scaling_equilib_sym_SOURCES = examples/Fortran/scaling/equilib_sym.f90
examples_Fortran_scaling_equilib_unsym_SOURCES = examples/Fortran/scaling/equilib_unsym.f90
examples_Fortran_scaling_hungarian_sym_SOURCES = examples/Fortran/scaling/hungarian_sym.f90
examples_Fortran_scaling_hungarian_unsym_SOURCES = examples/Fortran/scaling/hungarian_unsym.f90
examples_C_scaling_auction_sym_SOURCES = examples/C/scaling/auction_sym.c
#examples_C_scaling_auction_unsym_SOURCES = examples/C/scaling/auction_unsym.c
examples_C_scaling_equilib_sym_SOURCES = examples/C/scaling/equilib_sym.c
#examples_C_scaling_equilib_unsym_SOURCES = examples/C/scaling/equilib_unsym.c
examples_C_scaling_hungarian_sym_SOURCES = examples/C/scaling/hungarian_sym.c
#examples_C_scaling_hungarian_unsym_SOURCES = examples/C/scaling/hungarian_unsym.c
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
TESTS = scaling_test

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
if HAVE_NVCC
include_HEADERS += include/spral_ssids.h
libspral_a_SOURCES += \
	src/ssids/alloc.f90 \
	src/ssids/analyse.f90 \
	src/ssids/assemble_kernels.cu \
	src/ssids/cuda_datatypes.f90 \
	src/ssids/cuda_interfaces.f90 \
	src/ssids/datatypes.f90 \
	src/ssids/dense_factor_gpu.f90 \
	src/ssids/dense_factor_kernels.cu \
	src/ssids/dtrsv.h \
	src/ssids/factor_gpu.f90 \
	src/ssids/kernel_datatypes.h \
	src/ssids/node_solve_kernels.cu \
	src/ssids/node_solve_kernels.h \
	src/ssids/reorder_kernels.cu \
	src/ssids/solve_cpu.f90 \
	src/ssids/solve_gpu.f90 \
	src/ssids/solve_kernels.cu \
	src/ssids/ssids.f90 \
	src/ssids/syrk_kernels.cu \
	interfaces/C/ssids.f90
bin_PROGRAMS = spral_ssids
spral_ssids_SOURCES = \
	driver/spral_ssids.f90
spral_ssids_LINK = $(NVCCLINK)
spral_ssids_LDADD = $(LDADD) $(CUDA_LIBS)
check_PROGRAMS += \
	ssids_test \
	examples/Fortran/ssids \
	examples/C/ssids
ssids_test_SOURCES = tests/ssids.f90
ssids_test_LINK = $(NVCCLINK)
ssids_test_LDADD = $(LDADD) $(CUDA_LIBS)
examples_Fortran_ssids_SOURCES = examples/Fortran/ssids.f90
examples_Fortran_ssids_LDADD = $(LDADD) $(CUDA_LIBS)
examples_Fortran_ssids_LINK = $(NVCCLINK)
examples/Fortran/ssids.$(OBJEXT): libspral.a
examples_C_ssids_SOURCES = examples/C/ssids.c
examples_C_ssids_LDADD = $(LDADD) $(CUDA_LIBS)
examples_C_ssids_LINK = $(NVCCLINK)
tests/ssids.$(OBJEXT): libspral.a
examples/C/ssids.$(OBJEXT): libspral.a
TESTS += ssids_test
endif

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
src/ssids/ssids.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
				     	  		   src/match_order.$(OBJEXT) \
									src/matrix_util.$(OBJEXT) \
									src/metis$(METIS_VERSION)_wrapper.$(OBJEXT) \
									src/scaling.$(OBJEXT) \
									src/ssids/alloc.$(OBJEXT) \
									src/ssids/analyse.$(OBJEXT) \
									src/ssids/datatypes.$(OBJEXT) \
									src/ssids/factor_gpu.$(OBJEXT) \
									src/ssids/solve_cpu.$(OBJEXT) \
									src/ssids/solve_gpu.$(OBJEXT)
src/ssids/alloc.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
									src/ssids/datatypes.$(OBJEXT)
src/ssids/analyse.$(OBJEXT): src/ssids/datatypes.$(OBJEXT) \
									  src/core_analyse.$(OBJEXT) \
									  src/cuda/cuda.$(OBJEXT) \
									  src/pgm.$(OBJEXT)
src/ssids/cuda_interfaces.$(OBJEXT): src/ssids/cuda_datatypes.$(OBJEXT)
src/ssids/datatypes.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
										 src/scaling.$(OBJEXT)  \
										 src/ssids/cuda_datatypes.$(OBJEXT)
src/ssids/dense_factor_gpu.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
												  src/ssids/alloc.$(OBJEXT) \
												  src/ssids/cuda_datatypes.$(OBJEXT) \
												  src/ssids/cuda_interfaces.$(OBJEXT) \
												  src/ssids/datatypes.$(OBJEXT)
src/ssids/factor_gpu.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
							 	  		  src/ssids/alloc.$(OBJEXT) \
							 	  		  src/ssids/cuda_datatypes.$(OBJEXT) \
							 	  	 	  src/ssids/cuda_interfaces.$(OBJEXT) \
								  		  src/ssids/datatypes.$(OBJEXT) \
							 	  		  src/ssids/dense_factor_gpu.$(OBJEXT) \
										  src/ssids/solve_gpu.$(OBJEXT)
src/ssids/solve_cpu.$(OBJEXT): src/ssids/datatypes.$(OBJEXT)
src/ssids/solve_gpu.$(OBJEXT): src/cuda/cuda.$(OBJEXT) \
										 src/ssids/alloc.$(OBJEXT) \
								 	    src/ssids/cuda_datatypes.$(OBJEXT) \
								 	    src/ssids/cuda_interfaces.$(OBJEXT) \
								 		 src/ssids/datatypes.$(OBJEXT)
driver/spral_ssids.$(OBJEXT): libspral.a

# CUDA header deps
src/ssids/solve_kernels.$(OBJEXT): src/ssids/dtrsv.h

# C interface deps
interfaces/C/random.$(OBJEXT): src/random.$(OBJEXT)
interfaces/C/random_matrix.$(OBJEXT): src/random_matrix.$(OBJEXT)
interfaces/C/scaling.$(OBJEXT): src/scaling.$(OBJEXT)
interfaces/C/ssids.$(OBJEXT): src/ssids/ssids.$(OBJEXT)
