AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

# NVCC setup
SUFFIXES = .cu
PTX_FLAGS = -v
AM_NVCC_FLAGS = -Iinclude -arch=sm_20 -g -Xptxas="${PTX_FLAGS}"
AM_LD_FLAGS = -lcuda
NVCCLINK = \
	$(NVCC) $(NVCC_FLAGS) $(AM_NVCC_FLAGS) $(AM_LD_FLAGS) $(LDFLAGS) -o $@
.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_NVCC_FLAGS) -dc -o $@ $<

# FIXME: make below configurable
LDADD = -L. -lspral $(METIS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FCLIBS) -lrt
CUDA_LIBS = -lcublas

# Specify test codes to run
TESTS = \
	ssids_test \
	random_matrix_test

# Compile code to library
lib_LIBRARIES = libspral.a
libspral_a_SOURCES = \
	src/core_analyse.f90 \
	src/cuda/api_wrappers.cu \
	src/cuda/cuda.f90 \
	src/cuda/cuda_check.h \
	src/match_order.f90 \
	src/matrix_util.f90 \
	src/metis$(METIS_VERSION)_wrapper.f90 \
	src/pgm.f90 \
	src/random.f90 \
	src/random_matrix.f90 \
	src/rutherford_boeing.f90 \
  	src/scaling.f90 \
  	src/timer.f90 \
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
	src/ssids/syrk_kernels.cu

# Build a driver program
bin_PROGRAMS = \
	run_prob
run_prob_SOURCES = \
	driver/run_prob.f90
run_prob_LINK = $(NVCCLINK)
run_prob_LDADD = $(LDADD) $(CUDA_LIBS)

# Test programs
check_PROGRAMS = \
	ssids_test \
	random_matrix_test \
	examples/Fortran/random_matrix \
	examples/Fortran/ssids
ssids_test_SOURCES = tests/ssids.f90
ssids_test_LINK = $(NVCCLINK)
ssids_test_LDADD = $(LDADD) $(CUDA_LIBS)
random_matrix_test_SOURCES = tests/random_matrix.f90
examples_Fortran_random_matrix_SOURCES = examples/Fortran/random_matrix.f90
examples_Fortran_ssids_SOURCES = examples/Fortran/ssids.f90
examples_Fortran_ssids_LDADD = $(LDADD) $(CUDA_LIBS)
examples_Fortran_ssids_LINK = $(NVCCLINK)

# Automake doesn't know much about modules, so get rid of them ourself
clean-local:
	-rm *.mod

# Fortran 90 dependencies
driver/run_prob.$(OBJEXT): libspral.a
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
tests/ssids.$(OBJEXT): libspral.a
tests/random_matrix.$(OBJEXT): libspral.a

# CUDA header deps
src/ssids/solve_kernels.$(OBJEXT): src/ssids/dtrsv.h
