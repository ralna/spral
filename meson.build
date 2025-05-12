project(
  'SPRAL',
  'c', 'cpp', 'fortran',
  version: '2025.03.06',
  license: 'BSD-3',
  meson_version: '>= 0.63.0',
  default_options: [
    'buildtype=release',
    'libdir=lib',
    'default_library=shared',
    'warning_level=0',
    'c_std=c99',
    'cpp_std=c++11',
  ],
)

# Compilers
cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')
fc = meson.get_compiler('fortran')
nvcc_available = add_languages('cuda', required: false, native: false)

# Remove messages about deprecated Intel compilers
if cc.get_id() == 'intel'
  add_global_arguments('-diag-disable=10441', language : 'c')
  add_global_link_arguments('-diag-disable=10441', language : 'c')
endif
if cc.get_id() == 'intel-cl'
  add_global_arguments('/Qdiag-disable=10441', language : 'c')
  add_global_link_arguments('/Qdiag-disable=10441', language : 'c')
endif
if cxx.get_id() == 'intel'
  add_global_arguments('-diag-disable=10441', language : 'cpp')
  add_global_link_arguments('-diag-disable=10441', language : 'cpp')
endif
if cxx.get_id() == 'intel-cl'
  add_global_arguments('/Qdiag-disable=10441', language : 'cpp')
  add_global_link_arguments('/Qdiag-disable=10441', language : 'cpp')
endif

# Recognise old non-standard double complex intrinsics
if fc.get_id() == 'nagfor'
  add_global_arguments('-dcfuns', language : 'fortran')
endif

# Options
install_modules = get_option('modules')
build_gpu = get_option('gpu')
build_tests = get_option('tests')
build_examples = get_option('examples')
build_openmp = get_option('openmp')

libblas_name = get_option('libblas')
libblas_path = get_option('libblas_path')
libblas_include = include_directories(get_option('libblas_include'))

liblapack_name = get_option('liblapack')
liblapack_path = get_option('liblapack_path')

libhwloc_name = get_option('libhwloc')
libhwloc_path = get_option('libhwloc_path')
libhwloc_include = include_directories(get_option('libhwloc_include'))

libmetis_name = get_option('libmetis')
libmetis_path = get_option('libmetis_path')
libmetis_version = get_option('libmetis_version')
metis64 = get_option('metis64')

# Dependencies
if libblas_path == []
    libblas = dependency(libblas_name, required : false)
endif
if libblas_path != [] or not libblas.found()
    libblas = fc.find_library(libblas_name, dirs : libblas_path, required : true)
endif
if liblapack_path == []
  liblapack = dependency(liblapack_name, required : false)
endif
if liblapack_path != [] or not liblapack.found()
  liblapack = fc.find_library(liblapack_name, dirs : liblapack_path, required : true)
endif
libmetis = fc.find_library(libmetis_name, dirs : libmetis_path, required : true)
libhwloc = fc.find_library(libhwloc_name, dirs : libhwloc_path, required : false)
libcuda = dependency('cuda', version : '>=10', modules : ['cublas'], required : false)
lm = cc.find_library('m', required : false)
libspral_deps = [libblas, liblapack, libmetis, libhwloc, libcuda, lm]

# CBLAS
has_cblash = cc.has_header('cblas.h', include_directories : libblas_include)

# METIS
if metis64
  add_global_arguments('-DSPRAL_HAVE_METIS_H=1', language : 'fortran')
  add_global_arguments('-DSIZEOF_IDX_T=8', language : 'fortran')
endif

# HWLOC
has_hwloch = cc.has_header('hwloc.h', include_directories : libhwloc_include)
if libhwloc.found() and has_hwloch
  add_global_arguments('-DHAVE_HWLOC', language : 'cpp')
endif

# SCHED_GETCPU
if host_machine.system() == 'linux'
  add_global_arguments('-DHAVE_SCHED_GETCPU', language : 'cpp')
endif

# OpenMP
if build_openmp
  if fc.get_id() == 'nvidia_hpc'
    add_global_arguments('-mp', language : 'fortran')
  elif fc.get_id() == 'nagfor'
    add_global_arguments('-openmp', language : 'fortran')
  elif fc.get_id() == 'gcc' or fc.get_id() == 'llvm-flang'
    add_global_arguments('-fopenmp', language : 'fortran')
  elif fc.get_id() == 'intel' or fc.get_id() == 'intel-llvm'
    add_global_arguments('-qopenmp', language : 'fortran')
  elif fc.get_id() == 'intel-cl' or fc.get_id() == 'intel-llvm-cl'
    add_global_arguments('/Qopenmp', language : 'fortran')
  endif

  if cxx.get_id() == 'nvidia_hpc'
    add_global_arguments('-mp', language : 'cpp')
  elif cxx.get_id() == 'gcc' or cxx.get_id() == 'clang' or cxx.get_id() == 'clang-cl'
    add_global_arguments('-fopenmp', language : 'cpp')
  elif cxx.get_id() == 'intel' or cxx.get_id() == 'intel-llvm'
    add_global_arguments('-qopenmp', language : 'cpp')
  elif cxx.get_id() == 'intel-cl' or cxx.get_id() == 'intel-llvm-cl'
    add_global_arguments('/Qopenmp', language : 'cpp')
  endif

  lomp = '-lgomp'
  if cxx.get_id() == 'intel' or cxx.get_id() == 'intel-llvm'
    lomp = '-liomp5'
  endif
  if cxx.get_id() == 'nvidia_hpc'
    lomp = '-lomp'
  endif
  if host_machine.system() == 'darwin' or host_machine.system() == 'freebsd'
    if cxx.get_id() == 'clang'
      lomp = '-lomp'
    endif
  endif
  add_global_link_arguments(lomp, language : 'fortran')
  add_global_link_arguments(lomp, language : 'c')
  add_global_link_arguments(lomp, language : 'cpp')
endif

# Link flag for C++
lstdcpp = '-lstdc++'
if host_machine.system() != 'linux'
  if cxx.get_id() == 'clang' or cxx.get_id() == 'clang-cl'
    lstdcpp = '-lc++'
  endif
  if cxx.get_id() == 'intel' or cxx.get_id() == 'intel-llvm'
    lstdcpp = '-cxxlib'
  endif
  if cxx.get_id() == 'intel-cl' or cxx.get_id() == 'intel-llvm-cl'
    lstdcpp = '/Qcxxlib'
  endif
endif

binspral_src = []
libspral_src = []
libspral_cpp_src = []
libspral_nvcc_src = []

spral_examples = []
spral_c_examples = []

spral_tests = []
spral_c_tests = []
spral_cpp_tests = []

# Headers
spral_headers = []
libspral_include = []
libspral_include += include_directories('include', 'src')
libspral_include += libhwloc_include
libspral_include += libblas_include

# Sources
subdir('include')
subdir('interfaces/C')
subdir('src')
subdir('driver')
subdir('examples')
subdir('tests')

# Library
libspral = library('spral',
                   sources : libspral_src + libspral_cpp_src + libspral_nvcc_src,
                   dependencies : libspral_deps,
                   link_language : 'fortran',
                   link_args : lstdcpp,
                   include_directories: libspral_include,
                   install : true)

# Binary
executable('spral_ssids',
           sources : binspral_src,
           dependencies : libcuda,
           link_with : libspral,
           link_language : 'fortran',
           install : true)

# Headers
install_headers(spral_headers)

# Fortran modules
if install_modules
  script_modules = files('install_modules.py')
  meson.add_install_script(script_modules)
endif

# OMP-related environment variables
omp_env = ['OMP_CANCELLATION=true', 'OMP_PROC_BIND=true']

# Tests
if build_tests

  foreach test: spral_tests
    package = test[0]
    name = test[1]
    file = test[2]
    test(name,
         executable(name, file, link_with : libspral, dependencies : libspral_deps,
                    link_language : 'fortran', include_directories: libspral_include),
         suite : [package, 'tests', 'fortran'],
         timeout : 300, is_parallel : false, env: omp_env)
  endforeach

  foreach test: spral_c_tests
    package = test[0]
    name = test[1]
    file = test[2]
    test(name,
         executable(name, file, link_with : libspral, dependencies : libspral_deps,
                    link_language : 'c', include_directories : libspral_include),
         suite : [package, 'tests', 'C'],
         timeout : 300, is_parallel : false, env: omp_env)
  endforeach

  foreach test: spral_cpp_tests
    package = test[0]
    name = test[1]
    file = test[2]
    test(name,
         executable(name, file, link_with : libspral, dependencies : libspral_deps,
                    link_language : 'cpp', include_directories : libspral_include),
         suite : [package, 'tests', 'C++'],
         timeout : 300, is_parallel : false, env: omp_env)
  endforeach
endif

# Examples
if build_examples

  fortran_examples_folder = 'examples/Fortran'

  foreach example: spral_examples
    package = example[0]
    name = example[1]
    file = example[2]
    test(name,
         executable(name, file, link_with : libspral, dependencies : libspral_deps,
                    link_language : 'fortran', include_directories : libspral_include,
                    install : true, install_dir : fortran_examples_folder),
         suite : [package, 'examples', 'fortran'],
         timeout : 300, is_parallel : false, env: omp_env)
  endforeach

  c_examples_folder = 'examples/C'

  foreach example: spral_c_examples
    package = example[0]
    name = example[1]
    file = example[2]
    test(name,
         executable(name, file, link_with : libspral, dependencies : libspral_deps,
                    link_language : 'c', include_directories : libspral_include,
                    install : true, install_dir : c_examples_folder),
         suite : [package, 'examples', 'C'],
         timeout : 300, is_parallel : false, env: omp_env)
  endforeach
endif
