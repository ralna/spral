project(
  'SPRAL',
  'c', 'cpp', 'fortran',
  version: '2023.8.2',
  license: 'BSD-3',
  meson_version: '>= 0.63.0',
  default_options: [
    'buildtype=debug',
    'libdir=lib',
    'default_library=shared',
    'warning_level=0',
    'c_std=c99',
    'cpp_std=c++17',
  ],
)

cc = meson.get_compiler('c')
cxx = meson.get_compiler('cpp')
fc = meson.get_compiler('fortran')

# Remove messages about deprecated Intel compilers
if cc.get_id() == 'intel' or cc.get_id() == 'intel-cl'
  add_global_arguments('-diag-disable=10441', language : 'c')
  add_global_link_arguments('-diag-disable=10441', language : 'c')
endif
if cxx.get_id() == 'intel' or cxx.get_id() == 'intel-cl'
  add_global_arguments('-diag-disable=10441', language : 'cpp')
  add_global_link_arguments('-diag-disable=10441', language : 'cpp')
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

libblas_name = get_option('libblas')
libblas_path = get_option('libblas_path')

liblapack_name = get_option('liblapack')
liblapack_path = get_option('liblapack_path')

libmetis_name = get_option('libmetis')
libmetis_path = get_option('libmetis_path')
libmetis_version = get_option('libmetis_version')

# Dependencies
libblas = fc.find_library(libblas_name, dirs : libblas_path, required : true)
liblapack = fc.find_library(liblapack_name, dirs : liblapack_path, required : true)
libmetis = fc.find_library(libmetis_name, dirs : libmetis_path, required : true)
libhwloc = dependency('hwloc', required : true)
libcuda = dependency('cuda', version : '>=10', modules : ['cublas'], required : false)
# lm = cc.find_library('m', required : false)

libspral_deps = [libblas, liblapack, libmetis, libhwloc, libcuda]

# OpenMP
if fc.get_id() == 'nagfor'
  add_global_arguments('-openmp', language : 'fortran')
  add_global_link_arguments('-openmp', language : 'fortran')
elif fc.get_id() == 'gcc'
  add_global_arguments('-fopenmp', language : 'fortran')
  add_global_link_arguments('-fopenmp', language : 'fortran')
elif fc.get_id() == 'intel' or fc.get_id() == 'intel-llvm' or fc.get_id() == 'intel-cl' or fc.get_id() == 'intel-llvm-cl'
  add_global_arguments('-qopenmp', language : 'fortran')
  add_global_link_arguments('-qopenmp', language : 'fortran')
endif

if cxx.get_id() == 'gcc' or cxx.get_id() == 'clang'
  add_global_arguments('-fopenmp', language : 'cpp')
elif cxx.get_id() == 'intel' or cxx.get_id() == 'intel-cl' or cxx.get_id() == 'intel-llvm' or cxx.get_id() == 'intel-llvm-cl'
  add_global_arguments('-qopenmp', language : 'cpp')
endif


libspral_src = []
libspral_cpp_src = []
libspral_nvcc_src = []

spral_examples = []
spral_c_examples = []

spral_tests = []
spral_c_tests = []

# Headers
spral_headers = []
libspral_include = []
libspral_include += include_directories('src')

# Sources
subdir('include')
subdir('interfaces/C')
subdir('src')
subdir('driver')

# subdir('examples')
# subdir('tests')

# Library
libspral = library('spral',
                   sources : libspral_src + libspral_cpp_src,
                   dependencies : libspral_deps,
                   link_language : 'fortran',
                   link_args : '-lstdc++',
                   include_directories: libspral_include,
                   install : true)

# Headers
install_headers(spral_headers)

# Fortran modules
if install_modules
  script_modules = files('install_modules.py')
  meson.add_install_script(script_modules)
endif

# Tests
if build_tests

  fortran_tests_folder = 'tests/Fortran'

  foreach test: spral_tests
    name = test[0]
    file = test[1]
    test(name,
         executable(name, file, link_with : libspral, dependencies : libspral_deps, link_language : 'fortran',
                    include_directories: libspral_include , install : true, install_dir : fortran_tests_folder),
         is_parallel : false)
  endforeach

  c_tests_folder = 'tests/C'

  foreach test: spral_c_tests
    name = test[0]
    file = test[1]
    test(name,
         executable(name, file, link_with : libspral, dependencies : libspral_deps, link_language : 'c',
                    include_directories : libspral_include, install : true, install_dir : c_tests_folder),
         is_parallel : false)
  endforeach
endif

# Examples
if build_examples

  fortran_examples_folder = 'examples/Fortran'

  foreach example: spral_examples
    name = example[0]
    file = example[1]
    executable(name, file, link_with : libspral, dependencies : libspral_deps, link_language : 'fortran',
               include_directories : libspral_include, install : true, install_dir : fortran_examples_folder)
  endforeach

  c_examples_folder = 'examples/C'

  foreach example: spral_c_examples
    name = example[0]
    file = example[1]
    executable(name, file, link_with : libspral, dependencies : libspral_deps, link_language : 'c',
               include_directories : libspral_include, install : true, install_dir : c_examples_folder)
  endforeach
endif
