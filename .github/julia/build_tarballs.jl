using BinaryBuilder, Pkg

haskey(ENV, "SPRAL_RELEASE") || error("The environment variable SPRAL_RELEASE is not defined.")
haskey(ENV, "SPRAL_COMMIT") || error("The environment variable SPRAL_COMMIT is not defined.")

name = "SPRAL"
version = VersionNumber(ENV["SPRAL_RELEASE"])

# Collection of sources required to complete build
sources = [
    GitSource("https://github.com/ralna/spral.git", ENV["SPRAL_COMMIT"])
]

# Bash recipe for building across all platforms
script = raw"""
# Export dependencies
mkdir ${prefix}/deps
cd ${libdir}
for file in $(ls .); do
   if [[ -f $file ]]; then
      if [[ -z $(ls -la $file | grep 'artifacts') ]]; then
         cp -P ${file} ${prefix}/deps/${file}
      else
         cp -L ${file} ${prefix}/deps/${file}
      fi
   fi
done
cd ${prefix}
cp -rL share/licenses deps/licenses
chmod -R u=rwx deps
tar -czvf deps.tar.gz deps
rm -r deps

# Update Ninja
cp ${host_prefix}/bin/ninja /usr/bin/ninja

cd ${WORKSPACE}/srcdir/spral

if [[ "${target}" == *mingw* ]]; then
  HWLOC="hwloc-15"
else
  HWLOC="hwloc"
fi

meson setup builddir --cross-file=${MESON_TARGET_TOOLCHAIN%.*}_gcc.meson \
                     --prefix=$prefix \
                     -Dlibhwloc=$HWLOC \
                     -Dlibblas=openblas \
                     -Dliblapack=openblas \
                     -Dexamples=true \
                     -Dtests=true

meson compile -C builddir
meson install -C builddir
"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = supported_platforms()
platforms = expand_gfortran_versions(platforms)

# The products that we will ensure are always built
products = [
    LibraryProduct("libspral", :libspral)
]

# Dependencies that must be installed before this package can be built
dependencies = [
    HostBuildDependency(PackageSpec(name="Ninja_jll", uuid="76642167-d241-5cee-8c94-7a494e8cb7b7")),
    Dependency(PackageSpec(name="METIS_jll", uuid="d00139f3-1899-568f-a2f0-47f597d42d70")),
    Dependency(PackageSpec(name="CompilerSupportLibraries_jll", uuid="e66e0078-7015-5450-92f7-15fbd957f2ae")),
    Dependency(PackageSpec(name="OpenBLAS32_jll", uuid="656ef2d0-ae68-5445-9ca0-591084a874a2")),
    Dependency(PackageSpec(name="Hwloc_jll", uuid="e33a78d0-f292-5ffc-b300-72abe9b543c8")),
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; preferred_gcc_version=v"9.1", julia_compat="1.6")
