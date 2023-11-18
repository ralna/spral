# Platforms:
#   - aarch64-apple-darwin
#   - aarch64-linux-gnu
#   - aarch64-linux-musl
#   - armv6l-linux-gnueabihf
#   - armv6l-linux-musleabihf
#   - armv7l-linux-gnueabihf
#   - armv7l-linux-musleabihf
#   - i686-linux-gnu
#   - i686-linux-musl
#   - i686-w64-mingw32
#   - powerpc64le-linux-gnu
#   - x86_64-apple-darwin
#   - x86_64-linux-gnu
#   - x86_64-linux-musl
#   - x86_64-unknown-freebsd
#   - x86_64-w64-mingw32
#
# Targets:
#   - platform-libgfortran3
#   - platform-libgfortran4
#   - platform-libgfortran5
#
# Examples:
# julia --color=yes build_tarballs.jl x86_64-linux-gnu-libgfortran5 --verbose
# julia --color=yes build_tarballs.jl x86_64-w64-mingw32-libgfortran5 --verbose --debug
# julia --color=yes build_tarballs.jl x86_64-apple-darwin-libgfortran5 --verbose
# julia --color=yes build_tarballs.jl aarch64-apple-darwin-libgfortran5 --verbose
using BinaryBuilder, Pkg

name = "SPRAL"
version = v"2023.11.18"

# Collection of sources required to complete build
sources = [
    GitSource("https://github.com/ralna/spral.git", "e723071ce2e0e6181bb65e1b365dc47449e1a912")
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

# Install a version of Meson ≥ 0.63.0
python3 -m pip install --user --upgrade meson

cd ${WORKSPACE}/srcdir/spral

if [[ "${target}" == *mingw* ]]; then
  HWLOC="hwloc-15"
else
  HWLOC="hwloc"
fi

meson setup builddir --cross-file=${MESON_TARGET_TOOLCHAIN%.*}_gcc.meson \
                     --prefix=$prefix -Dlibhwloc=$HWLOC \
                     -Dlibblas=openblas -Dliblapack=openblas \
                     -Dexamples=true -Dtests=true

for i in {1..10}
do
    meson compile -C builddir || true
done
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
    Dependency(PackageSpec(name="METIS_jll", uuid="d00139f3-1899-568f-a2f0-47f597d42d70")),
    Dependency(PackageSpec(name="CompilerSupportLibraries_jll", uuid="e66e0078-7015-5450-92f7-15fbd957f2ae")),
    Dependency(PackageSpec(name="OpenBLAS32_jll", uuid="656ef2d0-ae68-5445-9ca0-591084a874a2")),
    Dependency(PackageSpec(name="Hwloc_jll", uuid="e33a78d0-f292-5ffc-b300-72abe9b543c8")),
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; julia_compat="1.6")
