# Kidscpp 3.1

Kidscpp is the TolTEC KIDs timestream-processing library. The
`v3.x_spack` branch is an ordinary CMake package with an owned decentralized
Spack recipe.

## Ownership boundary

Kidscpp owns:

- raw TolTEC NetCDF metadata inspection;
- file-slice to `KidsData<RawTimeStream>` conversion;
- timestream solver/model behavior; and
- the `kids::kids` installed CMake target.

Citlali selects observation files and sample slices, then calls the Kidscpp
reader and solver APIs. Citlali does not duplicate the NetCDF adapter. The old
multipurpose Kidscpp CLI and sweep fitter are not part of this focused library.

## Use

```cmake
find_package(kidscpp 3.1 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE kids::kids)
```

Its exported config discovers Tula and propagates the complete installed target
closure.

The Spack `+openmp` variant (enabled by default) selects the matching Tula
performance closure. `kidscpp~openmp` preserves the same Kidscpp APIs while
building its transitive Tula/GrPPI layer without an OpenMP runtime.

## Development and tests

From the workspace dev container:

```console
spack -e tula_cmake/environments/production/gcc14 \
  install --test=all --overwrite kidscpp
spack -e tula_cmake/environments/production/llvm20 \
  install --test=all --overwrite kidscpp
```

The supported lanes are GCC 14 and LLVM/Clang 20, both C++23. Seven tests cover
metadata and slice ingestion, invalid stride and observation-type behavior,
PSD construction, and the timestream solver. The real tests automatically use
the sibling `tolteca_test_data` checkout; a missing fixture is visible as a
skip rather than silently treated as coverage.

`tests/installed_consumer` independently verifies the installed `kids::kids`
package. Run the complete production chain with `just production` from the
workspace root.
