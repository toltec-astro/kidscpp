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

`spack_repo/develop.yaml` declares the local `kidscpp` development spec and
source path. Deployment environments compose this repository-owned metadata
with the equivalent declarations from TulaCMake, Tula, and Citlali.

The Spack `+openmp` variant (enabled by default) selects the matching Tula
performance closure. `kidscpp~openmp` preserves the same Kidscpp APIs while
building its transitive Tula/GrPPI layer without an OpenMP runtime.

## Development and tests

From an activated `tolteca_deploy` development location:

```console
cd ../toltec_astro_dev
source dotbashrc
just cpp-install
spack -e "$TOLTECA_CPP_ENV" install --test=root --overwrite kidscpp
```

Select `development/linux-gcc14` or `development/linux-llvm20` in
`location.yaml`; both use C++23. The complete development suite has seven tests covering
metadata and slice ingestion, invalid stride and observation-type behavior,
PSD construction, and the timestream solver. When the exact sibling
`tolteca_test_data` fixture exists, configuration enables the real reader
case. Source-release builds retain generated NetCDF coverage without reporting
the unavailable large-data case as a skipped success.

`tests/installed_consumer` independently verifies the installed `kids::kids`
package. TulaCMake's focused matrix recipes remain available for package-level
regression work; deployment and full-chain installation use the location.

The installed `<kids/version.h>` reports artifact identity: semantic/source
version, tree state, compiler, C++ standard, package spec, and DAG hash.
Deployment profile and lock identity are runtime properties and are not
embedded in a shareable Kidscpp installation. The current clean-source recipe
selects immutable commit
`06b3130ba7f6f96b509011c5855f635f76a25087`; final tag naming is intentionally
deferred.
