# kidscpp

`kidscpp` is the C++ TolTEC timestream-processing package. The v3 port keeps
the intentionally trimmed scope from the archived v3 rewrite: timestream
solving remains in C++, while sweep finding and fitting are owned by the
Python pipeline.

The package is a downstream consumer of `tula/3.1.0`; dependency acquisition
and CMake target normalization are provided by `tula_cmake`.

## TolTEC raw-data boundary

Kidscpp owns the narrow NetCDF-to-timestream boundary used by Citlali:

```cpp
#include <kids/toltec/timestream.h>

auto meta = kids::toltec::get_raw_timestream_meta(filepath);
auto raw = kids::toltec::read_raw_timestream_slice(
    filepath, kids::toltec::SampleSlice{start, stop, step});
auto solved = kids::TimeStreamSolver{config}(raw);
```

The reader validates `ObsType=1`, maps TolTEC metadata, reads sliced `Ts`,
`Is`, and `Qs`, and constructs the absolute tone-frequency/model axis needed
by `TimeStreamSolver`. It does not restore the former generic sweep/data
dispatcher, Kidscpp CLI, or sweep fitter.

Production compatibility details are explicit: calibration fit-report names
retain their zero-padded observation/sub-observation/scan pattern, raw files
select the same first tone/model block as v1, and early files may supply the
time axis as `Data.Toltec.Xs` when `Data.Toltec.Ts` is absent.

Citlali owns observation orchestration and chooses sample slices; it calls this
Kidscpp API rather than maintaining a second NetCDF parser.

## Build

```sh
./build
```

The launcher obtains the pinned `tula_cmake` CLI and runs the Conan install
plus generated CMake preset workflow. A configured TolTEC Conan remote
supplies `tula-cmake/3.1.0` and `tula/3.1.0`; Tula is not fetched as a CMake
subproject.

For this multi-repository development workspace:

```sh
TULA_CMAKE_DEV_PROJECT=../tula_cmake ./build
```

The package publishes `kids::kids`. Its Conan `test_package` compiles a
separate consumer after package creation.

When the sibling `tolteca_test_data` repository is available, the workspace
gate sets `TOLTECA_TEST_DATA_ROOT` and verifies the reader against the 2024
`toltec0_018230_111_0000` timestream fixture. Without that environment
variable, only the real-file cases are reported as skipped.
