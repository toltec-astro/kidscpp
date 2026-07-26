# kidscpp

`kidscpp` is the C++ TolTEC timestream-processing package. The v3 port keeps
the intentionally trimmed scope from the archived v3 rewrite: timestream
solving remains in C++, while sweep finding and fitting are owned by the
Python pipeline.

The package is a downstream consumer of `tula/3.1.0`; dependency acquisition
and CMake target normalization are provided by `tula_cmake`.

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
