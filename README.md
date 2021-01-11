# kidscpp

`kidscpp` processes TolTEC LeKID/ROACH-2 detector signal.

Author: Zhiyuan Ma (2019)

License: BSD-3-Clause


## Build

0. Install dependencies. The two tested platforms are macOS and Ubuntu 20.04

    1. CMake >= 3.15
        * On macOS:

            ```bash
            $ brew install cmake
            ```
        * On Ubuntu 20.04:

            Follow instructions here: https://apt.kitware.com

    2. NetCDF (TolTEC raw data file format)
        * On macOS:

            ```bash
            $ brew install netcdf
            ```
        * On Ubuntu 20.04:

            ```bash
            $ apt install libnetcdf-dev
            ```

    3. C++ compiler (tested LLVM >=10 with GCC >= 9)
        * On macOS:

            ```bash
            $ brew install llvm
            ```
        * On Ubuntu 20.04:

            ```bash
            $ apt install gcc-10 g++-10
            # make gcc-10 the default compiler
            $ update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 1000 --slave /usr/bin/g++ g++ /usr/bin/g++-10
            ```

    4. Python with matplotlib (TODO: make this optional).

1. Clone the repo. Note the second line is needed initialize the submodules.
    ```bash
    $ git clone git@github.com:toltec-astro/kidscpp.git
    $ cd kidscpp
    $ git submodule update --recursive --init
    ```

2. Build

    ```bash
    # in the kidscpp repo
    $ mkdir build
    $ cd build
    $ cmake .. -DCMAKE_BUILD_TYPE=Release
    $ make
    ```

    It is likely that one needs to specify additional commandline arguments to
    the cmake command above.

    For example, on macOS, the llvm clang compiler installed via Homebrew may
    not be in the PATH, therefore the following is needed:
    ```bash
    $ cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++
    ```

3. Run

    The build will create the executable named `kids` in `build/bin/`.
    To see the command line options, run `./bin/kids -h`.
