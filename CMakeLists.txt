cmake_minimum_required(VERSION 3.12)
project(kids LANGUAGES CXX C)
include("${CMAKE_CURRENT_SOURCE_DIR}/common_utils/cmake/CommonSettings.cmake")
# dependencies
find_package(OpenMP REQUIRED)
include(FileSystem)
include(SpdlogAndFmt)
include(Re2)
include(Grppi)
include(Clipp)
# include(Cli11)
# include(Uriparser)
# default to use installed Eigen3 and ceres-solver
option(USE_INSTALLED_NETCDFCXX4 "Use Installed NetCDFCXX4" ON)
option(USE_INSTALLED_EIGEN3 "Use Installed Eigen3" ON)
option(USE_INSTALLED_CERES "Use Installed ceres-solver" ON)
include(NetCDFCXX4)
include(Eigen3)
include(Ceres)
# plotting
find_package(Python REQUIRED COMPONENTS Development Interpreter)
find_package(NumPy REQUIRED)
if (VERBOSE)
    print_target_properties(spdlog::spdlog)
    print_target_properties(google::re2)
    print_target_properties(Eigen3::Eigen)
    print_target_properties(ceres)
endif()
# gitversion
include(gitversion)
GenVersionHeader(${CMAKE_CURRENT_BINARY_DIR})
# targets
add_library(kids_core STATIC)
target_sources(kids_core
    PRIVATE
        "src/core/io.cpp"
        "src/toltec/dirfile.cpp"
        "src/toltec/ncfile.cpp"
        "src/sweep/finder.cpp"
        "src/sweep/finderresult.cpp"
        "src/sweep/fitter.cpp"
        "src/sweep/fitresult.cpp"
        "src/timestream/solver.cpp"
        "src/timestream/solverresult.cpp"
        )

target_link_libraries(kids_core
    PUBLIC
        cmake_utils::spdlog_and_fmt
        google::re2
        grppi::grppi
        Eigen3::Eigen
        ceres::ceres
        cpp::filesystem
        NetCDF::NetCDFCXX4
        Python::Python
    )
target_include_directories(kids_core
    PUBLIC
        # ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${NUMPY_INCLUDE_DIRS}
    )
add_executable(kids)
target_sources(kids
    PRIVATE
        "src/cli_main.cpp"
        )
target_link_libraries(kids
    PRIVATE
        kids_core
        utils::gitversion
        clipp::clipp
    )
set_property(TARGET kids PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
add_custom_target(related_files SOURCES README.md)
# gui
option(KIDS_BUILD_GUI "Build GUI" OFF)
if (KIDS_BUILD_GUI)
    add_subdirectory(src/gui)
endif()
# testing
option(KIDS_BUILD_TESTS "Build tests" ON)
if (KIDS_BUILD_TESTS)
    include(GBenchAndGTest)
    add_subdirectory(tests)
endif()
# doc
option(KIDS_BUILD_DOCS "Build docs" ON)
if (KIDS_BUILD_DOCS)
    include(DoxygenDocTarget)
endif()
