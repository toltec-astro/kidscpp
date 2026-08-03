cmake_minimum_required(VERSION 3.25)

project(
    kidscpp
    VERSION 3.1.0
    DESCRIPTION "KIDs timestream processing for TolTEC"
    HOMEPAGE_URL "https://github.com/toltec-astro/kidscpp"
    LANGUAGES CXX
)

include(CTest)
include(GNUInstallDirs)

find_package(TulaCMake 3 CONFIG REQUIRED)
find_package(
    tula 3.1 CONFIG REQUIRED
    COMPONENTS
        core
        logging
        ecsv
        eigen
        netcdf
        enum
        grppi
)

add_library(kidscpp_config INTERFACE)
add_library(kids::config ALIAS kidscpp_config)
set_target_properties(kidscpp_config PROPERTIES EXPORT_NAME config)
tula_cmake_generate_config_header(
    TARGET kidscpp_config
    INPUT "${CMAKE_CURRENT_SOURCE_DIR}/include/kids/config.h.in"
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/kids/config.h"
    BUILD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated"
)

add_library(kids STATIC)
add_library(kids::kids ALIAS kids)
set_target_properties(kids PROPERTIES EXPORT_NAME kids)
tula_cmake_target_defaults(TARGET kids)
set(KIDSCPP_PACKAGE_SPEC "" CACHE STRING "Concrete Spack package spec")
set(KIDSCPP_DAG_HASH "" CACHE STRING "Concrete Spack dependency DAG hash")
set(TOLTECA_BUILD_PROFILE "" CACHE STRING "tolteca_deploy profile identity")
set(TOLTECA_LOCK_SHA256 "" CACHE STRING "SHA-256 of the deployment Spack lock")
tula_cmake_generate_version_header(
    TARGET kids
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/kids/version.h"
    BUILD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated"
    NAMESPACE "kids::build"
    VERSION "${PROJECT_VERSION}"
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    PACKAGE_SPEC "${KIDSCPP_PACKAGE_SPEC}"
    DAG_HASH "${KIDSCPP_DAG_HASH}"
    BUILD_PROFILE "${TOLTECA_BUILD_PROFILE}"
    LOCK_SHA256 "${TOLTECA_LOCK_SHA256}"
)
target_sources(
    kids
    PRIVATE
        src/kids/toltec/timestream.cpp
        src/kids/timestream/solver.cpp
        src/kids/timestream/solverresult.cpp
)
target_include_directories(
    kids
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(
    kids
    PUBLIC
        tula::core
        tula::logging
        tula::ecsv
        tula::eigen
        tula::netcdf
        tula::enum
        tula::grppi
        kids::config
)

install(
    DIRECTORY include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    PATTERN "config.h.in" EXCLUDE
)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/generated/kids/config.h"
        "${CMAKE_CURRENT_BINARY_DIR}/generated/kids/version.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/kids"
)
tula_cmake_install_package(
    PACKAGE kidscpp
    EXPORT kidscppTargets
    NAMESPACE kids::
    VERSION "${PROJECT_VERSION}"
    TARGETS kids kidscpp_config
    CONFIG_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/cmake/kidscppConfig.cmake.in"
)

option(
    KIDS_BUILD_TESTS
    "Build and run Kidscpp behavior tests"
    "${PROJECT_IS_TOP_LEVEL}"
)
if(KIDS_BUILD_TESTS)
    set(
        TOLTECA_TEST_DATA_ROOT
        ""
        CACHE PATH
        "Path to the tolteca_test_data checkout used by real NetCDF tests"
    )
    if(NOT TOLTECA_TEST_DATA_ROOT)
        set(_tolteca_test_data_candidate
            "${CMAKE_CURRENT_SOURCE_DIR}/../tolteca_test_data")
        if(IS_DIRECTORY "${_tolteca_test_data_candidate}")
            cmake_path(
                ABSOLUTE_PATH _tolteca_test_data_candidate
                NORMALIZE
                OUTPUT_VARIABLE TOLTECA_TEST_DATA_ROOT
            )
        endif()
    endif()
    add_subdirectory(tests)
endif()
