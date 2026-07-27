cmake_minimum_required(VERSION 3.25)

project(kidscpp
    VERSION 3.1.0
    DESCRIPTION "KIDs timestream processing for TolTEC"
    LANGUAGES CXX
)

include(TulaProject)
tula_resolve_features()
find_package(tula CONFIG REQUIRED)
tula_add_config_header(
    TARGET kidscpp_config
    ALIAS kids::config
    PREFIX KIDSCPP
    OUTPUT kids/config.h
)

add_library(kids STATIC)
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
        "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
    kids
    PUBLIC
        tula::headers
        tula::logging
        tula::yaml_cpp
        tula::csv_parser
        tula::netcdf_cxx4
        tula::bitmask
        tula::meta_enum
        tula::perflibs
        tula::eigen
        tula::grppi
        kids::config
)
target_compile_features(kids PUBLIC cxx_std_23)
add_library(kids::kids ALIAS kids)

option(KIDS_BUILD_TESTS "Build kidscpp behavior tests" "${PROJECT_IS_TOP_LEVEL}")
if(KIDS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

install(
    TARGETS kids
    EXPORT kidscpp_targets
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/generated/kids/config.h"
    DESTINATION include/kids
)
