cmake_minimum_required(VERSION 3.20)
project(kidscpp
    DESCRIPTION
        "KIDs data processing for TolTEC"
    HOMEPAGE_URL
        "https://toltec-astro.github.io/kidscpp"
    )

# Fetch tula
include(FetchContent)
FetchContent_Declare(
  tula
  GIT_REPOSITORY https://github.com/toltec-astro/tula.git
  GIT_TAG        "main"
)
FetchContent_MakeAvailable(tula)

# Add the tula cmake modules
set(CMAKE_MODULE_PATH "${tula_SOURCE_DIR}/tula_cmake/cmake" ${CMAKE_MODULE_PATH})
include(common_settings)

# Check if this project is embedded as sub project.
if(NOT DEFINED KIDSCPP_STANDALONE)
    if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
        set(KIDSCPP_STANDALONE ON)
    else()
        set(KIDSCPP_STANDALONE OFF)
    endif()
endif()

# Add the needed tula modules
set(dep_libs tula::headers)

include(perflibs)
set(dep_libs ${dep_libs} tula::perflibs)

include(logging)
set(dep_libs ${dep_libs} tula::logging)

include(Enum)
set(dep_libs ${dep_libs} tula::Enum)

include(Clipp)
set(dep_libs ${dep_libs} tula::Clipp)

include(Yaml)
set(dep_libs ${dep_libs} tula::Yaml)

include(NetCDFCXX4)
set(dep_libs ${dep_libs} tula::NetCDFCXX4)

include(Eigen3)
set(dep_libs ${dep_libs} tula::Eigen3)

include(Ceres)
set(dep_libs ${dep_libs} tula::Ceres)

include(Grppi)
set(dep_libs ${dep_libs} tula::Grppi)

include(Csv)
set(dep_libs ${dep_libs} tula::Csv)

include(Re2)
set(dep_libs ${dep_libs} tula::Re2)

add_library(kids STATIC)
target_sources(kids
    PRIVATE
        # "src/kids/proc/proc.cpp"
        "src/kids/core/io.cpp"
        "src/kids/toltec/ncfile.cpp"
        "src/kids/sweep/finder.cpp"
        "src/kids/sweep/finderresult.cpp"
        "src/kids/sweep/fitter.cpp"
        "src/kids/sweep/fitresult.cpp"
        "src/kids/timestream/solver.cpp"
        "src/kids/timestream/solverresult.cpp"
        )
target_include_directories(kids PUBLIC "include")
target_link_libraries(kids PUBLIC ${dep_libs})

option(KIDS_WITH_PLOTTING "Enable plotting" OFF)
if (KIDS_WITH_PLOTTING)

    include(Matplotlibcpp)
    target_link_libraries(kids PUBLIC tula::Matplotlibcpp)
    target_compile_definitions(kids PUBLIC WITH_PLOTTING)
endif()

add_library(kids::kids ALIAS kids)

## optional targets
option(KIDS_BUILD_CLI "Build CLI" ON)
if (KIDS_BUILD_CLI)

    set(dep_libs kids::kids)

    include(Clipp)
    set(dep_libs ${dep_libs} tula::Clipp)

    add_executable(kids_cli)
    target_sources(kids_cli
        PRIVATE
            "src/kids/cli/main.cpp"
            )
    target_link_libraries(kids_cli PRIVATE ${dep_libs})
    target_include_directories(kids_cli PUBLIC "include")
    set_property(TARGET kids_cli PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
    set_property(TARGET kids_cli PROPERTY OUTPUT_NAME kids)
endif()

option(KIDS_BUILD_GUI "Build GUI" OFF)
if (KIDS_BUILD_GUI)
    add_subdirectory(src/kids/gui)
endif()

if (KIDSCPP_STANDALONE)
    # Tests, examples, and docs
    option(KIDSCPP_BUILD_TESTS "Build tests" ON)
    if (KIDSCPP_BUILD_TESTS)
        add_subdirectory(tests)
    endif()

    option(KIDSCPP_BUILD_EXAMPLES "Build example programs" OFF)
    if (KIDSCPP_BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()

    option(KIDSCPP_BUILD_DOC "Build docs" OFF)
    if (KIDSCPP_BUILD_DOC)
        include(docxgen_target)
    endif()
endif()

# Create concfig headers
generate_gitversion_header(${PROJECT_NAME})
generate_config_header(${PROJECT_NAME})
