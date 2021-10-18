project(kids_gui LANGUAGES CXX C)
if (NOT Qt5_PREFIX)
    # set some heuristic Qt5 paths
    # TODO detect Qt5 versions
    if (APPLE)
        set(Qt5_PREFIX "/Applications/Qt/5.12.2/clang_64")
        message("Use Qt5 from ${Qt5_PREFIX}")
    else()
        message(FATAL_ERROR "please set -DQt5_PREFIX=</path/to/Qt/prefix> to compile kids_gui")
    endif()
endif()
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${Qt5_PREFIX}")
message("module path ${CMAKE_MODULE_PATH}")
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

# Find the QtWidgets library
find_package(Qt5 COMPONENTS Widgets PrintSupport REQUIRED)
add_executable(kids_gui EXCLUDE_FROM_ALL)
target_link_libraries(kids_gui
    PRIVATE
        kids_core
        Qt5::Widgets
        Qt5::PrintSupport
        )
target_sources(kids_gui
    PRIVATE
        main.cpp
        kids_gui.cpp
        plot.cpp
        qcustomplot/qcustomplot.cpp
        # doublerangeslider.cpp
        # sciencespinbox.cpp
        # qcustomplot/qcustomplot.cpp
        # kids_gui.ui
    )
add_custom_target(ui_files SOURCES kids_gui.ui)
