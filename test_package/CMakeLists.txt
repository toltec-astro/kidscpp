cmake_minimum_required(VERSION 3.25)
project(kidscpp_package_test LANGUAGES CXX)

find_package(kidscpp CONFIG REQUIRED)

add_executable(kidscpp_package_test src/main.cpp)
target_link_libraries(kidscpp_package_test PRIVATE kids::kids)
target_compile_features(kidscpp_package_test PRIVATE cxx_std_23)
