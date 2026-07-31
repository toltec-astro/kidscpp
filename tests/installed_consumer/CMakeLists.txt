cmake_minimum_required(VERSION 3.25)
project(kidscpp_installed_consumer LANGUAGES CXX)

find_package(kidscpp CONFIG REQUIRED)

add_executable(kidscpp_installed_consumer main.cpp)
target_link_libraries(kidscpp_installed_consumer PRIVATE kids::kids)
target_compile_features(kidscpp_installed_consumer PRIVATE cxx_std_23)

enable_testing()
add_test(NAME kidscpp_installed_consumer COMMAND kidscpp_installed_consumer)
