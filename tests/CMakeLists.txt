find_package(GTest CONFIG REQUIRED)

include(GoogleTest)

foreach(_name welch solver)
    add_executable("kids_test_${_name}" "test_${_name}.cpp")
    target_link_libraries(
        "kids_test_${_name}"
        PRIVATE kids::kids GTest::gtest_main
    )
    target_compile_features("kids_test_${_name}" PRIVATE cxx_std_23)
    gtest_discover_tests(
        "kids_test_${_name}"
        TEST_PREFIX "kidscpp::"
    )
endforeach()
