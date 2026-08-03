find_package(GTest CONFIG REQUIRED)

include(GoogleTest)

foreach(_name toltec_timestream welch solver)
    add_executable("kids_test_${_name}" "test_${_name}.cpp")
    target_link_libraries(
        "kids_test_${_name}"
        PRIVATE kids::kids GTest::gtest_main
    )
    target_compile_features("kids_test_${_name}" PRIVATE cxx_std_23)
    gtest_discover_tests(
        "kids_test_${_name}"
        TEST_PREFIX "kidscpp::"
        PROPERTIES
            ENVIRONMENT
                "TOLTECA_TEST_DATA_ROOT=${TOLTECA_TEST_DATA_ROOT}"
    )
endforeach()

target_compile_definitions(
    kids_test_toltec_timestream
    PRIVATE KIDS_HAS_REAL_TEST_DATA=${KIDS_HAS_REAL_TEST_DATA}
)
