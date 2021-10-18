project (kids_test LANGUAGES CXX C)

include(testing)

add_executable(kids_test EXCLUDE_FROM_ALL)
set_target_properties(kids_test
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
target_sources(kids_test
    PRIVATE
        main.cpp
        # test_utils.cpp
        test_toltec_ncfile.cpp
        test_welch.cpp
    )
target_link_libraries(kids_test
    PRIVATE
        kids::kids
        tula::testing
    )
add_dependencies(check kids_test)
gtest_discover_tests(
    kids_test
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        TEST_PREFIX "kids::"
        DISCOVERY_MODE PRE_TEST
    )
