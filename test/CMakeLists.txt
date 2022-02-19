cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 17)

add_executable(elements_test elements_test.cpp )
add_executable(tiny_pga_test tiny_pga_test.cpp)
add_executable(auotdf_pga_test auotdf_pga_test.cpp)
add_executable(pga3d_test pga3d_test.cpp)
add_executable(autodf_pga3d_test autodf_pga3d_test.cpp)

target_link_libraries(elements_test gtest gtest_main pthread)
target_link_libraries(tiny_pga_test gtest gtest_main pthread)
target_link_libraries(auotdf_pga_test gtest gtest_main pthread)
target_link_libraries(pga3d_test gtest gtest_main pthread)
target_link_libraries(autodf_pga3d_test gtest gtest_main pthread)

add_test(NAME elements_test COMMAND elements_test)
add_test(NAME tiny_pga_test COMMAND tiny_pga_test)
add_test(NAME auotdf_pga_test COMMAND auotdf_pga_test)
add_test(NAME autodf_pga3d_test COMMAND autodf_pga3d_test)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")

        set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
        include(CodeCoverage)
        append_coverage_compiler_flags()

        setup_target_for_coverage_gcovr_html(
                NAME tinypga_test_coverage
                EXECUTABLE ctest --verbose
                DEPENDENCIES elements_test tiny_pga_test auotdf_pga_test autodf_pga3d_test
                BASE_DIRECTORY "../"
                EXCLUDE "test" "tiny_autodf/test"
        )

endif()
