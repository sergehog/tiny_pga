cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 14)

find_package(GTest CONFIG REQUIRED)

enable_testing()

add_executable(elements_test elements_test.cpp )
add_executable(pga_test pga_test.cpp)
add_executable(autodf_test autodf_test.cpp)
add_executable(auotdf_pga_test auotdf_pga_test.cpp)
add_executable(pga3d_test pga3d_test.cpp)
add_executable(autodf_pga3d_test autodf_pga3d_test.cpp)

target_link_libraries(elements_test gtest gtest_main pthread)
target_link_libraries(pga_test gtest gtest_main pthread)
target_link_libraries(autodf_test gtest gtest_main pthread)
target_link_libraries(auotdf_pga_test gtest gtest_main pthread)
target_link_libraries(pga3d_test gtest gtest_main pthread)
target_link_libraries(autodf_pga3d_test gtest gtest_main pthread)

add_test(NAME elements_test COMMAND elements_test)
add_test(NAME pga_test COMMAND pga_test)
add_test(NAME autodf_test COMMAND autodf_test)
add_test(NAME auotdf_pga_test COMMAND auotdf_pga_test)
add_test(NAME autodf_pga3d_test COMMAND autodf_pga3d_test)

