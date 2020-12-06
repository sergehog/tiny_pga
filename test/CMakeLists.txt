cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 14)

find_package(GTest CONFIG REQUIRED)

add_executable(pga_test pga_test.cpp elements_test.cpp)

target_link_libraries(pga_test gtest gtest_main pthread)

add_test(NAME pga_test COMMAND pga_test)

