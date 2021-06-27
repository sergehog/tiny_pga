cmake_minimum_required(VERSION 3.15)
project(tiny_pga)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_subdirectory(tiny_autodf)

enable_testing()
add_subdirectory(test)