cmake_minimum_required(VERSION 3.15)
project(tiny_pga)

set(CMAKE_CXX_STANDARD 14)

add_executable(tiny_pga main.cpp tiny_pga.h)

add_subdirectory(test)