cmake_minimum_required(VERSION 3.1)
project(metnum_tp1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -pedantic -std=c++11 -Wall -ggdb")

set(SOURCE_FILES src/main.cpp)
add_executable(metnum_tp1 ${SOURCE_FILES})