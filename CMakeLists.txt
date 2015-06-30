cmake_minimum_required(VERSION 2.7)
project(metnum_tp1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -pedantic -std=c++11 -Wall -ggdb")

set(SOURCE_FILES src/main.cpp)
add_executable(tp ${SOURCE_FILES})
