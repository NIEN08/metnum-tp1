cmake_minimum_required(VERSION 3.1)
project(TP1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Isrc/")

set(SOURCE_FILES src/main.cpp src/Problem.cpp)
add_executable(TP1 ${SOURCE_FILES})