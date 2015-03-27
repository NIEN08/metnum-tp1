cmake_minimum_required(VERSION 3.1)
project(TP1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Isrc/")

set(SOURCE_FILES src/main.cpp src/Matrix.cpp)
add_executable(TP1 ${SOURCE_FILES})