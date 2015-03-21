cmake_minimum_required(VERSION 3.1)
project(TP1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(TP1 ${SOURCE_FILES} src/Matrix.cpp src/Matrix.h)