cmake_minimum_required(VERSION 3.1)
project(TP1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ggdb -Isrc/")
add_executable(TP1 src/main.cpp)