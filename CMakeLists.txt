cmake_minimum_required(VERSION 3.20)
project(ParallelProgramming)

set(CMAKE_CXX_STANDARD 11)

set(OMP_COMPILER_FLAG "-fopenmp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OMP_COMPILER_FLAG}")

add_executable(ParallelProgramming calculatePI.cpp)

