cmake_minimum_required(VERSION 3.19)

project(mesh-refinement CXX)

add_subdirectory(algo)
add_subdirectory(ui)

include(CTest)
enable_testing()
add_subdirectory(tests)
