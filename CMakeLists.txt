cmake_minimum_required (VERSION 3.3)

project (UltrasonicLevitator)

set (CMAKE_CXX_STANDARD 14)

set( CMAKE_FLAGS_RELEASE -O3 -Wall )

file( GLOB_RECURSE src "src/*.cpp" )

include_directories (eigen/)
include_directories (~/stator/)

# # for a c executable
# add_executable (main main.cpp ${src})
# target_link_libraries(main "/usr/local/lib/libnlopt.so")

add_subdirectory(pybind11)
pybind11_add_module(ultrasonics src/pybind/pybinder.cpp ${src})
target_link_libraries(ultrasonics PRIVATE "/usr/local/lib/libnlopt.so")