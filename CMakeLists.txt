cmake_minimum_required (VERSION 3.3)

project (UltrasonicLevitator)

set (CMAKE_CXX_STANDARD 14)

set( CMAKE_CXX_FLAGS "-O3 -Wall -DHAVE_CSTDDEF" )

file( GLOB_RECURSE src "src/*.cpp" )
file( GLOB_RECURSE remove_pybinders "src/pybind/*.cpp" )
list( REMOVE_ITEM src ${remove_pybinders} )
list( APPEND ${src} "src/pybind/pybinder.cpp" )

include_directories (eigen/)
include_directories (LBFGSpp/include/)
include_directories (Ipopt/include/)
link_directories (Ipopt/lib/)
link_directories (/usr/lib/)
link_directories (/usr/local/lib/)
# include_directories (~/stator/)


# # for a c executable
# add_executable (main main.cpp ${src})
# target_link_libraries(main "/usr/local/lib/libnlopt.so")

add_subdirectory(pybind11)
pybind11_add_module(ultrasonics src/pybind/pybinder.cpp ${src})
target_link_libraries(ultrasonics PRIVATE "libnlopt.so")
target_link_libraries(ultrasonics PRIVATE "libipopt.so")
