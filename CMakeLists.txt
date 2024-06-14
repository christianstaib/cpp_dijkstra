cmake_minimum_required(VERSION 3.15)
project(ParallelDijkstra)

# Conan integration
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

add_executable(ParallelDijkstra main.cpp)
target_link_libraries(ParallelDijkstra ${CONAN_LIBS})

