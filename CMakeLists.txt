cmake_minimum_required(VERSION 3.15)
project(simulation CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(glm REQUIRED)
find_package(OpenMP REQUIRED)
find_package(CLI11 REQUIRED)
find_package(indicators REQUIRED)
find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
add_executable(${PROJECT_NAME} src/main.cpp src/space.cpp src/octree.cpp src/naive_calculations.cpp src/csv_writer.cpp src/cli.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")


target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX glm::glm CLI11::CLI11 indicators::indicators MPI::MPI_CXX)

