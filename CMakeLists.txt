cmake_minimum_required(VERSION 3.15)
project(paths CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(glm REQUIRED)
find_package(OpenMP REQUIRED)
find_package(CLI11 REQUIRED)
find_package(indicators REQUIRED)

include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
add_executable(${PROJECT_NAME} src/main.cpp src/space.cpp src/octree.cpp)

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX glm::glm CLI11::CLI11 indicators::indicators)

