cmake_minimum_required(VERSION 3.15)
project(paths CXX)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(MPI REQUIRED)
find_package(zpp_bits REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
add_executable(${PROJECT_NAME} src/main.cpp src/graph.cpp)

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX nlohmann_json::nlohmann_json zpp_bits::zpp_bits MPI::MPI_CXX)

