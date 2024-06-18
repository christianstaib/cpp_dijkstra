# cmake_minimum_required(VERSION 3.15)
# project(ParallelDijkstra)
#
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#
# # # Conan integration
# # include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
# # conan_basic_setup()
#
# add_executable(ParallelDijkstra main.cpp)
# target_link_libraries(ParallelDijkstra ${CONAN_LIBS})

cmake_minimum_required(VERSION 3.15)
project(compressor C)

find_package(nlohmann_json REQUIRED)

add_executable(${PROJECT_NAME} src/main.c)
target_link_libraries(${PROJECT_NAME} nlohmann_json::nlohmann_json)

