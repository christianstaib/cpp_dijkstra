cmake_minimum_required(VERSION 3.15)
project(compressor CXX)

find_package(OpenMP REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(${PROJECT_NAME} src/main.cpp src/graph.cpp)

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX nlohmann_json::nlohmann_json)


