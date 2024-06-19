cmake_minimum_required(VERSION 3.15)
project(compressor C)

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(${PROJECT_NAME} src/main.c)

target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)


