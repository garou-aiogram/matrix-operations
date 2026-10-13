cmake_minimum_required(VERSION 3.16)
project(MatrixOperations LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(matrix_core src/Matrix.cpp)
target_include_directories(matrix_core PUBLIC src)

add_executable(matrix_demo src/main.cpp)
target_link_libraries(matrix_demo PRIVATE matrix_core)
