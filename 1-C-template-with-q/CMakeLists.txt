cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME "template-with-q")
project(${PROJECT_NAME})

# Add source to this project's executable.
add_executable (${PROJECT_NAME} main.cpp)