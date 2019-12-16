cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project("4-A-big-integer")
add_executable (${PROJECT_NAME} main.cpp biginteger.h)