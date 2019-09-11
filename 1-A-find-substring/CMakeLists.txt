# CMakeList.txt : CMake project for 1-A-find-substring, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project("find-substring")

# Add source to this project's executable.
add_executable ("find-substring" main.cpp)

# TODO: Add tests and install targets if needed.
