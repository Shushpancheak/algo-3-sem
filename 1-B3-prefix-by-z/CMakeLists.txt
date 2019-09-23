cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project("prefix-by-z")

# Add source to this project's executable.
add_executable ("prefix-by-z" main.cpp ConvertZFuncToPrefixFunc.hpp ConvertZFuncToPrefixFunc.cpp)

################################
# GTest
################################
enable_testing()
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
include_directories(${gtest_SOURCE_DIR})

################################
# Unit Tests
################################
add_executable(runUnitTests tests.cpp ConvertZFuncToPrefixFunc.cpp)
# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest_main)
add_test(runUnitTests runUnitTests)