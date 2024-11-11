# Minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name
project(MyProgram)

# Define project directory for use in the program
add_definitions(-DPROJECT_DIR="${CMAKE_SOURCE_DIR}")

# Set C/C++ standards
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall -Wextra)

set(SRC_FILES src/main.cpp src/binomial_expectation.cpp)
add_executable(my_program ${SRC_FILES})

####################### TESTING STUFF STARTS HERE ###########################################################

enable_testing()
add_subdirectory(tests)