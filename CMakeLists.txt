# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.17)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  chenx
  VERSION 0.0.1
  LANGUAGES CXX)

option(BUILD_TEST "Enable testing" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT TARGET phenx AND BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

if(BUILD_TEST)
  enable_testing()
  add_subdirectory(extern/googletest)
  add_subdirectory(tests)
endif()
