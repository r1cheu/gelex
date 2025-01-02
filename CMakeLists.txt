# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.18)
# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  chenx
  VERSION 0.0.1
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "DEBUG")

option(USE_MKL "Use the MKL" OFF)
option(BUILD_TEST "Enable testing" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)

if(USE_MKL)
  find_package(MKL CONFIG REQUIRED)
else()
  set(BLA_VENDOR "OpenBLAS")
  find_package(BLAS REQUIRED)
endif()
find_package(OpenMP)
find_package(spdlog REQUIRED)

set(ARMA_DEFINITIONS ARMA_DONT_USE_WRAPPER ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
                     ARMA_OPENMP_THREADS=64)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND ARMA_DEFINITIONS ARMA_DEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  list(APPEND ARMA_DEFINITIONS ARMA_NO_DEBUG)
endif()
add_compile_definitions(${ARMA_DEFINITIONS})

file(GLOB_RECURSE chenx_SOURCES src/*.cpp)
add_library(chenx SHARED ${chenx_SOURCES})
target_include_directories(chenx PUBLIC include)
target_include_directories(chenx PUBLIC extern/armadillo/include)
target_link_libraries(chenx PUBLIC OpenMP::OpenMP_CXX ${BLAS_LIBRARIES}
                                   spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)

if(NOT TARGET phenx AND BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

if(BUILD_TEST)
  add_subdirectory(tests)
endif()
