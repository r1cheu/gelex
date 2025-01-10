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

if(NOT CMAKE_BUILD_TYPE)
  if(DEFINED ENV{CMAKE_BUILD_TYPE})
    message(
      STATUS
        "Using environment variable CMAKE_BUILD_TYPE: $ENV{CMAKE_BUILD_TYPE}")
    set(CMAKE_BUILD_TYPE
        $ENV{CMAKE_BUILD_TYPE}
        CACHE STRING "Build type" FORCE)
  else()
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING "Build type" FORCE)
  endif()
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(USE_MKL "Use the MKL" OFF)
option(BUILD_TEST "Enable testing" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)

if(USE_MKL)
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS REQUIRED)
else()
  set(BLA_VENDOR "OpenBLAS")
  find_package(BLAS REQUIRED)
endif()

find_package(OpenMP)
find_package(spdlog REQUIRED)
add_library(armadillo_config INTERFACE)

# Set public definitions
target_compile_definitions(
  armadillo_config
  INTERFACE ARMA_DONT_USE_WRAPPER ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
            ARMA_OPENMP_THREADS=64 $<$<CONFIG:Debug>:ARMA_DEBUG>
            $<$<CONFIG:Release>:ARMA_NO_DEBUG>)

file(GLOB_RECURSE chenx_SOURCES src/*.cpp)
add_library(chenx SHARED ${chenx_SOURCES})
target_include_directories(chenx PUBLIC include)
target_include_directories(chenx PUBLIC extern/armadillo/include)
target_link_libraries(
  chenx PUBLIC armadillo_config OpenMP::OpenMP_CXX ${BLAS_LIBRARIES}
               spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)

if(NOT TARGET phenx AND BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

if(BUILD_TEST)
  add_subdirectory(tests)
endif()

install(TARGETS chenx DESTINATION lib)
