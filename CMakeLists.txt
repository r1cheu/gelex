cmake_minimum_required(VERSION 3.20)

# --------------------------------------------------
# ------------ Project Info ------------------------
# --------------------------------------------------
project(
  gelex
  VERSION 0.7.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --------------------------------------------------
# ------------ Options -----------------------------
# --------------------------------------------------
option(GELEX_USE_MKL "Use Intel MKL for Linear Algebra" OFF)
option(GELEX_BUILD_TESTS "Build unit tests" OFF)
option(GELEX_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(GELEX_ENABLE_COVERAGE "Enable Clang source-based code coverage" OFF)
option(GELEX_NATIVE_OPTIMIZATION "Use -march=native" OFF)
option(GELEX_ENABLE_SANITIZERS "Enable ASAN/UBSAN" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# --------------------------------------------------
# ------------ Dependencies ------------------------
# --------------------------------------------------
# list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)
find_package(spdlog REQUIRED)

# BLAS/LAPACK Selection
if(GELEX_USE_MKL)
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
else()
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
endif()

# --------------------------------------------------
# ------------ Build Targets -----------------------
# --------------------------------------------------
configure_file("${PROJECT_SOURCE_DIR}/config.h.in"
               "${PROJECT_BINARY_DIR}/config.h")

add_subdirectory(src)

if(GELEX_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(GELEX_BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif()
