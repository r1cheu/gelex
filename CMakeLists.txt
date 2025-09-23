cmake_minimum_required(VERSION 3.18)
project(
  "Gelex"
  VERSION 0.4.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")

option(USE_MKL "Use the MKL" OFF)
option(BUILD_TEST "Enable testing" OFF)

find_package(OpenMP)
find_package(spdlog REQUIRED)
include(FetchContent)
FetchContent_Declare(
  argparse
  GIT_REPOSITORY https://github.com/p-ranav/argparse.git
  GIT_TAG d924b84eba1f0f0adf38b20b7b4829f6f65b6570)
FetchContent_MakeAvailable(argparse)

# FetchContent_Declare( fmt GIT_REPOSITORY https://github.com/fmtlib/fmt GIT_TAG
# 40626af88bd7df9a5fb80be7b25ac85b122d6c21) # 11.2.0
# FetchContent_MakeAvailable(fmt)
#
# core library
file(GLOB_RECURSE CORE_SOURCE "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM CORE_SOURCE "${PROJECT_SOURCE_DIR}/src/main.cpp")

add_library(core ${CORE_SOURCE})
set_target_properties(core PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(
  core
  PUBLIC ${PROJECT_SOURCE_DIR}/include
         ${PROJECT_SOURCE_DIR}/ext/armadillo/include
         ${PROJECT_SOURCE_DIR}/ext/eigen
  PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_compile_definitions(
  core
  PUBLIC BARKEEP_ENABLE_FMT_FORMAT ARMA_DONT_USE_WRAPPER
         ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS ARMA_OPENMP_THREADS=64
         $<$<CONFIG:Debug>:ARMA_DEBUG> $<$<CONFIG:Release>:ARMA_NO_DEBUG>)

if(USE_MKL)
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS REQUIRED)
  find_path(BLAS_INCLUDE_DIR mkl.h)
  target_compile_definitions(core INTERFACE EIGEN_USE_MKL_ALL USE_MKL)
else()
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  find_path(BLAS_INCLUDE_DIR cblas.h)
  target_compile_definitions(core INTERFACE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
endif()

target_link_libraries(core PUBLIC OpenMP::OpenMP_CXX BLAS::BLAS spdlog::spdlog
                                  $<$<BOOL:${MINGW}>:ws2_32>)

add_executable(gelex src/main.cpp)
target_link_libraries(gelex PRIVATE core argparse)

if(BUILD_TEST)
  add_subdirectory(tests/)
endif()
