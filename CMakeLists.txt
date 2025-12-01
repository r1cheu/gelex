cmake_minimum_required(VERSION 3.18)

# --------------------------------------------------
# ------------ Basic Project Setup -----------------
# --------------------------------------------------

project(
  "gelex"
  VERSION 0.6.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CONFIGURE_HEADER "${PROJECT_BINARY_DIR}/config.h")
configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${CONFIGURE_HEADER}")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      Debug
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# --------------------------------------------------
# ------------ Options and Dependencies ------------
# --------------------------------------------------

option(USE_MKL "Use the MKL" OFF)
option(BUILD_TEST "Enable testing" OFF)
option(ENABLE_COVERAGE "Enable Clang source-based code coverage" OFF)

if(ENABLE_COVERAGE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Enabled Source-based Code Coverage")
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING "Force Debug for coverage" FORCE)
    set(COVERAGE_FLAGS -fprofile-instr-generate -fcoverage-mapping)
  else()
    message(FATAL_ERROR "Code coverage requires Clang compiler.")
  endif()
endif()

find_package(OpenMP)
find_package(spdlog REQUIRED)
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
  PUBLIC SYSTEM ${PROJECT_SOURCE_DIR}/ext
  PUBLIC SYSTEM ${PROJECT_SOURCE_DIR}/ext/armadillo/include
  PUBLIC SYSTEM ${PROJECT_SOURCE_DIR}/ext/eigen
  PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_compile_definitions(
  core
  PUBLIC BARKEEP_ENABLE_FMT_FORMAT ARMA_DONT_USE_WRAPPER
         ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS ARMA_OPENMP_THREADS=64 ARMA_NO_DEBUG
         $<$<CONFIG:Release>:EIGEN_NO_DEBUG>)

target_compile_options(
  core PUBLIC $<$<CONFIG:Release>:-march=native -O3> -Wall -Wextra
              $<$<CONFIG:Debug>:-fsanitize=undefined>)

if(ENABLE_COVERAGE)
  target_compile_options(core PRIVATE ${COVERAGE_FLAGS})
  target_link_options(core PRIVATE ${COVERAGE_FLAGS})
endif()

if(USE_MKL)
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS REQUIRED)
  target_compile_definitions(core INTERFACE EIGEN_USE_MKL_ALL USE_MKL)
else()
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  target_compile_definitions(core INTERFACE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
endif()

target_link_libraries(core PUBLIC OpenMP::OpenMP_CXX BLAS::BLAS spdlog::spdlog
                                  $<$<BOOL:${MINGW}>:ws2_32>)

add_executable(gelex src/main.cpp)
target_include_directories(gelex PUBLIC "${PROJECT_BINARY_DIR}")

target_link_libraries(gelex PRIVATE core
                                    $<$<CONFIG:Debug>:-fsanitize=undefined>)

if(BUILD_TEST)
  add_subdirectory(tests/)
  if(ENABLE_COVERAGE)
    find_program(LLVM_PROFDATA llvm-profdata REQUIRED)
    find_program(LLVM_COV llvm-cov REQUIRED)

    set(COVERAGE_DIR "${CMAKE_BINARY_DIR}/coverage_report")
    set(PROF_RAW "${CMAKE_BINARY_DIR}/default.profraw")
    set(PROF_DATA "${CMAKE_BINARY_DIR}/coverage.profdata")

    add_custom_target(
      coverage
      COMMAND LLVM_PROFILE_FILE=${PROF_RAW} $<TARGET_FILE:test>
      COMMAND ${LLVM_PROFDATA} merge -sparse ${PROF_RAW} -o ${PROF_DATA}
      COMMAND
        ${LLVM_COV} report $<TARGET_FILE:test> -instr-profile=${PROF_DATA}
        -show-region-summary=false "${PROJECT_SOURCE_DIR}/src"
        "${PROJECT_SOURCE_DIR}/include/"
      COMMAND
        ${LLVM_COV} show $<TARGET_FILE:test> -instr-profile=${PROF_DATA}
        -format=html -output-dir=${COVERAGE_DIR} "${PROJECT_SOURCE_DIR}/src"
        "${PROJECT_SOURCE_DIR}/include/"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      DEPENDS test
      COMMENT "Running tests and generating coverage report... ")
  endif()
endif()

# Add benchmarks with test dependencies disabled
add_subdirectory(benchmark/)
