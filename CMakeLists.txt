cmake_minimum_required(VERSION 3.18)
project(${SKBUILD_PROJECT_NAME})
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "gelexy bindings Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")

option(USE_MKL "Use the MKL" OFF)
option(BUILD_TEST "Enable testing" OFF)

if(USE_MKL)
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS REQUIRED)
else()
  set(BLA_VENDOR "OpenBLAS")
  find_package(BLAS REQUIRED)
endif()

# ---- Dependencies ----
find_package(OpenMP)
find_package(spdlog REQUIRED)
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

add_library(common INTERFACE)
target_compile_definitions(
  common
  INTERFACE ARMA_DONT_USE_WRAPPER ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
            ARMA_OPENMP_THREADS=64 $<$<CONFIG:Debug>:ARMA_DEBUG>
            $<$<CONFIG:Release>:ARMA_NO_DEBUG>)

target_include_directories(common INTERFACE include)
target_include_directories(common INTERFACE ext/armadillo/include)
target_link_libraries(
  common INTERFACE OpenMP::OpenMP_CXX ${BLAS_LIBRARIES} spdlog::spdlog
                   $<$<BOOL:${MINGW}>:ws2_32>)
file(GLOB_RECURSE CORE_SOURCE src/*.cpp)
nanobind_add_module(_core NB_STATIC ${CORE_SOURCE})
target_link_libraries(_core PRIVATE common)

if(BUILD_TEST)
  nanobind_add_module(_test NB_STATIC tests/python/test_dense_caster.cpp)
  target_link_libraries(_test PRIVATE common)
  install(TARGETS _test LIBRARY DESTINATION gelexy)

  set(GELEX_TESTS_DIR "${PROJECT_SOURCE_DIR}/tests")
  add_compile_definitions(GELEX_TESTS_DIR="${GELEX_TESTS_DIR}")
  find_package(Catch2 3 REQUIRED)
  file(GLOB test_SOURCE ${PROJECT_SOURCE_DIR}/tests/cpp/*.cpp)
  list(FILTER CORE_SOURCE EXCLUDE REGEX "src/python/")
  add_executable(test ${test_SOURCE} ${CORE_SOURCE})
  target_link_libraries(test PRIVATE common Catch2::Catch2WithMain)
endif()

install(TARGETS _core LIBRARY DESTINATION gelexy)
