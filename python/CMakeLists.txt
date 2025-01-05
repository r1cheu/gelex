# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.18)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(CMAKE_BUILD_TYPE "Debug")

find_package(
  Python 3.8
  COMPONENTS Interpreter Development.Module
  REQUIRED)
if(NOT TARGET chenx)
  find_package(chenx CONFIG)
  if(NOT chenx_FOUND)
    message(STATUS "Using bundled chenx sources")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. _deps/chenx)
  endif()
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# the definition won't be passed from chenx due it's another project re-add it
# set(ARMA_DEFINITIONS ARMA_DONT_USE_WRAPPER ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
# ARMA_OPENMP_THREADS=64)
#
# if(CMAKE_BUILD_TYPE STREQUAL "Debug") list(APPEND ARMA_DEFINITIONS ARMA_DEBUG)
# elseif(CMAKE_BUILD_TYPE STREQUAL "Release") list(APPEND ARMA_DEFINITIONS
# ARMA_NO_DEBUG) endif() add_compile_definitions(${ARMA_DEFINITIONS})
#
nanobind_add_module(_chenx src/bind.cpp)
target_link_libraries(_chenx PRIVATE chenx)
install(TARGETS _chenx LIBRARY DESTINATION phenx)

set_target_properties(_chenx PROPERTIES BUILD_WITH_INSTALL_RPATH ON
                                        INSTALL_RPATH "$ORIGIN/../lib")
