# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.18)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

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

if(NOT CMAKE_BUILD_TYPE)
  if(DEFINED ENV{CMAKE_BUILD_TYPE})
    set(CMAKE_BUILD_TYPE $ENV{CMAKE_BUILD_TYPE})
  else()
    set(CMAKE_BUILD_TYPE Debug)
  endif()
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

execute_process( COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_chenx src/bind.cpp)
target_link_libraries(_chenx PRIVATE chenx)
install(TARGETS _chenx LIBRARY DESTINATION phenx)

nanobind_add_stub(
  _chenx_stub
  INSTALL_TIME
  MODULE
  _chenx
  OUTPUT
  phenx/_chenx.pyi
  PYTHON_PATH
  "phenx")

set_target_properties(_chenx PROPERTIES BUILD_WITH_INSTALL_RPATH ON
                                        INSTALL_RPATH "$ORIGIN/../lib")
