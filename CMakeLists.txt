cmake_minimum_required(VERSION 3.20)

# --------------------------------------------------
# ------------ Project Info ------------------------
# --------------------------------------------------
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" GELEX_VERSION)
project(
  gelex
  VERSION ${GELEX_VERSION}
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_CXX_COMPILER_LAUNCHER)
  find_program(SCCACHE_PROGRAM sccache)
  find_program(CCACHE_PROGRAM ccache)
  if(SCCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${SCCACHE_PROGRAM}")
    message(STATUS "Using sccache: ${SCCACHE_PROGRAM}")
  elseif(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
  endif()
endif()

# --------------------------------------------------
# ------------ Options -----------------------------
# --------------------------------------------------
option(GELEX_BUILD_TESTS "Build unit tests" OFF)
option(GELEX_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(GELEX_BUILD_CLI "Build the command-line application" OFF)
option(GELEX_BUILD_PYTHON "Build the Python bindings" OFF)
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
find_package(OpenMP REQUIRED)
find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)

# --------------------------------------------------
# ------------ Build Targets -----------------------
# --------------------------------------------------
add_subdirectory(src)

if(GELEX_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(GELEX_BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif()

if(GELEX_BUILD_CLI)
  add_subdirectory(apps)
endif()

if(GELEX_BUILD_PYTHON)
  add_subdirectory(bindings)
endif()

# --------------------------------------------------
# ------------ Install / Export --------------------
# --------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(
  TARGETS gelex
  EXPORT gelexTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/gelex
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  EXPORT gelexTargets
  NAMESPACE gelex::
  FILE gelexTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gelex)

configure_package_config_file(
  cmake/gelexConfig.cmake.in "${PROJECT_BINARY_DIR}/gelexConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gelex)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/gelexConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(
  FILES "${PROJECT_BINARY_DIR}/gelexConfig.cmake"
        "${PROJECT_BINARY_DIR}/gelexConfigVersion.cmake"
        "${PROJECT_SOURCE_DIR}/cmake/FindLAPACKE.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gelex)
