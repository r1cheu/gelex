# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(chenx)

if(NOT CMAKE_BUILD_TYPE)
  set(default_build_type "Debug")
  message(
    STATUS
      "Set the build type to `${default_build_type}` as none was specified.")
  set(CMAKE_BUILD_TYPE
      ${default_build_type}
      CACHE STRING "Chooce the build type." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "chenx Build Type: ${CMAKE_BUILD_TYPE}")

# Set the version
set(chenx_Version_Major 0)
set(chenx_Version_Minor 1)
set(chenx_Version_Patch 0)
set(chenx_Version_Status "-dev")
set(PROJECT_VERSION
    "${chenx_Version_Major}.${chenx_Version_Minor}.${chenx_Version_Patch}${chenx_Version_Status}"
)
message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
