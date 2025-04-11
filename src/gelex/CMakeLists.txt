cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(gelex)

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
message(STATUS "gelex Build Type: ${CMAKE_BUILD_TYPE}")

# Set the version
set(gelex_Version_Major 0)
set(gelex_Version_Minor 2)
set(gelex_Version_Patch 0)
set(gelex_Version_Status "-dev")
set(PROJECT_VERSION
    "${gelex_Version_Major}.${gelex_Version_Minor}.${gelex_Version_Patch}${gelex_Version_Status}"
)
message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(USE_MKL "Use the MKL" OFF)
option(BUILD_TEST "Enable testing" OFF)

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

file(GLOB_RECURSE gelex_SOURCES src/*.cpp)
add_library(gelex SHARED ${gelex_SOURCES})
target_include_directories(gelex PUBLIC include)
target_include_directories(gelex PUBLIC ext/armadillo/include)
target_link_libraries(
  gelex PUBLIC armadillo_config OpenMP::OpenMP_CXX ${BLAS_LIBRARIES}
               spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)

if(BUILD_TEST)
  add_subdirectory(tests)
endif()

install(TARGETS gelex DESTINATION lib)
