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

option(USE_MKL "Use the MKL" OFF)
option(BUILD_TEST "Enable testing" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)

if(USE_MKL)
  find_package(MKL CONFIG REQUIRED)
else()
  set(BLA_VENDOR "OpenBLAS")
  find_package(BLAS REQUIRED)
endif()
find_package(OpenMP)
find_package(spdlog REQUIRED)

file(GLOB_RECURSE chenx_SOURCES src/*.cpp)

add_compile_definitions(ARMA_DONT_USE_WRAPPER)

add_library(chenx SHARED ${chenx_SOURCES})
target_include_directories(chenx PUBLIC include)
target_include_directories(chenx PRIVATE extern/armadillo/include)
target_link_libraries(chenx PRIVATE OpenMP::OpenMP_CXX ${BLAS_LIBRARIES}
                                    spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)

if(NOT TARGET phenx AND BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

if(BUILD_TEST)
  enable_testing()
  add_subdirectory(extern/googletest)
  add_subdirectory(tests)
endif()
