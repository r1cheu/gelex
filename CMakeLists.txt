cmake_minimum_required(VERSION 3.18)
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
  find_path(BLAS_INCLUDE_DIR mkl.h)
else()
  set(BLA_VENDOR OpenBLAS)
  find_package(BLAS REQUIRED)
  find_path(BLAS_INCLUDE_DIR cblas.h)
endif()
find_package(OpenMP)
find_package(spdlog REQUIRED)
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Create common config target
add_library(common_config INTERFACE)

# Set common properties ONCE
target_compile_definitions(
  common_config
  INTERFACE 
    BARKEEP_ENABLE_FMT_FORMAT
    ARMA_DONT_USE_WRAPPER
    ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS
    ARMA_OPENMP_THREADS=64
    $<$<CONFIG:Debug>:ARMA_DEBUG>
    $<$<CONFIG:Release>:ARMA_NO_DEBUG>
)

target_include_directories(
  common_config
  INTERFACE 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/ext/armadillo/include
    ${PROJECT_SOURCE_DIR}/ext/barkeep/barkeep
)

target_link_libraries(
  common_config
  INTERFACE 
    OpenMP::OpenMP_CXX 
    BLAS::BLAS 
    spdlog::spdlog
    $<$<BOOL:${MINGW}>:ws2_32>
)

if(USE_MKL)
  target_compile_definitions(common_config INTERFACE USE_MKL)
endif()

# Update common target to use config
add_library(common INTERFACE)
target_link_libraries(common INTERFACE common_config)
# Core sources (excluding python bindings)
file(GLOB_RECURSE CORE_SOURCE src/*.cpp)
list(FILTER CORE_SOURCE EXCLUDE REGEX "src/python/")

# Python bindings
file(GLOB_RECURSE PYTHON_SOURCE src/python/*.cpp)
nanobind_add_module(_core NB_STATIC ${PYTHON_SOURCE})
target_link_libraries(_core PRIVATE common)
nanobind_add_stub(
  _core_stub
  MODULE
  _core
  OUTPUT
  _core.pyi
  PYTHON_PATH
  $<TARGET_FILE_DIR:_core>
  MARKER_FILE
  py.typed
  DEPENDS
  _core)

if(BUILD_TEST)
  # Get test sources without recompiling python bindings
  set(TEST_CORE_SOURCE ${CORE_SOURCE})
  
  # Add tests subdirectory
  add_subdirectory(tests/cpp)
endif()

install(TARGETS _core LIBRARY DESTINATION gelexy)
install(FILES ${CMAKE_BINARY_DIR}/_core.pyi ${CMAKE_BINARY_DIR}/py.typed
        DESTINATION gelexy)
