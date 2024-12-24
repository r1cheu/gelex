# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.18)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(Python3 REQUIRED COMPONENTS NumPy)
find_package(pybind11 CONFIG REQUIRED)
find_package(OpenMP)

option(USE_MKL "Use the MKL" OFF)

if(USE_MKL)
  find_package(MKL CONFIG REQUIRED)
else()
  set(BLA_VENDOR "OpenBLAS")
  find_package(BLAS REQUIRED)
endif()
add_compile_definitions(ARMA_DONT_USE_WRAPPER)
add_compile_definitions(FMT_HEADER_ONLY=1)

# Function to add a pybind11 module with common settings
function(add_pybind11_module module_name source_file)
  pybind11_add_module(${module_name} MODULE ${source_file})
  target_include_directories(
    ${module_name} PRIVATE ../extern/carma/include ../extern/fmt/include
                           ../include ../extern/armadillo/include)
  target_link_libraries(
    ${module_name} PRIVATE pybind11::module Python3::NumPy OpenMP::OpenMP_CXX
                           ${BLAS_LIBRARIES})
endfunction()

# Add _core module
add_pybind11_module(_core src/bind_core.cpp)
# Add _optim module set_target_properties(_core PROPERTIES INSTALL_RPATH
# "$ORIGIN/lib") The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION .)
