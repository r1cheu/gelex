# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15)

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

set(BLA_VENDOR "OpenBLAS")
find_package(BLAS REQUIRED)
find_package(fmt REQUIRED)
add_compile_definitions(ARMA_DONT_USE_WRAPPER)
# Add a library using FindPython's tooling (pybind11 also provides a helper like
# this)
pybind11_add_module(_core MODULE src/bind_core.cpp)
target_include_directories(
  _core
  PRIVATE ${PhenxCpp_SOURCE_DIR}/extern/carma/include
          ${PhenxCpp_SOURCE_DIR}/include
          ${PhenxCpp_SOURCE_DIR}/extern/armadillo/include)
target_link_libraries(_core PRIVATE fmt::fmt-header-only pybind11::module
                                    OpenMP::OpenMP_CXX openblas Python3::NumPy)
set_target_properties(_core PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION .)
