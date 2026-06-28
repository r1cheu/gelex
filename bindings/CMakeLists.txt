cmake_minimum_required(VERSION 3.20)
project(gelex_python LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(gelex CONFIG REQUIRED)
find_package(
  Python 3.12
  COMPONENTS Interpreter Development.Module
  REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

file(GLOB_RECURSE BINDING_SOURCES CONFIGURE_DEPENDS src/*.cpp)

nanobind_add_module(_gelex STABLE_ABI NB_STATIC ${BINDING_SOURCES})
target_include_directories(_gelex PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(_gelex PRIVATE gelex::gelex)

set(GELEX_PY_PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/python/gelex)
set_target_properties(_gelex PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                                        ${GELEX_PY_PACKAGE_DIR})

option(GELEX_PY_GENERATE_STUBS "Generate .pyi type stubs for the extension" ON)

if(GELEX_PY_GENERATE_STUBS)
  nanobind_add_stub(
    _gelex_stub
    MODULE
    _gelex
    OUTPUT
    ${GELEX_PY_PACKAGE_DIR}/_gelex.pyi
    PYTHON_PATH
    $<TARGET_FILE_DIR:_gelex>
    DEPENDS
    _gelex)
endif()

# Install the package into the environment's site-packages so that a conda
# package built from this project exposes `import gelex` natively -- no
# PYTHONPATH, no pip. The destination is resolved relative to the install prefix
# to stay correct inside the pixi-build sandbox.
execute_process(
  COMMAND
    "${Python_EXECUTABLE}" -c
    "import os,sysconfig;print(os.path.relpath(sysconfig.get_path('platlib'),sysconfig.get_path('data')))"
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE GELEX_PY_SITELIB)

install(TARGETS _gelex LIBRARY DESTINATION ${GELEX_PY_SITELIB}/gelex)
install(FILES ${GELEX_PY_PACKAGE_DIR}/__init__.py
        DESTINATION ${GELEX_PY_SITELIB}/gelex)

if(GELEX_PY_GENERATE_STUBS)
  install(FILES ${GELEX_PY_PACKAGE_DIR}/_gelex.pyi
          DESTINATION ${GELEX_PY_SITELIB}/gelex)
endif()
