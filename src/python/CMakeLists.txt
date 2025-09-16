cmake_minimum_required(VERSION 3.18)

find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core NB_STATIC ${PYTHON_BINDING})
target_link_libraries(_core PRIVATE core)
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
install(TARGETS _core LIBRARY DESTINATION gelexy)
install(FILES ${CMAKE_BINARY_DIR}/_core.pyi ${CMAKE_BINARY_DIR}/py.typed
        DESTINATION gelexy)
