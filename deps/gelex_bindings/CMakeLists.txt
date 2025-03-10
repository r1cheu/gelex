cmake_minimum_required(VERSION 3.18)

project(
  _gelex
  VERSION 0.2.0
  LANGUAGES CXX)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(
  Python 3.8
  COMPONENTS Interpreter Development.Module
  REQUIRED)

if(NOT TARGET gelex)
  find_package(gelex CONFIG)
  if(NOT gelex_FOUND)
    message(STATUS "Using bundled gelex sources")
    add_subdirectory(${PROJECT_SOURCE_DIR}/../gelex _deps)
  endif()
endif()

execute_process(
  COMMAND ${Python_EXECUTABLE} -c
          "import sysconfig; print(sysconfig.get_path('purelib'))"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
  OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(${PROJECT_NAME} src/bind.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE gelex)
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}/gelexy
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${BINDIR})
