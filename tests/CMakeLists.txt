cmake_minimum_required(VERSION 3.18)

add_compile_definitions(GELEX_TESTS_DIR="${PROJECT_SOURCE_DIR}/tests")

find_package(Catch2 3 REQUIRED)

file(GLOB test_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_executable(test ${test_SOURCE})

target_link_libraries(
  test PRIVATE -Wl,--whole-archive core -Wl,--no-whole-archive
               Catch2::Catch2WithMain $<$<CONFIG:Debug>:-fsanitize=undefined>)
if(ENABLE_COVERAGE)
  target_compile_options(test PRIVATE ${COVERAGE_FLAGS})
  target_link_options(test PRIVATE ${COVERAGE_FLAGS})
endif()

include(CTest)
include(Catch)
catch_discover_tests(test)
