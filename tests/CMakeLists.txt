cmake_minimum_required(VERSION 3.18)

add_compile_definitions(GELEX_TESTS_DIR="${PROJECT_SOURCE_DIR}/tests")

find_package(Catch2 3 REQUIRED)

file(GLOB test_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_executable(test ${test_SOURCE})

target_link_libraries(
  test PRIVATE core Catch2::Catch2WithMain
               $<$<CONFIG:Debug>:-fsanitize=undefined,address>)

include(CTest)
include(Catch)
catch_discover_tests(test)
