# tests/CMakeLists.txt Dataset Test
set(GELEX_TESTS_DIR "${PROJECT_SOURCE_DIR}/tests")
add_compile_definitions(GELEX_TESTS_DIR="${GELEX_TESTS_DIR}")

find_package(Catch2 3 REQUIRED)
file(GLOB test_SOURCE ${PROJECT_SOURCE_DIR}/tests/*.cpp)
add_executable(test ${test_SOURCE})
target_link_libraries(test PRIVATE gelex Catch2::Catch2WithMain)
