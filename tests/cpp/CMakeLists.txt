cmake_minimum_required(VERSION 3.18)

# Set tests directory definition
add_compile_definitions(GELEX_TESTS_DIR="${PROJECT_SOURCE_DIR}/tests")

# Find test-specific dependency
find_package(Catch2 3 REQUIRED)

# Add test executable
file(GLOB test_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_executable(test ${test_SOURCE} ${TEST_CORE_SOURCE})

# Use common configuration
target_link_libraries(test PRIVATE common_config Catch2::Catch2WithMain)
