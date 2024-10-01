# tests/CMakeLists.txt Dataset Test

add_executable(all_tests test_dataset.cpp) # test_optim.cpp)
target_include_directories(all_tests PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/extern/armadillo/include)
target_link_libraries(all_tests openblas GTest::gtest_main)

# Optim Test Discover and enable tests
include(GoogleTest)
gtest_discover_tests(all_tests)
