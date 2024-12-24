# tests/CMakeLists.txt Dataset Test


if(USE_MKL)
  find_package(MKL CONFIG REQUIRED)
else()
  set(BLA_VENDOR "OpenBLAS")
  find_package(BLAS REQUIRED)
endif()

find_package(OpenMP)
add_compile_definitions(ARMA_DONT_USE_WRAPPER)

add_executable(all_tests test_dataset.cpp) # test_optim.cpp)
target_include_directories(all_tests PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/extern/armadillo/include)
target_link_libraries(all_tests ${BLAS_LIBRARIES} OpenMP::OpenMP_CXX GTest::gtest_main)

# Optim Test Discover and enable tests
include(GoogleTest)
gtest_discover_tests(all_tests)
