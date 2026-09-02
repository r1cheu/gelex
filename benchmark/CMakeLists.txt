if(NOT TARGET gelex::gelex)
  find_package(gelex CONFIG REQUIRED)
endif()

find_package(Catch2 3 REQUIRED)

add_library(nanobench STATIC nanobench.cpp)
target_include_directories(nanobench SYSTEM
                           PUBLIC ${PROJECT_SOURCE_DIR}/ext/nanobench)
target_link_libraries(nanobench PUBLIC Eigen3::Eigen)

file(GLOB BENCHMARK_SOURCES CONFIGURE_DEPENDS "benchmark_*.cpp")

add_executable(gelex_benchmarks ${BENCHMARK_SOURCES})
target_include_directories(gelex_benchmarks PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(gelex_benchmarks PRIVATE gelex::gelex nanobench
                                               Catch2::Catch2WithMain)

if(GELEX_NATIVE_OPTIMIZATION)
  target_compile_options(gelex_benchmarks PRIVATE -march=native -O3)
endif()
