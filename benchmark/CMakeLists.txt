find_package(Catch2 3 REQUIRED)

add_library(nanobench STATIC nanobench.cpp)
target_include_directories(
  nanobench PUBLIC ${PROJECT_SOURCE_DIR}/ext ${PROJECT_SOURCE_DIR}/include
                   ${PROJECT_SOURCE_DIR}/ext/eigen)

file(GLOB BENCHMARK_SOURCES CONFIGURE_DEPENDS "benchmark_*.cpp")

add_executable(gelex_benchmarks ${BENCHMARK_SOURCES})
target_link_libraries(gelex_benchmarks PRIVATE gelex::core nanobench
                                               Catch2::Catch2WithMain)

if(GELEX_NATIVE_OPTIMIZATION)
  target_compile_options(gelex_benchmarks PRIVATE -march=native -O3)
endif()
