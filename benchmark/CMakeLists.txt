include(FetchContent)
FetchContent_Declare(
  googlebenchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.8.3)
set(BENCHMARK_ENABLE_TESTING
    OFF
    CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL
    OFF
    CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googlebenchmark)

add_executable(bench ./benchmark_gwas_writer.cpp)

target_link_libraries(bench PRIVATE gelex::core benchmark::benchmark)

if(GELEX_NATIVE_OPTIMIZATION)
  target_compile_options(bench PRIVATE -march=native -O3)
endif()
