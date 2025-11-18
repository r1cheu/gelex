# Benchmark configuration for gelex

# Add Google Benchmark as a subdirectory
add_subdirectory(benchmark)

# Create benchmark executables
add_executable(bench_sample benchmark_sample.cpp)
add_executable(bench_parser parser_benchmark.cpp)

# Link against our core library and benchmark
target_compile_options(bench_sample PUBLIC -march=native -O3)

target_link_libraries(bench_sample PRIVATE benchmark::benchmark)
