add_library(nanobench STATIC nanobench.cpp)
target_include_directories(
  nanobench PUBLIC ${PROJECT_SOURCE_DIR}/ext ${PROJECT_SOURCE_DIR}/include
                   ${PROJECT_SOURCE_DIR}/ext/eigen)

add_executable(bench_sample ./benchmark_sample.cpp)
target_link_libraries(bench_sample PRIVATE nanobench)

add_executable(bench_tsv ./benchmark_tsv_approaches.cpp)
target_link_libraries(bench_tsv PRIVATE gelex::core nanobench)

add_executable(bench_zvz ./benchmark_zvz.cpp)
target_link_libraries(bench_zvz PRIVATE gelex::core nanobench)

add_executable(bench_gwas ./benchmark_gwas_writer.cpp)
target_link_libraries(bench_gwas PRIVATE gelex::core nanobench)

add_executable(bench_bim_access ./benchmark_bim_access.cpp)
target_link_libraries(bench_bim_access PRIVATE gelex::core nanobench)

add_executable(bench_reml_proj ./benchmark_reml_proj.cpp)
target_link_libraries(bench_reml_proj PRIVATE gelex::core nanobench)

if(GELEX_NATIVE_OPTIMIZATION)
  foreach(target bench_sample bench_tsv bench_zvz bench_gwas bench_bim_access
                 bench_reml_proj)
    target_compile_options(${target} PRIVATE -march=native -O3)
  endforeach()
endif()
