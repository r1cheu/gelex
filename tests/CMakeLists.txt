find_package(Catch2 3 REQUIRED)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "*.cpp")

add_executable(gelex_tests ${TEST_SOURCES})

target_compile_definitions(
  gelex_tests PRIVATE GELEX_TESTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(
  gelex_tests PRIVATE "${PROJECT_SOURCE_DIR}/src")

target_link_libraries(
  gelex_tests PRIVATE -Wl,--whole-archive gelex::core -Wl,--no-whole-archive
                      Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(gelex_tests)

# --------------------------------------------------
# ------------ Coverage Report ---------------------
# --------------------------------------------------
if(GELEX_ENABLE_COVERAGE)
  find_program(LLVM_PROFDATA llvm-profdata REQUIRED)
  find_program(LLVM_COV llvm-cov REQUIRED)

  set(COVERAGE_DIR "${PROJECT_BINARY_DIR}/coverage_report")
  set(PROF_RAW "${PROJECT_BINARY_DIR}/default.profraw")
  set(PROF_DATA "${PROJECT_BINARY_DIR}/coverage.profdata")

  add_custom_target(
    coverage
    COMMAND LLVM_PROFILE_FILE=${PROF_RAW} $<TARGET_FILE:gelex_tests>
    COMMAND ${LLVM_PROFDATA} merge -sparse ${PROF_RAW} -o ${PROF_DATA}
    COMMAND
      ${LLVM_COV} show $<TARGET_FILE:gelex_tests> -instr-profile=${PROF_DATA}
      -format=html -output-dir=${COVERAGE_DIR}
      -ignore-filename-regex="tests/*|ext/*" "${PROJECT_SOURCE_DIR}/src"
      "${PROJECT_SOURCE_DIR}/include"
    COMMAND
      ${LLVM_COV} report $<TARGET_FILE:gelex_tests> -instr-profile=${PROF_DATA}
      -ignore-filename-regex="tests/*|ext/*" "${PROJECT_SOURCE_DIR}/src"
      "${PROJECT_SOURCE_DIR}/include"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS gelex_tests
    COMMENT "Generating HTML coverage report in ${COVERAGE_DIR} and terminal summary")
endif()
