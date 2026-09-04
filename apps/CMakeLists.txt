cmake_minimum_required(VERSION 3.20)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(gelex_cli LANGUAGES CXX)
  set(CMAKE_CXX_STANDARD 23)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

if(NOT TARGET gelex::gelex)
  find_package(gelex CONFIG REQUIRED)
endif()
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(CLI11 CONFIG REQUIRED)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/version.h")

file(GLOB_RECURSE CLI_SOURCES CONFIGURE_DEPENDS cli/*.cpp)
list(APPEND CLI_SOURCES main.cpp)

add_executable(gelex_cli ${CLI_SOURCES})
set_target_properties(gelex_cli PROPERTIES OUTPUT_NAME "gelex")

target_include_directories(gelex_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(gelex_cli PRIVATE gelex::gelex fmt::fmt CLI11::CLI11
                                        spdlog::spdlog)

install(TARGETS gelex_cli RUNTIME DESTINATION bin COMPONENT cli)

option(GELEX_CLI_BUILD_TESTS "Build CLI tests" OFF)

if(GELEX_CLI_BUILD_TESTS OR GELEX_BUILD_TESTS)
  enable_testing()
  find_package(Catch2 3 REQUIRED)

  add_executable(
    gelex_cli_tests
    tests/test_color_formatter.cpp
    tests/test_predict_computation.cpp
    tests/test_predict_io.cpp
    tests/test_predict_engine.cpp
    tests/test_table.cpp
    tests/test_completion_generator.cpp
    tests/test_mcmc_command.cpp
    tests/test_mcmc_data.cpp
    tests/test_mcmc_recipe.cpp
    tests/test_random_design_data.cpp
    tests/test_summary.cpp
    cli/color_formatter.cpp
    cli/theme.cpp
    cli/completion/completion_generator.cpp
    cli/completion/choice_descriptions.cpp
    cli/formatter.cpp
    cli/logging.cpp
    cli/mcmc/command.cpp
    cli/mcmc/config.cpp
    cli/mcmc/data.cpp
    cli/mcmc/progress.cpp
    cli/progress.cpp
    cli/random_design_data.cpp
    cli/reml_data.cpp
    cli/report_printer.cpp
    cli/summary.cpp
    cli/table.cpp
    cli/predict/compute.cpp
    cli/predict/io.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/bed_fixture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/file_fixture.cpp)
  target_include_directories(
    gelex_cli_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}/../tests)
  target_link_libraries(
    gelex_cli_tests PRIVATE gelex::gelex fmt::fmt CLI11::CLI11 spdlog::spdlog
                            Catch2::Catch2WithMain)

  include(Catch)
  catch_discover_tests(gelex_cli_tests)
endif()
