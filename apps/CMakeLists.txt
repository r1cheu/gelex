set(CLI_SOURCES
    main.cpp
    cli/assoc/assoc_args.cpp
    cli/assoc/assoc_config.cpp
    cli/assoc/assoc_command.cpp
    cli/assoc/gwas_runner.cpp
    cli/fit/fit_args.cpp
    cli/fit/fit_config.cpp
    cli/fit/fit_command.cpp
    cli/fit/fit_reporter.cpp
    cli/grm/grm_args.cpp
    cli/grm/grm_config.cpp
    cli/grm/grm_command.cpp
    cli/grm/grm_reporter.cpp
    cli/predict/predict_args.cpp
    cli/predict/predict_config.cpp
    cli/predict/predict_command.cpp
    cli/simulate/simulate_args.cpp
    cli/simulate/simulate_config.cpp
    cli/simulate/simulate_command.cpp
    cli/simulate/simulator_reporter.cpp
    cli/post/post_args.cpp
    cli/post/post_config.cpp
    cli/post/post_command.cpp
    cli/post/post_reporter.cpp
    cli/data_pipe_reporter.cpp
    cli/data_pipe_config.cpp
    cli/cli_helper.cpp)

add_executable(gelex_cli ${CLI_SOURCES})
set_target_properties(gelex_cli PROPERTIES OUTPUT_NAME "gelex")

target_include_directories(gelex_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(gelex_cli PRIVATE gelex::core)

install(TARGETS gelex_cli RUNTIME DESTINATION bin)
