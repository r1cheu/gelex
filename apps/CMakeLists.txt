set(CLI_SOURCES
    main.cpp
    cli/assoc/args.cpp
    cli/assoc/config.cpp
    cli/assoc/command.cpp
    cli/assoc/reporter.cpp
    cli/mcmc/args.cpp
    cli/mcmc/config.cpp
    cli/mcmc/command.cpp
    cli/mcmc/reporter.cpp
    cli/grm/args.cpp
    cli/grm/config.cpp
    cli/grm/command.cpp
    cli/grm/reporter.cpp
    cli/predict/args.cpp
    cli/predict/config.cpp
    cli/predict/command.cpp
    cli/predict/reporter.cpp
    cli/simulate/args.cpp
    cli/simulate/config.cpp
    cli/simulate/command.cpp
    cli/simulate/reporter.cpp
    cli/post/args.cpp
    cli/post/config.cpp
    cli/post/command.cpp
    cli/post/reporter.cpp
    cli/report_printer.cpp
    cli/fit_reporter.cpp
    cli/reml_reporter.cpp
    cli/dataset_reporter.cpp
    cli/pheno_reporter.cpp
    cli/geno_reporter.cpp
    cli/grm_pipe_reporter.cpp
    cli/data_pipe_config.cpp
    cli/cli_helper.cpp
    cli/bayes_recipe_options.cpp)

add_executable(gelex_cli ${CLI_SOURCES})
set_target_properties(gelex_cli PROPERTIES OUTPUT_NAME "gelex")

target_include_directories(gelex_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(gelex_cli PRIVATE gelex::core)

install(TARGETS gelex_cli RUNTIME DESTINATION bin)
