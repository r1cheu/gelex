set(CLI_SOURCES
    main.cpp
    cli/assoc/args.cpp
    cli/assoc/command.cpp
    cli/assoc/reporter.cpp
    cli/mcmc/args.cpp
    cli/mcmc/command.cpp
    cli/mcmc/reporter.cpp
    cli/grm/args.cpp
    cli/grm/command.cpp
    cli/grm/reporter.cpp
    cli/predict/args.cpp
    cli/predict/command.cpp
    cli/predict/reporter.cpp
    cli/simulate/args.cpp
    cli/simulate/command.cpp
    cli/simulate/reporter.cpp
    cli/post/args.cpp
    cli/post/command.cpp
    cli/post/reporter.cpp
    cli/reml/args.cpp
    cli/reml/command.cpp
    cli/reml/reporter.cpp
    cli/report_printer.cpp
    cli/formatter.cpp
    cli/timer.cpp
    cli/fit_reporter.cpp
    cli/reml_reporter.cpp
    cli/geno_reporter.cpp
    cli/cli_helper.cpp
    cli/color_formatter.cpp
    cli/bayes_recipe_options.cpp)

add_executable(gelex_cli ${CLI_SOURCES})
set_target_properties(gelex_cli PROPERTIES OUTPUT_NAME "gelex")

target_include_directories(gelex_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(gelex_cli PRIVATE gelex::core CLI11::CLI11)

install(TARGETS gelex_cli RUNTIME DESTINATION bin)
