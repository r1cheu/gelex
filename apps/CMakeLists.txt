set(CLI_SOURCES
    main.cpp
    cli/assoc_args.cpp
    cli/assoc_command.cpp
    cli/fit_args.cpp
    cli/fit_command.cpp
    cli/grm_args.cpp
    cli/grm_command.cpp
    cli/gwas_runner.cpp
    cli/predict_args.cpp
    cli/predict_command.cpp
    cli/simulate_args.cpp
    cli/simulate_command.cpp
    cli/cli_helper.cpp)

add_executable(gelex_cli ${CLI_SOURCES})
set_target_properties(gelex_cli PROPERTIES OUTPUT_NAME "gelex")

target_include_directories(gelex_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                             ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(gelex_cli PRIVATE gelex::core)

install(TARGETS gelex_cli RUNTIME DESTINATION bin)
