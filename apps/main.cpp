/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <iostream>
#include <string_view>

#include <CLI/CLI.hpp>

#include "version.h"

#include "cli/assoc/args.h"
#include "cli/cli_helper.h"
#include "cli/color_formatter.h"
#include "cli/grm/args.h"
#include "cli/mcmc/args.h"
#include "cli/post/args.h"
#include "cli/predict/args.h"
#include "cli/reml/args.h"
#include "cli/simulate/args.h"

namespace
{
constexpr std::string_view ERROR_MARKER = "[\033[31merror\033[0m] ";
}  // namespace

auto main(int argc, char* argv[]) -> int
{
    CLI::App program{
        "High-Performance Genomic Prediction with Bayesian and Frequentist "
        "Models",
        PROJECT_NAME};
    program.formatter(cli::make_cli_formatter());
    program.set_version_flag("-v,--version", PROJECT_VERSION);
    program.require_subcommand(1);

    int exit_code = 0;
    setup_mcmc_command(program, exit_code);
    setup_simulate_command(program, exit_code);
    setup_predict_command(program, exit_code);
    setup_grm_command(program, exit_code);
    setup_assoc_command(program, exit_code);
    setup_post_command(program, exit_code);
    setup_reml_command(program, exit_code);

    if (argc <= 1)
    {
        cli::print_gelex_banner_message(PROJECT_VERSION);
        std::cerr << "\n" << program.help();
        return 1;
    }

    try
    {
        program.parse(argc, argv);
    }
    catch (const CLI::ParseError& err)
    {
        if (err.get_exit_code() == 0)
        {
            if (argc > 1
                && (std::string_view{argv[1]} == "-h"
                    || std::string_view{argv[1]} == "--help"))
            {
                cli::print_gelex_banner_message(PROJECT_VERSION);
                std::cout << "\n";
            }
            return program.exit(err);
        }

        std::cerr << ERROR_MARKER << err.what() << "\n";

        if (argc > 1)
        {
            auto* subcommand = program.get_subcommand_no_throw(argv[1]);
            if (subcommand != nullptr)
            {
                std::cerr << subcommand->help();
                return 1;
            }
        }
        std::cerr << program.help();
        return 1;
    }

    return exit_code;
}
