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

#include <array>
#include <chrono>
#include <string_view>

#include <argparse.h>

#include "config.h"

#include "cli/assoc/assoc_args.h"
#include "cli/assoc/assoc_command.h"
#include "cli/cli_helper.h"
#include "cli/fit/fit_args.h"
#include "cli/fit/fit_command.h"
#include "cli/grm/grm_args.h"
#include "cli/grm/grm_command.h"
#include "cli/post/post_args.h"
#include "cli/post/post_command.h"
#include "cli/predict/predict_args.h"
#include "cli/predict/predict_command.h"
#include "cli/simulate/simulate_args.h"
#include "cli/simulate/simulate_command.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/utils/formatter.h"

namespace
{
constexpr std::string_view kErrorMarker = "[\033[31merror\033[0m] ";
}  // namespace

struct CommandDescriptor
{
    std::string_view name;
    argparse::ArgumentParser* parser;
    void (*setup_fn)(argparse::ArgumentParser&);
    int (*execute_fn)(argparse::ArgumentParser&);
};

auto execute_command(
    argparse::ArgumentParser& parser,
    int (*execute_fn)(argparse::ArgumentParser&)) -> int
{
    try
    {
        gelex::logging::initialize(parser.get("--out"));
        auto start = std::chrono::steady_clock::now();
        auto result = execute_fn(parser);
        auto elapsed = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - start)
                           .count();
        if (auto logger = gelex::logging::get())
        {
            logger->info("");
            logger->info("{}", gelex::done_message(elapsed));
        }
        return result;
    }
    catch (const std::exception& e)
    {
        auto logger = gelex::logging::get();
        if (logger)
        {
            logger->error("{}", e.what());
        }
        else
        {
            std::cerr << kErrorMarker << e.what() << "\n";
        }
        return 1;
    }
    catch (...)
    {
        auto logger = gelex::logging::get();
        if (logger)
        {
            logger->error("unknown exception");
        }
        else
        {
            std::cerr << kErrorMarker << "unknown exception\n";
        }
        return 1;
    }
}

auto main(int argc, char* argv[]) -> int
{
    argparse::ArgumentParser program(PROJECT_NAME, PROJECT_VERSION);
    argparse::ArgumentParser fit("fit");
    argparse::ArgumentParser simulate("simulate");
    argparse::ArgumentParser predict("predict");
    argparse::ArgumentParser grm("grm");
    argparse::ArgumentParser assoc("assoc");
    argparse::ArgumentParser post("post");

    const std::array commands
        = {CommandDescriptor{"fit", &fit, setup_fit_args, fit_execute},
           CommandDescriptor{
               "simulate", &simulate, setup_simulate_args, simulate_execute},
           CommandDescriptor{
               "predict", &predict, setup_predict_args, predict_execute},
           CommandDescriptor{"grm", &grm, setup_grm_args, grm_execute},
           CommandDescriptor{"assoc", &assoc, setup_assoc_args, assoc_execute},
           CommandDescriptor{"post", &post, setup_post_args, post_execute}};

    for (const auto& cmd : commands)
    {
        cmd.setup_fn(*cmd.parser);
        program.add_subparser(*cmd.parser);
    }

    if (argc <= 1)
    {
        gelex::cli::print_gelex_banner_message(PROJECT_VERSION);
        std::cerr << "\n" << program;
        return 1;
    }

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err)
    {
        std::cerr << kErrorMarker << err.what() << "\n";

        bool found = false;
        for (const auto& cmd : commands)
        {
            if (program.is_subcommand_used(cmd.name))
            {
                std::cerr << *cmd.parser;
                found = true;
                break;
            }
        }
        if (!found)
        {
            std::cerr << program;
        }
        return 1;
    }

    for (const auto& cmd : commands)
    {
        if (program.is_subcommand_used(cmd.name))
        {
            return execute_command(*cmd.parser, cmd.execute_fn);
        }
    }

    std::cerr << program;
    return 1;
}
