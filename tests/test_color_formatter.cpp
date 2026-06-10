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

#include <memory>
#include <string>

#include <CLI/CLI.hpp>
#include <catch2/catch_test_macros.hpp>

#include "cli/color_formatter.h"

namespace cli
{

auto is_tty() -> bool
{
    return false;
}

}  // namespace cli

namespace
{

auto make_test_program() -> std::unique_ptr<CLI::App>
{
    auto program = std::make_unique<CLI::App>("Root description", "gelex");
    program->formatter(cli::make_cli_formatter());
    program->set_version_flag("-v,--version", "0.0.0");
    program->require_subcommand(1);

    auto* subcommand = program->add_subcommand("run");
    auto& command = *subcommand;
    command.description("Run the model");
    command.add_option("-p,--pheno", "Phenotype file")
        ->group("I/O")
        ->type_name("<PHENOTYPE>")
        ->required();
    command.add_option(
               "--grm", "GRM file prefix(es). Can specify multiple GRMs.")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->required();
    command.add_option("-o,--out", "Output file prefix")
        ->group("I/O")
        ->type_name("<OUT>")
        ->default_val(std::string{"gelex"});
    command.add_option("--pi", "Mixture proportions")
        ->group("Model")
        ->expected(1, -1)
        ->allow_extra_args();
    command.footer(
        "Docs:\n"
        "  https://example.invalid/run");

    return program;
}

}  // namespace

TEST_CASE("CLI color formatter preserves root help layout", "[cli][formatter]")
{
    auto program = make_test_program();

    REQUIRE(
        program->help()
        == R"(
Root description

Usage:
  gelex [OPTIONS] [COMMAND]

Commands:
  run               Run the model

Options:
  -h, --help        Print this help message and exit
  -v, --version     Display program version information and exit
)");
}

TEST_CASE(
    "CLI color formatter preserves subcommand help layout",
    "[cli][formatter]")
{
    auto program = make_test_program();

    REQUIRE(
        program->get_subcommand("run")->help()
        == R"(
Run the model

Usage:
  gelex run -p <PHENOTYPE> --grm <GRM> ... [OPTIONS]

I/O:
  -p, --pheno <PHENOTYPE>
                    Phenotype file
      --grm <GRM> ...
                    GRM file prefix(es). Can specify multiple GRMs.
  -o, --out <OUT>   Output file prefix [default: gelex]

Model:
      --pi ...      Mixture proportions

Options:
  -h, --help        Print this help message and exit

Docs:
  https://example.invalid/run
)");
}
