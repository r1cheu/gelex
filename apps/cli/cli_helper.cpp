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

#include "cli_helper.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/base.h>
#include <fmt/format.h>
#include <omp.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/logging.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "report_printer.h"
#include "version.h"

namespace cli
{

namespace
{

constexpr std::string_view ERROR_MARKER = "[\033[31merror\033[0m] ";

}  // namespace

auto parse_genotype_method(std::string_view value) -> gelex::GenotypeMethod
{
    for (const auto& [code, method] : gelex::GENOTYPE_METHOD_CODES)
    {
        if (value.size() != code.size())
        {
            continue;
        }
        const bool matches = std::equal(
            code.begin(),
            code.end(),
            value.begin(),
            [](unsigned char expected, unsigned char actual)
            { return std::tolower(expected) == std::tolower(actual); });
        if (matches)
        {
            return method;
        }
    }

    throw gelex::GelexException(
        fmt::format(
            "invalid genotype process method: \"{}\". Valid: "
            "SH, CH, OSH, OCH, S, C, OS, OC, NS, NC",
            value));
}

auto parse_genetic_modes(std::string_view sv) -> std::vector<gelex::GeneticMode>
{
    if (sv == "A")
    {
        return {gelex::GeneticMode::A};
    }
    if (sv == "D")
    {
        return {gelex::GeneticMode::D};
    }
    if (sv == "AD")
    {
        return {gelex::GeneticMode::A, gelex::GeneticMode::D};
    }
    throw gelex::GelexException(
        fmt::format("invalid --mode: \"{}\". Valid: A, D, AD", sv));
}

auto is_tty() -> bool
{
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

auto setup_parallelization(int num_threads) -> void
{
    if (num_threads > 0)
    {
        omp_set_num_threads(num_threads);
        Eigen::setNbThreads(num_threads);
    }
}

auto add_common_io_options(CLI::App& cmd, BaseDataConfig& config) -> void
{
    cmd.add_option(
           "-p,--pheno",
           config.pheno_path,
           "Phenotype TSV with FID, IID, trait columns")
        ->group("I/O")
        ->type_name("<PHENOTYPE>")
        ->check(CLI::ExistingFile)
        ->required();
    cmd.add_option(
           "--pheno-col",
           config.pheno_col,
           "0-based phenotype column after FID/IID; first trait=0")
        ->group("I/O")
        ->type_name("<COL>")
        ->check(cli::non_negative_number())
        ->capture_default_str();
    cmd.add_option(
           "--qcovar",
           config.qcovar_path,
           "Quantitative covariate TSV with FID, IID, numeric columns")
        ->group("I/O")
        ->type_name("<QCOVAR>")
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "--dcovar",
           config.dcovar_path,
           "Discrete covariate TSV with FID, IID, factor columns")
        ->group("I/O")
        ->type_name("<DCOVAR>")
        ->check(CLI::ExistingFile);
}

auto open_unit_interval() -> CLI::Validator
{
    return CLI::Validator{
        [](std::string& input)
        {
            double value{};
            if (!CLI::detail::lexical_cast(input, value) || value <= 0.0
                || value >= 1.0)
            {
                return std::string{
                    "expected a value in the open interval (0, 1)"};
            }
            return std::string{};
        },
        "FLOAT in (0 - 1)",
        "PROB"};
}

auto non_negative_number() -> CLI::Validator
{
    return CLI::Validator{
        [](std::string& input)
        {
            double value{};
            if (!CLI::detail::lexical_cast(input, value) || !(value >= 0.0))
            {
                return std::string{"expected a non-negative value"};
            }
            return std::string{};
        },
        "NONNEGATIVE",
        "NONNEGATIVE"};
}

auto genotype_method_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::GENOTYPE_METHOD_CODES.size());
    for (const auto& entry : gelex::GENOTYPE_METHOD_CODES)
    {
        names.emplace_back(entry.first);
    }
    return CLI::IsMember(names, CLI::ignore_case);
}

auto report_command_line(const CLI::App& cmd) -> void
{
    auto& p = cli::printer();
    p.block(gelex::section("Command line:"));

    std::vector<const CLI::Option*> printed;
    for (const auto* option : cmd.parse_order())
    {
        if (option->count() == 0)
        {
            continue;
        }
        if (std::ranges::find(printed, option) != printed.end())
        {
            continue;
        }
        printed.push_back(option);

        std::string name;
        const auto& long_names = option->get_lnames();
        const auto& short_names = option->get_snames();
        if (!long_names.empty())
        {
            name = "--" + long_names.front();
        }
        else if (!short_names.empty())
        {
            name = "-" + short_names.front();
        }
        if (name.empty())
        {
            continue;
        }

        std::vector<std::string> values;
        for (const auto& value : option->results())
        {
            if (!value.empty() && value != "true")
            {
                values.push_back(value);
            }
        }
        if (values.empty())
        {
            p.line("  {}", name);
        }
        else
        {
            p.line("  {} {}", name, fmt::join(values, " "));
        }
    }
}

auto execute_cli_command(
    const CLI::App& cmd,
    std::string_view banner_title,
    const std::function<int()>& execute_fn) -> int
{
    try
    {
        cli::logging::initialize(cmd.get_option("--out")->as<std::string>());
        cli::printer().block(
            gelex::command_banner(PROJECT_VERSION, banner_title));
        report_command_line(cmd);
        auto start = std::chrono::steady_clock::now();
        auto result = execute_fn();
        auto elapsed = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - start)
                           .count();
        if (cli::logging::get())
        {
            cli::printer().block(gelex::done_message(elapsed));
        }
        return result;
    }
    catch (const std::exception& e)
    {
        auto logger = cli::logging::get();
        if (logger)
        {
            logger->error("{}", e.what());
        }
        else
        {
            std::cerr << ERROR_MARKER << e.what() << "\n";
        }
        return 1;
    }
    catch (...)
    {
        auto logger = cli::logging::get();
        if (logger)
        {
            logger->error("unknown exception");
        }
        else
        {
            std::cerr << ERROR_MARKER << "unknown exception\n";
        }
        return 1;
    }
}

}  // namespace cli
