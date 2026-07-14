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

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fmt/base.h>
#include <fmt/format.h>
#include <functional>
#include <iostream>
#include <omp.h>
#include <string>
#include <string_view>
#include <vector>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include "gelex/bayes/recipe_options.h"
#include "gelex/data/genotype_method.h"
#include "gelex/infra/stats/rank_inverse_norm_transform.h"
#include "gelex/types/genetic_mode.h"

#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/logging.h"
#include "cli/reml_data.h"
#include "report_printer.h"
#include "theme.h"
#include "version.h"

namespace gelex
{

auto lexical_cast(const std::string& input, GenotypeMethod& output) -> bool
{
    for (const auto& [code, method] : GENOTYPE_METHOD_CODES)
    {
        if (input.size() == code.size()
            && std::equal(
                code.begin(),
                code.end(),
                input.begin(),
                [](unsigned char expected, unsigned char actual)
                { return std::tolower(expected) == std::tolower(actual); }))
        {
            output = method;
            return true;
        }
    }
    return false;
}

auto lexical_cast(const std::string& input, GeneticMode& output) -> bool
{
    for (const auto& [mode, name] : GENETIC_MODE_NAMES)
    {
        if (input == name)
        {
            output = mode;
            return true;
        }
    }
    return false;
}

auto lexical_cast(const std::string& input, GeneticModeSet& output) -> bool
{
    for (const auto& [set, name] : GENETIC_MODE_SET_NAMES)
    {
        if (input == name)
        {
            output = set;
            return true;
        }
    }
    return false;
}

auto lexical_cast(const std::string& input, RintType& output) -> bool
{
    for (const auto& [type, name] : RINT_TYPE_NAMES)
    {
        if (input == name)
        {
            output = type;
            return true;
        }
    }
    return false;
}

}  // namespace gelex

namespace cli
{

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

auto add_random_effect_options(CLI::App& cmd, RemlDataConfig& config) -> void
{
    cmd.add_option(
           "--grm",
           config.grm,
           "GRM prefix(es) without suffix; reads <prefix>.bin/.id")
        ->group("I/O")
        ->type_name("<GRM>")
        ->expected(1, -1)
        ->allow_extra_args();
    cmd.add_option(
           "--drand",
           config.drand_path,
           "Discrete random-effect TSV (FID\tIID\tFactor\t...); each factor "
           "column becomes a variance component via one-hot ZZ^T")
        ->group("I/O")
        ->type_name("<DRAND>")
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "--qrand",
           config.qrand_paths,
           "Quantitative random-effect matrix TSV (FID\tIID\tValue\t...); "
           "each file forms one linear-kernel component ZZ^T")
        ->group("I/O")
        ->type_name("<QRAND>")
        ->expected(1, -1)
        ->allow_extra_args()
        ->check(CLI::ExistingFile);
    cmd.add_option(
           "--interaction",
           config.interactions,
           "Interaction random effect '<a>:<b>', the rescaled Hadamard product "
           "of two kernels. Each operand is a loaded effect name or a GRM "
           "prefix")
        ->group("I/O")
        ->type_name("<A:B>")
        ->expected(1, -1)
        ->allow_extra_args();
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

auto bayes_recipe_scheme_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::bayes::BAYES_RECIPE_SCHEME_NAMES.size());
    for (const auto& [scheme, name] : gelex::bayes::BAYES_RECIPE_SCHEME_NAMES)
    {
        names.emplace_back(name);
    }
    return CLI::IsMember(names);
}

auto genetic_mode_set_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::GENETIC_MODE_SET_NAMES.size());
    for (const auto& [set, name] : gelex::GENETIC_MODE_SET_NAMES)
    {
        names.emplace_back(name);
    }
    return CLI::IsMember(names);
}

auto rint_type_validator() -> CLI::Validator
{
    std::vector<std::string> names;
    names.reserve(gelex::RINT_TYPE_NAMES.size());
    for (const auto& [type, name] : gelex::RINT_TYPE_NAMES)
    {
        names.emplace_back(name);
    }
    return CLI::IsMember(names);
}

auto report_command_line(const CLI::App& cmd) -> void
{
    auto& p = cli::printer();

    auto cwd = std::filesystem::current_path().string();
    if (const char* home = std::getenv("HOME");
        home != nullptr && cwd.starts_with(home))
    {
        cwd.replace(0, std::string_view{home}.size(), "~");
    }
    p.block(cli::section("Command line ({})", cwd));

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
            cli::command_banner(PROJECT_VERSION, banner_title));
        report_command_line(cmd);
        auto start = std::chrono::steady_clock::now();
        auto result = execute_fn();
        auto elapsed = std::chrono::duration<double>(
                           std::chrono::steady_clock::now() - start)
                           .count();
        if (cli::logging::get())
        {
            cli::printer().block(cli::done_message(elapsed));
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
            std::cerr << error_marker() << e.what() << "\n";
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
            std::cerr << error_marker() << "unknown exception\n";
        }
        return 1;
    }
}

}  // namespace cli
