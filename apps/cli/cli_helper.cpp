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

#include <iostream>
#include <string>
#include <string_view>

#include <argparse.h>
#include <fmt/base.h>
#include <fmt/format.h>
#include <omp.h>
#include <Eigen/Core>

#include "gelex/data/genotype/method.h"
#include "gelex/exception.h"

namespace gelex::cli
{

auto parse_genotype_method(std::string_view value) -> GenotypeMethod
{
    std::string lower(value);
    std::transform(
        lower.begin(),
        lower.end(),
        lower.begin(),
        [](unsigned char c) { return std::tolower(c); });

    static const std::unordered_map<std::string, GenotypeMethod> METHOD_MAP = {
        {"standardizehwe", GenotypeMethod::StandardizeHWE},
        {"sh", GenotypeMethod::StandardizeHWE},
        {"centerhwe", GenotypeMethod::CenterHWE},
        {"ch", GenotypeMethod::CenterHWE},
        {"orthstandardizehwe", GenotypeMethod::OrthStandardizeHWE},
        {"osh", GenotypeMethod::OrthStandardizeHWE},
        {"orthcenterhwe", GenotypeMethod::OrthCenterHWE},
        {"och", GenotypeMethod::OrthCenterHWE},
        {"standardize", GenotypeMethod::Standardize},
        {"s", GenotypeMethod::Standardize},
        {"center", GenotypeMethod::Center},
        {"c", GenotypeMethod::Center},
        {"orthstandardize", GenotypeMethod::OrthStandardize},
        {"os", GenotypeMethod::OrthStandardize},
        {"orthcenter", GenotypeMethod::OrthCenter},
        {"oc", GenotypeMethod::OrthCenter},
        {"noiastandardize", GenotypeMethod::NOIAStandardize},
        {"ns", GenotypeMethod::NOIAStandardize},
        {"noiacenter", GenotypeMethod::NOIACenter},
        {"nc", GenotypeMethod::NOIACenter},
    };

    auto it = METHOD_MAP.find(lower);
    if (it == METHOD_MAP.end())
    {
        throw gelex::GelexException(
            fmt::format(
                "Invalid genotype process method: \"{}\". Valid: "
                "StandardizeHWE(SH), CenterHWE(CH), "
                "OrthStandardizeHWE(OSH), OrthCenterHWE(OCH), "
                "Standardize(S), Center(C), OrthStandardize(OS), "
                "OrthCenter(OC), NOIAStandardize(NS), NOIACenter(NC)",
                value));
    }
    return it->second;
}

auto parse_genetic_modes(std::string_view sv) -> std::vector<GeneticMode>
{
    if (sv == "A")
    {
        return {GeneticMode::A};
    }
    if (sv == "D")
    {
        return {GeneticMode::D};
    }
    if (sv == "AD")
    {
        return {GeneticMode::A, GeneticMode::D};
    }
    throw gelex::GelexException(
        fmt::format("invalid --mode: \"{}\". Valid: A, D, AD", sv));
}

auto is_tty() -> bool
{
    return isatty(fileno(stdout)) != 0;
}

auto setup_parallelization(int num_threads) -> void
{
    if (num_threads > 0)
    {
        omp_set_num_threads(num_threads);
        Eigen::setNbThreads(num_threads);
    }
}

auto print_gelex_banner_message(std::string_view version) -> void
{
    std::cout << "Gelex [version " << version
              << "] - High-Performance Genomic Prediction with Bayesian and "
                 "Frequentist Models\n\n";
    std::cout
        << R"(Gelex is a specialized CLI tool designed for genomic selection and prediction in breeding
programs and quantitative genetics research. Built on modern C++23 with memory-mapped I/O
and BLAS/LAPACK acceleration, Gelex offers seamless integration with PLINK binary formats
and efficient processing of large-scale genomic data.

Basic Usage:
    $ gelex mcmc --bfile genotypes --pheno phenotypes.tsv --method RR --out results
    $ gelex predict --bfile genotypes --effects results.snp_effects --out pred
    $ gelex grm --bfile genotypes --out grm_output
    $ gelex assoc --bfile genotypes --pheno phenotypes.tsv --out gwas_results

Found a Bug or Have a Feature Request?
    Open an issue at: https://github.com/r1cheu/gelex/issues

For more information, see the documentation at: https://github.com/r1cheu/gelex
)";
}

auto format_epilog(std::string_view text) -> std::string
{
    namespace c = argparse::colors;
    const bool enabled = c::enabled();
    std::string bg = enabled ? fmt::format("{}{}", c::BOLD, c::GREEN) : "";
    std::string bc = enabled ? fmt::format("{}{}", c::BOLD, c::CYAN) : "";
    std::string cy(enabled ? c::CYAN : "");
    std::string gy = enabled ? "\033[90m" : "";
    std::string rs(enabled ? c::RESET : "");

    return fmt::format(
        fmt::runtime(text),
        fmt::arg("bg", bg),
        fmt::arg("bc", bc),
        fmt::arg("cy", cy),
        fmt::arg("gy", gy),
        fmt::arg("rs", rs));
}

}  // namespace gelex::cli
