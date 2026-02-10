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

#include <barkeep.h>
#include <unistd.h>

#include <fmt/format.h>
#include <omp.h>
#include <Eigen/Core>
#include "config.h"

#include "gelex/logger.h"
#include "utils/formatter.h"

namespace gelex::cli
{
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

auto build_chr_groups(bool do_loco, const gelex::SnpEffects& snp_effects)
    -> std::vector<ChrGroup>
{
    std::vector<ChrGroup> groups;
    auto num_snps = static_cast<Eigen::Index>(snp_effects.size());

    if (do_loco)
    {
        std::string current_chr;
        Eigen::Index range_start = 0;

        for (Eigen::Index i = 0; i < num_snps; ++i)
        {
            if (snp_effects[i].chrom != current_chr)
            {
                if (!current_chr.empty())
                {
                    groups.push_back(
                        {current_chr, {{range_start, i}}, i - range_start});
                }
                current_chr = snp_effects[i].chrom;
                range_start = i;
            }
        }
        if (!current_chr.empty())
        {
            groups.push_back(
                {current_chr,
                 {{range_start, num_snps}},
                 num_snps - range_start});
        }
    }
    else
    {
        groups.push_back({"all", {{0, num_snps}}, num_snps});
    }
    return groups;
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
    $ gelex fit --bfile genotypes --pheno phenotypes.tsv --method RR --out results
    $ gelex predict --bfile genotypes --effects results.snp_effects --out pred
    $ gelex grm --bfile genotypes --out grm_output
    $ gelex assoc --bfile genotypes --pheno phenotypes.tsv --out gwas_results

Found a Bug or Have a Feature Request?
    Open an issue at: https://github.com/r1cheu/gelex/issues

For more information, see the documentation at: https://github.com/r1cheu/gelex
)";
}

auto print_fit_header(
    std::string_view model_name,
    bool has_dominance,
    int iters,
    int burn_in,
    int threads) -> void
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: Model Fitting (MCMC)", PROJECT_VERSION);
    std::string model_str = fmt::format(
        "Bayes{} ({})",
        model_name,
        has_dominance ? "Additive + Dominance" : "Additive");
    std::string chain_str = fmt::format(
        "{} iters ({} burn-in, {} sampling)", iters, burn_in, iters - burn_in);
    std::string comp_str = fmt::format("{}", threads);

    std::vector<std::pair<std::string, std::string>> items
        = {{"Model", model_str}, {"Chain", chain_str}, {"Threads", comp_str}};

    logger->info(gelex::header_box(title, items, 70));
    logger->info("");
}

auto print_grm_header(
    std::string_view method,
    bool do_additive,
    bool do_dominant,
    int chunk_size,
    int threads) -> void
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: GRM Computation", PROJECT_VERSION);

    std::string mode_str;
    if (do_additive && do_dominant)
    {
        mode_str = "Additive + Dominance";
    }
    else if (do_additive)
    {
        mode_str = "Additive";
    }
    else
    {
        mode_str = "Dominance";
    }

    std::string method_str = std::string(method);
    std::string chunk_str = fmt::format("{}", chunk_size);
    std::string comp_str = fmt::format("{}", threads);

    std::vector<std::pair<std::string, std::string>> items
        = {{"Method", method_str},
           {"Mode", mode_str},
           {"Chunk Size", chunk_str},
           {"Threads", comp_str}};

    logger->info(gelex::header_box(title, items, 70));
    logger->info("");
}

auto print_simulate_header(bool has_dominance) -> void
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: Phenotype Simulation", PROJECT_VERSION);
    std::string mode_str = has_dominance ? "Additive + Dominance" : "Additive";

    std::vector<std::pair<std::string, std::string>> items
        = {{"Mode", mode_str}};

    logger->info(gelex::header_box(title, items, 70));
    logger->info("");
}

auto print_assoc_header(int threads) -> void
{
    auto logger = gelex::logging::get();

    std::string title
        = fmt::format("gelex v{} :: GWAS Analysis", PROJECT_VERSION);
    std::vector<std::pair<std::string, std::string>> header_items
        = {{"Method", "AI-REML (Average Information)"},
           {"Threads", fmt::format("{}", threads)}};
    logger->info(gelex::header_box(title, header_items, 70));
    logger->info("");
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
