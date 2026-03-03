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

#ifndef GELEX_CLI_CLI_HELPER_H_
#define GELEX_CLI_CLI_HELPER_H_

#include <argparse.h>
#include <algorithm>
#include <atomic>
#include <cctype>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <barkeep.h>
#include <Eigen/Core>

#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/exception.h"
#include "gelex/infra/detail/indicator.h"
#include "gelex/types/snp_info.h"

namespace gelex::cli
{

struct ChrGroup
{
    std::string name;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
    Eigen::Index total_snps;
};

inline auto parse_genotype_process_method(std::string_view value)
    -> GenotypeProcessMethod
{
    std::string lower(value);
    std::transform(
        lower.begin(),
        lower.end(),
        lower.begin(),
        [](unsigned char c) { return std::tolower(c); });

    static const std::unordered_map<std::string, GenotypeProcessMethod>
        kMethodMap = {
            {"standardizehwe", GenotypeProcessMethod::StandardizeHWE},
            {"sh", GenotypeProcessMethod::StandardizeHWE},
            {"centerhwe", GenotypeProcessMethod::CenterHWE},
            {"ch", GenotypeProcessMethod::CenterHWE},
            {"orthstandardizehwe", GenotypeProcessMethod::OrthStandardizeHWE},
            {"osh", GenotypeProcessMethod::OrthStandardizeHWE},
            {"orthcenterhwe", GenotypeProcessMethod::OrthCenterHWE},
            {"och", GenotypeProcessMethod::OrthCenterHWE},
            {"standardize", GenotypeProcessMethod::Standardize},
            {"s", GenotypeProcessMethod::Standardize},
            {"center", GenotypeProcessMethod::Center},
            {"c", GenotypeProcessMethod::Center},
            {"orthstandardize", GenotypeProcessMethod::OrthStandardize},
            {"os", GenotypeProcessMethod::OrthStandardize},
            {"orthcenter", GenotypeProcessMethod::OrthCenter},
            {"oc", GenotypeProcessMethod::OrthCenter},
        };

    auto it = kMethodMap.find(lower);
    if (it == kMethodMap.end())
    {
        throw gelex::InvalidInputException(
            "Invalid genotype process method: \"" + std::string(value)
            + "\". Valid: StandardizeHWE(SH), CenterHWE(CH),"
              " OrthStandardizeHWE(OSH), OrthCenterHWE(OCH),"
              " Standardize(S), Center(C), OrthStandardize(OS),"
              " OrthCenter(OC)");
    }
    return it->second;
}

auto is_tty() -> bool;

auto setup_parallelization(int num_threads) -> void;

auto build_chr_groups(bool do_loco, const gelex::SnpEffects& snp_effects)
    -> std::vector<ChrGroup>;

inline auto create_progress_bar(
    std::atomic<size_t>& counter,
    size_t total,
    std::string_view format = "{bar}") -> gelex::detail::ProgressBar
{
    return gelex::detail::create_progress_bar(counter, total, format);
}

auto print_gelex_banner_message(std::string_view version) -> void;

auto print_assoc_header(int threads) -> void;

auto format_epilog(std::string_view text) -> std::string;

}  // namespace gelex::cli

#endif  // GELEX_CLI_CLI_HELPER_H_
