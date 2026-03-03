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

#ifndef GELEX_CLI_FIT_CONFIG_H_
#define GELEX_CLI_FIT_CONFIG_H_

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gelex/algo/infer/params.h"
#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/types/effects.h"

namespace argparse
{
class ArgumentParser;
}

struct FitConfig
{
    std::string method_name;
    gelex::BayesAlphabet method;
    bool use_dominance;
    int threads;
    int seed;
    std::string out_prefix;
    gelex::MCMCParams mcmc_params;
    std::filesystem::path phenotype_path;
    int phenotype_column;
    std::filesystem::path bed_path;
    bool use_mmap;
    int chunk_size;
    std::filesystem::path qcovar_path;
    std::filesystem::path dcovar_path;
    gelex::GenotypeProcessMethod genotype_method;
    std::optional<std::vector<double>> pi;
    std::optional<std::vector<double>> dpi;
    std::optional<std::vector<double>> scale;
    std::optional<std::vector<double>> dscale;

    static auto make(argparse::ArgumentParser& cmd) -> FitConfig;
    auto validate() const -> void;
};

#endif  // GELEX_CLI_FIT_CONFIG_H_
