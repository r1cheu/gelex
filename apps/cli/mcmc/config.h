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

#ifndef APPS_CLI_MCMC_CONFIG_H_
#define APPS_CLI_MCMC_CONFIG_H_

#include <algorithm>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

#include "cli/common_data.h"
#include "cli/random_design_data.h"

namespace cli
{

struct McmcConfig
{
    static constexpr auto option_modes
        = gelex::GeneticMode::A | gelex::GeneticMode::D;

    BaseDataConfig base_data;
    RandomDesignDataConfig random;
    std::string bfile;
    std::string out{"gelex"};
    gelex::GenotypeMethod geno_method{
        gelex::GenotypeMethod::OrthStandardizeHWE};
    gelex::BayesMethod method{gelex::BayesMethod::RR};
    gelex::GeneticModeSet mode{gelex::GeneticMode::A};
    gelex::HomogeneousModeValues<option_modes, std::optional<double>>
        genetic_variance_shares;
    gelex::JointModeValues<
        gelex::HomogeneousModeValues<option_modes, std::vector<double>>,
        std::vector<double>>
        mixture_probabilities;
    gelex::HomogeneousModeValues<option_modes, std::vector<double>>
        mixture_scales;
    std::optional<double> random_pve;
    int iters{5000};
    int burn_in{3000};
    int thin{1};
    int seed{42};
    int threads{
        std::max(1, static_cast<int>(std::thread::hardware_concurrency() / 2))};
};

auto validate_mcmc_config(const McmcConfig& config) -> void;

}  // namespace cli

#endif  // APPS_CLI_MCMC_CONFIG_H_
