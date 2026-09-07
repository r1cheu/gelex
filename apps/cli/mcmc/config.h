// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
    std::string manno;
    std::string out{"gelex"};
    gelex::GenotypeMethod geno_method{
        gelex::GenotypeMethod::OrthStandardizeHWE};
    gelex::BayesMethod method{gelex::BayesMethod::RR};
    gelex::GeneticModeSet mode{gelex::GeneticMode::A};
    gelex::HomogeneousModeValues<option_modes, std::optional<double>>
        genetic_variance_proportion;
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
