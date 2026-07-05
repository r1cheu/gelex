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

#include "cli/common_data.h"
#include "gelex/data/genotype_method.h"
#include "gelex/types/genetic_mode.h"

namespace cli
{

struct McmcConfig
{
    BaseDataConfig base_data;
    std::string bfile;
    std::string out{"gelex"};
    gelex::GenotypeMethod geno_method{
        gelex::GenotypeMethod::OrthStandardizeHWE};
    std::string method{"RR"};
    gelex::GeneticModeSet mode{gelex::GeneticMode::A};
    std::optional<double> h2;
    std::optional<double> d2;
    std::optional<double> dom_pos_prob;
    std::optional<double> random_pve;
    std::vector<double> pi;
    std::vector<double> dpi;
    std::vector<double> scale;
    std::vector<double> dscale;
    std::vector<double> jpi;
    std::optional<bool> sample_pi;
    std::optional<bool> sample_dpi;
    std::optional<bool> sample_jpi;
    int iters{5000};
    int burn_in{3000};
    int thin{1};
    int seed{42};
    int chunk_size{10000};
    int threads{
        std::max(1, static_cast<int>(std::thread::hardware_concurrency() / 2))};
    std::optional<int> checkpoint_step;
    std::optional<std::string> from_ckpt;
    bool mmap{false};
};

}  // namespace cli

#endif  // APPS_CLI_MCMC_CONFIG_H_
