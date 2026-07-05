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

#ifndef APPS_CLI_ASSOC_CONFIG_H_
#define APPS_CLI_ASSOC_CONFIG_H_

#include <algorithm>
#include <string>
#include <thread>
#include <vector>

#include "cli/common_data.h"
#include "gelex/data/genotype_method.h"
#include "gelex/types/genetic_mode.h"

namespace cli
{

struct AssocConfig
{
    BaseDataConfig base_data;
    std::string bfile;
    std::vector<std::string> grm;
    std::string out{"gelex"};
    bool write_cov{false};
    gelex::GeneticModeSet mode{gelex::GeneticMode::A};
    bool loco{false};
    gelex::GenotypeMethod geno_method{gelex::GenotypeMethod::OrthCenterHWE};
    int max_iter{100};
    double tolerance{1e-6};
    int chunk_size{10000};
    int threads{
        std::max(1, static_cast<int>(std::thread::hardware_concurrency() / 2))};
};

}  // namespace cli

#endif  // APPS_CLI_ASSOC_CONFIG_H_
