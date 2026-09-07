// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_ASSOC_CONFIG_H_
#define APPS_CLI_ASSOC_CONFIG_H_

#include <algorithm>
#include <string>
#include <thread>

#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

#include "cli/common_data.h"
#include "cli/reml_data.h"

namespace cli
{

struct AssocConfig
{
    BaseDataConfig base_data;
    RemlDataConfig random;
    std::string bfile;
    std::string out{"gelex"};
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
