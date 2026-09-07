// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_GRM_CONFIG_H_
#define APPS_CLI_GRM_CONFIG_H_

#include <string>
#include <thread>

#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

namespace cli
{

struct GrmConfig
{
    std::string bfile;
    std::string out{"grm"};
    gelex::GenotypeMethod geno_method{
        gelex::GenotypeMethod::OrthStandardizeHWE};
    gelex::GeneticModeSet mode{gelex::GeneticMode::A};
    bool loco{false};
    int chunk_size{10000};
    int threads{static_cast<int>(std::thread::hardware_concurrency() / 2)};
};

}  // namespace cli

#endif  // APPS_CLI_GRM_CONFIG_H_
