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

#ifndef APPS_CLI_GRM_CONFIG_H_
#define APPS_CLI_GRM_CONFIG_H_

#include <string>
#include <thread>

#include "gelex/data/genotype_method.h"
#include "gelex/types/genetic_mode.h"

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
