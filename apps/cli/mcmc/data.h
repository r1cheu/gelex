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

#ifndef APPS_CLI_MCMC_DATA_H_
#define APPS_CLI_MCMC_DATA_H_

#include <string>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"

#include "cli/random_design_data.h"

namespace cli
{

struct McmcDesignData
{
    gelex::Bed bed;
    std::vector<gelex::bayes::RandomDesign> random;
};

class McmcDataLoader
{
   public:
    McmcDataLoader(gelex::Bed bed, const RandomDesignDataConfig& random_config);

    auto load_indices(
        std::vector<const gelex::DataFrameIndex<std::string>*>& indices)
        -> void;
    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void;
    auto results() && -> McmcDesignData;

   private:
    gelex::Bed bed_;
    RandomDesignDataLoader random_loader_;
};

}  // namespace cli

#endif  // APPS_CLI_MCMC_DATA_H_
