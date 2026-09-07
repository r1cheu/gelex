// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_MCMC_DATA_H_
#define APPS_CLI_MCMC_DATA_H_

#include <string>
#include <vector>

#include "gelex/bayes/random_design.h"
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
