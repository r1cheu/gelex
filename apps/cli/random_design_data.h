// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_RANDOM_DESIGN_DATA_H_
#define APPS_CLI_RANDOM_DESIGN_DATA_H_

#include <optional>
#include <span>
#include <string>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"

namespace cli
{

struct RandomDesignDataConfig
{
    std::optional<std::string> drand_path;
    std::vector<std::string> qrand_paths;

    [[nodiscard]] auto has_random_design() const noexcept -> bool
    {
        return drand_path.has_value() || !qrand_paths.empty();
    }
};

struct QuantitativeRandomData
{
    std::string name;
    gelex::DataFrame<std::string> frame;
};

struct RandomDesignData
{
    std::optional<gelex::DataFrame<std::string>> discrete;
    std::vector<QuantitativeRandomData> quantitative;
};

class RandomDesignDataLoader
{
   public:
    explicit RandomDesignDataLoader(const RandomDesignDataConfig& config);

    auto load_indices(
        std::vector<const gelex::DataFrameIndex<std::string>*>& indices)
        -> void;
    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void;
    auto results() && -> RandomDesignData;

    [[nodiscard]] auto effect_names() const noexcept
        -> std::span<const std::string>
    {
        return effect_names_;
    }

   private:
    RandomDesignDataConfig config_;
    RandomDesignData data_;
    std::vector<std::string> effect_names_;
};

}  // namespace cli

#endif  // APPS_CLI_RANDOM_DESIGN_DATA_H_
