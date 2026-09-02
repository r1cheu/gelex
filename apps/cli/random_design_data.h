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
