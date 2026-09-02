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

#include "data.h"

#include <utility>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"

#include "cli/random_design_data.h"

namespace cli
{

McmcDataLoader::McmcDataLoader(
    gelex::Bed bed,
    const RandomDesignDataConfig& random_config)
    : bed_{std::move(bed)}, random_loader_{random_config}
{
}

auto McmcDataLoader::load_indices(
    std::vector<const gelex::DataFrameIndex<std::string>*>& indices) -> void
{
    indices.push_back(&bed_.sample_index());
    random_loader_.load_indices(indices);
}

auto McmcDataLoader::gather(
    const gelex::DataFrameIndex<std::string>& common_index) -> void
{
    bed_.gather(common_index);
    random_loader_.gather(common_index);
}

auto McmcDataLoader::results() && -> McmcDesignData
{
    auto random_data = std::move(random_loader_).results();
    std::vector<gelex::bayes::RandomDesign> random;
    if (random_data.discrete)
    {
        random = gelex::bayes::make_random_designs(*random_data.discrete);
    }
    random.reserve(random.size() + random_data.quantitative.size());
    for (auto& quantitative : random_data.quantitative)
    {
        random.push_back(
            gelex::bayes::make_quantitative_random_design(
                quantitative.frame, std::move(quantitative.name)));
    }
    return McmcDesignData{.bed = std::move(bed_), .random = std::move(random)};
}

}  // namespace cli
