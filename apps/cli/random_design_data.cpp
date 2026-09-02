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

#include "cli/random_design_data.h"

#include <filesystem>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"

namespace cli
{

RandomDesignDataLoader::RandomDesignDataLoader(
    const RandomDesignDataConfig& config)
    : config_{config}
{
}

auto RandomDesignDataLoader::load_indices(
    std::vector<const gelex::DataFrameIndex<std::string>*>& indices) -> void
{
    if (config_.drand_path)
    {
        data_.discrete = gelex::read_dcovar(*config_.drand_path);
        indices.push_back(&data_.discrete->index());
        auto names = data_.discrete->names();
        effect_names_.insert(
            effect_names_.end(),
            std::make_move_iterator(names.begin()),
            std::make_move_iterator(names.end()));
    }

    data_.quantitative.reserve(config_.qrand_paths.size());
    effect_names_.reserve(effect_names_.size() + config_.qrand_paths.size());
    for (const auto& path : config_.qrand_paths)
    {
        auto name = std::filesystem::path(path).stem().string();
        data_.quantitative.push_back(
            QuantitativeRandomData{
                .name = name, .frame = gelex::read_qcovar(path)});
        indices.push_back(&data_.quantitative.back().frame.index());
        effect_names_.push_back(std::move(name));
    }
}

auto RandomDesignDataLoader::gather(
    const gelex::DataFrameIndex<std::string>& common_index) -> void
{
    if (data_.discrete)
    {
        data_.discrete->gather(common_index);
    }
    for (auto& quantitative : data_.quantitative)
    {
        quantitative.frame.gather(common_index);
    }
}

auto RandomDesignDataLoader::results() && -> RandomDesignData
{
    return std::move(data_);
}

}  // namespace cli
