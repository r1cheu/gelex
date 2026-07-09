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

#include "reml_data.h"

#include <cstddef>
#include <iterator>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

#include "gelex/data/covariates.h"
#include "gelex/data/reader.h"

namespace cli
{

RemlDataLoader::RemlDataLoader(const RemlDataConfig& config) noexcept
    : config_(config)
{
}

auto RemlDataLoader::load_indices(
    std::vector<const gelex::DataFrameIndex<std::string>*>& indices) -> void
{
    drand_ = config_.drand_path
                 ? std::make_optional(gelex::read_dcovar(*config_.drand_path))
                 : std::nullopt;
    if (drand_)
    {
        indices.push_back(&drand_->index());
    }

    qrand_.reserve(config_.qrand_paths.size());
    for (const auto& path : config_.qrand_paths)
    {
        qrand_.emplace_back(gelex::read_qcovar(path));
        indices.push_back(&qrand_.back().index());
    }

    grm_indices_.reserve(config_.grm.size());
    for (const auto& path : config_.grm)
    {
        grm_indices_.emplace_back(gelex::read_grm_ids(path));
        indices.push_back(&grm_indices_.back());
    }
}

auto RemlDataLoader::gather(
    const gelex::DataFrameIndex<std::string>& common_index) -> void
{
    if (drand_)
    {
        drand_->gather(common_index);
        random_designs_ = gelex::make_random_designs(*drand_);
    }
    for (auto&& [i, frame] : std::views::enumerate(qrand_))
    {
        frame.gather(common_index);
        random_designs_.push_back(
            gelex::make_quantitative_random_design(
                frame, config_.qrand_paths[static_cast<std::size_t>(i)]));
    }
    auto grm_designs = gelex::make_grm_designs(config_.grm, common_index);
    random_designs_.insert(
        random_designs_.end(),
        std::make_move_iterator(grm_designs.begin()),
        std::make_move_iterator(grm_designs.end()));
}

auto RemlDataLoader::results() && -> std::vector<gelex::freq::RandomDesign>
{
    return std::move(random_designs_);
}

}  // namespace cli
