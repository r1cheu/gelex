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

#include "gelex/data/pipe/grm.h"

#include <cstddef>
#include <filesystem>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/freq/design.h"
#include "gelex/infra/logging/grm_pipe_event.h"
#include "gelex/infra/logging/notify.h"

namespace gelex
{

GrmPipe::GrmPipe(
    std::vector<std::filesystem::path> grm_paths,
    GrmPipeObserver observer)
    : grm_paths_(std::move(grm_paths)), observer_(std::move(observer))
{
    sample_indices_.reserve(grm_paths_.size());
    grm_types_.reserve(grm_paths_.size());
    for (const auto& grm_path : grm_paths_)
    {
        sample_indices_.push_back(read_grm_ids(grm_path.string()));

        auto type = GeneticMode::A;
        const auto path = grm_path.string();
        if (!path.contains("add") && path.contains("dom"))
        {
            type = GeneticMode::D;
        }
        grm_types_.push_back(type);

        notify(
            observer_,
            GrmLoadedEvent{
                .num_samples
                = static_cast<size_t>(sample_indices_.back().size()),
                .type = grm_types_.back()});
    }
}

GrmPipe::~GrmPipe() = default;
GrmPipe::GrmPipe(GrmPipe&&) noexcept = default;
GrmPipe& GrmPipe::operator=(GrmPipe&&) noexcept = default;

auto GrmPipe::sample_indices() const
    -> std::vector<const dataframe::Index<std::string>*>
{
    return sample_indices_
           | std::views::transform([](const auto& sample_index)
                                   { return &sample_index; })
           | std::ranges::to<std::vector>();
}

auto GrmPipe::load(const dataframe::Index<std::string>& sample_index) -> void
{
    grms_.clear();
    grms_.reserve(grm_paths_.size());
    for (auto&& [i, grm_path] : std::views::enumerate(grm_paths_))
    {
        grms_.push_back(
            freq::GeneticDesign(
                grm_types_[static_cast<size_t>(i)],
                read_grm(grm_path.string(), &sample_index)));
    }
}

}  // namespace gelex
