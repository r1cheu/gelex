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

#include "gelex/pipeline/grm_pipe.h"

#include <filesystem>
#include <ranges>
#include <utility>
#include <vector>

#include "gelex/data/grm/detail/grm_reader.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/notify.h"

namespace gelex
{

GrmPipe::GrmPipe(
    std::vector<std::filesystem::path> grm_paths,
    DataPipeObserver observer)
    : grm_paths_(std::move(grm_paths)), observer_(std::move(observer))
{
    for (const auto& grm_path : grm_paths_)
    {
        grm_readers_.emplace_back(grm_path);
        notify(
            observer_,
            GrmLoadedEvent{
                .num_samples
                = static_cast<size_t>(grm_readers_.back().num_samples()),
                .type = grm_readers_.back().type()});
    }
}

GrmPipe::~GrmPipe() = default;
GrmPipe::GrmPipe(GrmPipe&&) noexcept = default;
GrmPipe& GrmPipe::operator=(GrmPipe&&) noexcept = default;

auto GrmPipe::sample_indices() const
    -> std::vector<const dataframe::Index<std::string>*>
{
    return grm_readers_
           | std::views::transform([](const auto& r)
                                   { return &r.sample_index(); })
           | std::ranges::to<std::vector>();
}

auto GrmPipe::load(const dataframe::Index<std::string>& sample_index) -> void
{
    grms_ = grm_readers_
            | std::views::transform(
                [&](auto& r)
                { return freq::GeneticEffect(r.type(), r.load(sample_index)); })
            | std::ranges::to<std::vector>();
}

}  // namespace gelex
