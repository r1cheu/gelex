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
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/grm/grm_loader.h"
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
        grm_loaders_.emplace_back(grm_path);
        notify(
            observer_,
            GrmLoadedEvent{
                .num_samples
                = static_cast<size_t>(grm_loaders_.back().num_samples()),
                .type = grm_loaders_.back().type()});
    }
}

GrmPipe::~GrmPipe() = default;
GrmPipe::GrmPipe(GrmPipe&&) noexcept = default;
GrmPipe& GrmPipe::operator=(GrmPipe&&) noexcept = default;

auto GrmPipe::sample_id_sets() const
    -> std::vector<std::span<const std::string>>
{
    std::vector<std::span<const std::string>> result;
    result.reserve(grm_loaders_.size());
    for (const auto& loader : grm_loaders_)
    {
        result.emplace_back(loader.sample_ids());
    }
    return result;
}

auto GrmPipe::load(const std::shared_ptr<SampleManager>& sample_manager) -> void
{
    const auto& id_map = sample_manager->common_id_map();
    grms_.reserve(grm_loaders_.size());
    for (auto& loader : grm_loaders_)
    {
        grms_.emplace_back(loader.type(), loader.load(id_map));
    }
}

}  // namespace gelex
