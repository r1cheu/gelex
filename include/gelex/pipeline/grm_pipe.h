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

#ifndef GELEX_PIPELINE_GRM_PIPE_H_
#define GELEX_PIPELINE_GRM_PIPE_H_

#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "gelex/data/genotype/sample_manager.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/types/freq_effect.h"

namespace gelex
{

namespace detail
{
class GrmLoader;
}

class GrmPipe
{
   public:
    explicit GrmPipe(
        std::vector<std::filesystem::path> grm_paths,
        DataPipeObserver observer = {});
    GrmPipe(const GrmPipe&) = delete;
    GrmPipe(GrmPipe&&) noexcept;
    GrmPipe& operator=(const GrmPipe&) = delete;
    GrmPipe& operator=(GrmPipe&&) noexcept;
    ~GrmPipe();

    auto sample_id_sets() const -> std::vector<std::span<const std::string>>;

    auto load(const std::shared_ptr<SampleManager>& sample_manager) -> void;

    auto take_grms() && -> std::vector<freq::GeneticEffect>
    {
        return std::move(grms_);
    }

    auto grm_paths() const -> const std::vector<std::filesystem::path>&
    {
        return grm_paths_;
    }

   private:
    std::vector<std::filesystem::path> grm_paths_;
    std::vector<detail::GrmLoader> grm_loaders_;
    std::vector<freq::GeneticEffect> grms_;
    DataPipeObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_GRM_PIPE_H_
