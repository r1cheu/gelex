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

#ifndef GELEX_DATA_PIPE_GRM_H_
#define GELEX_DATA_PIPE_GRM_H_

#include <filesystem>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/freq/design.h"
#include "gelex/infra/logging/grm_pipe_event.h"

namespace gelex
{

class GrmPipe
{
   public:
    explicit GrmPipe(
        std::vector<std::filesystem::path> grm_paths,
        GrmPipeObserver observer = {});
    GrmPipe(const GrmPipe&) = delete;
    GrmPipe(GrmPipe&&) noexcept;
    GrmPipe& operator=(const GrmPipe&) = delete;
    GrmPipe& operator=(GrmPipe&&) noexcept;
    ~GrmPipe();

    auto sample_indices() const
        -> std::vector<const dataframe::Index<std::string>*>;

    auto load(const dataframe::Index<std::string>& sample_index) -> void;

    auto take_random_designs() && -> std::vector<freq::RandomDesign>
    {
        return std::move(grms_);
    }

    auto grm_paths() const -> const std::vector<std::filesystem::path>&
    {
        return grm_paths_;
    }

   private:
    std::vector<std::filesystem::path> grm_paths_;
    std::vector<dataframe::Index<std::string>> sample_indices_;
    std::vector<freq::RandomDesign> grms_;
    GrmPipeObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_DATA_PIPE_GRM_H_
