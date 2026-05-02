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

#ifndef GELEX_PIPELINE_POST_POSTERIOR_ANALYSIS_ENGINE_H_
#define GELEX_PIPELINE_POST_POSTERIOR_ANALYSIS_ENGINE_H_

#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/infra/logging/post_event.h"
#include "gelex/io/binary_reader.h"

namespace gelex
{

class PosteriorAnalysisEngine
{
   public:
    explicit PosteriorAnalysisEngine(
        std::span<const std::string_view> samples_prefix,
        double hdpi_threshold = 0.94,
        std::optional<std::string> gfile = std::nullopt);

    auto run(const PostObserver& observer = {}) -> void;

   private:
    double hdpi_threshold_;
    std::vector<io::detail::BinaryReader> readers_;
    std::optional<std::string> gfile_;

    auto check_consistency() const -> bool;
    auto process_gebv_variance() -> std::vector<ParameterDiag>;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_POST_POSTERIOR_ANALYSIS_ENGINE_H_
