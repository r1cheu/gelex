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

#ifndef GELEX_ESTIMATOR_BAYES_PARAMETER_WRITER_H_
#define GELEX_ESTIMATOR_BAYES_PARAMETER_WRITER_H_

#include <filesystem>
#include <memory>
#include <span>
#include <string>

#include <Eigen/Core>

#include "gelex/types/mcmc_results.h"

namespace gelex::detail
{
class TextWriter;
}

namespace gelex
{

class ParameterWriter
{
   public:
    ParameterWriter(
        const MCMCResult& result,
        const std::filesystem::path& output_path);
    ~ParameterWriter();

    auto write() -> void;

   private:
    const MCMCResult* result_;
    std::unique_ptr<detail::TextWriter> writer_;

    auto write_fixed_effects() -> void;
    auto write_random_effects() -> void;
    auto write_residual_variance() -> void;

    auto write_genetic_effect(
        const std::string& variance_label,
        const std::string& heritability_label,
        const BaseMarkerSummary* effect) -> void;

    auto write_summary_statistics(
        std::span<const std::string> terms,
        const PosteriorSummary& stats,
        Eigen::Index n_params) -> void;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_PARAMETER_WRITER_H_
