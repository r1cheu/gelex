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
#include <fstream>
#include <functional>

#include <Eigen/Core>

#include "gelex/types/mcmc_results.h"

namespace gelex
{

class ParameterWriter
{
   public:
    explicit ParameterWriter(const MCMCResult& result);

    void write(const std::filesystem::path& path) const;

   private:
    const MCMCResult* result_;

    void write_fixed_effects(std::ofstream& stream) const;
    void write_random_effects(std::ofstream& stream) const;
    void write_residual_variance(std::ofstream& stream) const;
    void write_additive_effect(std::ofstream& stream) const;
    void write_dominant_effect(std::ofstream& stream) const;

    static void write_genetic_effect(
        std::ofstream& stream,
        const std::string& variance_label,
        const std::string& heritability_label,
        const std::function<const BaseMarkerSummary*()>& effect_getter);

    static void write_summary_statistics(
        std::span<const std::string> terms,
        std::ofstream& stream,
        const PosteriorSummary& stats,
        Eigen::Index n_params);
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_PARAMETER_WRITER_H_
