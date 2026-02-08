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

#ifndef GELEX_TYPES_FIXED_EFFECTS_H_
#define GELEX_TYPES_FIXED_EFFECTS_H_

#include <optional>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "covariates.h"

namespace gelex
{

struct FixedEffect
{
    std::vector<std::string> names;
    std::vector<std::optional<std::vector<std::string>>> levels;
    std::vector<std::optional<std::string>> reference_levels;
    Eigen::MatrixXd X;
    Eigen::VectorXd cols_norm;

    struct CovariateInfoView
    {
        std::string_view name;
        std::span<const std::string> levels;
        std::string_view reference_level;
    };

    CovariateInfoView operator[](size_t i)
    {
        return CovariateInfoView{
            .name = names[i],
            .levels = levels[i] ? std::span<const std::string>(
                                      levels[i]->data(), levels[i]->size())
                                : std::span<const std::string>{},
            .reference_level
            = reference_levels[i] ? *reference_levels[i] : std::string_view{}};
    }

    // for freq model
    static auto build(
        std::optional<QuantitativeCovariate> qcovariate,
        std::optional<DiscreteCovariate> dcovariate) -> FixedEffect;

    static auto build(Eigen::Index n_samples) -> FixedEffect;
};

}  // namespace gelex

#endif  // GELEX_TYPES_FIXED_EFFECTS_H_
