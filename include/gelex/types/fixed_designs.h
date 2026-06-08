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

#ifndef GELEX_TYPES_FIXED_DESIGNS_H_
#define GELEX_TYPES_FIXED_DESIGNS_H_

#include <optional>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

namespace infra
{
class FieldVisitor;
}

struct QuantitativeCovariate
{
    std::vector<std::string> names;
    Eigen::MatrixXd X;
};

struct DiscreteCovariate
{
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> levels;
    std::vector<std::string> reference_levels;
    Eigen::MatrixXd X;
};

struct FixedDesign
{
    std::vector<std::string> names;
    std::vector<std::optional<std::vector<std::string>>> levels;
    std::vector<std::optional<std::string>> reference_levels;
    Eigen::MatrixXd X;
    Eigen::VectorXd XtX_diag;

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

    static auto make(
        std::optional<QuantitativeCovariate> qcovariate,
        std::optional<DiscreteCovariate> dcovariate) -> FixedDesign;

    static auto make(Eigen::Index n_samples) -> FixedDesign;

    auto visit(infra::FieldVisitor& visitor) const -> void;
};

}  // namespace gelex

#endif  // GELEX_TYPES_FIXED_DESIGNS_H_
