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

#ifndef GELEX_DATA_FIXED_DESIGN_H_
#define GELEX_DATA_FIXED_DESIGN_H_

#include <Eigen/Core>
#include <span>
#include <string>
#include <vector>

#include "gelex/data/covariates.h"

namespace gelex
{

class FixedDesign
{
   public:
    [[nodiscard]] auto quantitative_names() const noexcept
        -> std::span<const std::string>
    {
        return quantitative_names_;
    }

    [[nodiscard]] auto discrete_terms() const noexcept
        -> std::span<const DiscreteCovariateTerm>
    {
        return discrete_terms_;
    }

    [[nodiscard]] auto column_names() const noexcept
        -> std::span<const std::string>
    {
        return column_names_;
    }

    [[nodiscard]] auto X() const noexcept -> const Eigen::MatrixXd&
    {
        return matrix_;
    }

    [[nodiscard]] auto xtx_diag() const noexcept -> const Eigen::VectorXd&
    {
        return xtx_diag_;
    }

    static auto make(Eigen::Index n_samples) -> FixedDesign;

    static auto make(QuantitativeCovariate quantitative) -> FixedDesign;

    static auto make(DiscreteCovariate discrete) -> FixedDesign;

    static auto make(
        QuantitativeCovariate quantitative,
        DiscreteCovariate discrete) -> FixedDesign;

   private:
    FixedDesign() = default;

    std::vector<std::string> quantitative_names_;
    std::vector<DiscreteCovariateTerm> discrete_terms_;
    std::vector<std::string> column_names_;
    Eigen::MatrixXd matrix_;
    Eigen::VectorXd xtx_diag_;
};

}  // namespace gelex

#endif  // GELEX_DATA_FIXED_DESIGN_H_
