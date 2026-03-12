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

#ifndef GELEX_TYPES_SIM_TYPES_H_
#define GELEX_TYPES_SIM_TYPES_H_

#include <Eigen/Core>

namespace gelex
{

struct EffectSizeClass
{
    double proportion;
    double variance;
};

struct CausalEffects
{
    Eigen::VectorXd additive;
    Eigen::VectorXd dominance;
    Eigen::VectorXi add_class;
    Eigen::VectorXi dom_class;

    auto resize(Eigen::Index n_snps) -> void
    {
        additive.resize(n_snps);
        dominance.resize(n_snps);
        add_class.resize(n_snps);
        dom_class.resize(n_snps);
    }

    [[nodiscard]] auto size() const -> Eigen::Index { return additive.size(); }
};

}  // namespace gelex

#endif  // GELEX_TYPES_SIM_TYPES_H_
