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

#ifndef GELEX_BAYES_GENETIC_DETAIL_FITTED_UPDATE_H_
#define GELEX_BAYES_GENETIC_DETAIL_FITTED_UPDATE_H_

#include <Eigen/Core>
#include <array>
#include <variant>

#include "gelex/bayes/genotype/operations.h"

namespace gelex::detail
{

// A negative component index denotes an absent contribution.
[[nodiscard]] inline auto make_fitted_update(
    Eigen::Ref<Eigen::MatrixXd> fitted_values,
    Eigen::Index old_component,
    Eigen::Index new_component,
    double old_coefficient,
    double new_coefficient) -> std::
    variant<std::monostate, bayes::AxpyTarget, std::array<bayes::AxpyTarget, 2>>
{
    const auto make_target = [&](Eigen::Index component, double delta)
    { return bayes::AxpyTarget{delta, fitted_values.col(component)}; };
    if (old_component == new_component)
    {
        const double delta = new_coefficient - old_coefficient;
        if (old_component < 0 || delta == 0.0)
        {
            return std::monostate{};
        }
        return make_target(old_component, delta);
    }
    const bool remove_old = old_component >= 0 && old_coefficient != 0.0;
    const bool add_new = new_component >= 0 && new_coefficient != 0.0;
    if (remove_old && add_new)
    {
        return std::array{
            make_target(old_component, -old_coefficient),
            make_target(new_component, new_coefficient)};
    }
    if (remove_old)
    {
        return make_target(old_component, -old_coefficient);
    }
    if (add_new)
    {
        return make_target(new_component, new_coefficient);
    }
    return std::monostate{};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_FITTED_UPDATE_H_
