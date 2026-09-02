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

#ifndef GELEX_BAYES_DETAIL_FITTED_VALUE_UPDATE_H_
#define GELEX_BAYES_DETAIL_FITTED_VALUE_UPDATE_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <span>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/component_layout.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

struct FittedValueTransition
{
    std::size_t old_class_index{};
    std::size_t new_class_index{};
    double old_coefficient{};
    double new_coefficient{};
};

template <GeneticMode Mode, ComponentLayout Layout, typename FittedMatrix>
auto apply_fitted_value_transition(
    const bayes::GeneticProjection& projection,
    Eigen::Index marker,
    FittedValueTransition transition,
    FittedMatrix& component_fitted_values,
    Eigen::VectorXd& adjusted_response) -> void
{
    std::array<bayes::AxpyTarget, 3> targets{};
    std::size_t target_count = 0;
    const double coefficient_difference
        = transition.old_coefficient - transition.new_coefficient;
    if (coefficient_difference != 0.0)
    {
        targets.at(target_count++)
            = bayes::AxpyTarget{coefficient_difference, adjusted_response};
    }

    const int old_component
        = Layout::component_index(Mode, transition.old_class_index);
    const int new_component
        = Layout::component_index(Mode, transition.new_class_index);
    if (old_component == new_component)
    {
        if (old_component != Layout::no_component
            && coefficient_difference != 0.0)
        {
            targets.at(target_count++) = bayes::AxpyTarget{
                -coefficient_difference,
                component_fitted_values.col(old_component)};
        }
    }
    else
    {
        if (old_component != Layout::no_component
            && transition.old_coefficient != 0.0)
        {
            targets.at(target_count++) = bayes::AxpyTarget{
                -transition.old_coefficient,
                component_fitted_values.col(old_component)};
        }
        if (new_component != Layout::no_component
            && transition.new_coefficient != 0.0)
        {
            targets.at(target_count++) = bayes::AxpyTarget{
                transition.new_coefficient,
                component_fitted_values.col(new_component)};
        }
    }
    projection.axpy(
        marker,
        std::span<const bayes::AxpyTarget>{targets.data(), target_count});
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_FITTED_VALUE_UPDATE_H_
