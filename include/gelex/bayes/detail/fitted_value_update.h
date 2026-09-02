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
#include "gelex/bayes/genetic/state.h"

namespace gelex::detail
{

struct FittedValueTransition
{
    std::size_t old_class_index{};
    std::size_t new_class_index{};
    double old_coefficient{};
    double new_coefficient{};
};

template <std::size_t ClassCount>
auto apply_fitted_value_transition(
    const bayes::GeneticProjection& projection,
    Eigen::Index marker,
    FittedValueTransition transition,
    ScaledMixtureState<ClassCount>& state,
    Eigen::VectorXd& adjusted_response) -> void
{
    std::array<bayes::AxpyTarget, 3> targets{};
    std::size_t target_count = 0;
    const double coefficient_delta
        = transition.new_coefficient - transition.old_coefficient;
    if (coefficient_delta != 0.0)
    {
        targets.at(target_count++)
            = bayes::AxpyTarget{-coefficient_delta, adjusted_response};
    }

    const auto append_component_target
        = [&](std::size_t class_index, double delta)
    {
        if (class_index == 0 || delta == 0.0)
        {
            return;
        }

        const auto component_index = static_cast<Eigen::Index>(class_index - 1);
        targets.at(target_count++) = bayes::AxpyTarget{
            delta, state.fitted_values.col(component_index)};
    };

    if (transition.old_class_index == transition.new_class_index)
    {
        append_component_target(transition.old_class_index, coefficient_delta);
    }
    else
    {
        append_component_target(
            transition.old_class_index, -transition.old_coefficient);
        append_component_target(
            transition.new_class_index, transition.new_coefficient);
    }

    if (target_count != 0)
    {
        projection.axpy(
            marker,
            std::span<const bayes::AxpyTarget>{targets.data(), target_count});
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_FITTED_VALUE_UPDATE_H_
