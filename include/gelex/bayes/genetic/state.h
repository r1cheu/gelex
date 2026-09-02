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

#ifndef GELEX_BAYES_GENETIC_STATE_H_
#define GELEX_BAYES_GENETIC_STATE_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "gelex/bayes/genetic_family.h"

namespace gelex
{

namespace detail
{
template <VarianceLayout Kind>
using marker_variance_state_t = std::
    conditional_t<Kind == VarianceLayout::Pooled, double, Eigen::VectorXd>;

}  // namespace detail

template <VarianceLayout Kind>
struct GaussianState
{
    detail::marker_variance_state_t<Kind> variance{};
    Eigen::VectorXd fitted_values;  // total
};

struct HalfNormalState
{
    double variance{};
    Eigen::Vector2d probit_coefficients = Eigen::Vector2d::Zero();
    Eigen::VectorXd fitted_values;  // total
};

template <VarianceLayout Kind>
struct SpikeSlabState
{
    detail::marker_variance_state_t<Kind> variance{};
    Eigen::VectorX<std::uint8_t> assignment;
    double probability{};
    Eigen::VectorXd fitted_values;  // total
};

template <std::size_t ClassCount>
    requires(
        ClassCount > 1
        && ClassCount <= static_cast<std::size_t>(
                             std::numeric_limits<std::uint8_t>::max())
                             + 1)
struct ScaledMixtureState
{
    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t component_count = class_count - 1;

    double variance{};
    Eigen::VectorX<std::uint8_t> assignment;
    std::array<double, class_count> probabilities{};
    Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(component_count)>
        fitted_values;
};

// Classes are NULL, A-only, D-only and AD; fitted_values holds one column per
// (mode, class) cell in which that mode is active, so every column carries a
// single mode and the two columns of a mode sum to that mode's total.
template <std::size_t ClassCount>
    requires(ClassCount == 4)
struct JointSpikeSlabState
{
    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t component_count = 4;
    static constexpr int no_component = -1;
    static constexpr std::array<int, class_count> additive_components{
        no_component,
        0,
        no_component,
        1};
    static constexpr std::array<int, class_count> dominance_components{
        no_component,
        no_component,
        2,
        3};

    Eigen::VectorX<std::uint8_t> assignment;
    std::array<double, class_count> probabilities{};
    Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(component_count)>
        fitted_values;
};

template <typename FamilyState>
struct GeneticModeState
{
    Eigen::VectorXd coefficients;
    FamilyState family_state;
};

template <VarianceLayout Kind>
[[nodiscard]] auto genetic_value(const GaussianState<Kind>& state)
    -> const Eigen::VectorXd&
{
    return state.fitted_values;
}

[[nodiscard]] inline auto genetic_value(const HalfNormalState& state)
    -> const Eigen::VectorXd&
{
    return state.fitted_values;
}

template <VarianceLayout Kind>
[[nodiscard]] auto genetic_value(const SpikeSlabState<Kind>& state)
    -> const Eigen::VectorXd&
{
    return state.fitted_values;
}

template <std::size_t ClassCount>
[[nodiscard]] auto genetic_value(const ScaledMixtureState<ClassCount>& state)
{
    return state.fitted_values.rowwise().sum();
}

template <typename FamilyState>
[[nodiscard]] auto genetic_value(const GeneticModeState<FamilyState>& state)
    -> decltype(auto)
{
    return genetic_value(state.family_state);
}

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_STATE_H_
