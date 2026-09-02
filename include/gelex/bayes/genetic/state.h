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

#include "gelex/bayes/semantic_method.h"

namespace gelex
{

namespace detail
{

template <VarianceLayout Kind>
struct MarkerVarianceStateType;

template <>
struct MarkerVarianceStateType<VarianceLayout::Pooled>
{
    using type = double;
};

template <>
struct MarkerVarianceStateType<VarianceLayout::Unpooled>
{
    using type = Eigen::VectorXd;
};

template <VarianceLayout Kind>
using marker_variance_state_t = typename MarkerVarianceStateType<Kind>::type;

}  // namespace detail

template <VarianceLayout Kind>
struct GaussianState
{
    detail::marker_variance_state_t<Kind> variance{};
    Eigen::VectorXd fitted_values;  // total
};

template <HalfNormalAsymmetry Axis>
struct HalfNormalState;

template <>
struct HalfNormalState<HalfNormalAsymmetry::Count>
{
    double variance{};
    Eigen::VectorX<std::uint8_t> assignment;
    double positive_probability{};
    Eigen::VectorXd fitted_values;  // total
};

template <>
struct HalfNormalState<HalfNormalAsymmetry::Magnitude>
{
    std::array<double, 2> variances{};  // neg, pos
    Eigen::VectorX<std::uint8_t> assignment;
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

template <std::size_t ClassCount>
    requires(
        ClassCount > 1
        && ClassCount <= static_cast<std::size_t>(
                             std::numeric_limits<std::uint8_t>::max())
                             + 1)
struct JointSpikeSlabState
{
    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t component_count = class_count - 1;

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

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_STATE_H_
