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

#include "gelex/bayes/genetic/component_layout.h"
#include "gelex/bayes/genetic/prior.h"
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
};

template <VarianceLayout Kind>
struct SpikeSlabState
{
    detail::marker_variance_state_t<Kind> variance{};
    Eigen::VectorX<std::uint8_t> assignment;
    double probability{};
};

struct ScaledMixtureState
{
    double variance{};
    Eigen::VectorX<std::uint8_t> assignment;
    std::array<double, ScaledMixturePrior<>::class_count> probabilities{};
};

struct JointSpikeSlabState
{
    Eigen::VectorX<std::uint8_t> assignment;
    std::array<double, JointSpikeSlabPrior<>::class_count> probabilities{};

    Eigen::VectorX<std::uint8_t> dominance_sign;
    double positive_probability{};
};

template <typename FamilyState, ComponentLayout Layout>
struct GeneticModeState
{
    static_assert(
        Layout::class_count
        <= static_cast<std::size_t>(std::numeric_limits<std::uint8_t>::max())
               + 1);

    using component_layout = Layout;
    using ComponentFittedMatrix = Eigen::Matrix<
        double,
        Eigen::Dynamic,
        static_cast<int>(Layout::component_count)>;

    Eigen::VectorXd coefficients;
    ComponentFittedMatrix component_fitted_values;
    FamilyState family_state;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_STATE_H_
