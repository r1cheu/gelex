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

#ifndef GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_H_
#define GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/genetic_spec.h"
#include "gelex/bayes/genetic/detail/draws_support.h"
#include "gelex/bayes/genetic/detail/prior_support.h"
#include "gelex/bayes/genetic/detail/result_support.h"
#include "gelex/bayes/genetic/detail/state_support.h"
#include "gelex/bayes/genetic/detail/summary_support.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
struct JointSpikeSlabFamily
{
};

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct HalfNormalPrior
{
    VarianceParameter variance;
};

template <
    std::size_t ClassCount,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
struct JointSpikeSlabPrior
{
    static constexpr std::size_t class_count = ClassCount;

    SimplexParameter<class_count, WeightUpdate> probabilities;
};

struct HalfNormalState
{
    double variance{};
    Eigen::Vector2d probit_coefficients = Eigen::Vector2d::Zero();
    Eigen::VectorXd fitted_values;  // total
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

[[nodiscard]] inline auto genetic_value(const HalfNormalState& state)
    -> const Eigen::VectorXd&
{
    return state.fitted_values;
}

struct HalfNormalDraws
{
    ScalarDraw variance;
    VectorDraw probit_coefficients;

    auto append(const HalfNormalState& state) -> void
    {
        variance.append(state.variance);
        probit_coefficients.append(state.probit_coefficients);
    }
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct JointSpikeSlabDraws
{
    CategoryDraw<ClassCount> assignment;
    detail::weight_draw_t<WeightUpdate, VectorDraw> probabilities;
    VectorDraw component_explained_variance;

    auto append(const JointSpikeSlabState<ClassCount>& state) -> void
    {
        assignment.append(state.assignment);
        probabilities.append(state.probabilities);
        component_explained_variance.append(
            matvar<0>(state.fitted_values, VarNormType::Population));
    }
};

struct HalfNormalPosteriorResult
{
    ScalarPosteriorResult variance;
    VectorPosteriorResult probit_coefficients;
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct JointSpikeSlabPosteriorResult
{
    detail::weight_result_t<WeightUpdate, VectorPosteriorResult> probabilities;
    VectorPosteriorResult component_explained_variance;
};

}  // namespace gelex

namespace gelex::detail
{

template <MixtureWeightUpdate WeightUpdate>
struct GeneticSpecFor<
    GeneticMode::A | GeneticMode::D,
    JointSpikeSlabFamily<WeightUpdate>>
{
    using type = JointModeValues<
        ModeValues<GeneticMode::A | GeneticMode::D, Gaussian, HalfNormal>,
        JointSpikeSlab>;
};

template <GeneticMode Mode>
    requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
constexpr auto initial_activity(const JointSpikeSlab& spec) -> double
{
    const auto& probabilities = spec.probabilities();
    if constexpr (Mode == GeneticMode::A)
    {
        return probabilities.at(1) + probabilities.at(3);
    }
    else
    {
        return probabilities.at(2) + probabilities.at(3);
    }
}

template <GeneticModeSet Modes, MixtureWeightUpdate WeightUpdate>
    requires(Modes == (GeneticMode::A | GeneticMode::D))
auto make_prior(
    JointSpikeSlabFamily<WeightUpdate> /*family*/,
    const genetic_spec_t<Modes, JointSpikeSlabFamily<WeightUpdate>>&
        genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    const auto& joint_spec = genetic_spec.joint();
    auto mode_priors = transform_mode_values(
        genetic_spec.mode_values(),
        [&]<GeneticMode Mode>(const auto&)
        {
            auto variance = calibrator.calibrate(
                Mode, initial_activity<Mode>(joint_spec));
            if constexpr (Mode == GeneticMode::A)
            {
                return GaussianPrior<VarianceLayout::Pooled>{
                    .variance = std::move(variance)};
            }
            else
            {
                return HalfNormalPrior{.variance = std::move(variance)};
            }
        });

    using JointPrior
        = JointSpikeSlabPrior<JointSpikeSlab::class_count, WeightUpdate>;
    return JointModeValues{
        std::move(mode_priors),
        JointPrior{
            .probabilities = make_parameter<WeightUpdate>(
                joint_spec.probabilities(),
                make_uniform_dirichlet_prior<JointPrior::class_count>())}};
}

inline auto make_state(
    const HalfNormalPrior& prior,
    GeneticStateDimensions dimensions) -> HalfNormalState
{
    return {
        .variance = prior.variance.initial,
        .probit_coefficients = Eigen::Vector2d::Zero(),
        .fitted_values = Eigen::VectorXd::Zero(dimensions.individual_count)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const JointSpikeSlabPrior<ClassCount, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> JointSpikeSlabState<ClassCount>
{
    auto state = JointSpikeSlabState<ClassCount>{
        .assignment
        = Eigen::VectorX<std::uint8_t>::Zero(dimensions.marker_count),
        .probabilities = prior.probabilities.initial,
        .fitted_values = {},
    };
    state.fitted_values.setZero(
        dimensions.individual_count,
        static_cast<Eigen::Index>(state.component_count));
    return state;
}

[[nodiscard]] inline auto make_draws(
    const HalfNormalPrior& /*prior*/,
    GeneticDrawsBuilder& builder) -> HalfNormalDraws
{
    return {
        .variance = builder.scalar("variance"),
        .probit_coefficients = builder.vector("probit_coefficients", 2)};
}

// Rows follow the fitted column layout of JointSpikeSlabState: A in A-only,
// A in AD, D in D-only, D in AD.
template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const JointSpikeSlabPrior<ClassCount, WeightUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder)
    -> JointSpikeSlabDraws<ClassCount, WeightUpdate>
{
    return {
        .assignment
        = builder.category<ClassCount>("assignment", builder.marker_count()),
        .probabilities
        = make_probabilities_draw<WeightUpdate, ClassCount>(builder),
        .component_explained_variance = make_component_explained_variance_draw<
            JointSpikeSlabState<ClassCount>::component_count>(builder)};
}

inline auto make_result(const HalfNormalDraws& draws)
    -> HalfNormalPosteriorResult
{
    return {
        .variance = make_result(draws.variance),
        .probit_coefficients = make_result(draws.probit_coefficients)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_result(const JointSpikeSlabDraws<ClassCount, WeightUpdate>& draws)
    -> JointSpikeSlabPosteriorResult<ClassCount, WeightUpdate>
{
    return {
        .probabilities = make_result(draws.probabilities),
        .component_explained_variance
        = make_result(draws.component_explained_variance)};
}

inline auto write_family_summary_rows(
    TextWriter& writer,
    const HalfNormalPosteriorResult& result) -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probit_coefficients);
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto write_family_summary_rows(
    TextWriter& writer,
    const JointSpikeSlabPosteriorResult<ClassCount, WeightUpdate>& result)
    -> void
{
    write_summary_rows(writer, result.probabilities);
    write_summary_rows(writer, result.component_explained_variance);
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_H_
