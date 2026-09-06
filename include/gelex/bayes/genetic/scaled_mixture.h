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

#ifndef GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
#define GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ranges>
#include <variant>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/genetic_spec.h"
#include "gelex/bayes/genetic/detail/draws_support.h"
#include "gelex/bayes/genetic/detail/fitted_update.h"
#include "gelex/bayes/genetic/detail/pip_support.h"
#include "gelex/bayes/genetic/detail/prior_support.h"
#include "gelex/bayes/genetic/detail/result_support.h"
#include "gelex/bayes/genetic/detail/state_support.h"
#include "gelex/bayes/genetic/detail/summary_support.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/genotype/operations.h"
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
struct ScaledMixtureFamily
{
};

template <
    std::size_t ClassCount,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct ScaledMixturePrior
{
    static constexpr std::size_t class_count = ClassCount;

    VarianceParameter variance;
    SimplexParameter<class_count, WeightUpdate> probabilities;
    std::array<double, class_count> scales{};
};

template <std::size_t ClassCount>
    requires(
        ClassCount > 1
        && ClassCount <= static_cast<std::size_t>(
                             std::numeric_limits<std::uint8_t>::max())
                             + 1)
class ScaledMixtureState;

namespace detail
{
template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const ScaledMixturePrior<ClassCount, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> ScaledMixtureState<ClassCount>;
}  // namespace detail

template <std::size_t ClassCount>
    requires(
        ClassCount > 1
        && ClassCount <= static_cast<std::size_t>(
                             std::numeric_limits<std::uint8_t>::max())
                             + 1)
class ScaledMixtureState
{
   public:
    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t component_count = class_count - 1;

    auto coefficients() const -> const Eigen::VectorXd&
    {
        return coefficients_;
    }
    auto assignments() const -> const Eigen::VectorX<std::uint8_t>&
    {
        return assignments_;
    }
    auto class_counts() const -> const std::array<std::size_t, class_count>&
    {
        return class_counts_;
    }
    auto fitted_values() const -> const Eigen::
        Matrix<double, Eigen::Dynamic, static_cast<int>(component_count)>&
    {
        return fitted_values_;
    }
    auto variance() const -> double { return variance_; }
    auto variance() -> double& { return variance_; }
    auto probabilities() const -> const std::array<double, class_count>&
    {
        return probabilities_;
    }
    auto probabilities() -> std::array<double, class_count>&
    {
        return probabilities_;
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] auto transition(
        Eigen::Index marker_index,
        double coefficient,
        std::uint8_t assignment)
        -> std::variant<
            std::monostate,
            bayes::AxpyTarget,
            std::array<bayes::AxpyTarget, 2>>
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        const double old_value = coefficients_(marker_index);
        const std::uint8_t old_assignment = assignments_(marker_index);
        const double new_value = assignment == 0 ? 0.0 : coefficient;
        if (old_assignment != assignment)
        {
            --class_counts_[old_assignment];
            ++class_counts_[assignment];
        }
        coefficients_(marker_index) = new_value;
        assignments_(marker_index) = assignment;

        return detail::make_fitted_update(
            fitted_values_,
            static_cast<Eigen::Index>(old_assignment) - 1,
            static_cast<Eigen::Index>(assignment) - 1,
            old_value,
            new_value);
    }

   private:
    template <std::size_t Count, MixtureWeightUpdate WeightUpdate>
    friend auto detail::make_state(
        const ScaledMixturePrior<Count, WeightUpdate>& prior,
        detail::GeneticStateDimensions dimensions) -> ScaledMixtureState<Count>;

    ScaledMixtureState(
        double variance,
        std::array<double, class_count> probabilities,
        Eigen::Index num_markers,
        Eigen::Index num_individuals)
        : coefficients_(Eigen::VectorXd::Zero(num_markers)),
          assignments_(Eigen::VectorX<std::uint8_t>::Zero(num_markers)),
          class_counts_{static_cast<std::size_t>(num_markers)},
          fitted_values_(
              Eigen::Matrix<
                  double,
                  Eigen::Dynamic,
                  static_cast<int>(component_count)>::
                  Zero(num_individuals, component_count)),
          variance_(variance),
          probabilities_(probabilities)
    {
    }

    Eigen::VectorXd coefficients_;
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, class_count> class_counts_;
    Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(component_count)>
        fitted_values_;
    double variance_;
    std::array<double, class_count> probabilities_;
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct ScaledMixtureDraws
{
    ScalarDraw variance;
    CategoryDraw<ClassCount> assignment;
    detail::weight_draw_t<WeightUpdate, VectorDraw> probabilities;
    VectorDraw component_explained_variance;

    auto append(const ScaledMixtureState<ClassCount>& state) -> void
    {
        variance.append(state.variance());
        assignment.append(state.assignments());
        probabilities.append(state.probabilities());
        component_explained_variance.append(
            matvar<0>(state.fitted_values(), VarNormType::Population));
    }
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct ScaledMixtureResult
{
    ScalarResult variance;
    detail::weight_result_t<WeightUpdate, VectorResult> probabilities;
    VectorResult component_explained_variance;
};

}  // namespace gelex

namespace gelex::detail
{

template <GeneticModeSet Modes, MixtureWeightUpdate WeightUpdate>
struct GeneticSpecFor<Modes, ScaledMixtureFamily<WeightUpdate>>
{
    using type = HomogeneousModeValues<Modes, ScaledMixture>;
};

constexpr auto initial_activity(const ScaledMixture& spec) -> double
{
    auto activity = 0.0;
    for (const auto [probability, scale] :
         std::views::zip(spec.probabilities(), spec.scales()))
    {
        activity += probability * scale;
    }
    return activity;
}

template <GeneticModeSet Modes, MixtureWeightUpdate WeightUpdate>
auto make_prior(
    ScaledMixtureFamily<WeightUpdate> /*family*/,
    const genetic_spec_t<Modes, ScaledMixtureFamily<WeightUpdate>>&
        genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    return transform_mode_values(
        genetic_spec,
        [&]<GeneticMode Mode>(const ScaledMixture& spec)
            -> ScaledMixturePrior<ScaledMixture::class_count, WeightUpdate>
        {
            return {
                .variance = calibrator.calibrate(Mode, initial_activity(spec)),
                .probabilities = make_parameter<WeightUpdate>(
                    spec.probabilities(),
                    make_uniform_dirichlet_prior<ScaledMixture::class_count>()),
                .scales = spec.scales()};
        });
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const ScaledMixturePrior<ClassCount, WeightUpdate>& prior,
    GeneticStateDimensions dimensions) -> ScaledMixtureState<ClassCount>
{
    return {
        prior.variance.initial,
        prior.probabilities.initial,
        dimensions.marker_count,
        dimensions.individual_count};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const ScaledMixturePrior<ClassCount, WeightUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder)
    -> ScaledMixtureDraws<ClassCount, WeightUpdate>
{
    return {
        .variance = builder.scalar("variance"),
        .assignment
        = builder.category<ClassCount>("assignment", builder.marker_count()),
        .probabilities
        = make_probabilities_draw<WeightUpdate, ClassCount>(builder),
        .component_explained_variance = make_component_explained_variance_draw<
            ScaledMixtureState<ClassCount>::component_count>(builder)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_result(const ScaledMixtureDraws<ClassCount, WeightUpdate>& draws)
    -> ScaledMixtureResult<ClassCount, WeightUpdate>
{
    return {
        .variance = make_result(draws.variance),
        .probabilities = make_result(draws.probabilities),
        .component_explained_variance
        = make_result(draws.component_explained_variance)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_pip(
    const ScaledMixtureDraws<ClassCount, WeightUpdate>& draws)
    -> MarkerPipResult
{
    return MarkerPipResult{
        draws.assignment.probability_of(is_non_null_category)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto write_family_summary_rows(
    TextWriter& writer,
    const ScaledMixtureResult<ClassCount, WeightUpdate>& result) -> void
{
    write_summary_rows(writer, result.variance);
    write_summary_rows(writer, result.probabilities);
    write_summary_rows(writer, result.component_explained_variance);
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
