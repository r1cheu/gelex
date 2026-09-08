// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
#define GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <ranges>
#include <string_view>
#include <utility>
#include <variant>

#include "gelex/bayes/genetic/detail/fitted_update.h"
#include "gelex/bayes/genetic/detail/marker_variance.h"
#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/operations.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct ScaledMixturePrior
{
    static constexpr std::size_t class_count = ScaledMixtureSpec<>::class_count;

    VarianceParameter variance;
    SimplexParameter<class_count, WeightUpdate> probabilities;
    std::array<double, class_count> scales{};
};

GELEX_NAMESPACE_BEGIN(detail)
template <MixtureWeightUpdate WeightUpdate>
constexpr auto initial_activity(const ScaledMixtureSpec<WeightUpdate>& spec)
    -> double
{
    auto activity = 0.0;
    for (const auto [probability, scale] :
         std::views::zip(spec.probabilities(), spec.scales()))
    {
        activity += probability * scale;
    }
    return activity;
}

GELEX_NAMESPACE_END(detail)

template <GeneticMode Mode, MixtureWeightUpdate WeightUpdate>
auto make_prior(
    const ScaledMixtureSpec<WeightUpdate>& spec,
    const MarkerVarianceCalibrator& calibrator)
    -> ScaledMixturePrior<WeightUpdate>
{
    return {
        .variance = calibrator.calibrate(Mode, detail::initial_activity(spec)),
        .probabilities = detail::make_parameter<WeightUpdate>(
            spec.probabilities(),
            make_uniform_dirichlet_prior<ScaledMixtureSpec<>::class_count>()),
        .scales = spec.scales()};
}

template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
class ScaledMixtureState
{
   public:
    static constexpr std::size_t class_count = ScaledMixtureSpec<>::class_count;
    static constexpr std::size_t component_count = class_count - 1;

    using fitted_values_type = Eigen::
        Matrix<double, Eigen::Dynamic, static_cast<int>(component_count)>;

    ScaledMixtureState(
        double variance,
        std::array<double, class_count> probabilities,
        GeneticDimensions dimensions)
        : coefficients_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          assignments_(
              Eigen::VectorX<std::uint8_t>::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          class_counts_{dimensions.marker},
          fitted_values_(
              fitted_values_type::Zero(
                  static_cast<Eigen::Index>(dimensions.individual),
                  component_count)),
          variance_(variance),
          probabilities_(probabilities)
    {
        detail::validate_probability_simplex(
            probabilities_, "scaled-mixture probabilities");
    }

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
    auto fitted_values() const -> const fitted_values_type&
    {
        return fitted_values_;
    }
    auto variance() const -> double { return variance_; }
    auto variance() -> double& { return variance_; }
    auto probabilities() const -> const std::array<double, class_count>&
    {
        return probabilities_;
    }
    auto set_probabilities(std::array<double, class_count> probabilities)
        -> void
        requires(WeightUpdate == MixtureWeightUpdate::Enabled)
    {
        detail::validate_probability_simplex(
            probabilities, "scaled-mixture probabilities");
        probabilities_ = probabilities;
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
    Eigen::VectorXd coefficients_;
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, class_count> class_counts_;
    fitted_values_type fitted_values_;
    double variance_;
    std::array<double, class_count> probabilities_;
};

template <MixtureWeightUpdate WeightUpdate>
auto make_state(
    const ScaledMixturePrior<WeightUpdate>& prior,
    GeneticDimensions dimensions) -> ScaledMixtureState<WeightUpdate>
{
    return {prior.variance.initial, prior.probabilities.initial, dimensions};
}

template <MixtureWeightUpdate WeightUpdate>
class ScaledMixtureDraws
{
   public:
    using probability_writer_type = probability_writer_t<WeightUpdate>;

    explicit ScaledMixtureDraws(
        PayloadWriter<double> variance,
        PayloadWriter<float> coefficients,
        PayloadWriter<std::uint8_t> assignments,
        probability_writer_type probabilities,
        PayloadWriter<double> component_explained_variance)
        : variance_{std::move(variance)},
          coefficients_{std::move(coefficients)},
          assignments_{std::move(assignments)},
          probabilities_{std::move(probabilities)},
          component_explained_variance_{std::move(component_explained_variance)}
    {
    }

    auto append(const ScaledMixtureState<WeightUpdate>& state) -> void
    {
        variance_.append(state.variance());
        coefficients_.append(
            state.coefficients().template cast<float>().eval());
        assignments_.append(state.assignments());
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            probabilities_.append(state.probabilities());
        }
        component_explained_variance_.append(
            matvar<0>(state.fitted_values(), VarNormType::Population));
    }

   private:
    PayloadWriter<double> variance_;
    PayloadWriter<float> coefficients_;
    PayloadWriter<std::uint8_t> assignments_;
    [[no_unique_address]] probability_writer_type probabilities_;
    PayloadWriter<double> component_explained_variance_;
};

template <MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const ScaledMixtureState<WeightUpdate>& state,
    BinaryWriter& writer,
    std::string_view prefix,
    std::size_t draw_count) -> ScaledMixtureDraws<WeightUpdate>
{
    const auto marker_count
        = static_cast<std::size_t>(state.coefficients().size());
    auto variance = writer.reserve<double>(
        fmt::format("{}/{}", prefix, variance_id), BinaryShape{1, draw_count});
    auto coefficients = writer.reserve<float>(
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{marker_count, draw_count});
    auto assignments = writer.reserve<std::uint8_t>(
        fmt::format("{}/{}", prefix, assignment_id),
        BinaryShape{marker_count, draw_count});
    auto probabilities = [&]() -> probability_writer_t<WeightUpdate>
    {
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            return writer.reserve<double>(
                fmt::format("{}/{}", prefix, probabilities_id),
                BinaryShape{
                    ScaledMixtureState<WeightUpdate>::class_count, draw_count});
        }
        else
        {
            return {};
        }
    }();
    auto component_explained_variance = writer.reserve<double>(
        fmt::format("{}/{}", prefix, component_explained_variance_id),
        BinaryShape{
            ScaledMixtureState<WeightUpdate>::component_count, draw_count});

    return ScaledMixtureDraws<WeightUpdate>{
        std::move(variance),
        std::move(coefficients),
        std::move(assignments),
        std::move(probabilities),
        std::move(component_explained_variance)};
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
