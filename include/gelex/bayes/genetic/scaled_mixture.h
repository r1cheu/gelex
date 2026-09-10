// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
#define GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <ranges>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/detail/marker_variance.h"
#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/dense_writer.h"
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

    auto transition(
        Eigen::Index marker,
        double coefficient,
        std::uint8_t assignment) -> void
    {
        assert(marker >= 0 && marker < coefficients_.size());
        assert(assignment < class_count);
        const auto old_assignment = assignments_(marker);
        if (old_assignment != assignment)
        {
            --class_counts_[old_assignment];
            ++class_counts_[assignment];
        }
        coefficients_(marker) = assignment == 0 ? 0.0 : coefficient;
        assignments_(marker) = assignment;
    }

   private:
    Eigen::VectorXd coefficients_;
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, class_count> class_counts_;
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
        DenseStream<double> variance,
        CscStream<double> coefficients,
        assignment_writer_t assignments,
        probability_writer_type probabilities)
        : variance_{std::move(variance)},
          coefficients_{std::move(coefficients)},
          assignments_{std::move(assignments)},
          probabilities_{std::move(probabilities)}
    {
    }

    auto operator<<(const ScaledMixtureState<WeightUpdate>& state)
        -> ScaledMixtureDraws&
    {
        variance_ << state.variance();
        coefficients_ << state.coefficients();
        assignments_ << state.assignments();
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            probabilities_ << state.probabilities();
        }
        return *this;
    }

   private:
    DenseStream<double> variance_;
    CscStream<double> coefficients_;
    assignment_writer_t assignments_;
    [[no_unique_address]] probability_writer_type probabilities_;
};

template <MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const ScaledMixtureState<WeightUpdate>& state,
    DrawWriters writers,
    std::string_view prefix,
    std::size_t draw_count) -> ScaledMixtureDraws<WeightUpdate>
{
    const auto marker_count
        = static_cast<std::size_t>(state.coefficients().size());
    auto variance = writers.dense.reserve<double>(
        fmt::format("{}/{}", prefix, variance_id), BinaryShape{1, draw_count});
    auto coefficients = writers.sparse.reserve<double>(
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{marker_count, draw_count});
    auto assignments = writers.sparse.reserve<std::uint8_t>(
        fmt::format("{}/{}", prefix, assignment_id),
        BinaryShape{marker_count, draw_count});
    auto probabilities = [&]() -> probability_writer_t<WeightUpdate>
    {
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            return writers.dense.reserve<double>(
                fmt::format("{}/{}", prefix, probabilities_id),
                BinaryShape{
                    ScaledMixtureState<WeightUpdate>::class_count, draw_count});
        }
        else
        {
            return {};
        }
    }();

    return ScaledMixtureDraws<WeightUpdate>{
        std::move(variance),
        std::move(coefficients),
        std::move(assignments),
        std::move(probabilities)};
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_SCALED_MIXTURE_H_
