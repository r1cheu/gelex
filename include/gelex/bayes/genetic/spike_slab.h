// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_SPIKE_SLAB_H_
#define GELEX_BAYES_GENETIC_SPIKE_SLAB_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/detail/marker_variance.h"
#include "gelex/bayes/genetic/draw_schema.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct SpikeSlabPrior
{
    VarianceParameter variance;
    ProbabilityParameter<WeightUpdate> probability;
};

template <
    GeneticMode Mode,
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate>
auto make_prior(
    const SpikeSlabSpec<Kind, WeightUpdate>& spec,
    const MarkerVarianceCalibrator& calibrator)
    -> SpikeSlabPrior<Kind, WeightUpdate>
{
    return {
        .variance = calibrator.calibrate(Mode, spec.probability()),
        .probability = detail::make_parameter<WeightUpdate>(
            spec.probability(), make_beta_prior(1.0, 1.0))};
}

template <VarianceLayout Kind>
class SpikeSlabState
{
   public:
    using variance_type = detail::marker_variance_state_t<Kind>;
    SpikeSlabState(
        variance_type variance,
        double probability,
        GeneticDimensions dimensions)
        : coefficients_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          assignments_(
              Eigen::VectorX<std::uint8_t>::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          class_counts_{dimensions.marker, 0},
          fitted_values_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.individual))),
          variance_(std::move(variance)),
          probability_(probability)
    {
    }

    auto coefficients() const -> const Eigen::VectorXd&
    {
        return coefficients_;
    }
    auto assignments() const -> const Eigen::VectorX<std::uint8_t>&
    {
        return assignments_;
    }
    auto class_counts() const -> const std::array<std::size_t, 2>&
    {
        return class_counts_;
    }
    auto fitted_values() const -> const Eigen::VectorXd&
    {
        return fitted_values_;
    }
    auto variance() const -> const variance_type& { return variance_; }
    auto variance() -> variance_type& { return variance_; }
    auto probability() const -> double { return probability_; }
    auto probability() -> double& { return probability_; }

    auto transition(Eigen::Index marker, double coefficient, bool active)
        -> void
    {
        const auto old_assignment = assignments_(marker);
        const auto new_assignment = static_cast<std::uint8_t>(active);
        if (old_assignment != new_assignment)
        {
            --class_counts_[old_assignment];
            ++class_counts_[new_assignment];
        }
        assignments_(marker) = new_assignment;
        coefficients_(marker) = active ? coefficient : 0.0;
    }

    auto transition(const Eigen::Ref<const Eigen::VectorXd>& delta) -> void
    {
        fitted_values_.noalias() += delta;
    }

   private:
    Eigen::VectorXd coefficients_;
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, 2> class_counts_;
    Eigen::VectorXd fitted_values_;
    variance_type variance_;
    double probability_;
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const SpikeSlabPrior<Kind, WeightUpdate>& prior,
    GeneticDimensions dimensions) -> SpikeSlabState<Kind>
{
    return {
        detail::initial_marker_variance<Kind>(
            prior.variance, static_cast<Eigen::Index>(dimensions.marker)),
        prior.probability.initial,
        dimensions};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
class SpikeSlabDraws
{
   public:
    using variance_writer_type = marker_variance_writer_t<Kind>;
    using probability_writer_type = probability_writer_t<WeightUpdate>;

    explicit SpikeSlabDraws(
        variance_writer_type variances,
        PayloadWriter<float> coefficients,
        PayloadWriter<std::uint8_t> assignments,
        probability_writer_type probability)
        : variances_{std::move(variances)},
          coefficients_{std::move(coefficients)},
          assignments_{std::move(assignments)},
          probability_{std::move(probability)}
    {
    }

    auto append(const SpikeSlabState<Kind>& state) -> void
    {
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variances_.append(state.variance());
        }
        else
        {
            variances_.append(state.variance().template cast<float>().eval());
        }
        coefficients_.append(
            state.coefficients().template cast<float>().eval());
        assignments_.append(state.assignments());
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            probability_.append(state.probability());
        }
    }

   private:
    variance_writer_type variances_;
    PayloadWriter<float> coefficients_;
    PayloadWriter<std::uint8_t> assignments_;
    [[no_unique_address]] probability_writer_type probability_;
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const SpikeSlabPrior<Kind, WeightUpdate>& /*prior*/,
    BinaryWriter& writer,
    std::string_view prefix,
    std::size_t draw_count,
    GeneticDimensions dimensions) -> SpikeSlabDraws<Kind, WeightUpdate>
{
    const std::size_t variance_size
        = (Kind == VarianceLayout::Pooled) ? 1 : dimensions.marker;

    auto variances = writer.reserve<marker_variance_dtype_t<Kind>>(
        fmt::format("{}/{}", prefix, variance_id),
        BinaryShape{variance_size, draw_count});
    auto coefficients = writer.reserve<float>(
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{dimensions.marker, draw_count});
    auto assignments = writer.reserve<std::uint8_t>(
        fmt::format("{}/{}", prefix, assignment_id),
        BinaryShape{dimensions.marker, draw_count});
    auto probability = [&]() -> probability_writer_t<WeightUpdate>
    {
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            return writer.reserve<double>(
                fmt::format("{}/{}", prefix, probability_id),
                BinaryShape{1, draw_count});
        }
        else
        {
            return {};
        }
    }();

    return SpikeSlabDraws<Kind, WeightUpdate>{
        std::move(variances),
        std::move(coefficients),
        std::move(assignments),
        std::move(probability)};
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_SPIKE_SLAB_H_
