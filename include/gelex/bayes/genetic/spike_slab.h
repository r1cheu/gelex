// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_SPIKE_SLAB_H_
#define GELEX_BAYES_GENETIC_SPIKE_SLAB_H_

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/detail/marker_variance.h"
#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/marker_effect_traits.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/calibration.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/dense_writer.h"
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

template <
    VarianceLayout Kind,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
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
          variance_(std::move(variance)),
          probability_(probability)
    {
        validate_probability(probability_);
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
    auto variance() const -> const variance_type& { return variance_; }
    auto variance() -> variance_type& { return variance_; }
    auto probability() const -> double { return probability_; }
    auto set_probability(double probability) -> void
        requires(WeightUpdate == MixtureWeightUpdate::Enabled)
    {
        validate_probability(probability);
        probability_ = probability;
    }

    auto transition(Eigen::Index marker, double coefficient, bool active)
        -> void
    {
        assert(marker >= 0 && marker < coefficients_.size());
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

   private:
    static auto validate_probability(double probability) -> void
    {
        if (!std::isfinite(probability) || probability <= 0.0
            || probability >= 1.0)
        {
            throw GelexException(
                fmt::format(
                    "spike-slab inclusion probability must lie in the open "
                    "interval (0, 1), got {}",
                    probability));
        }
    }

    Eigen::VectorXd coefficients_;
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, 2> class_counts_;
    variance_type variance_;
    double probability_;
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_state(
    const SpikeSlabPrior<Kind, WeightUpdate>& prior,
    GeneticDimensions dimensions) -> SpikeSlabState<Kind, WeightUpdate>
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
    static constexpr CoefficientLayout coefficient_layout
        = CoefficientLayout::Sparse;

    explicit SpikeSlabDraws(
        variance_writer_type variances,
        CscStream<double> coefficients,
        assignment_writer_t assignments,
        probability_writer_type probability)
        : variances_{std::move(variances)},
          coefficients_{std::move(coefficients)},
          assignments_{std::move(assignments)},
          probability_{std::move(probability)}
    {
    }

    auto operator<<(const SpikeSlabState<Kind, WeightUpdate>& state)
        -> SpikeSlabDraws&
    {
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variances_ << state.variance();
        }
        coefficients_ << state.coefficients();
        assignments_ << state.assignments();
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            probability_ << state.probability();
        }
        return *this;
    }

   private:
    [[no_unique_address]] variance_writer_type variances_;
    CscStream<double> coefficients_;
    assignment_writer_t assignments_;
    [[no_unique_address]] probability_writer_type probability_;
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const SpikeSlabState<Kind, WeightUpdate>& state,
    DrawWriters writers,
    std::string_view prefix,
    std::size_t draw_count) -> SpikeSlabDraws<Kind, WeightUpdate>
{
    const auto marker_count
        = static_cast<std::size_t>(state.coefficients().size());

    auto variances = reserve_marker_variance<Kind>(
        writers, fmt::format("{}/{}", prefix, variance_id), draw_count);
    auto coefficients = writers.sparse.reserve<double>(
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{marker_count, draw_count});
    auto assignments = writers.sparse.reserve<std::uint8_t>(
        fmt::format("{}/{}", prefix, assignment_id),
        BinaryShape{marker_count, draw_count});
    auto probability = [&]() -> probability_writer_t<WeightUpdate>
    {
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            return writers.dense.reserve<double>(
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

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
struct SpikeSlabDiagnostics
{
    [[no_unique_address]] marker_variance_diagnostics_t<Kind> variance;
    [[no_unique_address]] probability_diagnostics_t<WeightUpdate> probability;
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_diagnostics(
    std::type_identity<SpikeSlabDraws<Kind, WeightUpdate>> /*draws*/,
    DrawReaders readers,
    std::string_view prefix,
    double prob) -> SpikeSlabDiagnostics<Kind, WeightUpdate>
{
    return {
        .variance = diagnose_marker_variance<Kind>(
            readers, fmt::format("{}/{}", prefix, variance_id), prob),
        .probability = diagnose_probability<WeightUpdate>(
            readers, fmt::format("{}/{}", prefix, probability_id), prob)};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto append_diagnostic_entries(
    std::vector<DiagnosticEntry>& out,
    const SpikeSlabDiagnostics<Kind, WeightUpdate>& diagnostics,
    std::string_view prefix) -> void
{
    append_entry(
        out, fmt::format("{}/{}", prefix, variance_id), diagnostics.variance);
    append_entry(
        out,
        fmt::format("{}/{}", prefix, probability_id),
        diagnostics.probability);
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_marker_effects(
    std::type_identity<SpikeSlabDraws<Kind, WeightUpdate>> /*draws*/,
    DrawReaders readers,
    std::string_view prefix) -> MixtureMarkerEffects
{
    return {
        .coefficients = summarize_coefficients<CoefficientLayout::Sparse>(
            readers, fmt::format("{}/{}", prefix, coefficients_id)),
        .pip = inclusion_probability(
            readers.sparse,
            fmt::format("{}/{}", prefix, assignment_id),
            [](std::uint8_t assignment) { return assignment != 0; })};
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_SPIKE_SLAB_H_
