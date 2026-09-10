// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_H_
#define GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_H_

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/mode_values.h"
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

// NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
struct HalfNormalPrior
{
    VarianceParameter variance;
};

template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
struct JointSpikeSlabPrior
{
    static constexpr std::size_t class_count
        = JointSpikeSlabSpec<>::class_count;

    SimplexParameter<class_count, WeightUpdate> probabilities;
};

GELEX_NAMESPACE_BEGIN(detail)
template <GeneticMode Mode, MixtureWeightUpdate WeightUpdate>
    requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
constexpr auto initial_activity(const JointSpikeSlabSpec<WeightUpdate>& spec)
    -> double
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

GELEX_NAMESPACE_END(detail)

template <GeneticModeSet Modes, MixtureWeightUpdate WeightUpdate>
    requires(Modes == (GeneticMode::A | GeneticMode::D))
auto make_prior(
    const JointModeValues<
        ModeValues<Modes, GaussianSpec<>, HalfNormalSpec>,
        JointSpikeSlabSpec<WeightUpdate>>& genetic_spec,
    const MarkerVarianceCalibrator& calibrator)
{
    const auto& joint_spec = genetic_spec.joint();
    auto mode_priors = transform_mode_values(
        genetic_spec.mode_values(),
        [&]<GeneticMode Mode>(const auto&)
        {
            auto variance = calibrator.calibrate(
                Mode, detail::initial_activity<Mode>(joint_spec));
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

    using JointPrior = JointSpikeSlabPrior<WeightUpdate>;
    return JointModeValues{
        std::move(mode_priors),
        JointPrior{
            .probabilities = detail::make_parameter<WeightUpdate>(
                joint_spec.probabilities(),
                make_uniform_dirichlet_prior<JointPrior::class_count>())}};
}

class HalfNormalState
{
   public:
    HalfNormalState(double variance, GeneticDimensions dimensions)
        : coefficients_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          variance_(variance),
          annotation_coefficients_(Eigen::Vector2d::Zero())
    {
    }

    auto coefficients() const -> const Eigen::VectorXd&
    {
        return coefficients_;
    }
    auto variance() const -> double { return variance_; }
    auto variance() -> double& { return variance_; }
    auto annotation_coefficients() const -> const Eigen::Vector2d&
    {
        return annotation_coefficients_;
    }
    auto annotation_coefficients() -> Eigen::Vector2d&
    {
        return annotation_coefficients_;
    }

    auto transition(Eigen::Index marker, double coefficient) -> void
    {
        assert(marker >= 0 && marker < coefficients_.size());
        coefficients_(marker) = coefficient;
    }

   private:
    Eigen::VectorXd coefficients_;
    double variance_;
    Eigen::Vector2d annotation_coefficients_;
};

// Classes are NULL, A-only, D-only and AD.
template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
class JointSpikeSlabState
{
   public:
    static constexpr std::size_t class_count
        = JointSpikeSlabSpec<>::class_count;

    JointSpikeSlabState(
        std::array<double, class_count> probabilities,
        GeneticDimensions dimensions)
        : assignments_(
              Eigen::VectorX<std::uint8_t>::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          class_counts_{dimensions.marker},
          probabilities_(probabilities)
    {
        detail::validate_probability_simplex(
            probabilities_, "joint spike-slab probabilities");
    }

    auto assignments() const -> const Eigen::VectorX<std::uint8_t>&
    {
        return assignments_;
    }
    auto class_counts() const -> const std::array<std::size_t, class_count>&
    {
        return class_counts_;
    }
    auto probabilities() const -> const std::array<double, class_count>&
    {
        return probabilities_;
    }
    auto set_probabilities(std::array<double, class_count> probabilities)
        -> void
        requires(WeightUpdate == MixtureWeightUpdate::Enabled)
    {
        detail::validate_probability_simplex(
            probabilities, "joint spike-slab probabilities");
        probabilities_ = probabilities;
    }

    auto transition(Eigen::Index marker, std::uint8_t assignment) -> void
    {
        assert(marker >= 0 && marker < assignments_.size());
        assert(assignment < class_count);
        const auto old_assignment = assignments_(marker);
        if (old_assignment != assignment)
        {
            --class_counts_[old_assignment];
            ++class_counts_[assignment];
        }
        assignments_(marker) = assignment;
    }

   private:
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, class_count> class_counts_;
    std::array<double, class_count> probabilities_;
};

inline auto make_state(
    const HalfNormalPrior& prior,
    GeneticDimensions dimensions) -> HalfNormalState
{
    return {prior.variance.initial, dimensions};
}

template <MixtureWeightUpdate WeightUpdate>
auto make_state(
    const JointSpikeSlabPrior<WeightUpdate>& prior,
    GeneticDimensions dimensions) -> JointSpikeSlabState<WeightUpdate>
{
    return {prior.probabilities.initial, dimensions};
}

template <CoefficientLayout Layout>
class HalfNormalDraws
{
   public:
    using coefficients_writer_type = coefficients_writer_t<Layout>;

    explicit HalfNormalDraws(
        DenseStream<double> variance,
        coefficients_writer_type coefficients,
        DenseStream<double> annotation_coefficients)
        : variance_{std::move(variance)},
          coefficients_{std::move(coefficients)},
          annotation_coefficients_{std::move(annotation_coefficients)}
    {
    }

    auto operator<<(const HalfNormalState& state) -> HalfNormalDraws&
    {
        variance_ << state.variance();
        coefficients_ << state.coefficients();
        annotation_coefficients_ << state.annotation_coefficients();
        return *this;
    }

   private:
    DenseStream<double> variance_;
    coefficients_writer_type coefficients_;
    DenseStream<double> annotation_coefficients_;
};

template <MixtureWeightUpdate WeightUpdate>
class JointSpikeSlabDraws
{
   public:
    using probability_writer_type = probability_writer_t<WeightUpdate>;

    explicit JointSpikeSlabDraws(
        assignment_writer_t assignments,
        probability_writer_type probabilities)
        : assignments_{std::move(assignments)},
          probabilities_{std::move(probabilities)}
    {
    }

    auto operator<<(const JointSpikeSlabState<WeightUpdate>& state)
        -> JointSpikeSlabDraws&
    {
        assignments_ << state.assignments();
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            probabilities_ << state.probabilities();
        }
        return *this;
    }

   private:
    assignment_writer_t assignments_;
    [[no_unique_address]] probability_writer_type probabilities_;
};

template <CoefficientLayout Layout = CoefficientLayout::Dense>
[[nodiscard]] auto make_draws(
    const HalfNormalState& state,
    DrawWriters writers,
    std::string_view prefix,
    std::size_t draw_count) -> HalfNormalDraws<Layout>
{
    const auto marker_count
        = static_cast<std::size_t>(state.coefficients().size());
    auto variance = writers.dense.reserve<double>(
        fmt::format("{}/{}", prefix, variance_id), BinaryShape{1, draw_count});
    auto coefficients = reserve_coefficients<Layout>(
        writers,
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{marker_count, draw_count});
    auto annotation_coefficients = writers.dense.reserve<double>(
        fmt::format("{}/{}", prefix, annotation_coefficients_id),
        BinaryShape{
            static_cast<std::size_t>(state.annotation_coefficients().size()),
            draw_count});

    return HalfNormalDraws<Layout>{
        std::move(variance),
        std::move(coefficients),
        std::move(annotation_coefficients)};
}

template <MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const JointSpikeSlabState<WeightUpdate>& state,
    DrawWriters writers,
    std::string_view prefix,
    std::size_t draw_count) -> JointSpikeSlabDraws<WeightUpdate>
{
    const auto marker_count
        = static_cast<std::size_t>(state.assignments().size());
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
                    JointSpikeSlabState<WeightUpdate>::class_count,
                    draw_count});
        }
        else
        {
            return {};
        }
    }();

    return JointSpikeSlabDraws<WeightUpdate>{
        std::move(assignments), std::move(probabilities)};
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_H_
