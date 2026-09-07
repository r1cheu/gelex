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

#include "gelex/bayes/genetic/detail/fitted_update.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"
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

    using JointPrior = JointSpikeSlabPrior<WeightUpdate>;
    return JointModeValues{
        std::move(mode_priors),
        JointPrior{
            .probabilities = make_parameter<WeightUpdate>(
                joint_spec.probabilities(),
                make_uniform_dirichlet_prior<JointPrior::class_count>())}};
}
GELEX_NAMESPACE_END(detail)

class HalfNormalState
{
   public:
    HalfNormalState(double variance, GeneticDimensions dimensions)
        : coefficients_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          fitted_values_(
              Eigen::VectorXd::Zero(
                  static_cast<Eigen::Index>(dimensions.individual))),
          variance_(variance),
          annotation_coefficients_(Eigen::Vector2d::Zero())
    {
    }

    auto coefficients() const -> const Eigen::VectorXd&
    {
        return coefficients_;
    }
    auto fitted_values() const -> const Eigen::VectorXd&
    {
        return fitted_values_;
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
        coefficients_(marker) = coefficient;
    }
    auto transition(const Eigen::Ref<const Eigen::VectorXd>& delta) -> void
    {
        fitted_values_.noalias() += delta;
    }

   private:
    Eigen::VectorXd coefficients_;
    Eigen::VectorXd fitted_values_;
    double variance_;
    Eigen::Vector2d annotation_coefficients_;
};

// Classes are NULL, A-only, D-only and AD; fitted_values holds one column per
// (mode, class) cell in which that mode is active, so every column carries a
// single mode and the two columns of a mode sum to that mode's total.
class JointSpikeSlabState
{
   public:
    static constexpr std::size_t class_count
        = JointSpikeSlabSpec<>::class_count;
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

    using FittedValues = Eigen::
        Matrix<double, Eigen::Dynamic, static_cast<int>(component_count)>;
    using ModeCoefficients
        = HomogeneousModeValues<GeneticMode::A | GeneticMode::D, double>;

    JointSpikeSlabState(
        std::array<double, class_count> probabilities,
        GeneticDimensions dimensions)
        : assignments_(
              Eigen::VectorX<std::uint8_t>::Zero(
                  static_cast<Eigen::Index>(dimensions.marker))),
          class_counts_{dimensions.marker},
          probabilities_(probabilities),
          fitted_values_(
              FittedValues::Zero(
                  static_cast<Eigen::Index>(dimensions.individual),
                  component_count))
    {
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
    auto probabilities() -> std::array<double, class_count>&
    {
        return probabilities_;
    }
    auto fitted_values() const -> const FittedValues& { return fitted_values_; }

    template <GeneticMode Mode>
        requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
    [[nodiscard]] static constexpr auto fitted_component_index(
        std::size_t class_index) noexcept -> int
    {
        assert(class_index < class_count);
        if constexpr (Mode == GeneticMode::A)
        {
            return additive_components[class_index];
        }
        else
        {
            return dominance_components[class_index];
        }
    }

    [[nodiscard]] auto transition(
        Eigen::Index marker,
        std::uint8_t assignment,
        const ModeCoefficients& old_coefficients,
        const ModeCoefficients& new_coefficients)
    {
        const auto old_assignment = assignments_(marker);
        auto updates = generate_mode_values<GeneticMode::A | GeneticMode::D>(
            [&]<GeneticMode Mode>()
            {
                return detail::make_fitted_update(
                    fitted_values_,
                    fitted_component_index<Mode>(old_assignment),
                    fitted_component_index<Mode>(assignment),
                    old_coefficients.template get<Mode>(),
                    new_coefficients.template get<Mode>());
            });
        if (old_assignment != assignment)
        {
            --class_counts_[old_assignment];
            ++class_counts_[assignment];
        }
        assignments_(marker) = assignment;
        return updates;
    }

   private:
    Eigen::VectorX<std::uint8_t> assignments_;
    std::array<std::size_t, class_count> class_counts_;
    std::array<double, class_count> probabilities_;
    FittedValues fitted_values_;
};

GELEX_NAMESPACE_BEGIN(detail)
inline auto make_state(
    const HalfNormalPrior& prior,
    GeneticDimensions dimensions) -> HalfNormalState
{
    return {prior.variance.initial, dimensions};
}

template <MixtureWeightUpdate WeightUpdate>
auto make_state(
    const JointSpikeSlabPrior<WeightUpdate>& prior,
    GeneticDimensions dimensions) -> JointSpikeSlabState
{
    return {prior.probabilities.initial, dimensions};
}
GELEX_NAMESPACE_END(detail)

class HalfNormalDraws
{
   public:
    explicit HalfNormalDraws(
        PayloadWriter<double> variance,
        PayloadWriter<float> coefficients,
        PayloadWriter<float> annotation_coefficients)
        : variance_{std::move(variance)},
          coefficients_{std::move(coefficients)},
          annotation_coefficients_{std::move(annotation_coefficients)}
    {
    }

    auto append(const HalfNormalState& state) -> void
    {
        variance_.append(state.variance());
        coefficients_.append(state.coefficients().cast<float>().eval());
        annotation_coefficients_.append(
            state.annotation_coefficients().cast<float>().eval());
    }

   private:
    PayloadWriter<double> variance_;
    PayloadWriter<float> coefficients_;
    PayloadWriter<float> annotation_coefficients_;
};

template <MixtureWeightUpdate WeightUpdate>
class JointSpikeSlabDraws
{
   public:
    using probability_writer_type = probability_writer_t<WeightUpdate>;

    explicit JointSpikeSlabDraws(
        PayloadWriter<std::uint8_t> assignments,
        probability_writer_type probabilities,
        PayloadWriter<double> component_explained_variance)
        : assignments_{std::move(assignments)},
          probabilities_{std::move(probabilities)},
          component_explained_variance_{std::move(component_explained_variance)}
    {
    }

    auto append(const JointSpikeSlabState& state) -> void
    {
        assignments_.append(state.assignments());
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            probabilities_.append(state.probabilities());
        }
        component_explained_variance_.append(
            matvar<0>(state.fitted_values(), VarNormType::Population));
    }

   private:
    PayloadWriter<std::uint8_t> assignments_;
    [[no_unique_address]] probability_writer_type probabilities_;
    PayloadWriter<double> component_explained_variance_;
};

GELEX_NAMESPACE_BEGIN(detail)
[[nodiscard]] inline auto make_draws(
    const HalfNormalPrior& /*prior*/,
    BinaryWriter& writer,
    std::string_view prefix,
    std::size_t draw_count,
    GeneticDimensions dimensions) -> HalfNormalDraws
{
    auto variance = writer.reserve<double>(
        fmt::format("{}/{}", prefix, variance_id), BinaryShape{1, draw_count});
    auto coefficients = writer.reserve<float>(
        fmt::format("{}/{}", prefix, coefficients_id),
        BinaryShape{dimensions.marker, draw_count});
    auto annotation_coefficients = writer.reserve<float>(
        fmt::format("{}/{}", prefix, annotation_coefficients_id),
        BinaryShape{2, draw_count});

    return HalfNormalDraws{
        std::move(variance),
        std::move(coefficients),
        std::move(annotation_coefficients)};
}

// Rows follow the fitted column layout of JointSpikeSlabState: A in A-only,
// A in AD, D in D-only, D in AD.
template <MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_draws(
    const JointSpikeSlabPrior<WeightUpdate>& /*prior*/,
    BinaryWriter& writer,
    std::string_view prefix,
    std::size_t draw_count,
    GeneticDimensions dimensions) -> JointSpikeSlabDraws<WeightUpdate>
{
    auto assignments = writer.reserve<std::uint8_t>(
        fmt::format("{}/{}", prefix, assignment_id),
        BinaryShape{dimensions.marker, draw_count});
    auto probabilities = [&]() -> probability_writer_t<WeightUpdate>
    {
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            return writer.reserve<double>(
                fmt::format("{}/{}", prefix, probabilities_id),
                BinaryShape{JointSpikeSlabState::class_count, draw_count});
        }
        else
        {
            return {};
        }
    }();
    auto component_explained_variance = writer.reserve<double>(
        fmt::format("{}/{}", prefix, component_explained_variance_id),
        BinaryShape{JointSpikeSlabState::component_count, draw_count});

    return JointSpikeSlabDraws<WeightUpdate>{
        std::move(assignments),
        std::move(probabilities),
        std::move(component_explained_variance)};
}
GELEX_NAMESPACE_END(detail)

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_H_
