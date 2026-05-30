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

#ifndef GELEX_MODEL_BAYES_GENETIC_PRIORS_GAUSSIAN_H_
#define GELEX_MODEL_BAYES_GENETIC_PRIORS_GAUSSIAN_H_

#include <string_view>

#include <Eigen/Core>

#include "gelex/model/bayes/fields.h"
#include "gelex/model/bayes/genetic_prior_states/gaussian.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class SingleSharedGaussianPrior final : public Variance<SharedMarkerVariance>
{
   public:
    static constexpr std::string_view name = "shared_gaussian";

    SingleSharedGaussianPrior(GeneticMode mode, SharedMarkerVariance variance);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleSharedGaussianState;

   private:
    GeneticMode mode_;
};

class SinglePerMarkerGaussianPrior final : public Variance<PerMarkerVariance>
{
   public:
    static constexpr std::string_view name = "per_marker_gaussian";

    SinglePerMarkerGaussianPrior(GeneticMode mode, PerMarkerVariance variance);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SinglePerMarkerGaussianState;

   private:
    GeneticMode mode_;
};

class SingleFixedSharedSpikeSlabGaussianPrior final
    : public Variance<SharedMarkerVariance>
    , public Proportion<FixedMixtureProportion>
{
   public:
    static constexpr std::string_view name = "fixed_shared_spike_slab_gaussian";

    SingleFixedSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        FixedMixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleFixedSharedSpikeSlabGaussianState;

   private:
    GeneticMode mode_;
};

class SingleSampledSharedSpikeSlabGaussianPrior final
    : public Variance<SharedMarkerVariance>
    , public Proportion<SampledMixtureProportion>
{
   public:
    static constexpr std::string_view name
        = "sampled_shared_spike_slab_gaussian";

    SingleSampledSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        SampledMixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleSampledSharedSpikeSlabGaussianState;

   private:
    GeneticMode mode_;
};

class SingleFixedPerMarkerSpikeSlabGaussianPrior final
    : public Variance<PerMarkerVariance>
    , public Proportion<FixedMixtureProportion>
{
   public:
    static constexpr std::string_view name
        = "fixed_per_marker_spike_slab_gaussian";

    SingleFixedPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        FixedMixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleFixedPerMarkerSpikeSlabGaussianState;

   private:
    GeneticMode mode_;
};

class SingleSampledPerMarkerSpikeSlabGaussianPrior final
    : public Variance<PerMarkerVariance>
    , public Proportion<SampledMixtureProportion>
{
   public:
    static constexpr std::string_view name
        = "sampled_per_marker_spike_slab_gaussian";

    SingleSampledPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        SampledMixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleSampledPerMarkerSpikeSlabGaussianState;

   private:
    GeneticMode mode_;
};

class SingleFixedScaledMixtureGaussianPrior final
    : public Variance<SharedMarkerVariance>
    , public Proportion<FixedMixtureProportion>
    , public Multiplier<Eigen::VectorXd>
{
   public:
    static constexpr std::string_view name = "fixed_scaled_mixture_gaussian";

    SingleFixedScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        FixedMixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleFixedScaledMixtureGaussianState;

   private:
    GeneticMode mode_;
};

class SingleSampledScaledMixtureGaussianPrior final
    : public Variance<SharedMarkerVariance>
    , public Proportion<SampledMixtureProportion>
    , public Multiplier<Eigen::VectorXd>
{
   public:
    static constexpr std::string_view name = "sampled_scaled_mixture_gaussian";

    SingleSampledScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        SampledMixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleSampledScaledMixtureGaussianState;

   private:
    GeneticMode mode_;
};

class JointFixedGaussianMixturePrior final
    : public JointVariancesField<JointSharedMarkerVariance>
    , public Proportion<FixedMixtureProportion>
{
   public:
    static constexpr std::string_view name = "joint_fixed_mixture_gaussian";

    JointFixedGaussianMixturePrior(
        JointSharedMarkerVariance variance,
        FixedMixtureProportion proportion);

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> JointFixedGaussianMixtureState;
};

class JointSampledGaussianMixturePrior final
    : public JointVariancesField<JointSharedMarkerVariance>
    , public Proportion<SampledMixtureProportion>
{
   public:
    static constexpr std::string_view name = "joint_sampled_mixture_gaussian";

    JointSampledGaussianMixturePrior(
        JointSharedMarkerVariance variance,
        SampledMixtureProportion proportion);

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> JointSampledGaussianMixtureState;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIORS_GAUSSIAN_H_
