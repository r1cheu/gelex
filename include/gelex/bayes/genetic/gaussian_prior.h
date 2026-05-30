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

#ifndef GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_H_
#define GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_H_

#include <string_view>

#include <Eigen/Core>

#include "gelex/bayes/genetic/fields.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
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

class SingleSharedSpikeSlabGaussianPrior final
    : public Variance<SharedMarkerVariance>
    , public Proportion<MixtureProportion>
{
   public:
    static constexpr std::string_view name = "shared_spike_slab_gaussian";

    SingleSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleSharedSpikeSlabGaussianState;

   private:
    GeneticMode mode_;
};

class SinglePerMarkerSpikeSlabGaussianPrior final
    : public Variance<PerMarkerVariance>
    , public Proportion<MixtureProportion>
{
   public:
    static constexpr std::string_view name = "per_marker_spike_slab_gaussian";

    SinglePerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SinglePerMarkerSpikeSlabGaussianState;

   private:
    GeneticMode mode_;
};

class SingleScaledMixtureGaussianPrior final
    : public Variance<SharedMarkerVariance>
    , public Proportion<MixtureProportion>
    , public Multiplier<Eigen::VectorXd>
{
   public:
    static constexpr std::string_view name = "scaled_mixture_gaussian";

    SingleScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleScaledMixtureGaussianState;

   private:
    GeneticMode mode_;
};

class JointGaussianMixturePrior final
    : public JointVariancesField<JointSharedMarkerVariance>
    , public Proportion<MixtureProportion>
{
   public:
    static constexpr std::string_view name = "joint_mixture_gaussian";

    JointGaussianMixturePrior(
        JointSharedMarkerVariance variance,
        MixtureProportion proportion);

    auto visit(infra::FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> JointGaussianMixtureState;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_H_
