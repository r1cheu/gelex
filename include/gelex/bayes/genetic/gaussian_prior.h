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

#include <Eigen/Core>
#include <string_view>

#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes
{

class SingleSharedGaussianPrior final
{
   public:
    static constexpr std::string_view name = "shared_gaussian";

    SingleSharedGaussianPrior(GeneticMode mode, SharedMarkerVariance variance);

    auto mode() const -> GeneticMode { return mode_; }
    auto variance() -> SharedMarkerVariance& { return variance_; }
    auto variance() const -> const SharedMarkerVariance& { return variance_; }

    auto visit(FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleSharedGaussianState;

   private:
    SharedMarkerVariance variance_;
    GeneticMode mode_;
};

class SinglePerMarkerGaussianPrior final
{
   public:
    static constexpr std::string_view name = "per_marker_gaussian";

    SinglePerMarkerGaussianPrior(GeneticMode mode, PerMarkerVariance variance);

    auto mode() const -> GeneticMode { return mode_; }
    auto variance() -> PerMarkerVariance& { return variance_; }
    auto variance() const -> const PerMarkerVariance& { return variance_; }

    auto visit(FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SinglePerMarkerGaussianState;

   private:
    PerMarkerVariance variance_;
    GeneticMode mode_;
};

class SingleSharedSpikeSlabGaussianPrior final
{
   public:
    static constexpr std::string_view name = "shared_spike_slab_gaussian";

    SingleSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }
    auto variance() -> SharedMarkerVariance& { return variance_; }
    auto variance() const -> const SharedMarkerVariance& { return variance_; }
    auto proportion() -> MixtureProportion& { return proportion_; }
    auto proportion() const -> const MixtureProportion& { return proportion_; }

    auto visit(FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleSharedSpikeSlabGaussianState;

   private:
    SharedMarkerVariance variance_;
    MixtureProportion proportion_;
    GeneticMode mode_;
};

class SinglePerMarkerSpikeSlabGaussianPrior final
{
   public:
    static constexpr std::string_view name = "per_marker_spike_slab_gaussian";

    SinglePerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }
    auto variance() -> PerMarkerVariance& { return variance_; }
    auto variance() const -> const PerMarkerVariance& { return variance_; }
    auto proportion() -> MixtureProportion& { return proportion_; }
    auto proportion() const -> const MixtureProportion& { return proportion_; }

    auto visit(FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SinglePerMarkerSpikeSlabGaussianState;

   private:
    PerMarkerVariance variance_;
    MixtureProportion proportion_;
    GeneticMode mode_;
};

class SingleScaledMixtureGaussianPrior final
{
   public:
    static constexpr std::string_view name = "scaled_mixture_gaussian";

    SingleScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode { return mode_; }
    auto variance() -> SharedMarkerVariance& { return variance_; }
    auto variance() const -> const SharedMarkerVariance& { return variance_; }
    auto proportion() -> MixtureProportion& { return proportion_; }
    auto proportion() const -> const MixtureProportion& { return proportion_; }
    auto multiplier() -> Eigen::VectorXd& { return multiplier_; }
    auto multiplier() const -> const Eigen::VectorXd& { return multiplier_; }

    auto visit(FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> SingleScaledMixtureGaussianState;

   private:
    SharedMarkerVariance variance_;
    MixtureProportion proportion_;
    Eigen::VectorXd multiplier_;
    GeneticMode mode_;
};

class JointGaussianMixturePrior final
{
   public:
    static constexpr std::string_view name = "joint_mixture_gaussian";

    JointGaussianMixturePrior(
        JointSharedMarkerVariance variance,
        MixtureProportion proportion);

    auto variance(GeneticMode mode) -> SharedMarkerVariance&;
    auto variance(GeneticMode mode) const -> const SharedMarkerVariance&;
    auto variances() -> JointSharedMarkerVariance& { return variances_; }
    auto variances() const -> const JointSharedMarkerVariance&
    {
        return variances_;
    }
    auto proportion() -> MixtureProportion& { return proportion_; }
    auto proportion() const -> const MixtureProportion& { return proportion_; }

    auto visit(FieldVisitor& visitor) -> void;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> JointGaussianMixtureState;

   private:
    JointSharedMarkerVariance variances_;
    MixtureProportion proportion_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_H_
