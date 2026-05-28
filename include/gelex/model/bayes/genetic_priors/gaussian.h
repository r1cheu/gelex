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

#include <array>
#include <memory>
#include <string_view>

#include <Eigen/Core>

#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class SingleSharedGaussianPrior final
    : public SingleGeneticPrior
    , public SingleSharedMarkerVarianceCap
{
   public:
    static constexpr std::string_view name = "shared_gaussian";

    SingleSharedGaussianPrior(GeneticMode mode, SharedMarkerVariance variance);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> SharedMarkerVariance& override { return variance_; }
    auto variance() const -> const SharedMarkerVariance& override
    {
        return variance_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    SharedMarkerVariance variance_;
};

class SinglePerMarkerGaussianPrior final
    : public SingleGeneticPrior
    , public SinglePerMarkerVarianceCap
{
   public:
    static constexpr std::string_view name = "per_marker_gaussian";

    SinglePerMarkerGaussianPrior(GeneticMode mode, PerMarkerVariance variance);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> PerMarkerVariance& override { return variance_; }
    auto variance() const -> const PerMarkerVariance& override
    {
        return variance_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    PerMarkerVariance variance_;
};

class SingleSharedSpikeSlabGaussianPrior final
    : public SingleGeneticPrior
    , public SingleSharedMarkerVarianceCap
    , public SingleMixtureProportionCap
{
   public:
    static constexpr std::string_view name = "shared_spike_slab_gaussian";

    SingleSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> SharedMarkerVariance& override { return variance_; }
    auto variance() const -> const SharedMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> MixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const MixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    SharedMarkerVariance variance_;
    MixtureProportion mixture_proportion_;
};

class SinglePerMarkerSpikeSlabGaussianPrior final
    : public SingleGeneticPrior
    , public SinglePerMarkerVarianceCap
    , public SingleMixtureProportionCap
{
   public:
    static constexpr std::string_view name = "per_marker_spike_slab_gaussian";

    SinglePerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> PerMarkerVariance& override { return variance_; }
    auto variance() const -> const PerMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> MixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const MixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    PerMarkerVariance variance_;
    MixtureProportion mixture_proportion_;
};

class SingleScaledMixtureGaussianPrior final
    : public SingleGeneticPrior
    , public SingleSharedMarkerVarianceCap
    , public SingleMixtureProportionCap
    , public SingleMultiplierCap
{
   public:
    static constexpr std::string_view name = "scaled_mixture_gaussian";

    SingleScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        MixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> SharedMarkerVariance& override { return variance_; }
    auto variance() const -> const SharedMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> MixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const MixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto multiplier() -> Eigen::VectorXd& override { return multiplier_; }
    auto multiplier() const -> const Eigen::VectorXd& override
    {
        return multiplier_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    SharedMarkerVariance variance_;
    Eigen::VectorXd multiplier_;
    MixtureProportion mixture_proportion_;
};

class JointGaussianMixturePrior final
    : public JointGeneticPrior
    , public JointMarkerVarianceCap
    , public JointMixtureProportionCap
{
   public:
    static constexpr std::string_view name = "joint_mixture_gaussian";

    JointGaussianMixturePrior(
        std::array<JointMarkerVariance, 2> variances,
        MixtureProportion proportion);

    auto variance(GeneticMode mode) -> JointMarkerVariance& override;
    auto variance(GeneticMode mode) const
        -> const JointMarkerVariance& override;
    auto proportion() -> MixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const MixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<JointGeneticPriorState> override;

   private:
    std::array<JointMarkerVariance, 2> marker_variances_;
    MixtureProportion mixture_proportion_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIORS_GAUSSIAN_H_
