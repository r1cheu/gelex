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

class SingleFixedSharedSpikeSlabGaussianPrior final
    : public SingleGeneticPrior
    , public SingleSharedMarkerVarianceCap
    , public SingleFixedMixtureProportionCap
{
   public:
    static constexpr std::string_view name = "fixed_shared_spike_slab_gaussian";

    SingleFixedSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        FixedMixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> SharedMarkerVariance& override { return variance_; }
    auto variance() const -> const SharedMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> FixedMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const FixedMixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    SharedMarkerVariance variance_;
    FixedMixtureProportion mixture_proportion_;
};

class SingleSampledSharedSpikeSlabGaussianPrior final
    : public SingleGeneticPrior
    , public SingleSharedMarkerVarianceCap
    , public SingleSampledMixtureProportionCap
{
   public:
    static constexpr std::string_view name
        = "sampled_shared_spike_slab_gaussian";

    SingleSampledSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        SampledMixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> SharedMarkerVariance& override { return variance_; }
    auto variance() const -> const SharedMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> SampledMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const SampledMixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    SharedMarkerVariance variance_;
    SampledMixtureProportion mixture_proportion_;
};

class SingleFixedPerMarkerSpikeSlabGaussianPrior final
    : public SingleGeneticPrior
    , public SinglePerMarkerVarianceCap
    , public SingleFixedMixtureProportionCap
{
   public:
    static constexpr std::string_view name
        = "fixed_per_marker_spike_slab_gaussian";

    SingleFixedPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        FixedMixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> PerMarkerVariance& override { return variance_; }
    auto variance() const -> const PerMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> FixedMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const FixedMixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    PerMarkerVariance variance_;
    FixedMixtureProportion mixture_proportion_;
};

class SingleSampledPerMarkerSpikeSlabGaussianPrior final
    : public SingleGeneticPrior
    , public SinglePerMarkerVarianceCap
    , public SingleSampledMixtureProportionCap
{
   public:
    static constexpr std::string_view name
        = "sampled_per_marker_spike_slab_gaussian";

    SingleSampledPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        SampledMixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> PerMarkerVariance& override { return variance_; }
    auto variance() const -> const PerMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> SampledMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const SampledMixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<SingleGeneticPriorState> override;

   private:
    GeneticMode mode_;
    PerMarkerVariance variance_;
    SampledMixtureProportion mixture_proportion_;
};

class SingleFixedScaledMixtureGaussianPrior final
    : public SingleGeneticPrior
    , public SingleSharedMarkerVarianceCap
    , public SingleFixedMixtureProportionCap
    , public SingleMultiplierCap
{
   public:
    static constexpr std::string_view name = "fixed_scaled_mixture_gaussian";

    SingleFixedScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        FixedMixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> SharedMarkerVariance& override { return variance_; }
    auto variance() const -> const SharedMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> FixedMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const FixedMixtureProportion& override
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
    FixedMixtureProportion mixture_proportion_;
};

class SingleSampledScaledMixtureGaussianPrior final
    : public SingleGeneticPrior
    , public SingleSharedMarkerVarianceCap
    , public SingleSampledMixtureProportionCap
    , public SingleMultiplierCap
{
   public:
    static constexpr std::string_view name = "sampled_scaled_mixture_gaussian";

    SingleSampledScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        SampledMixtureProportion proportion);

    auto mode() const -> GeneticMode override { return mode_; }
    auto variance() -> SharedMarkerVariance& override { return variance_; }
    auto variance() const -> const SharedMarkerVariance& override
    {
        return variance_;
    }
    auto proportion() -> SampledMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const SampledMixtureProportion& override
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
    SampledMixtureProportion mixture_proportion_;
};

class JointFixedGaussianMixturePrior final
    : public JointGeneticPrior
    , public JointSharedMarkerVarianceCap
    , public JointFixedMixtureProportionCap
{
   public:
    static constexpr std::string_view name = "joint_fixed_mixture_gaussian";

    JointFixedGaussianMixturePrior(
        JointSharedMarkerVariance variance,
        FixedMixtureProportion proportion);

    auto variance(GeneticMode mode) -> SharedMarkerVariance& override;
    auto variance(GeneticMode mode) const
        -> const SharedMarkerVariance& override;
    auto proportion() -> FixedMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const FixedMixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<JointGeneticPriorState> override;

   private:
    JointSharedMarkerVariance marker_variance_;
    FixedMixtureProportion mixture_proportion_;
};

class JointSampledGaussianMixturePrior final
    : public JointGeneticPrior
    , public JointSharedMarkerVarianceCap
    , public JointSampledMixtureProportionCap
{
   public:
    static constexpr std::string_view name = "joint_sampled_mixture_gaussian";

    JointSampledGaussianMixturePrior(
        JointSharedMarkerVariance variance,
        SampledMixtureProportion proportion);

    auto variance(GeneticMode mode) -> SharedMarkerVariance& override;
    auto variance(GeneticMode mode) const
        -> const SharedMarkerVariance& override;
    auto proportion() -> SampledMixtureProportion& override
    {
        return mixture_proportion_;
    }
    auto proportion() const -> const SampledMixtureProportion& override
    {
        return mixture_proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<JointGeneticPriorState> override;

   private:
    JointSharedMarkerVariance marker_variance_;
    SampledMixtureProportion mixture_proportion_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIORS_GAUSSIAN_H_
