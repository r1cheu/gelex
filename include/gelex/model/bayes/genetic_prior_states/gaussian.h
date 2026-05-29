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

#ifndef GELEX_MODEL_BAYES_GENETIC_PRIOR_STATES_GAUSSIAN_H_
#define GELEX_MODEL_BAYES_GENETIC_PRIOR_STATES_GAUSSIAN_H_

#include <array>
#include <string_view>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class SingleSharedGaussianState final
    : public SingleGeneticPriorState
    , public SingleSharedVarianceStateCap
{
   public:
    static constexpr std::string_view name = "shared_gaussian";

    explicit SingleSharedGaussianState(double variance);

    auto variance() -> double& override { return variance_; }
    auto variance() const -> const double& override { return variance_; }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    double variance_{0.0};
};

class SinglePerMarkerGaussianState final
    : public SingleGeneticPriorState
    , public SinglePerMarkerVarianceStateCap
{
   public:
    static constexpr std::string_view name = "per_marker_gaussian";

    explicit SinglePerMarkerGaussianState(Eigen::VectorXd variance);

    auto variance() -> Eigen::VectorXd& override { return variance_; }
    auto variance() const -> const Eigen::VectorXd& override
    {
        return variance_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    Eigen::VectorXd variance_;
};

class SingleFixedSharedSpikeSlabGaussianState final
    : public SingleGeneticPriorState
    , public SingleSharedVarianceStateCap
    , public SingleMixtureAssignmentStateCap
{
   public:
    static constexpr std::string_view name = "fixed_shared_spike_slab_gaussian";

    SingleFixedSharedSpikeSlabGaussianState(
        double variance,
        Eigen::Index num_components,
        Eigen::Index num_markers);

    auto variance() -> double& override { return variance_; }
    auto variance() const -> const double& override { return variance_; }

    auto assignment() -> MixtureAssignmentState& override
    {
        return assignment_;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return assignment_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    double variance_{0.0};
    MixtureAssignmentState assignment_;
};

class SingleSampledSharedSpikeSlabGaussianState final
    : public SingleGeneticPriorState
    , public SingleSharedVarianceStateCap
    , public SingleMixtureAssignmentStateCap
    , public SingleSampledMixtureProportionStateCap
{
   public:
    static constexpr std::string_view name
        = "sampled_shared_spike_slab_gaussian";

    SingleSampledSharedSpikeSlabGaussianState(
        double variance,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers);

    auto variance() -> double& override { return variance_; }
    auto variance() const -> const double& override { return variance_; }

    auto assignment() -> MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }

    auto proportion() -> SampledMixtureProportionState& override
    {
        return proportion_;
    }
    auto proportion() const -> const SampledMixtureProportionState& override
    {
        return proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    double variance_{0.0};
    SampledMixtureProportionState proportion_;
};

class SingleFixedPerMarkerSpikeSlabGaussianState final
    : public SingleGeneticPriorState
    , public SinglePerMarkerVarianceStateCap
    , public SingleMixtureAssignmentStateCap
{
   public:
    static constexpr std::string_view name
        = "fixed_per_marker_spike_slab_gaussian";

    SingleFixedPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        Eigen::Index num_components,
        Eigen::Index num_markers);

    auto variance() -> Eigen::VectorXd& override { return variance_; }
    auto variance() const -> const Eigen::VectorXd& override
    {
        return variance_;
    }

    auto assignment() -> MixtureAssignmentState& override
    {
        return assignment_;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return assignment_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    Eigen::VectorXd variance_;
    MixtureAssignmentState assignment_;
};

class SingleSampledPerMarkerSpikeSlabGaussianState final
    : public SingleGeneticPriorState
    , public SinglePerMarkerVarianceStateCap
    , public SingleMixtureAssignmentStateCap
    , public SingleSampledMixtureProportionStateCap
{
   public:
    static constexpr std::string_view name
        = "sampled_per_marker_spike_slab_gaussian";

    SingleSampledPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers);

    auto variance() -> Eigen::VectorXd& override { return variance_; }
    auto variance() const -> const Eigen::VectorXd& override
    {
        return variance_;
    }

    auto assignment() -> MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }

    auto proportion() -> SampledMixtureProportionState& override
    {
        return proportion_;
    }
    auto proportion() const -> const SampledMixtureProportionState& override
    {
        return proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    Eigen::VectorXd variance_;
    SampledMixtureProportionState proportion_;
};

class SingleFixedScaledMixtureGaussianState final
    : public SingleGeneticPriorState
    , public SingleSharedVarianceStateCap
    , public SingleComponentStateCap
    , public SingleMixtureAssignmentStateCap
{
   public:
    static constexpr std::string_view name = "fixed_scaled_mixture_gaussian";

    SingleFixedScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance() -> double& override { return variance_; }
    auto variance() const -> const double& override { return variance_; }

    auto component() -> ComponentState& override { return component_; }
    auto component() const -> const ComponentState& override
    {
        return component_;
    }

    auto assignment() -> MixtureAssignmentState& override
    {
        return assignment_;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return assignment_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    double variance_{0.0};
    ComponentState component_;
    MixtureAssignmentState assignment_;
};

class SingleSampledScaledMixtureGaussianState final
    : public SingleGeneticPriorState
    , public SingleSharedVarianceStateCap
    , public SingleComponentStateCap
    , public SingleMixtureAssignmentStateCap
    , public SingleSampledMixtureProportionStateCap
{
   public:
    static constexpr std::string_view name = "sampled_scaled_mixture_gaussian";

    SingleSampledScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance() -> double& override { return variance_; }
    auto variance() const -> const double& override { return variance_; }

    auto component() -> ComponentState& override { return component_; }
    auto component() const -> const ComponentState& override
    {
        return component_;
    }

    auto assignment() -> MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }

    auto proportion() -> SampledMixtureProportionState& override
    {
        return proportion_;
    }
    auto proportion() const -> const SampledMixtureProportionState& override
    {
        return proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    double variance_{0.0};
    ComponentState component_;
    SampledMixtureProportionState proportion_;
};

class JointFixedGaussianMixtureState final
    : public JointGeneticPriorState
    , public JointSharedVarianceStateCap
    , public JointComponentStateCap
    , public JointMixtureAssignmentStateCap
{
   public:
    static constexpr std::string_view name = "joint_fixed_mixture_gaussian";

    JointFixedGaussianMixtureState(
        std::array<double, 2> variances,
        Eigen::Index num_components,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance(GeneticMode mode) -> double& override;
    auto variance(GeneticMode mode) const -> const double& override;

    auto component() -> ComponentState& override { return component_; }
    auto component() const -> const ComponentState& override
    {
        return component_;
    }

    auto assignment() -> MixtureAssignmentState& override
    {
        return assignment_;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return assignment_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    std::array<double, 2> variances_;
    ComponentState component_;
    MixtureAssignmentState assignment_;
};

class JointSampledGaussianMixtureState final
    : public JointGeneticPriorState
    , public JointSharedVarianceStateCap
    , public JointComponentStateCap
    , public JointMixtureAssignmentStateCap
    , public JointSampledMixtureProportionStateCap
{
   public:
    static constexpr std::string_view name = "joint_sampled_mixture_gaussian";

    JointSampledGaussianMixtureState(
        std::array<double, 2> variances,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance(GeneticMode mode) -> double& override;
    auto variance(GeneticMode mode) const -> const double& override;

    auto component() -> ComponentState& override { return component_; }
    auto component() const -> const ComponentState& override
    {
        return component_;
    }

    auto assignment() -> MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }
    auto assignment() const -> const MixtureAssignmentState& override
    {
        return proportion_.assignment;
    }

    auto proportion() -> SampledMixtureProportionState& override
    {
        return proportion_;
    }
    auto proportion() const -> const SampledMixtureProportionState& override
    {
        return proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    std::array<double, 2> variances_;
    ComponentState component_;
    SampledMixtureProportionState proportion_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIOR_STATES_GAUSSIAN_H_
