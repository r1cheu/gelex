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

#include "gelex/model/bayes/fields.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state_values.h"

namespace gelex::bayes
{

class SingleSharedGaussianState final : public Variance<double>
{
   public:
    static constexpr std::string_view name = "shared_gaussian";

    explicit SingleSharedGaussianState(double variance);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SinglePerMarkerGaussianState final : public Variance<Eigen::VectorXd>
{
   public:
    static constexpr std::string_view name = "per_marker_gaussian";

    explicit SinglePerMarkerGaussianState(Eigen::VectorXd variance);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SingleFixedSharedSpikeSlabGaussianState final
    : public Variance<double>
    , public AssignmentField<MixtureAssignmentState>
{
   public:
    static constexpr std::string_view name = "fixed_shared_spike_slab_gaussian";

    SingleFixedSharedSpikeSlabGaussianState(
        double variance,
        Eigen::Index num_components,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SingleSampledSharedSpikeSlabGaussianState final
    : public Variance<double>
    , public SampledProportionField<SampledMixtureProportionState>
{
   public:
    static constexpr std::string_view name
        = "sampled_shared_spike_slab_gaussian";

    SingleSampledSharedSpikeSlabGaussianState(
        double variance,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SingleFixedPerMarkerSpikeSlabGaussianState final
    : public Variance<Eigen::VectorXd>
    , public AssignmentField<MixtureAssignmentState>
{
   public:
    static constexpr std::string_view name
        = "fixed_per_marker_spike_slab_gaussian";

    SingleFixedPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        Eigen::Index num_components,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SingleSampledPerMarkerSpikeSlabGaussianState final
    : public Variance<Eigen::VectorXd>
    , public SampledProportionField<SampledMixtureProportionState>
{
   public:
    static constexpr std::string_view name
        = "sampled_per_marker_spike_slab_gaussian";

    SingleSampledPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SingleFixedScaledMixtureGaussianState final
    : public Variance<double>
    , public ComponentField<ComponentState>
    , public AssignmentField<MixtureAssignmentState>
{
   public:
    static constexpr std::string_view name = "fixed_scaled_mixture_gaussian";

    SingleFixedScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SingleSampledScaledMixtureGaussianState final
    : public Variance<double>
    , public ComponentField<ComponentState>
    , public SampledProportionField<SampledMixtureProportionState>
{
   public:
    static constexpr std::string_view name = "sampled_scaled_mixture_gaussian";

    SingleSampledScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class JointFixedGaussianMixtureState final
    : public Variances<double>
    , public ComponentField<ComponentState>
    , public AssignmentField<MixtureAssignmentState>
{
   public:
    static constexpr std::string_view name = "joint_fixed_mixture_gaussian";

    JointFixedGaussianMixtureState(
        std::array<double, 2> variances,
        Eigen::Index num_components,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class JointSampledGaussianMixtureState final
    : public Variances<double>
    , public ComponentField<ComponentState>
    , public SampledProportionField<SampledMixtureProportionState>
{
   public:
    static constexpr std::string_view name = "joint_sampled_mixture_gaussian";

    JointSampledGaussianMixtureState(
        std::array<double, 2> variances,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIOR_STATES_GAUSSIAN_H_
