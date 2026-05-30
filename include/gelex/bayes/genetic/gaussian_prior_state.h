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

#ifndef GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_STATE_H_
#define GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_STATE_H_

#include <array>
#include <string_view>

#include <Eigen/Core>

#include "gelex/bayes/genetic/fields.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior_state_values.h"

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

class SingleSharedSpikeSlabGaussianState final
    : public Variance<double>
    , public MixtureField<MixtureState>
{
   public:
    static constexpr std::string_view name = "shared_spike_slab_gaussian";

    SingleSharedSpikeSlabGaussianState(
        double variance,
        const MixtureProportion& proportion,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SinglePerMarkerSpikeSlabGaussianState final
    : public Variance<Eigen::VectorXd>
    , public MixtureField<MixtureState>
{
   public:
    static constexpr std::string_view name = "per_marker_spike_slab_gaussian";

    SinglePerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        const MixtureProportion& proportion,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class SingleScaledMixtureGaussianState final
    : public Variance<double>
    , public ComponentField<ComponentState>
    , public MixtureField<MixtureState>
{
   public:
    static constexpr std::string_view name = "scaled_mixture_gaussian";

    SingleScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        const MixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

class JointGaussianMixtureState final
    : public Variances<double>
    , public ComponentField<ComponentState>
    , public MixtureField<MixtureState>
{
   public:
    static constexpr std::string_view name = "joint_mixture_gaussian";

    JointGaussianMixtureState(
        std::array<double, 2> variances,
        const MixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_STATE_H_
