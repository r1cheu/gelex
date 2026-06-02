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

#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior_state_values.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class SingleSharedGaussianState final
{
   public:
    static constexpr std::string_view name = "shared_gaussian";

    explicit SingleSharedGaussianState(double variance);

    auto variance() -> double& { return variance_; }
    auto variance() const -> const double& { return variance_; }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    double variance_;
};

class SinglePerMarkerGaussianState final
{
   public:
    static constexpr std::string_view name = "per_marker_gaussian";

    explicit SinglePerMarkerGaussianState(Eigen::VectorXd variance);

    auto variance() -> Eigen::VectorXd& { return variance_; }
    auto variance() const -> const Eigen::VectorXd& { return variance_; }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    Eigen::VectorXd variance_;
};

class SingleSharedSpikeSlabGaussianState final
{
   public:
    static constexpr std::string_view name = "shared_spike_slab_gaussian";

    SingleSharedSpikeSlabGaussianState(
        double variance,
        const MixtureProportion& proportion,
        Eigen::Index num_markers);

    auto variance() -> double& { return variance_; }
    auto variance() const -> const double& { return variance_; }
    auto assignment() -> Eigen::VectorXi& { return mixture_.assignment; }
    auto assignment() const -> const Eigen::VectorXi&
    {
        return mixture_.assignment;
    }
    auto proportion() -> Eigen::VectorXd& { return mixture_.proportion; }
    auto proportion() const -> const Eigen::VectorXd&
    {
        return mixture_.proportion;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    double variance_;
    MixtureState mixture_;
};

class SinglePerMarkerSpikeSlabGaussianState final
{
   public:
    static constexpr std::string_view name = "per_marker_spike_slab_gaussian";

    SinglePerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        const MixtureProportion& proportion,
        Eigen::Index num_markers);

    auto variance() -> Eigen::VectorXd& { return variance_; }
    auto variance() const -> const Eigen::VectorXd& { return variance_; }
    auto assignment() -> Eigen::VectorXi& { return mixture_.assignment; }
    auto assignment() const -> const Eigen::VectorXi&
    {
        return mixture_.assignment;
    }
    auto proportion() -> Eigen::VectorXd& { return mixture_.proportion; }
    auto proportion() const -> const Eigen::VectorXd&
    {
        return mixture_.proportion;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    Eigen::VectorXd variance_;
    MixtureState mixture_;
};

class SingleScaledMixtureGaussianState final
{
   public:
    static constexpr std::string_view name = "scaled_mixture_gaussian";

    SingleScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        const MixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance() -> double& { return variance_; }
    auto variance() const -> const double& { return variance_; }
    auto component() -> ComponentState& { return component_; }
    auto component() const -> const ComponentState& { return component_; }
    auto assignment() -> Eigen::VectorXi& { return mixture_.assignment; }
    auto assignment() const -> const Eigen::VectorXi&
    {
        return mixture_.assignment;
    }
    auto proportion() -> Eigen::VectorXd& { return mixture_.proportion; }
    auto proportion() const -> const Eigen::VectorXd&
    {
        return mixture_.proportion;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    double variance_;
    ComponentState component_;
    MixtureState mixture_;
};

class JointGaussianMixtureState final
{
   public:
    static constexpr std::string_view name = "joint_mixture_gaussian";

    JointGaussianMixtureState(
        std::array<double, 2> variances,
        const MixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance(GeneticMode mode) -> double&;
    auto variance(GeneticMode mode) const -> const double&;
    auto component() -> ComponentState& { return component_; }
    auto component() const -> const ComponentState& { return component_; }
    auto assignment() -> Eigen::VectorXi& { return mixture_.assignment; }
    auto assignment() const -> const Eigen::VectorXi&
    {
        return mixture_.assignment;
    }
    auto proportion() -> Eigen::VectorXd& { return mixture_.proportion; }
    auto proportion() const -> const Eigen::VectorXd&
    {
        return mixture_.proportion;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    std::array<double, 2> variances_;
    ComponentState component_;
    MixtureState mixture_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_PRIOR_STATE_H_
