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

class SingleGaussianState final
    : public SingleGeneticPriorState
    , public SingleVarianceStateCap
{
   public:
    static constexpr std::string_view name = "gaussian";

    explicit SingleGaussianState(Eigen::VectorXd variance);

    auto variance() -> Eigen::VectorXd& override { return variance_; }
    auto variance() const -> const Eigen::VectorXd& override
    {
        return variance_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    Eigen::VectorXd variance_;
};

class SingleSpikeSlabGaussianState final
    : public SingleGeneticPriorState
    , public SingleVarianceStateCap
    , public SingleProportionStateCap
{
   public:
    static constexpr std::string_view name = "spike_slab_gaussian";

    SingleSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        const MixtureProportion& proportion,
        Eigen::Index num_markers);

    auto variance() -> Eigen::VectorXd& override { return variance_; }
    auto variance() const -> const Eigen::VectorXd& override
    {
        return variance_;
    }

    auto proportion() -> ProportionState& override { return proportion_; }
    auto proportion() const -> const ProportionState& override
    {
        return proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    Eigen::VectorXd variance_;
    ProportionState proportion_;
};

class SingleScaledMixtureGaussianState final
    : public SingleGeneticPriorState
    , public SingleVarianceStateCap
    , public SingleComponentStateCap
    , public SingleProportionStateCap
{
   public:
    static constexpr std::string_view name = "scaled_mixture_gaussian";

    SingleScaledMixtureGaussianState(
        Eigen::VectorXd variance,
        const Eigen::VectorXd& multiplier,
        const MixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance() -> Eigen::VectorXd& override { return variance_; }
    auto variance() const -> const Eigen::VectorXd& override
    {
        return variance_;
    }

    auto component() -> ComponentState& override { return component_; }
    auto component() const -> const ComponentState& override
    {
        return component_;
    }

    auto proportion() -> ProportionState& override { return proportion_; }
    auto proportion() const -> const ProportionState& override
    {
        return proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    Eigen::VectorXd variance_;
    ComponentState component_;
    ProportionState proportion_;
};

class JointGaussianMixtureState final
    : public JointGeneticPriorState
    , public JointVarianceStateCap
    , public JointComponentStateCap
    , public JointProportionStateCap
{
   public:
    static constexpr std::string_view name = "joint_mixture_gaussian";

    JointGaussianMixtureState(
        std::array<Eigen::VectorXd, 2> variances,
        const MixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance(GeneticMode mode) -> Eigen::VectorXd& override;
    auto variance(GeneticMode mode) const -> const Eigen::VectorXd& override;

    auto component() -> ComponentState& override { return component_; }
    auto component() const -> const ComponentState& override
    {
        return component_;
    }

    auto proportion() -> ProportionState& override { return proportion_; }
    auto proportion() const -> const ProportionState& override
    {
        return proportion_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;

   private:
    std::array<Eigen::VectorXd, 2> variances_;
    ComponentState component_;
    ProportionState proportion_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIOR_STATES_GAUSSIAN_H_
