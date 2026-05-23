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

#ifndef GELEX_MODEL_BAYES_GAUSSIAN_PRIOR_STATE_H_
#define GELEX_MODEL_BAYES_GAUSSIAN_PRIOR_STATE_H_

#include <array>
#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_specs.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::bayes
{

class GaussianState final
    : public GeneticPriorState
    , public VarianceStateCap
{
   public:
    explicit GaussianState(std::vector<Eigen::VectorXd> variances);

    auto variance() -> std::span<Eigen::VectorXd> override
    {
        return variances_;
    }
    auto variance() const -> std::span<const Eigen::VectorXd> override
    {
        return variances_;
    }

    auto visit_sample_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::MutableRecordVisitor& visitor)
        -> void override;

   private:
    std::vector<Eigen::VectorXd> variances_;
};

class SpikeSlabGaussianState final
    : public GeneticPriorState
    , public VarianceStateCap
    , public ProportionStateCap
{
   public:
    SpikeSlabGaussianState(
        std::vector<Eigen::VectorXd> variances,
        std::span<const ProportionSpec> proportion_specs,
        Eigen::Index num_markers);

    auto variance() -> std::span<Eigen::VectorXd> override
    {
        return variances_;
    }
    auto variance() const -> std::span<const Eigen::VectorXd> override
    {
        return variances_;
    }

    auto proportion() -> std::span<ProportionState> override
    {
        return proportions_;
    }
    auto proportion() const -> std::span<const ProportionState> override
    {
        return proportions_;
    }

    auto visit_sample_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::MutableRecordVisitor& visitor)
        -> void override;

   private:
    std::vector<Eigen::VectorXd> variances_;
    std::vector<ProportionState> proportions_;
};

class ScaledMixtureGaussianState final
    : public GeneticPriorState
    , public VarianceStateCap
    , public ComponentStateCap
    , public ProportionStateCap
{
   public:
    ScaledMixtureGaussianState(
        std::vector<Eigen::VectorXd> base_variances,
        std::span<const Eigen::VectorXd> multipliers,
        std::span<const ProportionSpec> proportion_specs,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance() -> std::span<Eigen::VectorXd> override
    {
        return variances_;
    }
    auto variance() const -> std::span<const Eigen::VectorXd> override
    {
        return variances_;
    }

    auto component() -> std::span<ComponentState> override
    {
        return components_;
    }
    auto component() const -> std::span<const ComponentState> override
    {
        return components_;
    }

    auto proportion() -> std::span<ProportionState> override
    {
        return proportions_;
    }
    auto proportion() const -> std::span<const ProportionState> override
    {
        return proportions_;
    }

    auto visit_sample_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::MutableRecordVisitor& visitor)
        -> void override;

   private:
    std::vector<Eigen::VectorXd> variances_;
    std::vector<ComponentState> components_;
    std::vector<ProportionState> proportions_;
};

class JointMixtureGaussianState final
    : public GeneticPriorState
    , public VarianceStateCap
    , public ComponentStateCap
    , public ProportionStateCap
{
   public:
    JointMixtureGaussianState(
        std::array<Eigen::VectorXd, 2> variances,
        const ProportionSpec& proportion_spec,
        Eigen::Index num_markers,
        Eigen::Index num_individuals);

    auto variance() -> std::span<Eigen::VectorXd> override
    {
        return variances_;
    }
    auto variance() const -> std::span<const Eigen::VectorXd> override
    {
        return variances_;
    }

    auto component() -> std::span<ComponentState> override
    {
        return components_;
    }
    auto component() const -> std::span<const ComponentState> override
    {
        return components_;
    }

    auto proportion() -> std::span<ProportionState> override
    {
        return proportions_;
    }
    auto proportion() const -> std::span<const ProportionState> override
    {
        return proportions_;
    }

    auto visit_sample_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::RecordVisitor& visitor) const
        -> void override;
    auto visit_checkpoint_records(infra::MutableRecordVisitor& visitor)
        -> void override;

   private:
    std::array<Eigen::VectorXd, 2> variances_;
    std::array<ComponentState, 1> components_;
    std::array<ProportionState, 1> proportions_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GAUSSIAN_PRIOR_STATE_H_
