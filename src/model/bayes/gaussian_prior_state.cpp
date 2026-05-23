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

#include "gelex/model/bayes/gaussian_prior_state.h"

#include <span>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_specs.h"
#include "gelex/model/bayes/prior_state.h"

namespace gelex::bayes
{

GaussianState::GaussianState(std::vector<Eigen::VectorXd> variances)
    : variances_(std::move(variances))
{
}

auto GaussianState::visit_sample_records(infra::RecordVisitor& visitor) const
    -> void
{
    VarianceStateCap::visit_sample_records(visitor);
}

auto GaussianState::visit_checkpoint_records(
    infra::RecordVisitor& visitor) const -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
}

auto GaussianState::visit_checkpoint_records(
    infra::MutableRecordVisitor& visitor) -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
}

SpikeSlabGaussianState::SpikeSlabGaussianState(
    std::vector<Eigen::VectorXd> variances,
    std::span<const ProportionSpec> proportion_specs,
    Eigen::Index num_markers)
    : variances_(std::move(variances))
{
    proportions_.reserve(proportion_specs.size());
    for (const auto& spec : proportion_specs)
    {
        proportions_.emplace_back(spec, num_markers);
    }
}

auto SpikeSlabGaussianState::visit_sample_records(
    infra::RecordVisitor& visitor) const -> void
{
    VarianceStateCap::visit_sample_records(visitor);
    ProportionStateCap::visit_sample_records(visitor);
}

auto SpikeSlabGaussianState::visit_checkpoint_records(
    infra::RecordVisitor& visitor) const -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
    ProportionStateCap::visit_checkpoint_records(visitor);
}

auto SpikeSlabGaussianState::visit_checkpoint_records(
    infra::MutableRecordVisitor& visitor) -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
    ProportionStateCap::visit_checkpoint_records(visitor);
}

ScaledMixtureGaussianState::ScaledMixtureGaussianState(
    std::vector<Eigen::VectorXd> base_variances,
    std::span<const Eigen::VectorXd> multipliers,
    std::span<const ProportionSpec> proportion_specs,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : variances_(std::move(base_variances))
{
    components_.reserve(multipliers.size());
    for (const auto& mult : multipliers)
    {
        components_.emplace_back(mult.size() - 1, num_individuals);
    }
    proportions_.reserve(proportion_specs.size());
    for (const auto& spec : proportion_specs)
    {
        proportions_.emplace_back(spec, num_markers);
    }
}

auto ScaledMixtureGaussianState::visit_sample_records(
    infra::RecordVisitor& visitor) const -> void
{
    VarianceStateCap::visit_sample_records(visitor);
    ComponentStateCap::visit_sample_records(visitor);
    ProportionStateCap::visit_sample_records(visitor);
}

auto ScaledMixtureGaussianState::visit_checkpoint_records(
    infra::RecordVisitor& visitor) const -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
    ComponentStateCap::visit_checkpoint_records(visitor);
    ProportionStateCap::visit_checkpoint_records(visitor);
}

auto ScaledMixtureGaussianState::visit_checkpoint_records(
    infra::MutableRecordVisitor& visitor) -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
    ComponentStateCap::visit_checkpoint_records(visitor);
    ProportionStateCap::visit_checkpoint_records(visitor);
}

JointMixtureGaussianState::JointMixtureGaussianState(
    std::array<Eigen::VectorXd, 2> variances,
    const ProportionSpec& proportion_spec,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : variances_(std::move(variances)),
      components_{ComponentState{1, num_individuals}},
      proportions_{ProportionState{proportion_spec, num_markers}}
{
}

auto JointMixtureGaussianState::visit_sample_records(
    infra::RecordVisitor& visitor) const -> void
{
    VarianceStateCap::visit_sample_records(visitor);
    ComponentStateCap::visit_sample_records(visitor);
    ProportionStateCap::visit_sample_records(visitor);
}

auto JointMixtureGaussianState::visit_checkpoint_records(
    infra::RecordVisitor& visitor) const -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
    ComponentStateCap::visit_checkpoint_records(visitor);
    ProportionStateCap::visit_checkpoint_records(visitor);
}

auto JointMixtureGaussianState::visit_checkpoint_records(
    infra::MutableRecordVisitor& visitor) -> void
{
    VarianceStateCap::visit_checkpoint_records(visitor);
    ComponentStateCap::visit_checkpoint_records(visitor);
    ProportionStateCap::visit_checkpoint_records(visitor);
}

}  // namespace gelex::bayes
