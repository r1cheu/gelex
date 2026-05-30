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

#include "gelex/model/bayes/genetic_prior_states/gaussian.h"

#include <array>
#include <utility>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/genetic_prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

SingleSharedGaussianState::SingleSharedGaussianState(double variance)
    : Variance<double>(variance)
{
}

auto SingleSharedGaussianState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
}

SinglePerMarkerGaussianState::SinglePerMarkerGaussianState(
    Eigen::VectorXd variance)
    : Variance<Eigen::VectorXd>(std::move(variance))
{
}

auto SinglePerMarkerGaussianState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
}

SingleFixedSharedSpikeSlabGaussianState::
    SingleFixedSharedSpikeSlabGaussianState(
        double variance,
        Eigen::Index num_components,
        Eigen::Index num_markers)
    : Variance<double>(variance),
      AssignmentField<AssignmentState>(
          AssignmentState{num_components, num_markers})
{
}

auto SingleFixedSharedSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    assignment().visit(visitor);
}

SingleSampledSharedSpikeSlabGaussianState::
    SingleSampledSharedSpikeSlabGaussianState(
        double variance,
        const SampledProportion& proportion,
        Eigen::Index num_markers)
    : Variance<double>(variance),
      SampledProportionField<SampledProportionState>(
          SampledProportionState{proportion, num_markers})
{
}

auto SingleSampledSharedSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    proportion().visit(visitor);
}

SingleFixedPerMarkerSpikeSlabGaussianState::
    SingleFixedPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        Eigen::Index num_components,
        Eigen::Index num_markers)
    : Variance<Eigen::VectorXd>(std::move(variance)),
      AssignmentField<AssignmentState>(
          AssignmentState{num_components, num_markers})
{
}

auto SingleFixedPerMarkerSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    assignment().visit(visitor);
}

SingleSampledPerMarkerSpikeSlabGaussianState::
    SingleSampledPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        const SampledProportion& proportion,
        Eigen::Index num_markers)
    : Variance<Eigen::VectorXd>(std::move(variance)),
      SampledProportionField<SampledProportionState>(
          SampledProportionState{proportion, num_markers})
{
}

auto SingleSampledPerMarkerSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    proportion().visit(visitor);
}

SingleFixedScaledMixtureGaussianState::SingleFixedScaledMixtureGaussianState(
    double variance,
    const Eigen::VectorXd& multiplier,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : Variance<double>(variance),
      ComponentField<ComponentState>(
          ComponentState{multiplier.size() - 1, num_individuals}),
      AssignmentField<AssignmentState>(
          AssignmentState{multiplier.size(), num_markers})
{
}

auto SingleFixedScaledMixtureGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    component().visit(visitor);
    assignment().visit(visitor);
}

SingleSampledScaledMixtureGaussianState::
    SingleSampledScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        const SampledProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals)
    : Variance<double>(variance),
      ComponentField<ComponentState>(
          ComponentState{multiplier.size() - 1, num_individuals}),
      SampledProportionField<SampledProportionState>(
          SampledProportionState{proportion, num_markers})
{
}

auto SingleSampledScaledMixtureGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    component().visit(visitor);
    proportion().visit(visitor);
}

JointFixedGaussianMixtureState::JointFixedGaussianMixtureState(
    std::array<double, 2> variances,
    Eigen::Index num_components,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : Variances<double>(std::move(variances)),
      ComponentField<ComponentState>(ComponentState{1, num_individuals}),
      AssignmentField<AssignmentState>(
          AssignmentState{num_components, num_markers})
{
}

auto JointFixedGaussianMixtureState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (const auto mode : modes)
    {
        auto mode_scope = visitor.scope(fmt::format("{}", mode));
        visitor.on(
            "variance",
            variance(mode),
            FieldFlag::checkpoint | FieldFlag::trace);
    }
    component().visit(visitor);
    assignment().visit(visitor);
}

JointSampledGaussianMixtureState::JointSampledGaussianMixtureState(
    std::array<double, 2> variances,
    const SampledProportion& proportion,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : Variances<double>(std::move(variances)),
      ComponentField<ComponentState>(ComponentState{1, num_individuals}),
      SampledProportionField<SampledProportionState>(
          SampledProportionState{proportion, num_markers})
{
}

auto JointSampledGaussianMixtureState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (const auto mode : modes)
    {
        auto mode_scope = visitor.scope(fmt::format("{}", mode));
        visitor.on(
            "variance",
            variance(mode),
            FieldFlag::checkpoint | FieldFlag::trace);
    }
    component().visit(visitor);
    proportion().visit(visitor);
}

}  // namespace gelex::bayes
