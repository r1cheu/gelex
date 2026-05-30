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

#include "gelex/bayes/genetic/gaussian_prior_state.h"

#include <array>
#include <utility>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior_state.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
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

SingleSharedSpikeSlabGaussianState::SingleSharedSpikeSlabGaussianState(
    double variance,
    const MixtureProportion& proportion,
    Eigen::Index num_markers)
    : Variance<double>(variance),
      MixtureField<MixtureState>(MixtureState{proportion, num_markers})
{
}

auto SingleSharedSpikeSlabGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    mixture().visit(visitor);
}

SinglePerMarkerSpikeSlabGaussianState::SinglePerMarkerSpikeSlabGaussianState(
    Eigen::VectorXd variance,
    const MixtureProportion& proportion,
    Eigen::Index num_markers)
    : Variance<Eigen::VectorXd>(std::move(variance)),
      MixtureField<MixtureState>(MixtureState{proportion, num_markers})
{
}

auto SinglePerMarkerSpikeSlabGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    mixture().visit(visitor);
}

SingleScaledMixtureGaussianState::SingleScaledMixtureGaussianState(
    double variance,
    const Eigen::VectorXd& multiplier,
    const MixtureProportion& proportion,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : Variance<double>(variance),
      ComponentField<ComponentState>(
          ComponentState{multiplier.size() - 1, num_individuals}),
      MixtureField<MixtureState>(MixtureState{proportion, num_markers})
{
}

auto SingleScaledMixtureGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on(
        "variance", variance(), FieldFlag::checkpoint | FieldFlag::trace);
    component().visit(visitor);
    mixture().visit(visitor);
}

JointGaussianMixtureState::JointGaussianMixtureState(
    std::array<double, 2> variances,
    const MixtureProportion& proportion,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : Variances<double>(std::move(variances)),
      ComponentField<ComponentState>(ComponentState{1, num_individuals}),
      MixtureField<MixtureState>(MixtureState{proportion, num_markers})
{
}

auto JointGaussianMixtureState::visit(infra::FieldVisitor& visitor) -> void
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
    mixture().visit(visitor);
}

}  // namespace gelex::bayes
