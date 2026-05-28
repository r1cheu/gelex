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
#include <ranges>
#include <utility>
#include <variant>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

SingleSharedGaussianState::SingleSharedGaussianState(double variance)
    : variance_(variance)
{
}

auto SingleSharedGaussianState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
}

SinglePerMarkerGaussianState::SinglePerMarkerGaussianState(
    Eigen::VectorXd variance)
    : variance_(std::move(variance))
{
}

auto SinglePerMarkerGaussianState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
}

SingleSharedSpikeSlabGaussianState::SingleSharedSpikeSlabGaussianState(
    double variance,
    const MixtureProportion& proportion,
    Eigen::Index num_markers)
    : variance_(variance), proportion_(proportion, num_markers)
{
}

auto SingleSharedSpikeSlabGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    proportion_.visit(visitor);
}

SinglePerMarkerSpikeSlabGaussianState::SinglePerMarkerSpikeSlabGaussianState(
    Eigen::VectorXd variance,
    const MixtureProportion& proportion,
    Eigen::Index num_markers)
    : variance_(std::move(variance)), proportion_(proportion, num_markers)
{
}

auto SinglePerMarkerSpikeSlabGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    proportion_.visit(visitor);
}

SingleScaledMixtureGaussianState::SingleScaledMixtureGaussianState(
    double variance,
    const Eigen::VectorXd& multiplier,
    const MixtureProportion& proportion,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : variance_(variance),
      component_(multiplier.size() - 1, num_individuals),
      proportion_(proportion, num_markers)
{
}

auto SingleScaledMixtureGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    component_.visit(visitor);
    proportion_.visit(visitor);
}

JointGaussianMixtureState::JointGaussianMixtureState(
    std::array<JointMarkerVarianceState, 2> variances,
    const MixtureProportion& proportion,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : variances_(std::move(variances)),
      component_(1, num_individuals),
      proportion_(proportion, num_markers)
{
}

auto JointGaussianMixtureState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto [i, mode] : std::views::enumerate(modes))
    {
        auto mode_scope = visitor.scope(fmt::format("{}", mode));
        std::visit(
            [&visitor](auto& variance)
            {
                visitor.on(
                    "variance",
                    variance,
                    FieldFlag::checkpoint | FieldFlag::trace);
            },
            variances_[i]);
    }
    component_.visit(visitor);
    proportion_.visit(visitor);
}

auto JointGaussianMixtureState::variance(GeneticMode mode)
    -> JointMarkerVarianceState&
{
    return variances_[std::to_underlying(mode)];
}

auto JointGaussianMixtureState::variance(GeneticMode mode) const
    -> const JointMarkerVarianceState&
{
    return variances_[std::to_underlying(mode)];
}

}  // namespace gelex::bayes
