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

SingleFixedSharedSpikeSlabGaussianState::
    SingleFixedSharedSpikeSlabGaussianState(
        double variance,
        Eigen::Index num_components,
        Eigen::Index num_markers)
    : variance_(variance), assignment_(num_components, num_markers)
{
}

auto SingleFixedSharedSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    assignment_.visit(visitor);
}

SingleSampledSharedSpikeSlabGaussianState::
    SingleSampledSharedSpikeSlabGaussianState(
        double variance,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers)
    : variance_(variance), proportion_(proportion, num_markers)
{
}

auto SingleSampledSharedSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    proportion_.visit(visitor);
}

SingleFixedPerMarkerSpikeSlabGaussianState::
    SingleFixedPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        Eigen::Index num_components,
        Eigen::Index num_markers)
    : variance_(std::move(variance)), assignment_(num_components, num_markers)
{
}

auto SingleFixedPerMarkerSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    assignment_.visit(visitor);
}

SingleSampledPerMarkerSpikeSlabGaussianState::
    SingleSampledPerMarkerSpikeSlabGaussianState(
        Eigen::VectorXd variance,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers)
    : variance_(std::move(variance)), proportion_(proportion, num_markers)
{
}

auto SingleSampledPerMarkerSpikeSlabGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    proportion_.visit(visitor);
}

SingleFixedScaledMixtureGaussianState::SingleFixedScaledMixtureGaussianState(
    double variance,
    const Eigen::VectorXd& multiplier,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : variance_(variance),
      component_(multiplier.size() - 1, num_individuals),
      assignment_(multiplier.size(), num_markers)
{
}

auto SingleFixedScaledMixtureGaussianState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    component_.visit(visitor);
    assignment_.visit(visitor);
}

SingleSampledScaledMixtureGaussianState::
    SingleSampledScaledMixtureGaussianState(
        double variance,
        const Eigen::VectorXd& multiplier,
        const SampledMixtureProportion& proportion,
        Eigen::Index num_markers,
        Eigen::Index num_individuals)
    : variance_(variance),
      component_(multiplier.size() - 1, num_individuals),
      proportion_(proportion, num_markers)
{
}

auto SingleSampledScaledMixtureGaussianState::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("variance", variance_, FieldFlag::checkpoint | FieldFlag::trace);
    component_.visit(visitor);
    proportion_.visit(visitor);
}

JointFixedGaussianMixtureState::JointFixedGaussianMixtureState(
    std::array<double, 2> variances,
    Eigen::Index num_components,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : variances_(std::move(variances)),
      component_(1, num_individuals),
      assignment_(num_components, num_markers)
{
}

auto JointFixedGaussianMixtureState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto [i, mode] : std::views::enumerate(modes))
    {
        auto mode_scope = visitor.scope(fmt::format("{}", mode));
        visitor.on(
            "variance",
            variances_[i],
            FieldFlag::checkpoint | FieldFlag::trace);
    }
    component_.visit(visitor);
    assignment_.visit(visitor);
}

auto JointFixedGaussianMixtureState::variance(GeneticMode mode) -> double&
{
    return variances_[std::to_underlying(mode)];
}

auto JointFixedGaussianMixtureState::variance(GeneticMode mode) const -> const
    double&
{
    return variances_[std::to_underlying(mode)];
}

JointSampledGaussianMixtureState::JointSampledGaussianMixtureState(
    std::array<double, 2> variances,
    const SampledMixtureProportion& proportion,
    Eigen::Index num_markers,
    Eigen::Index num_individuals)
    : variances_(std::move(variances)),
      component_(1, num_individuals),
      proportion_(proportion, num_markers)
{
}

auto JointSampledGaussianMixtureState::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto [i, mode] : std::views::enumerate(modes))
    {
        auto mode_scope = visitor.scope(fmt::format("{}", mode));
        visitor.on(
            "variance",
            variances_[i],
            FieldFlag::checkpoint | FieldFlag::trace);
    }
    component_.visit(visitor);
    proportion_.visit(visitor);
}

auto JointSampledGaussianMixtureState::variance(GeneticMode mode) -> double&
{
    return variances_[std::to_underlying(mode)];
}

auto JointSampledGaussianMixtureState::variance(GeneticMode mode) const -> const
    double&
{
    return variances_[std::to_underlying(mode)];
}

}  // namespace gelex::bayes
