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

#include "gelex/model/bayes/gaussian_prior.h"

#include <array>
#include <memory>
#include <ranges>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/gaussian_prior_state.h"

namespace gelex::bayes
{

namespace
{

auto make_variance_values(
    std::span<const MarkerVariance> marker_variances,
    Eigen::Index num_markers) -> std::vector<Eigen::VectorXd>
{
    std::vector<Eigen::VectorXd> values;
    values.reserve(marker_variances.size());
    for (const auto& marker_variance : marker_variances)
    {
        values.emplace_back(
            Eigen::VectorXd::Constant(
                marker_variance.marker_variance_size(num_markers),
                marker_variance.parameter().initial_value()));
    }
    return values;
}

}  // namespace

GaussianPrior::GaussianPrior(GeneticMode mode, MarkerVariance variance)
    : modes_{mode}, marker_variances_{variance}
{
}

auto GaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", modes_[0], FieldFlag::checkpoint | FieldFlag::report);
    marker_variances_[0].visit(visitor);
}

auto GaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<GaussianState>(
        make_variance_values(marker_variances_, num_markers));
}

SpikeSlabGaussianPrior::SpikeSlabGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance,
    MixtureProportion proportion)
    : modes_{mode},
      marker_variances_{variance},
      mixture_proportions_{std::move(proportion)}
{
    if (mixture_proportions_[0].size() != 2)
    {
        throw GelexException(
            "SpikeSlabGaussianPrior: proportion must have size 2");
    }
}

auto SpikeSlabGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", modes_[0], FieldFlag::checkpoint | FieldFlag::report);
    marker_variances_[0].visit(visitor);
    mixture_proportions_[0].visit(visitor);
}

auto SpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<SpikeSlabGaussianState>(
        make_variance_values(marker_variances_, num_markers),
        mixture_proportions_,
        num_markers);
}

ScaledMixtureGaussianPrior::ScaledMixtureGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance,
    Eigen::VectorXd multiplier,
    MixtureProportion proportion)
    : modes_{mode},
      marker_variances_{variance},
      multipliers_{std::move(multiplier)},
      mixture_proportions_{std::move(proportion)}
{
    if (multipliers_[0].size() != mixture_proportions_[0].size())
    {
        throw GelexException(
            "ScaledMixtureGaussianPrior: multiplier and proportion sizes "
            "differ");
    }
    if (multipliers_[0](0) != 0.0)
    {
        throw GelexException(
            "ScaledMixtureGaussianPrior: multiplier(0) must equal 0");
    }
}

auto ScaledMixtureGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", modes_[0], FieldFlag::checkpoint | FieldFlag::report);
    marker_variances_[0].visit(visitor);
    visitor.on(
        "multiplier",
        multipliers_[0],
        FieldFlag::checkpoint | FieldFlag::report);
    mixture_proportions_[0].visit(visitor);
}

auto ScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<ScaledMixtureGaussianState>(
        make_variance_values(marker_variances_, num_markers),
        multipliers_,
        mixture_proportions_,
        num_markers,
        num_individuals);
}

JointMixtureGaussianPrior::JointMixtureGaussianPrior(
    std::array<GeneticMode, 2> modes,
    std::array<MarkerVariance, 2> variances,
    MixtureProportion proportion)
    : modes_(modes),
      marker_variances_(variances),
      mixture_proportions_{std::move(proportion)}
{
    if (modes_[0] == modes_[1])
    {
        throw GelexException(
            "JointMixtureGaussianPrior: modes must be distinct");
    }
}

auto JointMixtureGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    for (auto [i, mode] : std::views::enumerate(modes_))
    {
        auto slot_scope = visitor.scope(std::to_string(i));
        visitor.on("mode", mode, FieldFlag::checkpoint | FieldFlag::report);
        marker_variances_[i].visit(visitor);
    }
    mixture_proportions_[0].visit(visitor);
}

auto JointMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> std::unique_ptr<GeneticPriorState>
{
    auto values = make_variance_values(marker_variances_, num_markers);
    return std::make_unique<JointMixtureGaussianState>(
        std::array<Eigen::VectorXd, 2>{
            std::move(values[0]), std::move(values[1])},
        mixture_proportions_[0],
        num_markers,
        num_individuals);
}

}  // namespace gelex::bayes
