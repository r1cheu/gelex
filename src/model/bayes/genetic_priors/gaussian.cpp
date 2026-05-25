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

#include "gelex/model/bayes/genetic_priors/gaussian.h"

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
#include "gelex/model/bayes/prior_parameters.h"

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

SingleGaussianPrior::SingleGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance)
    : mode_(mode), marker_variances_{variance}
{
}

auto SingleGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::source | FieldFlag::metadata);
    marker_variances_[0].visit(visitor);
}

auto SingleGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<GaussianState>(
        make_variance_values(marker_variances_, num_markers));
}

SingleSpikeSlabGaussianPrior::SingleSpikeSlabGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance,
    MixtureProportion proportion)
    : mode_(mode),
      marker_variances_{variance},
      mixture_proportions_{std::move(proportion)}
{
    if (mixture_proportions_[0].size() != 2)
    {
        throw GelexException(
            "SingleSpikeSlabGaussianPrior: proportion must have size 2");
    }
}

auto SingleSpikeSlabGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::source | FieldFlag::metadata);
    marker_variances_[0].visit(visitor);
    mixture_proportions_[0].visit(visitor);
}

auto SingleSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<SpikeSlabGaussianState>(
        make_variance_values(marker_variances_, num_markers),
        mixture_proportions_,
        num_markers);
}

SingleScaledMixtureGaussianPrior::SingleScaledMixtureGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance,
    Eigen::VectorXd multiplier,
    MixtureProportion proportion)
    : mode_(mode),
      marker_variances_{variance},
      multipliers_{std::move(multiplier)},
      mixture_proportions_{std::move(proportion)}
{
    if (multipliers_[0].size() != mixture_proportions_[0].size())
    {
        throw GelexException(
            "SingleScaledMixtureGaussianPrior: multiplier and proportion sizes "
            "differ");
    }
    if (multipliers_[0](0) != 0.0)
    {
        throw GelexException(
            "SingleScaledMixtureGaussianPrior: multiplier(0) must equal 0");
    }
}

auto SingleScaledMixtureGaussianPrior::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::source | FieldFlag::metadata);
    marker_variances_[0].visit(visitor);
    visitor.on("multiplier", multipliers_[0], FieldFlag::source);
    mixture_proportions_[0].visit(visitor);
}

auto SingleScaledMixtureGaussianPrior::make_state(
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

JointGaussianMixturePrior::JointGaussianMixturePrior(
    std::array<MarkerVariance, 2> variances,
    MixtureProportion proportion)
    : marker_variances_(variances), mixture_proportions_{std::move(proportion)}
{
}

auto JointGaussianMixturePrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    auto modes = std::array{GeneticMode::A, GeneticMode::D};
    for (auto [i, mode] : std::views::enumerate(modes))
    {
        auto slot_scope = visitor.scope(std::to_string(i));
        visitor.on("mode", mode, FieldFlag::source | FieldFlag::metadata);
        marker_variances_[i].visit(visitor);
    }
    mixture_proportions_[0].visit(visitor);
}

auto JointGaussianMixturePrior::make_state(
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
