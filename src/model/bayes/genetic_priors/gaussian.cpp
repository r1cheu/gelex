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

auto make_variance_value(
    const MarkerVariance& marker_variance,
    Eigen::Index num_markers) -> Eigen::VectorXd
{
    return Eigen::VectorXd::Constant(
        marker_variance.marker_variance_size(num_markers),
        marker_variance.parameter().initial_value());
}

}  // namespace

SingleGaussianPrior::SingleGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance)
    : mode_(mode), marker_variance_(variance)
{
}

auto SingleGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::source | FieldFlag::metadata);
    marker_variance_.visit(visitor);
}

auto SingleGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<GaussianState>(std::vector<Eigen::VectorXd>{
        make_variance_value(marker_variance_, num_markers)});
}

SingleSpikeSlabGaussianPrior::SingleSpikeSlabGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance,
    MixtureProportion proportion)
    : mode_(mode),
      marker_variance_(variance),
      mixture_proportion_(std::move(proportion))
{
    if (mixture_proportion_.size() != 2)
    {
        throw GelexException(
            "SingleSpikeSlabGaussianPrior: proportion must have size 2");
    }
}

auto SingleSpikeSlabGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::source | FieldFlag::metadata);
    marker_variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto SingleSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<SpikeSlabGaussianState>(
        std::vector<Eigen::VectorXd>{
            make_variance_value(marker_variance_, num_markers)},
        std::span<const MixtureProportion>{&mixture_proportion_, 1},
        num_markers);
}

SingleScaledMixtureGaussianPrior::SingleScaledMixtureGaussianPrior(
    GeneticMode mode,
    MarkerVariance variance,
    Eigen::VectorXd multiplier,
    MixtureProportion proportion)
    : mode_(mode),
      marker_variance_(variance),
      multiplier_(std::move(multiplier)),
      mixture_proportion_(std::move(proportion))
{
    if (multiplier_.size() != mixture_proportion_.size())
    {
        throw GelexException(
            "SingleScaledMixtureGaussianPrior: multiplier and proportion sizes "
            "differ");
    }
    if (multiplier_(0) != 0.0)
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
    marker_variance_.visit(visitor);
    visitor.on("multiplier", multiplier_, FieldFlag::source);
    mixture_proportion_.visit(visitor);
}

auto SingleScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<ScaledMixtureGaussianState>(
        std::vector<Eigen::VectorXd>{
            make_variance_value(marker_variance_, num_markers)},
        std::span<const Eigen::VectorXd>{&multiplier_, 1},
        std::span<const MixtureProportion>{&mixture_proportion_, 1},
        num_markers,
        num_individuals);
}

JointGaussianMixturePrior::JointGaussianMixturePrior(
    std::array<MarkerVariance, 2> variances,
    MixtureProportion proportion)
    : marker_variances_(variances), mixture_proportion_(std::move(proportion))
{
}

auto JointGaussianMixturePrior::variance(GeneticMode mode) -> MarkerVariance&
{
    return marker_variances_[std::to_underlying(mode)];
}

auto JointGaussianMixturePrior::variance(GeneticMode mode) const
    -> const MarkerVariance&
{
    return marker_variances_[std::to_underlying(mode)];
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
    mixture_proportion_.visit(visitor);
}

auto JointGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<JointMixtureGaussianState>(
        std::array<Eigen::VectorXd, 2>{
            make_variance_value(marker_variances_[0], num_markers),
            make_variance_value(marker_variances_[1], num_markers)},
        mixture_proportion_,
        num_markers,
        num_individuals);
}

}  // namespace gelex::bayes
