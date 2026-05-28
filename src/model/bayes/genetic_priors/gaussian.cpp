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
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/genetic_prior_states/gaussian.h"
#include "gelex/model/bayes/prior_parameters.h"

namespace gelex::bayes
{

SingleSharedGaussianPrior::SingleSharedGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance)
    : mode_(mode), variance_(std::move(variance))
{
}

auto SingleSharedGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
}

auto SingleSharedGaussianPrior::make_state(
    Eigen::Index /*num_markers*/,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleSharedGaussianState>(
        variance_.parameter().initial_value());
}

SinglePerMarkerGaussianPrior::SinglePerMarkerGaussianPrior(
    GeneticMode mode,
    PerMarkerVariance variance)
    : mode_(mode), variance_(std::move(variance))
{
}

auto SinglePerMarkerGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
}

auto SinglePerMarkerGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SinglePerMarkerGaussianState>(
        Eigen::VectorXd::Constant(
            num_markers, variance_.parameter().initial_value()));
}

SingleSharedSpikeSlabGaussianPrior::SingleSharedSpikeSlabGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance,
    MixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
    if (mixture_proportion_.size() != 2)
    {
        throw GelexException(
            "SingleSharedSpikeSlabGaussianPrior: proportion must have size 2");
    }
}

auto SingleSharedSpikeSlabGaussianPrior::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto SingleSharedSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleSharedSpikeSlabGaussianState>(
        variance_.parameter().initial_value(),
        mixture_proportion_,
        num_markers);
}

SinglePerMarkerSpikeSlabGaussianPrior::SinglePerMarkerSpikeSlabGaussianPrior(
    GeneticMode mode,
    PerMarkerVariance variance,
    MixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
    if (mixture_proportion_.size() != 2)
    {
        throw GelexException(
            "SinglePerMarkerSpikeSlabGaussianPrior: proportion must have size "
            "2");
    }
}

auto SinglePerMarkerSpikeSlabGaussianPrior::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto SinglePerMarkerSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SinglePerMarkerSpikeSlabGaussianState>(
        Eigen::VectorXd::Constant(
            num_markers, variance_.parameter().initial_value()),
        mixture_proportion_,
        num_markers);
}

SingleScaledMixtureGaussianPrior::SingleScaledMixtureGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance,
    Eigen::VectorXd multiplier,
    MixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
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
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    visitor.on(
        "multiplier", multiplier_, FieldFlag::checkpoint | FieldFlag::report);
    mixture_proportion_.visit(visitor);
}

auto SingleScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleScaledMixtureGaussianState>(
        variance_.parameter().initial_value(),
        multiplier_,
        mixture_proportion_,
        num_markers,
        num_individuals);
}

JointGaussianMixturePrior::JointGaussianMixturePrior(
    std::array<JointMarkerVariance, 2> variances,
    MixtureProportion proportion)
    : marker_variances_(std::move(variances)),
      mixture_proportion_(std::move(proportion))
{
}

auto JointGaussianMixturePrior::variance(GeneticMode mode)
    -> JointMarkerVariance&
{
    return marker_variances_[std::to_underlying(mode)];
}

auto JointGaussianMixturePrior::variance(GeneticMode mode) const
    -> const JointMarkerVariance&
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
        visitor.on("mode", mode, FieldFlag::checkpoint | FieldFlag::report);
        std::visit(
            [&visitor](auto& marker_variance)
            { marker_variance.visit(visitor); },
            marker_variances_[i]);
    }
    mixture_proportion_.visit(visitor);
}

auto JointGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> std::unique_ptr<JointGeneticPriorState>
{
    return std::make_unique<JointGaussianMixtureState>(
        std::array<JointMarkerVarianceState, 2>{
            std::visit(
                [num_markers](
                    const auto& marker_variance) -> JointMarkerVarianceState
                {
                    using MarkerVarianceType
                        = std::decay_t<decltype(marker_variance)>;
                    if constexpr (
                        std::
                            is_same_v<MarkerVarianceType, SharedMarkerVariance>)
                    {
                        return marker_variance.parameter().initial_value();
                    }
                    else
                    {
                        return Eigen::VectorXd::Constant(
                            marker_variance.marker_variance_size(num_markers),
                            marker_variance.parameter().initial_value());
                    }
                },
                marker_variances_[0]),
            std::visit(
                [num_markers](
                    const auto& marker_variance) -> JointMarkerVarianceState
                {
                    using MarkerVarianceType
                        = std::decay_t<decltype(marker_variance)>;
                    if constexpr (
                        std::
                            is_same_v<MarkerVarianceType, SharedMarkerVariance>)
                    {
                        return marker_variance.parameter().initial_value();
                    }
                    else
                    {
                        return Eigen::VectorXd::Constant(
                            marker_variance.marker_variance_size(num_markers),
                            marker_variance.parameter().initial_value());
                    }
                },
                marker_variances_[1])},
        mixture_proportion_,
        num_markers,
        num_individuals);
}

}  // namespace gelex::bayes
