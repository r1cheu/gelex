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
#include <utility>

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
    JointSharedMarkerVariance variance,
    MixtureProportion proportion)
    : marker_variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
}

auto JointGaussianMixturePrior::variance(GeneticMode mode)
    -> SharedMarkerVariance&
{
    return marker_variance_.variance(mode);
}

auto JointGaussianMixturePrior::variance(GeneticMode mode) const
    -> const SharedMarkerVariance&
{
    return marker_variance_.variance(mode);
}

auto JointGaussianMixturePrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    marker_variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto JointGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> std::unique_ptr<JointGeneticPriorState>
{
    return std::make_unique<JointGaussianMixtureState>(
        std::array{
            marker_variance_.variance(GeneticMode::A)
                .parameter()
                .initial_value(),
            marker_variance_.variance(GeneticMode::D)
                .parameter()
                .initial_value()},
        mixture_proportion_,
        num_markers,
        num_individuals);
}

}  // namespace gelex::bayes
