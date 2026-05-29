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
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

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

SingleFixedSharedSpikeSlabGaussianPrior::
    SingleFixedSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        FixedMixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
    if (mixture_proportion_.size() != 2)
    {
        throw GelexException(
            "SingleFixedSharedSpikeSlabGaussianPrior: proportion must have "
            "size 2");
    }
}

auto SingleFixedSharedSpikeSlabGaussianPrior::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto SingleFixedSharedSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleFixedSharedSpikeSlabGaussianState>(
        variance_.parameter().initial_value(),
        mixture_proportion_.size(),
        num_markers);
}

SingleSampledSharedSpikeSlabGaussianPrior::
    SingleSampledSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        SampledMixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
    if (mixture_proportion_.size() != 2)
    {
        throw GelexException(
            "SingleSampledSharedSpikeSlabGaussianPrior: proportion must have "
            "size 2");
    }
}

auto SingleSampledSharedSpikeSlabGaussianPrior::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto SingleSampledSharedSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleSampledSharedSpikeSlabGaussianState>(
        variance_.parameter().initial_value(),
        mixture_proportion_,
        num_markers);
}

SingleFixedPerMarkerSpikeSlabGaussianPrior::
    SingleFixedPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        FixedMixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
    if (mixture_proportion_.size() != 2)
    {
        throw GelexException(
            "SingleFixedPerMarkerSpikeSlabGaussianPrior: proportion must have "
            "size 2");
    }
}

auto SingleFixedPerMarkerSpikeSlabGaussianPrior::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto SingleFixedPerMarkerSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleFixedPerMarkerSpikeSlabGaussianState>(
        Eigen::VectorXd::Constant(
            num_markers, variance_.parameter().initial_value()),
        mixture_proportion_.size(),
        num_markers);
}

SingleSampledPerMarkerSpikeSlabGaussianPrior::
    SingleSampledPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        SampledMixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
    if (mixture_proportion_.size() != 2)
    {
        throw GelexException(
            "SingleSampledPerMarkerSpikeSlabGaussianPrior: proportion must "
            "have size 2");
    }
}

auto SingleSampledPerMarkerSpikeSlabGaussianPrior::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto SingleSampledPerMarkerSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleSampledPerMarkerSpikeSlabGaussianState>(
        Eigen::VectorXd::Constant(
            num_markers, variance_.parameter().initial_value()),
        mixture_proportion_,
        num_markers);
}

SingleFixedScaledMixtureGaussianPrior::SingleFixedScaledMixtureGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance,
    Eigen::VectorXd multiplier,
    FixedMixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      multiplier_(std::move(multiplier)),
      mixture_proportion_(std::move(proportion))
{
    if (multiplier_.size() != mixture_proportion_.size())
    {
        throw GelexException(
            "SingleFixedScaledMixtureGaussianPrior: multiplier and proportion "
            "sizes differ");
    }
    if (multiplier_(0) != 0.0)
    {
        throw GelexException(
            "SingleFixedScaledMixtureGaussianPrior: multiplier(0) must equal "
            "0");
    }
}

auto SingleFixedScaledMixtureGaussianPrior::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    visitor.on(
        "multiplier", multiplier_, FieldFlag::checkpoint | FieldFlag::report);
    mixture_proportion_.visit(visitor);
}

auto SingleFixedScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleFixedScaledMixtureGaussianState>(
        variance_.parameter().initial_value(),
        multiplier_,
        num_markers,
        num_individuals);
}

SingleSampledScaledMixtureGaussianPrior::
    SingleSampledScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        SampledMixtureProportion proportion)
    : mode_(mode),
      variance_(std::move(variance)),
      multiplier_(std::move(multiplier)),
      mixture_proportion_(std::move(proportion))
{
    if (multiplier_.size() != mixture_proportion_.size())
    {
        throw GelexException(
            "SingleSampledScaledMixtureGaussianPrior: multiplier and "
            "proportion sizes differ");
    }
    if (multiplier_(0) != 0.0)
    {
        throw GelexException(
            "SingleSampledScaledMixtureGaussianPrior: multiplier(0) must "
            "equal 0");
    }
}

auto SingleSampledScaledMixtureGaussianPrior::visit(
    infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance_.visit(visitor);
    visitor.on(
        "multiplier", multiplier_, FieldFlag::checkpoint | FieldFlag::report);
    mixture_proportion_.visit(visitor);
}

auto SingleSampledScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> std::unique_ptr<SingleGeneticPriorState>
{
    return std::make_unique<SingleSampledScaledMixtureGaussianState>(
        variance_.parameter().initial_value(),
        multiplier_,
        mixture_proportion_,
        num_markers,
        num_individuals);
}

JointFixedGaussianMixturePrior::JointFixedGaussianMixturePrior(
    JointSharedMarkerVariance variance,
    FixedMixtureProportion proportion)
    : marker_variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
}

auto JointFixedGaussianMixturePrior::variance(GeneticMode mode)
    -> SharedMarkerVariance&
{
    return marker_variance_.variance(mode);
}

auto JointFixedGaussianMixturePrior::variance(GeneticMode mode) const
    -> const SharedMarkerVariance&
{
    return marker_variance_.variance(mode);
}

auto JointFixedGaussianMixturePrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    marker_variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto JointFixedGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> std::unique_ptr<JointGeneticPriorState>
{
    return std::make_unique<JointFixedGaussianMixtureState>(
        std::array{
            marker_variance_.variance(GeneticMode::A)
                .parameter()
                .initial_value(),
            marker_variance_.variance(GeneticMode::D)
                .parameter()
                .initial_value()},
        mixture_proportion_.size(),
        num_markers,
        num_individuals);
}

JointSampledGaussianMixturePrior::JointSampledGaussianMixturePrior(
    JointSharedMarkerVariance variance,
    SampledMixtureProportion proportion)
    : marker_variance_(std::move(variance)),
      mixture_proportion_(std::move(proportion))
{
}

auto JointSampledGaussianMixturePrior::variance(GeneticMode mode)
    -> SharedMarkerVariance&
{
    return marker_variance_.variance(mode);
}

auto JointSampledGaussianMixturePrior::variance(GeneticMode mode) const
    -> const SharedMarkerVariance&
{
    return marker_variance_.variance(mode);
}

auto JointSampledGaussianMixturePrior::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    marker_variance_.visit(visitor);
    mixture_proportion_.visit(visitor);
}

auto JointSampledGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> std::unique_ptr<JointGeneticPriorState>
{
    return std::make_unique<JointSampledGaussianMixtureState>(
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
