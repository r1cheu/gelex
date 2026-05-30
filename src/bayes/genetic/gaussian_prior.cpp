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

#include "gelex/bayes/genetic/gaussian_prior.h"

#include <array>
#include <utility>

#include <Eigen/Core>

#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

SingleSharedGaussianPrior::SingleSharedGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance)
    : Variance<SharedMarkerVariance>(variance), mode_(mode)
{
}

auto SingleSharedGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance().visit(visitor);
}

auto SingleSharedGaussianPrior::make_state(
    Eigen::Index /*num_markers*/,
    Eigen::Index /*num_individuals*/) const -> SingleSharedGaussianState
{
    return SingleSharedGaussianState{variance().parameter().initial_value()};
}

SinglePerMarkerGaussianPrior::SinglePerMarkerGaussianPrior(
    GeneticMode mode,
    PerMarkerVariance variance)
    : Variance<PerMarkerVariance>(variance), mode_(mode)
{
}

auto SinglePerMarkerGaussianPrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance().visit(visitor);
}

auto SinglePerMarkerGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const -> SinglePerMarkerGaussianState
{
    return SinglePerMarkerGaussianState{Eigen::VectorXd::Constant(
        num_markers, variance().parameter().initial_value())};
}

SingleFixedSharedSpikeSlabGaussianPrior::
    SingleFixedSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        FixedProportion proportion)
    : Variance<SharedMarkerVariance>(variance),
      Proportion<FixedProportion>(std::move(proportion)),
      mode_(mode)
{
    if (this->proportion().size() != 2)
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
    variance().visit(visitor);
    proportion().visit(visitor);
}

auto SingleFixedSharedSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> SingleFixedSharedSpikeSlabGaussianState
{
    return SingleFixedSharedSpikeSlabGaussianState{
        variance().parameter().initial_value(),
        proportion().size(),
        num_markers};
}

SingleSampledSharedSpikeSlabGaussianPrior::
    SingleSampledSharedSpikeSlabGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        SampledProportion proportion)
    : Variance<SharedMarkerVariance>(variance),
      Proportion<SampledProportion>(std::move(proportion)),
      mode_(mode)
{
    if (this->proportion().size() != 2)
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
    variance().visit(visitor);
    proportion().visit(visitor);
}

auto SingleSampledSharedSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> SingleSampledSharedSpikeSlabGaussianState
{
    return SingleSampledSharedSpikeSlabGaussianState{
        variance().parameter().initial_value(), proportion(), num_markers};
}

SingleFixedPerMarkerSpikeSlabGaussianPrior::
    SingleFixedPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        FixedProportion proportion)
    : Variance<PerMarkerVariance>(variance),
      Proportion<FixedProportion>(std::move(proportion)),
      mode_(mode)
{
    if (this->proportion().size() != 2)
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
    variance().visit(visitor);
    proportion().visit(visitor);
}

auto SingleFixedPerMarkerSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> SingleFixedPerMarkerSpikeSlabGaussianState
{
    return SingleFixedPerMarkerSpikeSlabGaussianState{
        Eigen::VectorXd::Constant(
            num_markers, variance().parameter().initial_value()),
        proportion().size(),
        num_markers};
}

SingleSampledPerMarkerSpikeSlabGaussianPrior::
    SingleSampledPerMarkerSpikeSlabGaussianPrior(
        GeneticMode mode,
        PerMarkerVariance variance,
        SampledProportion proportion)
    : Variance<PerMarkerVariance>(variance),
      Proportion<SampledProportion>(std::move(proportion)),
      mode_(mode)
{
    if (this->proportion().size() != 2)
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
    variance().visit(visitor);
    proportion().visit(visitor);
}

auto SingleSampledPerMarkerSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> SingleSampledPerMarkerSpikeSlabGaussianState
{
    return SingleSampledPerMarkerSpikeSlabGaussianState{
        Eigen::VectorXd::Constant(
            num_markers, variance().parameter().initial_value()),
        proportion(),
        num_markers};
}

SingleFixedScaledMixtureGaussianPrior::SingleFixedScaledMixtureGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance,
    Eigen::VectorXd multiplier,
    FixedProportion proportion)
    : Variance<SharedMarkerVariance>(variance),
      Proportion<FixedProportion>(std::move(proportion)),
      Multiplier<Eigen::VectorXd>(std::move(multiplier)),
      mode_(mode)
{
    if (this->multiplier().size() != this->proportion().size())
    {
        throw GelexException(
            "SingleFixedScaledMixtureGaussianPrior: multiplier and proportion "
            "sizes differ");
    }
    if (this->multiplier()(0) != 0.0)
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
    variance().visit(visitor);
    visitor.on(
        "multiplier",
        this->multiplier(),
        FieldFlag::checkpoint | FieldFlag::report);
    proportion().visit(visitor);
}

auto SingleFixedScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> SingleFixedScaledMixtureGaussianState
{
    return SingleFixedScaledMixtureGaussianState{
        variance().parameter().initial_value(),
        multiplier(),
        num_markers,
        num_individuals};
}

SingleSampledScaledMixtureGaussianPrior::
    SingleSampledScaledMixtureGaussianPrior(
        GeneticMode mode,
        SharedMarkerVariance variance,
        Eigen::VectorXd multiplier,
        SampledProportion proportion)
    : Variance<SharedMarkerVariance>(variance),
      Proportion<SampledProportion>(std::move(proportion)),
      Multiplier<Eigen::VectorXd>(std::move(multiplier)),
      mode_(mode)
{
    if (this->multiplier().size() != this->proportion().size())
    {
        throw GelexException(
            "SingleSampledScaledMixtureGaussianPrior: multiplier and "
            "proportion sizes differ");
    }
    if (this->multiplier()(0) != 0.0)
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
    variance().visit(visitor);
    visitor.on(
        "multiplier",
        this->multiplier(),
        FieldFlag::checkpoint | FieldFlag::report);
    proportion().visit(visitor);
}

auto SingleSampledScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const
    -> SingleSampledScaledMixtureGaussianState
{
    return SingleSampledScaledMixtureGaussianState{
        variance().parameter().initial_value(),
        multiplier(),
        proportion(),
        num_markers,
        num_individuals};
}

JointFixedGaussianMixturePrior::JointFixedGaussianMixturePrior(
    JointSharedMarkerVariance variance,
    FixedProportion proportion)
    : JointVariancesField<JointSharedMarkerVariance>(variance),
      Proportion<FixedProportion>(std::move(proportion))
{
}

auto JointFixedGaussianMixturePrior::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    variances().visit(visitor);
    proportion().visit(visitor);
}

auto JointFixedGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> JointFixedGaussianMixtureState
{
    return JointFixedGaussianMixtureState{
        std::array{
            variance(GeneticMode::A).parameter().initial_value(),
            variance(GeneticMode::D).parameter().initial_value()},
        proportion().size(),
        num_markers,
        num_individuals};
}

JointSampledGaussianMixturePrior::JointSampledGaussianMixturePrior(
    JointSharedMarkerVariance variance,
    SampledProportion proportion)
    : JointVariancesField<JointSharedMarkerVariance>(variance),
      Proportion<SampledProportion>(std::move(proportion))
{
}

auto JointSampledGaussianMixturePrior::visit(infra::FieldVisitor& visitor)
    -> void
{
    auto scope = visitor.scope(name);
    variances().visit(visitor);
    proportion().visit(visitor);
}

auto JointSampledGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> JointSampledGaussianMixtureState
{
    return JointSampledGaussianMixtureState{
        std::array{
            variance(GeneticMode::A).parameter().initial_value(),
            variance(GeneticMode::D).parameter().initial_value()},
        proportion(),
        num_markers,
        num_individuals};
}

}  // namespace gelex::bayes
