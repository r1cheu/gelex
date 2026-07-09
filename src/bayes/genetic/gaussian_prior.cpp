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

#include <Eigen/Core>
#include <array>
#include <utility>

#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes
{

SingleSharedGaussianPrior::SingleSharedGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance)
    : variance_(variance), mode_(mode)
{
}

auto SingleSharedGaussianPrior::visit(FieldVisitor& visitor) -> void
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
    : variance_(variance), mode_(mode)
{
}

auto SinglePerMarkerGaussianPrior::visit(FieldVisitor& visitor) -> void
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

SingleSharedSpikeSlabGaussianPrior::SingleSharedSpikeSlabGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance,
    MixtureProportion proportion)
    : variance_(variance), proportion_(std::move(proportion)), mode_(mode)
{
    if (this->proportion().size() != 2)
    {
        throw GelexException(
            "SingleSharedSpikeSlabGaussianPrior: proportion must have "
            "size 2");
    }
}

auto SingleSharedSpikeSlabGaussianPrior::visit(FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance().visit(visitor);
    proportion().visit(visitor);
}

auto SingleSharedSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> SingleSharedSpikeSlabGaussianState
{
    return SingleSharedSpikeSlabGaussianState{
        variance().parameter().initial_value(), proportion(), num_markers};
}

SinglePerMarkerSpikeSlabGaussianPrior::SinglePerMarkerSpikeSlabGaussianPrior(
    GeneticMode mode,
    PerMarkerVariance variance,
    MixtureProportion proportion)
    : variance_(variance), proportion_(std::move(proportion)), mode_(mode)
{
    if (this->proportion().size() != 2)
    {
        throw GelexException(
            "SinglePerMarkerSpikeSlabGaussianPrior: proportion must have "
            "size 2");
    }
}

auto SinglePerMarkerSpikeSlabGaussianPrior::visit(FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("mode", mode_, FieldFlag::checkpoint | FieldFlag::report);
    variance().visit(visitor);
    proportion().visit(visitor);
}

auto SinglePerMarkerSpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> SinglePerMarkerSpikeSlabGaussianState
{
    return SinglePerMarkerSpikeSlabGaussianState{
        Eigen::VectorXd::Constant(
            num_markers, variance().parameter().initial_value()),
        proportion(),
        num_markers};
}

SingleScaledMixtureGaussianPrior::SingleScaledMixtureGaussianPrior(
    GeneticMode mode,
    SharedMarkerVariance variance,
    Eigen::VectorXd multiplier,
    MixtureProportion proportion)
    : variance_(variance),
      proportion_(std::move(proportion)),
      multiplier_(std::move(multiplier)),
      mode_(mode)
{
    if (this->multiplier().size() != this->proportion().size())
    {
        throw GelexException(
            "SingleScaledMixtureGaussianPrior: multiplier and proportion "
            "sizes differ");
    }
    if (this->multiplier()(0) != 0.0)
    {
        throw GelexException(
            "SingleScaledMixtureGaussianPrior: multiplier(0) must equal "
            "0");
    }
}

auto SingleScaledMixtureGaussianPrior::visit(FieldVisitor& visitor) -> void
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

auto SingleScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> SingleScaledMixtureGaussianState
{
    return SingleScaledMixtureGaussianState{
        variance().parameter().initial_value(),
        multiplier(),
        proportion(),
        num_markers,
        num_individuals};
}

JointGaussianMixturePrior::JointGaussianMixturePrior(
    JointSharedMarkerVariance variance,
    MixtureProportion proportion)
    : variances_(variance), proportion_(std::move(proportion))
{
}

auto JointGaussianMixturePrior::variance(GeneticMode mode)
    -> SharedMarkerVariance&
{
    return variances_.variance(mode);
}

auto JointGaussianMixturePrior::variance(GeneticMode mode) const
    -> const SharedMarkerVariance&
{
    return variances_.variance(mode);
}

auto JointGaussianMixturePrior::visit(FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    variances().visit(visitor);
    proportion().visit(visitor);
}

auto JointGaussianMixturePrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> JointGaussianMixtureState
{
    return JointGaussianMixtureState{
        std::array{
            variance(GeneticMode::A).parameter().initial_value(),
            variance(GeneticMode::D).parameter().initial_value()},
        proportion(),
        num_markers,
        num_individuals};
}

}  // namespace gelex::bayes
