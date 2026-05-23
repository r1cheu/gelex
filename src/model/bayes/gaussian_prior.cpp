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
#include <utility>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/gaussian_prior_state.h"

namespace gelex::bayes
{

GaussianPrior::GaussianPrior(GeneticMode mode, MarkerVarianceSpec variance)
    : modes_{mode}, variance_specs_{variance}
{
}

auto GaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<GaussianState>(
        make_variance_values(variance_specs_, num_markers));
}

SpikeSlabGaussianPrior::SpikeSlabGaussianPrior(
    GeneticMode mode,
    MarkerVarianceSpec variance,
    ProportionSpec proportion)
    : modes_{mode},
      variance_specs_{variance},
      proportion_specs_{std::move(proportion)}
{
    if (proportion_specs_[0].size() != 2)
    {
        throw GelexException(
            "SpikeSlabGaussianPrior: proportion must have size 2");
    }
}

auto SpikeSlabGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index /*num_individuals*/) const
    -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<SpikeSlabGaussianState>(
        make_variance_values(variance_specs_, num_markers),
        proportion_specs_,
        num_markers);
}

ScaledMixtureGaussianPrior::ScaledMixtureGaussianPrior(
    GeneticMode mode,
    MarkerVarianceSpec variance,
    Eigen::VectorXd multiplier,
    ProportionSpec proportion)
    : modes_{mode},
      variance_specs_{variance},
      multipliers_{std::move(multiplier)},
      proportion_specs_{std::move(proportion)}
{
    if (multipliers_[0].size() != proportion_specs_[0].size())
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

auto ScaledMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> std::unique_ptr<GeneticPriorState>
{
    return std::make_unique<ScaledMixtureGaussianState>(
        make_variance_values(variance_specs_, num_markers),
        multipliers_,
        proportion_specs_,
        num_markers,
        num_individuals);
}

JointMixtureGaussianPrior::JointMixtureGaussianPrior(
    std::array<GeneticMode, 2> modes,
    std::array<MarkerVarianceSpec, 2> variances,
    ProportionSpec proportion)
    : modes_(modes),
      variance_specs_(variances),
      proportion_specs_{std::move(proportion)}
{
    if (modes_[0] == modes_[1])
    {
        throw GelexException(
            "JointMixtureGaussianPrior: modes must be distinct");
    }
}

auto JointMixtureGaussianPrior::make_state(
    Eigen::Index num_markers,
    Eigen::Index num_individuals) const -> std::unique_ptr<GeneticPriorState>
{
    auto values = make_variance_values(variance_specs_, num_markers);
    return std::make_unique<JointMixtureGaussianState>(
        std::array<Eigen::VectorXd, 2>{
            std::move(values[0]), std::move(values[1])},
        proportion_specs_[0],
        num_markers,
        num_individuals);
}

}  // namespace gelex::bayes
