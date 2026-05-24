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

GaussianPrior::GaussianPrior(GeneticMode mode, MarkerVariance variance)
    : modes_{mode}, marker_variances_{variance}
{
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
