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

#include "gelex/model/bayes/genetic_priors/mixture_gaussian.h"

#include <cassert>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/genetic_prior_runtime_state.h"
#include "gelex/model/bayes/prior_specs.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

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

auto SpikeSlabGaussianPrior::modes() const -> std::span<const GeneticMode>
{
    return modes_;
}

auto SpikeSlabGaussianPrior::variance_specs() const
    -> std::span<const MarkerVarianceSpec>
{
    return variance_specs_;
}

auto SpikeSlabGaussianPrior::proportion_specs() const
    -> std::span<const ProportionSpec>
{
    return proportion_specs_;
}

auto SpikeSlabGaussianPrior::make_state(const GeneticPriorRuntimeInit& init)
    const -> std::unique_ptr<GeneticPriorRuntimeState>
{
    assert(init.effects.size() == modes_.size());
    assert(init.effects[0].mode == modes_[0]);
    const auto num_markers = init.effects[0].num_markers;
    assert(num_markers > 0);

    const auto& var_spec = variance_specs_[0];
    std::vector<Eigen::VectorXd> variances;
    variances.emplace_back(
        Eigen::VectorXd::Constant(
            var_spec.marker_variance_size(num_markers),
            var_spec.variance().initial_value()));
    return std::make_unique<MixtureGaussianRuntimeState>(
        MarkerVarianceRuntimeState(std::move(variances)),
        MixtureRuntimeState(
            Eigen::VectorXi::Zero(num_markers),
            proportion_specs_[0].initial_value().to_mat()));
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
    if (multipliers_[0].size()
        != static_cast<Eigen::Index>(proportion_specs_[0].size()))
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

auto ScaledMixtureGaussianPrior::modes() const -> std::span<const GeneticMode>
{
    return modes_;
}

auto ScaledMixtureGaussianPrior::variance_specs() const
    -> std::span<const MarkerVarianceSpec>
{
    return variance_specs_;
}

auto ScaledMixtureGaussianPrior::proportion_specs() const
    -> std::span<const ProportionSpec>
{
    return proportion_specs_;
}

auto ScaledMixtureGaussianPrior::multipliers() const
    -> std::span<const Eigen::VectorXd>
{
    return multipliers_;
}

auto ScaledMixtureGaussianPrior::make_state(const GeneticPriorRuntimeInit& init)
    const -> std::unique_ptr<GeneticPriorRuntimeState>
{
    assert(init.effects.size() == modes_.size());
    assert(init.effects[0].mode == modes_[0]);
    const auto num_markers = init.effects[0].num_markers;
    assert(num_markers > 0);

    const auto& var_spec = variance_specs_[0];
    std::vector<Eigen::VectorXd> variances;
    variances.emplace_back(
        Eigen::VectorXd::Constant(
            var_spec.marker_variance_size(num_markers),
            var_spec.variance().initial_value()));
    return std::make_unique<MixtureGaussianRuntimeState>(
        MarkerVarianceRuntimeState(std::move(variances)),
        MixtureRuntimeState(
            Eigen::VectorXi::Zero(num_markers),
            proportion_specs_[0].initial_value().to_mat()));
}

}  // namespace gelex::bayes
