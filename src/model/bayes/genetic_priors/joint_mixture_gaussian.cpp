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

#include "gelex/model/bayes/genetic_priors/joint_mixture_gaussian.h"

#include <cassert>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/runtime_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

JointMixtureGaussianPrior::JointMixtureGaussianPrior(
    std::array<GeneticMode, 2> modes,
    std::array<MarkerVarianceSpec, 2> variances,
    ProportionSpec proportion)
    : modes_(modes),
      variance_specs_(variances),
      proportion_spec_(std::move(proportion))
{
    if (modes_[0] == modes_[1])
    {
        throw GelexException(
            "JointMixtureGaussianPrior: modes must be distinct");
    }
}

auto JointMixtureGaussianPrior::modes() const -> std::span<const GeneticMode>
{
    return modes_;
}

auto JointMixtureGaussianPrior::variance_specs() const
    -> std::span<const MarkerVarianceSpec>
{
    return variance_specs_;
}

auto JointMixtureGaussianPrior::proportion_spec() const -> const ProportionSpec&
{
    return proportion_spec_;
}

auto JointMixtureGaussianPrior::make_state(const GeneticPriorRuntimeInit& init)
    const -> std::unique_ptr<GeneticPriorRuntimeState>
{
    assert(init.effects.size() == modes_.size());
    assert(init.effects[0].mode == modes_[0]);
    assert(init.effects[1].mode == modes_[1]);
    const auto num_markers = init.effects[0].num_markers;
    assert(num_markers > 0);
    assert(init.effects[1].num_markers == num_markers);

    std::vector<Eigen::VectorXd> variances;
    variances.reserve(variance_specs_.size());
    for (const auto& spec : variance_specs_)
    {
        variances.emplace_back(
            Eigen::VectorXd::Constant(
                spec.marker_variance_size(num_markers),
                spec.variance().initial_value()));
    }
    return std::make_unique<MixtureGaussianRuntimeState>(
        MarkerVarianceRuntimeState(std::move(variances)),
        MixtureRuntimeState(
            Eigen::VectorXi::Zero(num_markers),
            proportion_spec_.initial_value().to_mat()));
}

}  // namespace gelex::bayes
