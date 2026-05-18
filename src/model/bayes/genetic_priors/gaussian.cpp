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

#include <cassert>
#include <memory>
#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/genetic_prior_runtime_state.h"

namespace gelex::bayes
{

GaussianPrior::GaussianPrior(GeneticMode mode, MarkerVarianceSpec variance)
    : modes_{mode}, variance_specs_{variance}
{
}

auto GaussianPrior::modes() const -> std::span<const GeneticMode>
{
    return modes_;
}

auto GaussianPrior::variance_specs() const
    -> std::span<const MarkerVarianceSpec>
{
    return variance_specs_;
}

auto GaussianPrior::make_state(const GeneticPriorRuntimeInit& init) const
    -> std::unique_ptr<GeneticPriorRuntimeState>
{
    assert(init.effects.size() == modes_.size());
    assert(init.effects[0].mode == modes_[0]);
    const auto num_markers = init.effects[0].num_markers;
    assert(num_markers > 0);

    const auto& spec = variance_specs_[0];
    std::vector<Eigen::VectorXd> variances;
    variances.emplace_back(
        Eigen::VectorXd::Constant(
            spec.marker_variance_size(num_markers),
            spec.variance().initial_value()));
    return std::make_unique<GaussianRuntimeState>(
        MarkerVarianceRuntimeState(std::move(variances)));
}

}  // namespace gelex::bayes
