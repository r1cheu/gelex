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

#include "gelex/model/bayes/genetic_prior.h"

#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_specs.h"

namespace gelex::bayes
{

auto GeneticPrior::make_variance_values(
    std::span<const MarkerVarianceSpec> variance_specs,
    Eigen::Index num_markers) -> std::vector<Eigen::VectorXd>
{
    std::vector<Eigen::VectorXd> values;
    values.reserve(variance_specs.size());
    for (const auto& spec : variance_specs)
    {
        values.emplace_back(
            Eigen::VectorXd::Constant(
                spec.marker_variance_size(num_markers),
                spec.variance().initial_value()));
    }
    return values;
}

}  // namespace gelex::bayes
