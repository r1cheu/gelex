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

#include "gelex/model/bayes/samplers/detail/pi.h"

#include <random>

#include <Eigen/Core>

#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/model.h"

namespace gelex::detail
{
using Eigen::VectorXd;
using Eigen::VectorXi;

namespace AdditiveSampler
{
auto Pi::operator()(
    const BayesModel& /*model*/,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    // Check if the model has additive effects with pi estimation

    if (auto* state = states.additive(); state != nullptr)
    {
        VectorXi dirichlet_counts(state->pi.count.array() + 1);

        // Sample from Dirichlet distribution
        state->pi.prop = detail::dirichlet(dirichlet_counts, rng);
    }
}

}  // namespace AdditiveSampler
//
auto DominantSampler::Pi::operator()(
    const BayesModel& /*model*/,
    BayesState& states,
    std::mt19937_64& rng) const -> void
{
    // Check if the model has dominant effects with pi estimation

    if (auto* state = states.dominant(); state != nullptr)
    {
        VectorXi dirichlet_counts(state->pi.count.array() + 1);

        // Sample from Dirichlet distribution
        state->pi.prop = detail::dirichlet(dirichlet_counts, rng);
    }
}
}  // namespace gelex::detail
