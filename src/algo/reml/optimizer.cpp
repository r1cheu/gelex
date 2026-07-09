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

#include "gelex/algo/reml/optimizer.h"

#include <Eigen/Core>

#include "gelex/freq/model.h"

namespace gelex
{

auto collect_variance_components(const gelex::FreqState& state)
    -> Eigen::VectorXd
{
    auto n_random = static_cast<Eigen::Index>(state.random().size());
    Eigen::Index n_total = 1 + n_random;  // residual + random

    Eigen::VectorXd sigma(n_total);

    // residual variance first
    sigma(0) = state.residual().variance;

    // random effects
    Eigen::Index idx = 1;
    for (const auto& r : state.random())
    {
        sigma(idx++) = r.variance;
    }

    return sigma;
}

auto distribute_variance_components(
    gelex::FreqState& state,
    const Eigen::Ref<const Eigen::VectorXd>& sigma) -> void
{
    // residual variance first
    state.residual().variance = sigma(0);

    // random effects
    Eigen::Index idx = 1;
    for (auto& r : state.random())
    {
        r.variance = sigma(idx++);
    }
}

}  // namespace gelex
