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

#ifndef GELEX_MODEL_BAYES_SAMPLERS_DETAIL_STICK_BREAKING_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DETAIL_STICK_BREAKING_H_

#include <random>

#include <Eigen/Core>

namespace gelex::detail
{

// Sample stick-breaking probabilities from Beta posteriors.
// counts (K+1): counts[0] = null class, counts[1..K] = non-null classes
// q[k] ~ Beta(1 + sum(counts[k+1..K]), 1 + counts[k])
inline auto sample_stick_posteriors(
    const Eigen::Ref<const Eigen::VectorXi>& counts,
    std::mt19937_64& rng) -> Eigen::VectorXd
{
    const Eigen::Index K = counts.size() - 1;
    Eigen::VectorXd q(K);

    int remaining = counts.sum() - counts(0);
    for (Eigen::Index k = 0; k < K; ++k)
    {
        const double alpha = 1.0 + remaining;
        const double beta = 1.0 + counts(k);
        remaining -= counts(k + 1);

        std::gamma_distribution<double> gamma_a(alpha, 1.0);
        std::gamma_distribution<double> gamma_b(beta, 1.0);
        const double xa = gamma_a(rng);
        const double xb = gamma_b(rng);
        q(k) = xa / (xa + xb);
    }

    return q;
}

}  // namespace gelex::detail

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DETAIL_STICK_BREAKING_H_
