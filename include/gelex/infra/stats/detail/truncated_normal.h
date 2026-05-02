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

#ifndef GELEX_INFRA_STATS_DETAIL_TRUNCATED_NORMAL_H_
#define GELEX_INFRA_STATS_DETAIL_TRUNCATED_NORMAL_H_

#include <algorithm>
#include <cmath>
#include <random>

#include "gelex/infra/stats/detail/normal.h"

namespace gelex::stats::detail
{

// Robert (1995) accept-reject for TN(0,1,[alpha,+inf)) when alpha >= 3.5
// Uses shifted exponential proposal with optimal lambda
inline auto sample_left_truncated_standard_robert(
    double alpha,
    std::mt19937_64& rng) -> double
{
    const double lambda = (alpha + std::sqrt((alpha * alpha) + 4.0)) / 2.0;
    std::exponential_distribution<double> exp_dist(lambda);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    for (;;)
    {
        const double z = alpha + exp_dist(rng);
        const double rho = std::exp(-0.5 * (z - lambda) * (z - lambda));
        if (uniform(rng) <= rho)
        {
            return z;
        }
    }
}

// Sample from TN(0,1,[alpha,+inf))
// Uses inverse CDF when alpha < 3.5, Robert (1995) otherwise
inline auto sample_left_truncated_standard(double alpha, std::mt19937_64& rng)
    -> double
{
    constexpr double kThreshold = 3.5;

    if (alpha >= kThreshold)
    {
        return sample_left_truncated_standard_robert(alpha, rng);
    }

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double lo = norm_cdf(alpha);
    double u = lo + (uniform(rng) * (1.0 - lo));
    u = std::clamp(u, 1e-15, 1.0 - 1e-15);
    return norm_ppf(u);
}

// Sample from TN(mu, sigma, [a, +inf))
inline auto sample_left_truncated_normal(
    double mu,
    double sigma,
    double a,
    std::mt19937_64& rng) -> double
{
    const double alpha = (a - mu) / sigma;
    return mu + (sigma * sample_left_truncated_standard(alpha, rng));
}

// Sample from TN(mu, sigma, (-inf, b])
// By symmetry: -TN(-mu, sigma, [-b, +inf))
inline auto sample_right_truncated_normal(
    double mu,
    double sigma,
    double b,
    std::mt19937_64& rng) -> double
{
    return -sample_left_truncated_normal(-mu, sigma, -b, rng);
}

}  // namespace gelex::stats::detail

#endif  // GELEX_INFRA_STATS_DETAIL_TRUNCATED_NORMAL_H_
