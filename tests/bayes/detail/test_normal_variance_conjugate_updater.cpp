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

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <random>

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/stats/scaled_inv_chi2_distribution.h"
#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"

TEST_CASE(
    "NormalVarianceConjugateUpdater samples from the conjugate posterior",
    "[bayes][detail][normal_variance_conjugate_updater]")
{
    const auto prior = gelex::make_scaled_inv_chi2_prior(4.0, 0.5);
    gelex::detail::NormalVarianceConjugateUpdater updater{prior};
    double variance = 1.0;
    constexpr std::size_t count = 3;
    constexpr double sum_squares = 1.25;
    std::mt19937_64 rng{123};

    updater.update(variance, count, sum_squares, rng);

    std::mt19937_64 expected_rng{123};
    gelex::ScaledInvChi2Distribution<> distribution;
    const auto likelihood
        = gelex::make_normal_variance_likelihood(count, sum_squares);
    const auto posterior = prior + likelihood;
    const double expected
        = distribution(expected_rng, posterior.scaled_inv_chi2_parameters());
    REQUIRE(variance == expected);
    REQUIRE(rng() == expected_rng());
}
