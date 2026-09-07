// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
