// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <random>

#include "gelex/bayes/genetic/detail/probit_updater.h"
#include "gelex/bayes/stats/multi_quadratic_log_kernel.h"

TEST_CASE(
    "ProbitUpdater samples from supplied sufficient statistics",
    "[bayes][genetic][kernel][probit_updater]")
{
    Eigen::Vector2d annotation_coefficients{{0.2, -0.4}};
    const auto prior = gelex::make_multi_normal_prior(
        Eigen::Matrix2d{{2.0, 0.5}, {0.5, 1.0}});
    const gelex::MultiQuadraticLogKernel likelihood{
        Eigen::Matrix2d{{2.0, 0.25}, {0.25, 0.3125}},
        Eigen::Vector2d{{0.7, 0.275}},
        0.0};
    std::mt19937_64 rng{123};

    gelex::detail::ProbitUpdater updater{prior};
    updater.update(annotation_coefficients, likelihood, rng);

    std::mt19937_64 expected_rng{123};
    const auto posterior = prior + likelihood;
    const Eigen::LLT<Eigen::Matrix2d> precision_factor{posterior.quadratic()};
    const Eigen::Vector2d posterior_mean
        = precision_factor.solve(posterior.linear());
    std::normal_distribution<double> standard_normal_distribution;
    const Eigen::Vector2d standard_normal{
        {standard_normal_distribution(expected_rng),
         standard_normal_distribution(expected_rng)}};
    const Eigen::Vector2d expected
        = posterior_mean + precision_factor.matrixU().solve(standard_normal);

    REQUIRE(annotation_coefficients.isApprox(expected));
    REQUIRE(rng() == expected_rng());
}
