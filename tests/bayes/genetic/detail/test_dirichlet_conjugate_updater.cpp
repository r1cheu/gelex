// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <random>

#include "gelex/bayes/genetic/detail/dirichlet_conjugate_updater.h"
#include "gelex/bayes/stats/dirichlet_distribution.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"

using gelex::DirichletDistribution;
using gelex::detail::DirichletConjugateUpdater;

TEST_CASE(
    "DirichletConjugateUpdater samples from the conjugate posterior",
    "[bayes][genetic][kernel][dirichlet_conjugate_updater]")
{
    const auto prior = gelex::make_dirichlet_prior(std::array{0.5, 1.5, 2.0});
    DirichletConjugateUpdater<3> updater{prior};
    constexpr std::array<std::size_t, 3> counts{2, 3, 5};
    std::mt19937_64 rng{123};

    const auto weights = updater.draw(counts, rng);

    std::mt19937_64 expected_rng{123};
    DirichletDistribution<3> distribution;
    const auto posterior = prior + gelex::make_categorical_likelihood(counts);
    const auto expected
        = distribution(expected_rng, posterior.dirichlet_parameters());
    const Eigen::Map<const Eigen::Vector3d> actual_vector{weights.data()};
    const Eigen::Map<const Eigen::Vector3d> expected_vector{expected.data()};
    REQUIRE(actual_vector.isApprox(expected_vector));
    REQUIRE(rng == expected_rng);
}
