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

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <random>
#include <type_traits>

#include "gelex/bayes/genetic/kernel/dirichlet_conjugate_updater.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/stats/dirichlet_distribution.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"

using gelex::DirichletDistribution;
using gelex::DirichletLogKernel;
using gelex::FixedParameter;
using gelex::MixtureWeightUpdate;
using gelex::Parameter;
using gelex::detail::DirichletConjugateUpdater;
using gelex::detail::make_dirichlet_conjugate_updater;

static_assert(std::is_empty_v<
              DirichletConjugateUpdater<3, MixtureWeightUpdate::Disabled>>);

TEST_CASE(
    "DirichletConjugateUpdater samples from the conjugate posterior",
    "[bayes][genetic][kernel][dirichlet_conjugate_updater]")
{
    const auto prior = gelex::make_dirichlet_prior(std::array{0.5, 1.5, 2.0});
    const Parameter<std::array<double, 3>, DirichletLogKernel<3>> parameter{
        .initial = {0.2, 0.3, 0.5}, .prior = prior};
    auto updater = make_dirichlet_conjugate_updater<3>(parameter);
    std::array weights = parameter.initial;
    constexpr std::array<std::size_t, 3> counts{2, 3, 5};
    std::mt19937_64 rng{123};

    updater.update(weights, counts, rng);

    std::mt19937_64 expected_rng{123};
    DirichletDistribution<3> distribution;
    const auto posterior = prior + gelex::make_categorical_likelihood(counts);
    const auto expected
        = distribution(expected_rng, posterior.dirichlet_parameters());
    const Eigen::Map<const Eigen::Vector3d> actual_vector{weights.data()};
    const Eigen::Map<const Eigen::Vector3d> expected_vector{expected.data()};
    REQUIRE(actual_vector.isApprox(expected_vector));
}

TEST_CASE(
    "Disabled DirichletConjugateUpdater preserves weights and RNG state",
    "[bayes][genetic][kernel][dirichlet_conjugate_updater]")
{
    const FixedParameter<std::array<double, 3>> parameter{
        .initial = {0.2, 0.3, 0.5}};
    auto updater = make_dirichlet_conjugate_updater<3>(parameter);
    auto weights = parameter.initial;
    constexpr std::array<std::size_t, 3> counts{2, 3, 5};
    std::mt19937_64 rng{123};
    std::mt19937_64 expected_rng{123};

    updater.update(weights, counts, rng);

    REQUIRE(weights == parameter.initial);
    REQUIRE(rng() == expected_rng());
}
