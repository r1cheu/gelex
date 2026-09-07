// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>

#include "gelex/bayes/stats/dirichlet_log_kernel.h"

using gelex::DirichletLogKernel;
using gelex::make_beta_prior;
using gelex::make_categorical_likelihood;
using gelex::make_dirichlet_prior;
using gelex::make_uniform_dirichlet_prior;

template <std::size_t K>
concept ValidDirichletLogKernel = requires { typename DirichletLogKernel<K>; };

static_assert(ValidDirichletLogKernel<2>);
static_assert(!ValidDirichletLogKernel<1>);

TEST_CASE(
    "DirichletLogKernel combines a prior with categorical counts",
    "[bayes][stats][dirichlet_log_kernel]")
{
    const auto posterior
        = make_dirichlet_prior(std::array{0.5, 1.5, 2.0})
          + make_categorical_likelihood(std::array<std::size_t, 3>{2, 3, 5});
    const auto parameter_values
        = posterior.dirichlet_parameters().concentrations();
    const Eigen::Map<const Eigen::Vector3d> parameters{parameter_values.data()};
    const Eigen::Vector3d expected{{2.5, 4.5, 7.0}};

    REQUIRE(parameters.isApprox(expected));
}

TEST_CASE(
    "Uniform Dirichlet prior has unit concentrations",
    "[bayes][stats][dirichlet_log_kernel]")
{
    const auto concentrations = make_uniform_dirichlet_prior<3>()
                                    .dirichlet_parameters()
                                    .concentrations();
    REQUIRE(concentrations == std::array<double, 3>{1.0, 1.0, 1.0});
}

TEST_CASE(
    "Beta prior maps failure and success to categorical order",
    "[bayes][stats][dirichlet_log_kernel]")
{
    const auto posterior
        = make_beta_prior(2.0, 3.0)
          + make_categorical_likelihood(std::array<std::size_t, 2>{5, 7});
    const auto parameter_values
        = posterior.dirichlet_parameters().concentrations();
    const Eigen::Map<const Eigen::Vector2d> parameters{parameter_values.data()};
    const Eigen::Vector2d expected{{8.0, 9.0}};

    REQUIRE(parameters.isApprox(expected));
}
