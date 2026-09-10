// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/bayes/variance/heritability.h"
#include "gelex/exception.h"

TEST_CASE(
    "Heritability divides by total genetic plus residual variance",
    "[bayes][variance][heritability]")
{
    const Eigen::RowVectorXd explained{{1.0, 2.0}};
    const Eigen::RowVectorXd total{{2.0, 2.0}};
    const Eigen::RowVectorXd residual{{2.0, 6.0}};
    REQUIRE(
        gelex::heritability_draws(explained, total, residual)
            .isApprox(Eigen::RowVectorXd{{0.25, 0.25}}));
    REQUIRE_THROWS_AS(
        gelex::heritability_draws(explained, total, Eigen::RowVectorXd{{1.0}}),
        gelex::GelexException);
}
