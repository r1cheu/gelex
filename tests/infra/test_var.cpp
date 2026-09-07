// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/infra/var.h"

namespace gelex
{

TEST_CASE("vecvar computes vector variance", "[stats][var]")
{
    const Eigen::VectorXd values{{1.0, 2.0, 4.0, 7.0}};

    REQUIRE(vecvar(values, VarNormType::Sample) == 7.0);
    REQUIRE(vecvar(values, VarNormType::Population) == 5.25);
}

TEST_CASE("matvar computes axis-wise matrix variance", "[stats][var]")
{
    const Eigen::MatrixXd values{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 8.0},
    };
    const Eigen::RowVectorXd expected_colwise{{4.0, 28.0 / 3.0}};
    const Eigen::VectorXd expected_rowwise{{0.5, 0.5, 4.5}};

    REQUIRE(matvar<0>(values, VarNormType::Sample).isApprox(expected_colwise));
    REQUIRE(matvar<1>(values, VarNormType::Sample).isApprox(expected_rowwise));
}

}  // namespace gelex
