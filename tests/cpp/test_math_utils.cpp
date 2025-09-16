#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include "../src/data/math_utils.h"

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using namespace gelex;

TEST_CASE("centralize works correctly", "[utils]")
{
    MatrixXd x{{1, 4}, {2, 5}, {3, 6}};
    RowVectorXd expected_means{{2, 5}};
    MatrixXd expected_x{{-1, -1}, {0, 0}, {1, 1}};

    RowVectorXd actual_means = detail::centralize(x);

    REQUIRE(actual_means.isApprox(expected_means, 1e-7));
    REQUIRE(x.isApprox(expected_x, 1e-7));
}

TEST_CASE("standardize works correctly", "[utils]")
{
    MatrixXd x{{1, 4}, {2, 5}, {3, 6}};
    RowVectorXd expected_means{{2, 5}};
    RowVectorXd expected_stddevs{{1.0, 1.0}};  // After standardization
    MatrixXd expected_x{{-1, -1}, {0, 0}, {1, 1}};

    auto [actual_means, actual_stddevs] = detail::standardize(x);

    REQUIRE(actual_means.isApprox(expected_means, 1e-7));
    REQUIRE(actual_stddevs.isApprox(expected_stddevs, 1e-7));
    REQUIRE(x.isApprox(expected_x, 1e-7));
}

TEST_CASE("standardize handles constant columns", "[utils]")
{
    MatrixXd x{{1, 5}, {1, 5}, {1, 5}};
    RowVectorXd expected_means{{1, 5}};
    RowVectorXd expected_stddevs{{0, 0}};
    MatrixXd expected_x = MatrixXd::Zero(3, 2);  // All zeros

    auto [means, stddevs] = detail::standardize(x);

    REQUIRE(means.isApprox(expected_means, 1e-7));
    REQUIRE(stddevs.isApprox(expected_stddevs, 1e-7));
    REQUIRE(x.isApprox(expected_x, 1e-7));
}
