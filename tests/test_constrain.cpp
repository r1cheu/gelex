#include <cmath>

#include <catch2/catch_test_macros.hpp>

#include "gelex/optim/constrain.h"

namespace gelex
{
TEST_CASE("Constrain Tests", "[constrain]")
{
    SECTION("No negative values - no change")
    {
        Eigen::VectorXd varcmp(5);
        varcmp << 1.0, 2.0, 3.0, 4.0, 5.0;
        Eigen::VectorXd original = varcmp;

        constrain(varcmp, 100.0);

        REQUIRE(varcmp.isApprox(original));
    }

    SECTION("Zero values - not negative")
    {
        Eigen::VectorXd varcmp(3);
        varcmp << 0.0, 1.0, 2.0;
        Eigen::VectorXd original = varcmp;

        constrain(varcmp, 100.0);

        REQUIRE(varcmp.isApprox(original));
    }

    SECTION("Negative values are constrained")
    {
        Eigen::VectorXd varcmp(5);
        varcmp << 10.0, 20.0, -0.5, 40.0, 50.0;
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[2] == limit);
        REQUIRE(varcmp[0] < 10.0);
        REQUIRE(varcmp[1] < 20.0);
        REQUIRE(varcmp[3] < 40.0);
        REQUIRE(varcmp[4] < 50.0);
    }

    SECTION("Multiple negative values")
    {
        Eigen::VectorXd varcmp(5);
        varcmp << -0.5, 20.0, -0.3, 40.0, 50.0;
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[0] == limit);
        REQUIRE(varcmp[2] == limit);
        REQUIRE(varcmp[1] < 20.0);
        REQUIRE(varcmp[3] < 40.0);
        REQUIRE(varcmp[4] < 50.0);
    }

    SECTION("All negative values")
    {
        Eigen::VectorXd varcmp(3);
        varcmp << -1.0, -2.0, -3.0;
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[0] == limit);
        REQUIRE(varcmp[1] == limit);
        REQUIRE(varcmp[2] == limit);
    }

    SECTION("Half negative values - boundary")
    {
        Eigen::VectorXd varcmp(4);
        varcmp << -1.0, -2.0, 30.0, 40.0;
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[0] == limit);
        REQUIRE(varcmp[1] == limit);
        REQUIRE(varcmp[2] < 30.0);
        REQUIRE(varcmp[3] < 40.0);
    }

    SECTION("More than half negative values")
    {
        Eigen::VectorXd varcmp(4);
        varcmp << -1.0, -2.0, -3.0, 400.0;
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[0] == limit);
        REQUIRE(varcmp[1] == limit);
        REQUIRE(varcmp[2] == limit);
        REQUIRE(varcmp[3] < 400.0);
    }

    SECTION("Small positive values - not adjusted below zero")
    {
        Eigen::VectorXd varcmp(3);
        varcmp << 0.0001, 0.0001, -0.5;
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[2] == limit);
        REQUIRE(varcmp[0] >= 0.0);
        REQUIRE(varcmp[1] >= 0.0);
    }

    SECTION("Single element tests")
    {
        Eigen::VectorXd positive(1);
        positive << 5.0;
        constrain(positive, 100.0);
        REQUIRE(positive[0] == 5.0);

        Eigen::VectorXd negative(1);
        negative << -5.0;
        constrain(negative, 100.0);
        REQUIRE(negative[0] == 100.0 * 1e-6);
    }

    SECTION("Large vector with some negatives")
    {
        Eigen::VectorXd varcmp(10);
        varcmp << 10.0, 20.0, 30.0, 40.0, 50.0, -10.0, -20.0, 60.0, 70.0, 80.0;
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[5] == limit);
        REQUIRE(varcmp[6] == limit);
        for (int i = 0; i < 10; ++i)
        {
            if (i != 5 && i != 6)
            {
                REQUIRE(varcmp[i] > 0);
            }
        }
    }

    SECTION("All positive large vector")
    {
        Eigen::VectorXd varcmp = Eigen::VectorXd::Ones(100);
        Eigen::VectorXd original = varcmp;

        constrain(varcmp, 100.0);

        REQUIRE(varcmp.isApprox(original));
    }
}

TEST_CASE("Constrain Sum Preservation", "[constrain]")
{
    SECTION("Sum preserved after constraining")
    {
        Eigen::VectorXd varcmp(4);
        varcmp << -0.5, 20.0, -0.3, 40.0;
        double y_ssq = 100.0;
        double original_sum = varcmp.sum();

        constrain(varcmp, y_ssq);

        double new_sum = varcmp.sum();
        REQUIRE(std::abs(new_sum - original_sum) < 1e-10);
    }

    SECTION("Sum preserved with single negative")
    {
        Eigen::VectorXd varcmp(3);
        varcmp << 10.0, 20.0, -0.1;
        double y_ssq = 100.0;
        double original_sum = varcmp.sum();

        constrain(varcmp, y_ssq);

        double new_sum = varcmp.sum();
        REQUIRE(std::abs(new_sum - original_sum) < 1e-10);
    }

    SECTION("Sum with multiple negatives")
    {
        Eigen::VectorXd varcmp(4);
        varcmp << -0.1, -1.0, 500.0, 1000.0;
        double y_ssq = 100.0;
        double original_sum = varcmp.sum();

        constrain(varcmp, y_ssq);

        double new_sum = varcmp.sum();
        REQUIRE(std::abs(new_sum - original_sum) < 1e-9);
        REQUIRE(varcmp[0] == y_ssq * 1e-6);
        REQUIRE(varcmp[1] == y_ssq * 1e-6);
    }
}

TEST_CASE("Constrain Different y_ssq Values", "[constrain]")
{
    SECTION("Various y_ssq values produce appropriate limits")
    {
        Eigen::VectorXd varcmp1(3);
        varcmp1 << 100.0, -50.0, 200.0;
        constrain(varcmp1, 1e6);
        REQUIRE(varcmp1[1] == 1e6 * 1e-6);
        REQUIRE(varcmp1[0] < 100.0);
        REQUIRE(varcmp1[2] < 200.0);

        Eigen::VectorXd varcmp2(3);
        varcmp2 << 0.001, -0.00001, 0.002;
        constrain(varcmp2, 1e-6);
        REQUIRE(varcmp2[1] == 1e-6 * 1e-6);
        REQUIRE(varcmp2[0] < 0.001);
        REQUIRE(varcmp2[2] < 0.002);

        Eigen::VectorXd varcmp3(3);
        varcmp3 << 1.0, -0.5, 2.0;
        constrain(varcmp3, 1.0);
        REQUIRE(varcmp3[1] == 1e-6);
        REQUIRE(varcmp3[0] < 1.0);
        REQUIRE(varcmp3[2] < 2.0);
    }
}

TEST_CASE("Constrain Edge Cases", "[constrain]")
{
    SECTION("All values at limit")
    {
        Eigen::VectorXd varcmp(3);
        double y_ssq = 100.0;
        double limit = y_ssq * 1e-6;
        varcmp << limit, limit, limit;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[0] == limit);
        REQUIRE(varcmp[1] == limit);
        REQUIRE(varcmp[2] == limit);
    }

    SECTION("Values very close to zero positive")
    {
        Eigen::VectorXd varcmp(3);
        varcmp << 1e-15, 1e-15, -1e-10;
        double y_ssq = 100.0;

        constrain(varcmp, y_ssq);

        REQUIRE(varcmp[2] == y_ssq * 1e-6);
        REQUIRE(varcmp[0] >= 0.0);
        REQUIRE(varcmp[1] >= 0.0);
    }
}
}  // namespace gelex
