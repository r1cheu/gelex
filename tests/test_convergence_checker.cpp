#include <catch2/catch_test_macros.hpp>

#include "gelex/optim/convergence_checker.h"

namespace gelex
{
TEST_CASE("ConvergenceChecker Constructor", "[convergence]")
{
    SECTION("Default and custom tolerance")
    {
        ConvergenceChecker default_checker;
        ConvergenceChecker tight_checker(1e-12);
        ConvergenceChecker loose_checker(1e-2);

        Eigen::VectorXd sigma1 = Eigen::VectorXd::Random(5);
        Eigen::VectorXd sigma2 = sigma1 * 1.001;

        REQUIRE(default_checker.is_converged(sigma1, 0.0) == false);
        REQUIRE(tight_checker.is_converged(sigma1, 0.0) == false);
        REQUIRE(loose_checker.is_converged(sigma1, 1.0) == false);
        REQUIRE(loose_checker.is_converged(sigma2, 1.00001) == true);
    }
}

TEST_CASE("ConvergenceChecker First Call", "[convergence]")
{
    SECTION("Always returns false regardless of size or loglike")
    {
        ConvergenceChecker checker;

        REQUIRE(checker.is_converged(Eigen::VectorXd::Random(1), 0.0) == false);
        REQUIRE(checker.is_converged(Eigen::VectorXd::Random(2), 0.0) == false);
        REQUIRE(checker.is_converged(Eigen::VectorXd::Random(5), 0.0) == false);
        REQUIRE(
            checker.is_converged(Eigen::VectorXd::Random(10), 100.0) == false);
        REQUIRE(
            checker.is_converged(Eigen::VectorXd::Constant(5, 1.0), -100.0)
            == false);
    }
}

TEST_CASE("ConvergenceChecker Convergence Conditions", "[convergence]")
{
    SECTION("Converged with small differences")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.0 + 1e-12);

        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2, 100.00001) == true);
        REQUIRE(checker.is_converged(sigma2, 99.999) == true);
        REQUIRE(checker.is_converged(sigma2, 100.0) == true);
    }

    SECTION("Not converged with large differences")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.1);

        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2, 100.00001) == false);

        Eigen::VectorXd sigma3 = Eigen::VectorXd::Constant(5, 1.0 + 1e-12);
        REQUIRE(checker.is_converged(sigma3, 100.1) == false);
        REQUIRE(checker.is_converged(sigma3, 98.0) == false);
    }
}

TEST_CASE("ConvergenceChecker State", "[convergence]")
{
    SECTION("Stays converged after first convergence")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.0 + 1e-12);
        Eigen::VectorXd sigma3 = Eigen::VectorXd::Constant(5, 1.5);

        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2, 100.00001) == true);
        REQUIRE(checker.is_converged(sigma3, 100.0) == true);
        REQUIRE(checker.is_converged(sigma1 * 10, 1000.0) == true);
    }

    SECTION("Clear resets converged state")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.0 + 1e-12);

        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2, 100.00001) == true);

        checker.clear();

        Eigen::VectorXd sigma3 = Eigen::VectorXd::Constant(5, 2.0);
        Eigen::VectorXd sigma4 = Eigen::VectorXd::Constant(5, 2.0 + 1e-12);

        REQUIRE(checker.is_converged(sigma3, 100.0) == false);
        REQUIRE(checker.is_converged(sigma3, 200.0) == false);
        REQUIRE(checker.is_converged(sigma4, 200.00001) == true);
    }
}

TEST_CASE("ConvergenceChecker Vector Size Handling", "[convergence]")
{
    SECTION("Various vector sizes converge correctly")
    {
        ConvergenceChecker checker(1e-8);

        Eigen::VectorXd sigma1(1);
        sigma1 << 1.0;
        Eigen::VectorXd sigma1b(1);
        sigma1b << 1.0 + 1e-12;
        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma1b, 100.00001) == true);

        Eigen::VectorXd sigma2(2);
        sigma2 << 1.0, 2.0;
        Eigen::VectorXd sigma2b(2);
        sigma2b << 1.0 + 1e-12, 2.0 + 1e-12;
        checker.clear();
        REQUIRE(checker.is_converged(sigma2, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2b, 100.00001) == true);

        Eigen::VectorXd sigma3 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma3b = Eigen::VectorXd::Constant(5, 1.0 + 1e-12);
        checker.clear();
        REQUIRE(checker.is_converged(sigma3, 100.0) == false);
        REQUIRE(checker.is_converged(sigma3b, 100.00001) == true);

        Eigen::VectorXd sigma4 = Eigen::VectorXd::Constant(10, 1.0);
        Eigen::VectorXd sigma4b = Eigen::VectorXd::Constant(10, 1.0 + 1e-12);
        checker.clear();
        REQUIRE(checker.is_converged(sigma4, 100.0) == false);
        REQUIRE(checker.is_converged(sigma4b, 100.00001) == true);
    }

    SECTION("Size change between iterations prevents convergence")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(3, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.0 + 1e-12);

        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2, 100.00001) == false);
    }

    SECTION("Zero vector handling")
    {
        ConvergenceChecker checker;
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Zero(5);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Zero(5);

        REQUIRE(checker.is_converged(sigma1, 0.0) == false);
        REQUIRE(checker.is_converged(sigma2, 0.0) == false);
    }

    SECTION("Very small and very large values")
    {
        ConvergenceChecker checker(1e-8);

        Eigen::VectorXd small1 = Eigen::VectorXd::Constant(5, 1e-10);
        Eigen::VectorXd small2 = Eigen::VectorXd::Constant(5, 1e-10 + 1e-20);
        checker.is_converged(small1, 0.0);
        REQUIRE(checker.is_converged(small2, 0.0) == true);

        checker.clear();
        Eigen::VectorXd large1 = Eigen::VectorXd::Constant(5, 1e10);
        Eigen::VectorXd large2 = Eigen::VectorXd::Constant(5, 1e10 + 1e4);
        checker.is_converged(large1, 0.0);
        REQUIRE(checker.is_converged(large2, 0.0) == false);
    }
}

TEST_CASE("ConvergenceChecker Boundary Conditions", "[convergence]")
{
    SECTION("Boundary at tolerance")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.0 + 1.1e-8);

        checker.is_converged(sigma1, 100.0);
        REQUIRE(checker.is_converged(sigma2, 100.0) == false);
    }

    SECTION("Boundary loglike diff")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.0 + 1e-12);

        checker.is_converged(sigma1, 100.0);
        REQUIRE(checker.is_converged(sigma2, 100.0001) == false);
        REQUIRE(checker.is_converged(sigma2, 99.99) == false);
        REQUIRE(checker.is_converged(sigma2, 99.991) == false);
        REQUIRE(checker.is_converged(sigma2, 99.9909) == true);
    }
}

TEST_CASE("ConvergenceChecker Multiple Iterations", "[convergence]")
{
    SECTION("Gradual convergence over iterations")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma = Eigen::VectorXd::Constant(5, 1.0);
        double loglike = 100.0;

        REQUIRE(checker.is_converged(sigma, loglike) == false);

        sigma = Eigen::VectorXd::Constant(5, 1.0 + 1e-9);
        loglike = 100.00000001;
        REQUIRE(checker.is_converged(sigma, loglike) == true);
    }

    SECTION("Oscillating values never converge")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma1 = Eigen::VectorXd::Constant(5, 1.0);
        Eigen::VectorXd sigma2 = Eigen::VectorXd::Constant(5, 1.1);

        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2, 100.1) == false);
        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
        REQUIRE(checker.is_converged(sigma2, 100.1) == false);
        REQUIRE(checker.is_converged(sigma1, 100.0) == false);
    }

    SECTION("Multiple iterations with small improvements")
    {
        ConvergenceChecker checker(1e-8);
        Eigen::VectorXd sigma = Eigen::VectorXd::Constant(5, 1.0);
        double loglike = 100.0;

        for (int i = 0; i < 5; ++i)
        {
            sigma += Eigen::VectorXd::Constant(5, 1e-10);
            loglike += 1e-5;
            checker.is_converged(sigma, loglike);
        }

        Eigen::VectorXd final_sigma = sigma;
        REQUIRE(checker.is_converged(final_sigma, loglike) == true);
    }
}
}  // namespace gelex
