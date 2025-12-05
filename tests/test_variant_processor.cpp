#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gelex/data/variant_processor.h"
#include "gelex/exception.h"

using namespace gelex;
using Catch::Matchers::EndsWith;
using Catch::Matchers::MessageMatches;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("StandardizingProcessor - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0;
        Eigen::VectorXd original = variant;

        auto stats = StandardizingProcessor::process_variant(variant);

        // Check statistics
        REQUIRE_THAT(stats.mean, WithinRel(0.8, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(0.8366600265340756, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Check standardization
        Eigen::VectorXd expected(5);
        expected << -0.9561828874675147, 0.23904572186687866, 1.434274331201319,
            0.23904572186687866, -0.9561828874675147;

        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(expected(i), 1e-10));
        }
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 2.0, 2.0, 2.0, 2.0, 2.0;

        auto stats = StandardizingProcessor::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);

        // Variant should remain unchanged for monomorphic case
        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(2.0, 1e-10));
        }
    }

    SECTION("Exception path - variant size too small")
    {
        Eigen::VectorXd variant(1);
        variant << 0.0;

        REQUIRE_THROWS_MATCHES(
            StandardizingProcessor::process_variant(variant),
            InvalidInputException,
            MessageMatches(
                EndsWith("variant size 1 too small for processing")));
    }
}

TEST_CASE("RawProcessor - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0;
        Eigen::VectorXd original = variant;

        auto stats = RawProcessor::process_variant(variant);

        // Check statistics
        REQUIRE_THAT(stats.mean, WithinRel(0.8, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(0.8366600265340756, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Variant should remain unchanged
        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(original(i), 1e-10));
        }
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 2.0, 2.0, 2.0, 2.0, 2.0;

        auto stats = RawProcessor::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);
    }

    SECTION("Exception path - variant size too small")
    {
        Eigen::VectorXd variant(1);
        variant << 0.0;

        REQUIRE_THROWS_MATCHES(
            RawProcessor::process_variant(variant),
            InvalidInputException,
            MessageMatches(
                EndsWith("variant size 1 too small for processing")));
    }
}

TEST_CASE("HardWenbergProcessor - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0;
        Eigen::VectorXd original = variant;

        auto stats = HardWenbergProcessor::process_variant(variant);

        // Check statistics
        REQUIRE_THAT(stats.mean, WithinRel(0.8, 1e-10));
        double expected_stddev
            = std::sqrt(stats.mean * (1.0 - 0.5 * stats.mean));
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Check standardization with Hardy-Weinberg stddev
        // stddev = sqrt(mean * (1 - 0.5*mean)) = sqrt(0.8 * 0.6) = sqrt(0.48) =
        // 0.6928203230275509
        Eigen::VectorXd expected(5);
        expected << -1.1547005383792517, 0.28867513459481293,
            1.7320508075688772, 0.28867513459481293, -1.1547005383792517;

        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(expected(i), 1e-10));
        }
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 2.0, 2.0, 2.0, 2.0, 2.0;

        auto stats = HardWenbergProcessor::process_variant(variant);

        REQUIRE_THAT(stats.mean, WithinRel(2.0, 1e-10));
        double expected_stddev
            = std::sqrt(stats.mean * (1.0 - 0.5 * stats.mean));
        REQUIRE_THAT(stats.stddev, WithinRel(expected_stddev, 1e-10));
        REQUIRE(stats.is_monomorphic);

        // Variant should remain unchanged for monomorphic case
        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(2.0, 1e-10));
        }
    }

    SECTION("Exception path - variant size too small")
    {
        Eigen::VectorXd variant(1);
        variant << 0.0;

        REQUIRE_THROWS_MATCHES(
            HardWenbergProcessor::process_variant(variant),
            InvalidInputException,
            MessageMatches(
                EndsWith("variant size 1 too small for processing")));
    }
}

TEST_CASE("DominantStandardizingProcessor - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant with heterozygotes")
    {
        Eigen::VectorXd variant(6);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0, 2.0;

        auto stats = DominantStandardizingProcessor::process_variant(variant);

        // After converting 2.0 to 0.0, we have: 0.0, 1.0, 0.0, 1.0, 0.0, 0.0
        // Mean = (0+1+0+1+0+0)/6 = 2/6 = 0.333333
        // Variance = [(0-0.333)^2*4 + (1-0.333)^2*2]/5 = [0.111111*4 +
        // 0.444444*2]/5 = [0.444444 + 0.888888]/5 = 1.333332/5 = 0.2666664
        // Stddev = sqrt(0.2666664) = 0.5163978

        REQUIRE_THAT(stats.mean, WithinRel(0.3333333333333333, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(0.5163977794943222, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Check that 2.0 values were converted to 0.0 and then standardized
        // Original value 2.0 -> 0.0, then (0.0 - 0.333333)/0.5163978 =
        // -0.645497
        REQUIRE_THAT(variant(2), WithinRel(-0.6454972243679028, 1e-10));
        REQUIRE_THAT(variant(5), WithinRel(-0.6454972243679028, 1e-10));
    }

    SECTION("Happy path - variant with no heterozygotes")
    {
        Eigen::VectorXd variant(4);
        variant << 0.0, 2.0, 0.0, 2.0;

        auto stats = DominantStandardizingProcessor::process_variant(variant);

        // After converting 2.0 to 0.0, we have all zeros
        REQUIRE_THAT(stats.mean, WithinRel(0.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);
    }

    SECTION("Exception path - variant size too small")
    {
        Eigen::VectorXd variant(1);
        variant << 0.0;

        REQUIRE_THROWS_MATCHES(
            DominantStandardizingProcessor::process_variant(variant),
            InvalidInputException,
            MessageMatches(
                EndsWith("variant size 1 too small for processing")));
    }
}

TEST_CASE("DominantRawProcessor - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(6);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0, 2.0;
        Eigen::VectorXd original = variant;

        auto stats = DominantRawProcessor::process_variant(variant);

        // Statistics should be computed after converting 2.0 to 0.0
        REQUIRE_THAT(stats.mean, WithinRel(0.3333333333333333, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(0.5163977794943222, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Check that 2.0 values were converted to 0.0 (no standardization in
        // RawProcessor)
        REQUIRE_THAT(variant(2), WithinRel(0.0, 1e-10));
        REQUIRE_THAT(variant(5), WithinRel(0.0, 1e-10));

        // Other values should remain unchanged
        REQUIRE_THAT(variant(0), WithinRel(0.0, 1e-10));
        REQUIRE_THAT(variant(1), WithinRel(1.0, 1e-10));
        REQUIRE_THAT(variant(3), WithinRel(1.0, 1e-10));
        REQUIRE_THAT(variant(4), WithinRel(0.0, 1e-10));
    }

    SECTION("Happy path - monomorphic variant after conversion")
    {
        Eigen::VectorXd variant(4);
        variant << 2.0, 2.0, 2.0, 2.0;

        auto stats = DominantRawProcessor::process_variant(variant);

        // After converting all 2.0 to 0.0, we have all zeros
        REQUIRE_THAT(stats.mean, WithinRel(0.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);

        // All values should be 0.0
        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(0.0, 1e-10));
        }
    }
}

TEST_CASE("DominantOrthogonalHWEProcessor - Basic functionality", "[data]")
{
    SECTION("Happy path - polymorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 0.0, 1.0, 2.0, 1.0, 0.0;

        auto stats = DominantOrthogonalHWEProcessor::process_variant(variant);

        // p_freq = mean/2 = 0.8/2 = 0.4
        // stats.mean = 2 * p_freq^2 = 2 * 0.16 = 0.32
        // stats.stddev = 2 * p_freq * (1 - p_freq) = 2 * 0.4 * 0.6 = 0.48

        REQUIRE_THAT(stats.mean, WithinRel(0.32, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinRel(0.48, 1e-10));
        REQUIRE_FALSE(stats.is_monomorphic);

        // Check transformation
        // one_alt_encode = 2 * p_freq = 0.8
        // two_alt_encode = 4 * p_freq - 2 = 1.6 - 2 = -0.4
        // Then standardized: (x - 0.32) / 0.48
        // For 0.0: (0.0 - 0.32)/0.48 = -0.666667
        // For 1.0: (0.8 - 0.32)/0.48 = 1.0
        // For 2.0: (-0.4 - 0.32)/0.48 = -1.5

        Eigen::VectorXd expected(5);
        expected << -0.6666666666666666, 1.0, -1.5, 1.0, -0.6666666666666666;

        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(expected(i), 1e-10));
        }
    }

    SECTION("Happy path - monomorphic variant")
    {
        Eigen::VectorXd variant(5);
        variant << 2.0, 2.0, 2.0, 2.0, 2.0;

        auto stats = DominantOrthogonalHWEProcessor::process_variant(variant);

        // p_freq = 2.0/2 = 1.0
        // stats.mean = 2 * 1.0^2 = 2.0
        // stats.stddev = 2 * 1.0 * 0.0 = 0.0

        REQUIRE_THAT(stats.mean, WithinRel(2.0, 1e-10));
        REQUIRE_THAT(stats.stddev, WithinAbs(0.0, 1e-10));
        REQUIRE(stats.is_monomorphic);

        // Variant should remain unchanged for monomorphic case
        for (int i = 0; i < variant.size(); ++i)
        {
            REQUIRE_THAT(variant(i), WithinRel(2.0, 1e-10));
        }
    }
}
