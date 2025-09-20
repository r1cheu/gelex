#include <catch2/catch_test_macros.hpp>
#include "gelex/estimator/bayes/params.h"

TEST_CASE("MCMCParams constructor validates parameters", "[bayes][params]")
{
    SECTION("Valid parameters")
    {
        REQUIRE_NOTHROW(gelex::MCMCParams(1000, 100, 10, 4));
    }

    SECTION("Invalid burn-in")
    {
        REQUIRE_THROWS_AS(
            gelex::MCMCParams(1000, 1000, 10, 4), std::invalid_argument);
    }

    SECTION("Edge case: burn-in just below iterations")
    {
        REQUIRE_NOTHROW(gelex::MCMCParams(1000, 999, 10, 4));
    }
}
