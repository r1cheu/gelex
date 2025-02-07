
#include <armadillo>
#include <catch2/catch_test_macros.hpp>

#include "chenx/data/grm.h"

TEST_CASE("Genetic Relationship Matrix Test", "[Grm]")
{
    using arma::dmat;
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";

    SECTION("Initialization")
    {
        REQUIRE_NOTHROW(chenx::AddGrm(test_bed, 100));
    }

    SECTION("Small chunk size")
    {
        chenx::AddGrm grm_maker{test_bed, 2};
        dmat grm{grm_maker.Compute()};
        REQUIRE(!grm.is_zero());
    }
    SECTION("Big chunk size")
    {
        chenx::AddGrm grm_maker{test_bed, 10};
        dmat grm{grm_maker.Compute()};
        REQUIRE(!grm.is_zero());
    }

    SECTION("Result matches under different chunk size")
    {
        chenx::AddGrm small_maker{test_bed, 2};
        dmat small_grm{small_maker.Compute()};

        chenx::AddGrm big_maker{test_bed, 10};
        dmat big_grm{big_maker.Compute()};
        REQUIRE(arma::approx_equal(small_grm, big_grm, "absdiff", 1e-5));
    }
}
