#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "chenx/data/grm.h"

TEST_CASE("Genetic Relationship Matrix Test", "[Grm]")
{
    using arma::dmat;
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";
    using arma::approx_equal;

    SECTION("Initialization")
    {
        REQUIRE_NOTHROW(chenx::AddGrm(test_bed, {}, 100));
    }

    SECTION("Small chunk size")
    {
        chenx::AddGrm grm_maker{test_bed, {}, 2};
        dmat grm{grm_maker.Compute()};
        REQUIRE(!grm.is_zero());
    }
    SECTION("Big chunk size")
    {
        chenx::AddGrm grm_maker{test_bed, {}, 10};
        dmat grm{grm_maker.Compute()};
        REQUIRE(!grm.is_zero());
    }

    SECTION("Result matches under different chunk size")
    {
        chenx::AddGrm small_maker{test_bed, {}, 2};
        dmat small_grm{small_maker.Compute()};

        chenx::AddGrm big_maker{test_bed, {}, 10};
        dmat big_grm{big_maker.Compute()};
        REQUIRE(approx_equal(small_grm, big_grm, "absdiff", 1e-5));
    }

    SECTION("Check AddGrm Result")
    {
        chenx::AddGrm grm_maker{test_bed, {}, 100};
        dmat grm{grm_maker.Compute()};
        dmat expected_grm{
            {0.33333337, -0.33333331, 1.1589792e-08},
            {-0.33333331, 1.5, -1.1666666},
            {1.1589792e-08, -1.1666666, 1.1666666}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", 1e-5));
    }

    SECTION("Check DomGrm Result")
    {
        chenx::DomGrm grm_maker{test_bed, {}, 2};
        dmat grm{grm_maker.Compute()};
        dmat expected_grm{
            {0.88235295, 0.35294119, -0.52941173},
            {0.35294119, 0.88235295, -2.6490952e-08},
            {-0.52941173, -2.6490952e-08, 1.2352941}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", 1e-5));
    }
}
