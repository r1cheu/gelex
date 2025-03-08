#include <cstdint>

#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "chenx/data/grm.h"

TEST_CASE("Genetic Relationship Matrix Test", "[Grm]")
{
    static constexpr double ABS_DIFF = 1e-5;
    static constexpr uint64_t BIG_CHUNK_SIZE = 2;
    static constexpr uint64_t SMALL_CHUNK_SIZE = 10;

    using arma::dmat;
    using arma::rowvec;
    const std::string train_bed
        = std::string(CHENX_TESTS_DIR) + "/data/train.bed";
    using arma::approx_equal;

    SECTION("Initialization")
    {
        REQUIRE_NOTHROW(chenx::AddGrm(train_bed, BIG_CHUNK_SIZE));
    }

    SECTION("Small chunk size")
    {
        chenx::AddGrm grm_maker{train_bed, SMALL_CHUNK_SIZE};
        dmat grm{grm_maker.Compute()};
        REQUIRE(!grm.is_zero());
    }
    SECTION("Big chunk size")
    {
        chenx::AddGrm grm_maker{train_bed, BIG_CHUNK_SIZE};
        dmat grm{grm_maker.Compute()};
        REQUIRE(!grm.is_zero());
    }

    SECTION("Result matches under different chunk size")
    {
        chenx::AddGrm small_maker{train_bed, SMALL_CHUNK_SIZE};
        dmat small_grm{small_maker.Compute()};

        chenx::AddGrm big_maker{train_bed, BIG_CHUNK_SIZE};
        dmat big_grm{big_maker.Compute()};
        REQUIRE(approx_equal(small_grm, big_grm, "absdiff", ABS_DIFF));
    }

    SECTION("Check AddGrm Result")
    {
        chenx::AddGrm grm_maker{train_bed, BIG_CHUNK_SIZE};
        dmat grm{grm_maker.Compute()};
        dmat expected_grm{
            {0.33333337, -0.33333331, 1.1589792e-08},
            {-0.33333331, 1.5, -1.1666666},
            {1.1589792e-08, -1.1666666, 1.1666666}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", ABS_DIFF));
    }

    SECTION("Check AddGrm Center")
    {
        chenx::AddGrm grm_maker{train_bed, BIG_CHUNK_SIZE};
        dmat grm{grm_maker.Compute()};
        rowvec expected{1.0, 0.3333333, 1.3333333, 0.6666667};
        REQUIRE(
            approx_equal(grm_maker.center(), expected, "absdiff", ABS_DIFF));
    }

    SECTION("Check AddGrm Scale factor")
    {
        chenx::AddGrm grm_maker{train_bed, BIG_CHUNK_SIZE};
        dmat grm{grm_maker.Compute()};
        REQUIRE_THAT(
            grm_maker.scale_factor(), Catch::Matchers::WithinAbs(2, ABS_DIFF));
    }

    SECTION("Check DomGrm Result")
    {
        chenx::DomGrm grm_maker{train_bed, 2};
        dmat grm{grm_maker.Compute()};
        dmat expected_grm{
            {0.88235295, 0.35294119, -0.52941173},
            {0.35294119, 0.88235295, -2.6490952e-08},
            {-0.52941173, -2.6490952e-08, 1.2352941}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", ABS_DIFF));
    }

    SECTION("Check DomGrm Center")
    {
        chenx::DomGrm grm_maker{train_bed, BIG_CHUNK_SIZE};
        dmat grm{grm_maker.Compute()};
        rowvec expected{0.5, 0.2777778, 0.44444442, 0.44444442};
        REQUIRE(
            approx_equal(grm_maker.center(), expected, "absdiff", ABS_DIFF));
    }

    SECTION("Check DomGrm Scale factor")
    {
        chenx::DomGrm grm_maker{train_bed, 2};
        dmat grm{grm_maker.Compute()};
        REQUIRE_THAT(
            grm_maker.scale_factor(),
            Catch::Matchers::WithinAbs(0.9444445, ABS_DIFF));
    }
}
