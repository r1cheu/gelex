#include <string>
#include <vector>

#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include "chenx/data/bed_reader.h"
#include "chenx/data/cross_grm.h"

using arma::dmat;
using arma::rowvec;
using chenx::AddCrossGrm;
using chenx::DomCrossGrm;

TEST_CASE("CrossGrm computation", "[cross_grm]")
{
    // Mock file paths
    const std::string train_bed
        = std::string(CHENX_TESTS_DIR) + "/data/train.bed";
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";
    const std::string test_missmatch_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test_missmatch.bed";
    rowvec add_center{1.0, 0.3333333, 1.3333333, 0.6666667};
    double add_scale_factor = 2.0;

    rowvec dom_center{0.5, 0.2777778, 0.44444442, 0.44444442};
    double dom_scale_factor = 0.9444445;

    using arma::approx_equal;

    SECTION("AddCrossGrm computation")
    {
        // Setup
        //
        AddCrossGrm add_grm(train_bed, std::move(add_center), add_scale_factor);
        dmat grm = add_grm.Compute(test_bed);

        dmat expected_grm
            = {{3.3333334e-01, -3.3333334e-01, -3.3113690e-09},
               {1.6556845e-09, -1.6666663e-01, 1.6666661e-01},
               {-1.6666666e-01, -1.3333334e+00, 1.5000000e+00}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", 1e-7));
    }

    SECTION("AddCrossGrm chunk computation")
    {
        rowvec center{1.0, 0.3333333, 1.3333333, 0.6666667};
        double scale_factor = 2.0;

        AddCrossGrm add_grm(
            train_bed, std::move(add_center), add_scale_factor, 2);
        dmat grm = add_grm.Compute(test_bed);

        dmat expected_grm
            = {{3.3333334e-01, -3.3333334e-01, -3.3113690e-09},
               {1.6556845e-09, -1.6666663e-01, 1.6666661e-01},
               {-1.6666666e-01, -1.3333334e+00, 1.5000000e+00}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", 1e-7));
    }

    SECTION("DomCrossGrm computation")
    {
        DomCrossGrm dom_grm(train_bed, std::move(dom_center), dom_scale_factor);
        dmat grm = dom_grm.Compute(test_bed);
        dmat expected_grm
            = {{0.882353, 0.35294122, -0.5294118},
               {-1.0000001, -0.4705883, 0.76470584},
               {-0.2352941, 0.29411766, 0.47058827}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", 1e-5));
    }

    SECTION("DomCrossGrm Chunk computation")
    {
        DomCrossGrm dom_grm(
            train_bed, std::move(dom_center), dom_scale_factor, 2);
        dmat grm = dom_grm.Compute(test_bed);

        dmat expected_grm
            = {{0.882353, 0.35294122, -0.5294118},
               {-1.0000001, -0.4705883, 0.76470584},
               {-0.2352941, 0.29411766, 0.47058827}};

        REQUIRE(approx_equal(grm, expected_grm, "absdiff", 1e-5));
    }

    SECTION("CrossGrm with SNP mismatch throws exception")
    {
        AddCrossGrm add_grm(train_bed, std::move(add_center), add_scale_factor);
        REQUIRE_THROWS_WITH(
            add_grm.Compute(test_missmatch_bed),
            "SNPs in training and test sets do not match.");
    }

    SECTION("CrossGrm with excluded individuals")
    {
        std::vector<std::string> exclude_individuals = {"iid2"};
        rowvec center{0.5, 0.5, 1.5, 0.0};
        double scale_factor = 0.75;

        AddCrossGrm add_grm(
            train_bed,
            std::move(center),
            scale_factor,
            chenx::DEFAULT_CHUNK_SIZE,
            exclude_individuals);

        dmat grm = add_grm.Compute(test_bed);
        dmat expected_grm
            = {{1.0, -1.0}, {0.33333334, -0.33333334}, {-1.6666666, 1.6666666}};
        REQUIRE(approx_equal(grm, expected_grm, "absdiff", 1e-5));
    }
}
