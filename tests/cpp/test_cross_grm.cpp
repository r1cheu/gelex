#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include "gelex/data/bed_reader.h"
#include "gelex/data/grm.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

TEST_CASE("CrossGrm computation", "[cross_grm]")
{
    // Mock file paths
    const std::string train_bed
        = std::string(GELEX_TESTS_DIR) + "/data/train.bed";
    const std::string test_bed
        = std::string(GELEX_TESTS_DIR) + "/data/test.bed";
    const std::string test_missmatch_bed
        = std::string(GELEX_TESTS_DIR) + "/data/test_missmatch.bed";
    VectorXd add_center{{1.0, 0.3333333, 1.3333333, 0.6666667}};
    VectorXd p_major = add_center / 2;
    double add_scale_factor = 2.0;
    double dom_scale_factor = 0.9444445;

    SECTION("AddCrossGrm computation")
    {
        gelex::CrossGRM cgrm(train_bed, std::move(p_major), add_scale_factor);
        MatrixXd grm = cgrm.compute(test_bed);

        MatrixXd expected_grm{
            {3.3333334e-01, -3.3333334e-01, -3.3113690e-09},
            {1.6556845e-09, -1.6666663e-01, 1.6666661e-01},
            {-1.6666666e-01, -1.3333334e+00, 1.5000000e+00}};
        REQUIRE(grm.isApprox(expected_grm, 1e-5));
    }

    SECTION("AddCrossGrm chunk computation")
    {
        gelex::CrossGRM cgrm(train_bed, p_major, add_scale_factor, 2);
        MatrixXd grm = cgrm.compute(test_bed);

        MatrixXd expected_grm{
            {3.3333334e-01, -3.3333334e-01, -3.3113690e-09},
            {1.6556845e-09, -1.6666663e-01, 1.6666661e-01},
            {-1.6666666e-01, -1.3333334e+00, 1.5000000e+00}};
        REQUIRE(grm.isApprox(expected_grm, 1e-5));
    }

    SECTION("DomCrossGrm computation")
    {
        gelex::CrossGRM cgrm(train_bed, p_major, dom_scale_factor);
        MatrixXd grm = cgrm.compute(test_bed, false);
        MatrixXd expected_grm{
            {0.882353, 0.35294122, -0.5294118},
            {-1.0000001, -0.4705883, 0.76470584},
            {-0.2352941, 0.29411766, 0.47058827}};

        REQUIRE(grm.isApprox(expected_grm, 1e-5));
    }

    SECTION("DomCrossGrm Chunk computation")
    {
        gelex::CrossGRM dgrm(train_bed, p_major, dom_scale_factor, 2);
        MatrixXd grm = dgrm.compute(test_bed, false);

        MatrixXd expected_grm{
            {0.882353, 0.35294122, -0.5294118},
            {-1.0000001, -0.4705883, 0.76470584},
            {-0.2352941, 0.29411766, 0.47058827}};

        REQUIRE(grm.isApprox(expected_grm, 1e-5));
    }

    SECTION("CrossGrm with SNP mismatch throws exception")
    {
        gelex::CrossGRM add_grm(train_bed, p_major, add_scale_factor);
        REQUIRE_THROWS_WITH(
            add_grm.compute(test_missmatch_bed),
            "SNPs in training and test sets do not match.");
    }

    SECTION("CrossGrm with target individuals")
    {
        std::vector<std::string> target_individuals = {"iid1", "iid3"};
        VectorXd center{{0.5, 0.5, 1.5, 0.0}};
        p_major = center / 2;
        double scale_factor = 0.75;

        gelex::CrossGRM add_grm(
            train_bed,
            p_major,
            scale_factor,
            gelex::DEFAULT_CHUNK_SIZE,
            target_individuals);

        MatrixXd grm = add_grm.compute(test_bed);
        MatrixXd expected_grm{
            {1.0, -1.0}, {0.33333334, -0.33333334}, {-1.6666666, 1.6666666}};

        REQUIRE(grm.isApprox(expected_grm, 1e-7));
    }
}
