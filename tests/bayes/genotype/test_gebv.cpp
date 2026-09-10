// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <catch2/catch_test_macros.hpp>

#include "gelex/bayes/genotype/gebv.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/csc_reader.h"

#include "compact_genotype_fixture.h"

namespace
{

auto encoded_design(const gelex::bayes::GeneticProjection& projection)
    -> Eigen::MatrixXd
{
    Eigen::MatrixXd design
        = Eigen::MatrixXd::Zero(projection.rows(), projection.cols());
    for (Eigen::Index marker = 0; marker < projection.cols(); ++marker)
    {
        projection.multiply(marker, 1.0, design.col(marker));
    }
    return design;
}

}  // namespace

TEST_CASE("GEBV of one draw applies the projection", "[bayes][gebv]")
{
    const auto design = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0, 2.0}, {1.0, 0.0, 2.0}, {2.0, 1.0, 0.0}},
        gelex::GeneticModeSet{gelex::GeneticMode::A});
    const auto& projection = design.projection(gelex::GeneticMode::A);
    const Eigen::MatrixXd coefficients{
        {0.5, 0.0, -1.0, 0.0}, {0.0, 0.0, 2.0, 0.25}, {-0.75, 0.0, 0.0, 1.5}};
    const Eigen::MatrixXd expected = encoded_design(projection) * coefficients;
    Eigen::VectorXd target = Eigen::VectorXd::Constant(3, 9.0);

    SECTION("dense coefficients")
    {
        for (Eigen::Index draw = 0; draw < coefficients.cols(); ++draw)
        {
            gelex::gebv_draw(projection, coefficients, draw, target);
            REQUIRE(target.isApprox(expected.col(draw)));
        }
    }

    SECTION("sparse coefficients")
    {
        const gelex::CscReader::sparse_matrix_type<double> sparse
            = coefficients.sparseView();
        const gelex::CscReader::sparse_map_type<double> map{
            sparse.rows(),
            sparse.cols(),
            sparse.nonZeros(),
            sparse.outerIndexPtr(),
            sparse.innerIndexPtr(),
            sparse.valuePtr()};
        for (Eigen::Index draw = 0; draw < coefficients.cols(); ++draw)
        {
            gelex::gebv_draw(projection, map, draw, target);
            REQUIRE(target.isApprox(expected.col(draw)));
        }
    }

    SECTION("an all-zero draw clears the target")
    {
        gelex::gebv_draw(projection, coefficients, 1, target);
        REQUIRE(target.isZero());
    }
}
