// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/marker_effect_traits.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/dense_reader.h"
#include "gelex/io/dense_writer.h"

#include "file_fixture.h"

TEST_CASE(
    "Coefficient summaries agree between dense and sparse layouts",
    "[bayes][genetic][marker_effects]")
{
    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "coef.draws").string();
    const auto sparse_path = path + ".csc";
    const Eigen::MatrixXd coefficients{
        {0.5, 0.0, -1.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {-0.75, 0.5, 0.0, 1.5}};
    {
        gelex::DenseWriter dense{path};
        gelex::CscWriter sparse{sparse_path};
        auto dense_stream
            = dense.reserve<double>("dense", gelex::BinaryShape{3, 4});
        auto sparse_stream
            = sparse.reserve<double>("sparse", gelex::BinaryShape{3, 4});
        for (Eigen::Index draw = 0; draw < coefficients.cols(); ++draw)
        {
            const Eigen::VectorXd column = coefficients.col(draw);
            dense_stream << column;
            sparse_stream << std::span<const double>{
                column.data(), static_cast<std::size_t>(column.size())};
        }
        sparse.close();
        dense.close();
    }
    const gelex::DenseReader dense{path};
    const gelex::CscReader sparse{sparse_path};
    const gelex::DrawReaders readers{.dense = dense, .sparse = sparse};

    const auto from_dense
        = gelex::summarize_coefficients<gelex::CoefficientLayout::Dense>(
            readers, "dense");
    const auto from_sparse
        = gelex::summarize_coefficients<gelex::CoefficientLayout::Sparse>(
            readers, "sparse");

    const Eigen::VectorXd mean = coefficients.rowwise().mean();
    const Eigen::VectorXd sd
        = gelex::matvar<1>(coefficients, gelex::VarNormType::Sample)
              .cwiseSqrt();
    REQUIRE(from_dense.mean.isApprox(mean));
    REQUIRE(from_dense.sd.isApprox(sd));
    REQUIRE(from_sparse.mean.isApprox(mean));
    REQUIRE(from_sparse.sd.isApprox(sd));
    REQUIRE(from_sparse.mean(1) == 0.0);
    REQUIRE(from_sparse.sd(1) == 0.0);
}

TEST_CASE(
    "Inclusion probability counts stored classes per marker",
    "[bayes][genetic][marker_effects]")
{
    gelex::test::FileFixture fixture;
    const auto sparse_path
        = (fixture.get_test_dir() / "assign.draws.csc").string();
    const std::array<std::array<std::uint8_t, 3>, 4> draws{
        {{1, 0, 3}, {0, 0, 2}, {1, 0, 0}, {3, 0, 2}}};
    {
        gelex::CscWriter sparse{sparse_path};
        auto stream = sparse.reserve<std::uint8_t>(
            "assignment", gelex::BinaryShape{3, 4});
        for (const auto& draw : draws)
        {
            stream << std::span<const std::uint8_t>{draw};
        }
        sparse.close();
    }
    const gelex::CscReader sparse{sparse_path};

    const auto any = gelex::inclusion_probability(
        sparse, "assignment", [](std::uint8_t value) { return value != 0; });
    REQUIRE(any.isApprox(Eigen::VectorXd{{0.75, 0.0, 0.75}}));

    const auto class_two = gelex::inclusion_probability(
        sparse, "assignment", [](std::uint8_t value) { return value == 2; });
    REQUIRE(class_two.isApprox(Eigen::VectorXd{{0.0, 0.0, 0.5}}));
}

TEST_CASE(
    "MarkerEffectTable keeps uniquely named equal-length columns",
    "[bayes][genetic][marker_effects]")
{
    gelex::MarkerEffectTable table{2};
    table.add("BETA_A", Eigen::VectorXd{{1.0, 2.0}});
    REQUIRE(table.rows() == 2);
    REQUIRE(table.names().size() == 1);
    REQUIRE(table.column("BETA_A")(1) == 2.0);
    REQUIRE_THROWS_AS(
        table.add("SE_A", Eigen::VectorXd{{1.0}}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        table.add("BETA_A", Eigen::VectorXd{{1.0, 2.0}}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(table.column("PVE"), gelex::GelexException);
}

TEST_CASE(
    "Marker PVE scales squared means by column variance over the denominator",
    "[bayes][genetic][marker_effects]")
{
    const gelex::MarkerPveScale scale{Eigen::RowVectorXd{{0.5, 2.0}}, 4.0};
    REQUIRE(scale.pve(Eigen::VectorXd{{2.0, -1.0}})
                .isApprox(Eigen::VectorXd{{0.5, 0.5}}));
    REQUIRE_THROWS_AS(scale.pve(Eigen::VectorXd{{1.0}}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        gelex::MarkerPveScale(Eigen::RowVectorXd{{0.5}}, 0.0),
        gelex::GelexException);

    gelex::MarkerEffectTable table{2};
    gelex::append_marker_columns(
        table,
        gelex::MixtureMarkerEffects{
            .coefficients
            = {.mean = Eigen::VectorXd{{2.0, -1.0}},
               .sd = Eigen::VectorXd{{0.1, 0.2}}},
            .pip = Eigen::VectorXd{{1.0, 0.5}}},
        gelex::GeneticMode::D,
        scale);
    const std::vector<std::string> expected{"BETA_D", "SE_D", "PVE_D", "PIP_D"};
    REQUIRE(
        std::vector<std::string>(table.names().begin(), table.names().end())
        == expected);
    REQUIRE(table.column("PVE_D").isApprox(Eigen::VectorXd{{0.5, 0.5}}));
}
