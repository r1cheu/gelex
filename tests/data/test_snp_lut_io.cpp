// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>

#include "gelex/data/snp_lut.h"
#include "gelex/data/snp_lut_io.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/dense_writer.h"

#include "file_fixture.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::load_snp_luts;
using gelex::SnpLutMatrix;
using gelex::write_snp_luts;

TEST_CASE("SNP LUT round-trip preserves missing values", "[data][snp_lut][io]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr Eigen::Index num_snps = 150;
    SnpLutMatrix add = Eigen::VectorXd::LinSpaced(4 * num_snps, 0.1, 0.9)
                           .reshaped(4, num_snps)
                           .array();
    SnpLutMatrix dom = Eigen::VectorXd::LinSpaced(4 * num_snps, -0.8, 0.8)
                           .reshaped(4, num_snps)
                           .array();

    const auto path = dir / "test_ad.snplut";
    const gelex::ModeMap<SnpLutMatrix> expected{
        {GeneticMode::A, add}, {GeneticMode::D, dom}};
    write_snp_luts(path, expected);

    const auto actual = load_snp_luts(path);
    REQUIRE(actual.size() == 2);
    REQUIRE(actual.at(GeneticMode::A).isApprox(expected.at(GeneticMode::A)));
    REQUIRE(actual.at(GeneticMode::D).isApprox(expected.at(GeneticMode::D)));
    REQUIRE(actual.at(GeneticMode::A).row(1).isApprox(add.row(1)));
}

TEST_CASE("load_snp_luts rejects invalid LUT rows", "[data][snp_lut][io]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "invalid_rows.snplut";
    const Eigen::MatrixXd invalid = Eigen::MatrixXd::Zero(3, 2);

    {
        auto writer = gelex::open_dense_writer(path.string());
        writer.reserve<double>(
            "A/lut",
            gelex::BinaryShape{
                static_cast<std::uint64_t>(invalid.rows()),
                static_cast<std::uint64_t>(invalid.cols())})
            << invalid.reshaped();
        writer.close();
    }

    REQUIRE_THROWS_AS(load_snp_luts(path), GelexException);
}
