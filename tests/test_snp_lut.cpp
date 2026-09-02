/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <span>

#include "gelex/data/snp_lut.h"
#include "gelex/exception.h"
#include "gelex/io/binary_writer.h"
#include "gelex/io/snp_lut.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::load_snp_luts;
using gelex::SnpLutMatrix;
using gelex::write_snp_luts;

TEST_CASE("SNP LUT round-trip preserves missing values", "[snp_lut]")
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
    gelex::BinaryWriter writer(path.string());
    writer
        .reserve<double>(
            "A/lut",
            gelex::BinaryShape{
                static_cast<std::uint64_t>(add.rows()),
                static_cast<std::uint64_t>(add.cols())})
        .write(
            std::span<const double>{
                add.data(), static_cast<std::size_t>(add.size())});
    writer
        .reserve<double>(
            "D/lut",
            gelex::BinaryShape{
                static_cast<std::uint64_t>(dom.rows()),
                static_cast<std::uint64_t>(dom.cols())})
        .write(
            std::span<const double>{
                dom.data(), static_cast<std::size_t>(dom.size())});
    writer.close();

    const auto actual = load_snp_luts(path);
    REQUIRE(actual.size() == 2);
    REQUIRE(actual.at(GeneticMode::A).isApprox(add));
    REQUIRE(actual.at(GeneticMode::D).isApprox(dom));
    REQUIRE(actual.at(GeneticMode::A).row(1).isApprox(add.row(1)));
}

TEST_CASE("load_snp_luts rejects invalid LUT rows", "[snp_lut]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "invalid_rows.snplut";
    const Eigen::MatrixXd invalid = Eigen::MatrixXd::Zero(3, 2);

    gelex::BinaryWriter writer(path.string());
    writer
        .reserve<double>(
            "A/lut",
            gelex::BinaryShape{
                static_cast<std::uint64_t>(invalid.rows()),
                static_cast<std::uint64_t>(invalid.cols())})
        .write(
            std::span<const double>{
                invalid.data(), static_cast<std::size_t>(invalid.size())});
    writer.close();

    REQUIRE_THROWS_AS(load_snp_luts(path), GelexException);
}

TEST_CASE("SNP LUT writer persists a genetic design", "[snp_lut]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "genetic_design.snplut";
    auto design = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 2.0}, {2.0, 0.0}},
        GeneticMode::A | GeneticMode::D);

    write_snp_luts(path, design);

    const auto actual = load_snp_luts(path);
    REQUIRE(actual.size() == 2);
    REQUIRE(actual.at(GeneticMode::A)
                .isApprox(design.projection(GeneticMode::A).snp_luts()));
    REQUIRE(actual.at(GeneticMode::D)
                .isApprox(design.projection(GeneticMode::D).snp_luts()));
}
