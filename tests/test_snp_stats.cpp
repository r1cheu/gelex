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

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include <catch2/catch_test_macros.hpp>

#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_stats.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/snpstats.h"

#include "file_fixture.h"
#include "gelex/types/genetic_mode.h"

using gelex::BinaryReader;
using gelex::BinaryWriter;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::genotype_method_from_byte;
using gelex::GenotypeMethod;
using gelex::has_snp_stats;
using gelex::read_snp_stats;
using gelex::SnpStats;
using gelex::write_snp_stats;

TEST_CASE("genotype method byte round-trip", "[snpstats]")
{
    constexpr std::array METHODS{
        GenotypeMethod::StandardizeHWE,
        GenotypeMethod::CenterHWE,
        GenotypeMethod::Standardize,
        GenotypeMethod::Center,
        GenotypeMethod::OrthStandardizeHWE,
        GenotypeMethod::OrthCenterHWE,
        GenotypeMethod::OrthStandardize,
        GenotypeMethod::OrthCenter,
        GenotypeMethod::NOIAStandardize,
        GenotypeMethod::NOIACenter};

    for (auto method : METHODS)
    {
        REQUIRE(
            genotype_method_from_byte(std::to_underlying(method)) == method);
    }
}

TEST_CASE("invalid genotype method byte throws", "[snpstats]")
{
    constexpr std::array<uint8_t, 6> INVALID_BYTES{10, 11, 12, 13, 14, 15};

    for (uint8_t byte : INVALID_BYTES)
    {
        REQUIRE_THROWS_AS(genotype_method_from_byte(byte), GelexException);
    }
}

TEST_CASE("snpstats round-trip additive only", "[snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int NUM_SNPS = 200;
    SnpStats stats;
    stats.method = GenotypeMethod::StandardizeHWE;
    stats.mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.05, 0.95);
    stats.var = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.01, 0.50);
    stats.A1freq = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.03, 0.48);
    stats.valid_indices = {0, 1, 2, 99, 199};

    auto snpstats_path = dir / "test.snpstats";
    {
        BinaryWriter writer(snpstats_path.string());
        write_snp_stats(writer, GeneticMode::A, stats);
    }

    BinaryReader reader(snpstats_path.string());

    REQUIRE(has_snp_stats(reader, GeneticMode::A));
    REQUIRE_FALSE(has_snp_stats(reader, GeneticMode::D));

    auto data = read_snp_stats(reader, GeneticMode::A);

    REQUIRE(data.mean.size() == NUM_SNPS);
    REQUIRE(data.mean.isApprox(stats.mean));
    REQUIRE(data.var.isApprox(stats.var));
    REQUIRE(data.A1freq.isApprox(stats.A1freq));
    REQUIRE(data.valid_indices == stats.valid_indices);
    REQUIRE(data.method == stats.method);
}

TEST_CASE("snpstats round-trip additive and dominance", "[snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int NUM_SNPS = 150;
    SnpStats add;
    add.method = GenotypeMethod::StandardizeHWE;
    add.mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.1, 0.9);
    add.var = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.02, 0.30);
    add.A1freq = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.05, 0.45);
    add.valid_indices = {7, 88};

    SnpStats dom;
    dom.method = GenotypeMethod::OrthStandardizeHWE;
    dom.mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.2, 0.8);
    dom.var = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.05, 0.25);
    dom.A1freq = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.10, 0.40);

    auto snpstats_path = dir / "test_ad.snpstats";
    {
        BinaryWriter writer(snpstats_path.string());
        write_snp_stats(writer, GeneticMode::A, add);
        write_snp_stats(writer, GeneticMode::D, dom);
    }

    BinaryReader reader(snpstats_path.string());

    REQUIRE(has_snp_stats(reader, GeneticMode::A));
    REQUIRE(has_snp_stats(reader, GeneticMode::D));

    auto add_data = read_snp_stats(reader, GeneticMode::A);
    REQUIRE(add_data.mean.isApprox(add.mean));
    REQUIRE(add_data.var.isApprox(add.var));
    REQUIRE(add_data.A1freq.isApprox(add.A1freq));
    REQUIRE(add_data.valid_indices == add.valid_indices);
    REQUIRE(add_data.method == add.method);

    auto dom_data = read_snp_stats(reader, GeneticMode::D);
    REQUIRE(dom_data.mean.isApprox(dom.mean));
    REQUIRE(dom_data.var.isApprox(dom.var));
    REQUIRE(dom_data.valid_indices.empty());
    REQUIRE(dom_data.method == dom.method);
}

TEST_CASE("snpstats round-trip empty valid_indices", "[snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int NUM_SNPS = 100;
    SnpStats stats;
    stats.method = GenotypeMethod::CenterHWE;
    stats.mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.1, 0.9);
    stats.var = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.2, 0.6);
    stats.A1freq = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.05, 0.45);

    auto snpstats_path = dir / "test_center.snpstats";
    {
        BinaryWriter writer(snpstats_path.string());
        write_snp_stats(writer, GeneticMode::A, stats);
    }

    BinaryReader reader(snpstats_path.string());

    REQUIRE(has_snp_stats(reader, GeneticMode::A));

    auto data = read_snp_stats(reader, GeneticMode::A);

    REQUIRE(data.mean.size() == NUM_SNPS);
    REQUIRE(data.mean.isApprox(stats.mean));
    REQUIRE(data.var.isApprox(stats.var));
    REQUIRE(data.valid_indices.empty());
    REQUIRE(data.method == GenotypeMethod::CenterHWE);
}
