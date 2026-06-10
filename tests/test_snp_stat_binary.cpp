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

#include "gelex/data/genotype/method.h"
#include "gelex/exception.h"
#include "gelex/io/locistats/reader.h"
#include "gelex/io/locistats/writer.h"

#include "file_fixture.h"
#include "gelex/types/genetic_effect_type.h"

using gelex::EffectType;
using gelex::GenotypeMethod;
using gelex::LociStatsReader;
using gelex::LociStatsWriter;
using gelex::GelexException;
using gelex::genotype_method_from_byte;

TEST_CASE("genotype method byte round-trip", "[sbin][snpstats]")
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
        REQUIRE(genotype_method_from_byte(std::to_underlying(method)) == method);
    }
}

TEST_CASE("invalid genotype method byte throws", "[sbin][snpstats]")
{
    constexpr std::array<uint8_t, 6> INVALID_BYTES{10, 11, 12, 13, 14, 15};

    for (uint8_t byte : INVALID_BYTES)
    {
        REQUIRE_THROWS_AS(genotype_method_from_byte(byte), GelexException);
    }
}

TEST_CASE("sbin round-trip additive only", "[sbin][snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int NUM_SNPS = 200;
    Eigen::VectorXd mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.05, 0.95);
    Eigen::VectorXd stddev = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.01, 0.50);
    std::vector<int64_t> mono = {3, 42, 101};

    const auto METHOD = GenotypeMethod::StandardizeHWE;

    auto sbin_path = dir / "test.sbin";
    {
        LociStatsWriter writer(sbin_path.string());
        writer.write(
            EffectType::add(),
            std::to_underlying(METHOD),
            mean,
            &stddev,
            mono);
    }

    LociStatsReader reader(sbin_path.string());

    REQUIRE(reader.has(EffectType::add()));
    REQUIRE_FALSE(reader.has(EffectType::dom()));

    auto data = reader.read(EffectType::add());

    REQUIRE(data.mean.size() == NUM_SNPS);
    REQUIRE(data.stddev.has_value());
    REQUIRE(data.mean.isApprox(mean));
    REQUIRE(data.stddev->isApprox(stddev));
    REQUIRE(data.mono_indices == mono);
    REQUIRE(data.method == METHOD);
}

TEST_CASE("sbin round-trip additive and dominance", "[sbin][snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int NUM_SNPS = 150;
    Eigen::VectorXd add_mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.1, 0.9);
    Eigen::VectorXd add_stddev
        = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.02, 0.30);
    std::vector<int64_t> add_mono = {7, 88};

    Eigen::VectorXd dom_mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.2, 0.8);
    Eigen::VectorXd dom_stddev
        = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.05, 0.25);

    const auto ADD_METHOD = GenotypeMethod::StandardizeHWE;
    const auto DOM_METHOD = GenotypeMethod::OrthStandardizeHWE;

    auto sbin_path = dir / "test_ad.sbin";
    {
        LociStatsWriter writer(sbin_path.string());
        writer.write(
            EffectType::add(),
            std::to_underlying(ADD_METHOD),
            add_mean,
            &add_stddev,
            add_mono);
        writer.write(
            EffectType::dom(),
            std::to_underlying(DOM_METHOD),
            dom_mean,
            &dom_stddev);
    }

    LociStatsReader reader(sbin_path.string());

    REQUIRE(reader.has(EffectType::add()));
    REQUIRE(reader.has(EffectType::dom()));

    auto add_data = reader.read(EffectType::add());
    REQUIRE(add_data.mean.isApprox(add_mean));
    REQUIRE(add_data.stddev.has_value());
    REQUIRE(add_data.stddev->isApprox(add_stddev));
    REQUIRE(add_data.mono_indices == add_mono);
    REQUIRE(add_data.method == ADD_METHOD);

    auto dom_data = reader.read(EffectType::dom());
    REQUIRE(dom_data.mean.isApprox(dom_mean));
    REQUIRE(dom_data.stddev.has_value());
    REQUIRE(dom_data.stddev->isApprox(dom_stddev));
    REQUIRE(dom_data.mono_indices.empty());
    REQUIRE(dom_data.method == DOM_METHOD);
}

TEST_CASE("sbin round-trip center only (no stddev)", "[sbin][snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int NUM_SNPS = 100;
    Eigen::VectorXd mean = Eigen::VectorXd::LinSpaced(NUM_SNPS, 0.1, 0.9);
    std::vector<int64_t> mono = {5, 50};

    auto sbin_path = dir / "test_center.sbin";
    {
        LociStatsWriter writer(sbin_path.string());
        writer.write(
            EffectType::add(),
            std::to_underlying(GenotypeMethod::CenterHWE),
            mean,
            nullptr,
            mono);
    }

    LociStatsReader reader(sbin_path.string());

    REQUIRE(reader.has(EffectType::add()));

    auto data = reader.read(EffectType::add());

    REQUIRE(data.mean.size() == NUM_SNPS);
    REQUIRE_FALSE(data.stddev.has_value());
    REQUIRE(data.mean.isApprox(mean));
    REQUIRE(data.mono_indices == mono);
    REQUIRE(data.method == GenotypeMethod::CenterHWE);
}
