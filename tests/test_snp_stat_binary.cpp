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

#include <cstdint>
#include <vector>

#include <Eigen/Core>

#include <catch2/catch_test_macros.hpp>

#include "gelex/io/binary_format.h"
#include "gelex/io/sbin_reader.h"
#include "gelex/io/sbin_writer.h"

#include "file_fixture.h"

using gelex::SbinReader;
using gelex::SbinWriter;
using gelex::detail::EffectType;

TEST_CASE("sbin round-trip additive only", "[sbin][snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int kNumSnps = 200;
    Eigen::VectorXd mean = Eigen::VectorXd::LinSpaced(kNumSnps, 0.05, 0.95);
    Eigen::VectorXd stddev = Eigen::VectorXd::LinSpaced(kNumSnps, 0.01, 0.50);
    std::vector<int64_t> mono = {3, 42, 101};

    auto sbin_path = dir / "test.sbin";
    {
        SbinWriter writer(sbin_path.string());
        writer.write(EffectType::Add, mean, &stddev, mono);
        writer.finalize();
    }

    SbinReader reader(sbin_path.string());

    REQUIRE(reader.has(EffectType::Add));
    REQUIRE_FALSE(reader.has(EffectType::Dom));

    auto data = reader.read(EffectType::Add);

    REQUIRE(data.mean.size() == kNumSnps);
    REQUIRE(data.stddev.has_value());
    REQUIRE(data.mean.isApprox(mean));
    REQUIRE(data.stddev->isApprox(stddev));
    REQUIRE(data.mono_indices == mono);
}

TEST_CASE("sbin round-trip additive and dominance", "[sbin][snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int kNumSnps = 150;
    Eigen::VectorXd add_mean = Eigen::VectorXd::LinSpaced(kNumSnps, 0.1, 0.9);
    Eigen::VectorXd add_stddev
        = Eigen::VectorXd::LinSpaced(kNumSnps, 0.02, 0.30);
    std::vector<int64_t> add_mono = {7, 88};

    Eigen::VectorXd dom_mean = Eigen::VectorXd::LinSpaced(kNumSnps, 0.2, 0.8);
    Eigen::VectorXd dom_stddev
        = Eigen::VectorXd::LinSpaced(kNumSnps, 0.05, 0.25);

    auto sbin_path = dir / "test_ad.sbin";
    {
        SbinWriter writer(sbin_path.string());
        writer.write(EffectType::Add, add_mean, &add_stddev, add_mono);
        writer.write(EffectType::Dom, dom_mean, &dom_stddev);
        writer.finalize();
    }

    SbinReader reader(sbin_path.string());

    REQUIRE(reader.has(EffectType::Add));
    REQUIRE(reader.has(EffectType::Dom));

    auto add_data = reader.read(EffectType::Add);
    REQUIRE(add_data.mean.isApprox(add_mean));
    REQUIRE(add_data.stddev.has_value());
    REQUIRE(add_data.stddev->isApprox(add_stddev));
    REQUIRE(add_data.mono_indices == add_mono);

    auto dom_data = reader.read(EffectType::Dom);
    REQUIRE(dom_data.mean.isApprox(dom_mean));
    REQUIRE(dom_data.stddev.has_value());
    REQUIRE(dom_data.stddev->isApprox(dom_stddev));
    REQUIRE(dom_data.mono_indices.empty());
}

TEST_CASE("sbin round-trip center only (no stddev)", "[sbin][snpstats]")
{
    gelex::test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int kNumSnps = 100;
    Eigen::VectorXd mean = Eigen::VectorXd::LinSpaced(kNumSnps, 0.1, 0.9);
    std::vector<int64_t> mono = {5, 50};

    auto sbin_path = dir / "test_center.sbin";
    {
        SbinWriter writer(sbin_path.string());
        writer.write(EffectType::Add, mean, nullptr, mono);
        writer.finalize();
    }

    SbinReader reader(sbin_path.string());

    REQUIRE(reader.has(EffectType::Add));

    auto data = reader.read(EffectType::Add);

    REQUIRE(data.mean.size() == kNumSnps);
    REQUIRE_FALSE(data.stddev.has_value());
    REQUIRE(data.mean.isApprox(mean));
    REQUIRE(data.mono_indices == mono);
}
