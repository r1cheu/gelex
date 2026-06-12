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

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <type_traits>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "bed_fixture.h"
#include "gelex/data/detail/bed_source.h"
#include "gelex/data/detail/index_projection.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using gelex::GelexException;
using gelex::detail::BedSource;
using gelex::detail::open_bed_source;
using gelex::test::BedFixture;

static_assert(!std::is_copy_constructible_v<gelex::detail::BedSource>);
static_assert(!std::is_copy_assignable_v<gelex::detail::BedSource>);
static_assert(std::is_nothrow_move_constructible_v<gelex::detail::BedSource>);
static_assert(std::is_nothrow_move_assignable_v<gelex::detail::BedSource>);
static_assert(!std::is_copy_constructible_v<gelex::detail::IndexProjection>);
static_assert(!std::is_copy_assignable_v<gelex::detail::IndexProjection>);
static_assert(
    std::is_nothrow_move_constructible_v<gelex::detail::IndexProjection>);
static_assert(
    std::is_nothrow_move_assignable_v<gelex::detail::IndexProjection>);

TEST_CASE(
    "BedSource exposes contiguous BED variant bytes",
    "[data][bed_source]")
{
    BedFixture fixture;

    const auto num_samples = static_cast<std::size_t>(5);
    const auto num_variants = static_cast<std::size_t>(4);
    const auto bed_prefix = fixture
                                .create_bed_files(
                                    static_cast<Eigen::Index>(num_samples),
                                    static_cast<Eigen::Index>(num_variants),
                                    0.0)
                                .first;

    BedSource source{
        bed_prefix.string() + ".bed",
        num_variants,
        num_samples,
    };

    REQUIRE(source.size() == num_variants);
    REQUIRE(source.num_samples() == num_samples);
    REQUIRE(source.stride() == 2);
    REQUIRE(source.size_bytes() == num_variants * source.stride());
    REQUIRE(source.bytes().size() == source.size_bytes());
    REQUIRE(source.data() == source.bytes().data());

    REQUIRE(source[0].data() == source.data());
    REQUIRE(source[0].size() == source.stride());
    REQUIRE(source[1].data() == source.data() + source.stride());
    REQUIRE(source.at(2).data() == source.data() + (2 * source.stride()));

    const auto chunk = source.subspan(1, 2);
    REQUIRE(chunk.data() == source.data() + source.stride());
    REQUIRE(chunk.size() == 2 * source.stride());
}

TEST_CASE("BedSource rejects invalid shape and ranges", "[data][bed_source]")
{
    BedFixture fixture;
    const auto bed_prefix = fixture.create_bed_files(5, 4, 0.0).first;
    const auto bed_path = bed_prefix.string() + ".bed";

    REQUIRE_THROWS_AS(
        static_cast<void>(BedSource{bed_path, 0, 5}), GelexException);
    REQUIRE_THROWS_AS(
        static_cast<void>(BedSource{bed_path, 4, 0}), GelexException);

    BedSource source{bed_path, 4, 5};

    REQUIRE_THROWS_AS(source.at(4), GelexException);
    REQUIRE_THROWS_AS(source.subspan(4, 1), GelexException);
    REQUIRE_THROWS_AS(source.subspan(2, 3), GelexException);
}

TEST_CASE(
    "open_bed_source infers shape from PLINK sidecar files",
    "[data][bed_source]")
{
    BedFixture fixture;

    const auto num_samples = static_cast<std::size_t>(7);
    const auto num_variants = static_cast<std::size_t>(3);
    const auto bed_prefix = fixture
                                .create_bed_files(
                                    static_cast<Eigen::Index>(num_samples),
                                    static_cast<Eigen::Index>(num_variants),
                                    0.0)
                                .first;

    auto source = open_bed_source(bed_prefix.string());

    REQUIRE(source.size() == num_variants);
    REQUIRE(source.num_samples() == num_samples);
    REQUIRE(source.stride() == 2);
    REQUIRE(source.size_bytes() == num_variants * source.stride());
}

TEST_CASE("open_bed_source reports missing PLINK files", "[data][bed_source]")
{
    SECTION("missing fam")
    {
        BedFixture fixture;
        const auto bed_prefix = fixture.create_bed_files(5, 4, 0.0).first;

        fs::remove(bed_prefix.string() + ".fam");

        REQUIRE_THROWS_AS(open_bed_source(bed_prefix.string()), GelexException);
    }

    SECTION("missing bim")
    {
        BedFixture fixture;
        const auto bed_prefix = fixture.create_bed_files(5, 4, 0.0).first;

        fs::remove(bed_prefix.string() + ".bim");

        REQUIRE_THROWS_AS(open_bed_source(bed_prefix.string()), GelexException);
    }

    SECTION("missing bed")
    {
        BedFixture fixture;
        const auto bed_prefix = fixture.create_bed_files(5, 4, 0.0).first;

        fs::remove(bed_prefix.string() + ".bed");

        REQUIRE_THROWS_AS(open_bed_source(bed_prefix.string()), GelexException);
    }
}

TEST_CASE(
    "BedSource validates BED header and payload size",
    "[data][bed_source]")
{
    SECTION("invalid magic")
    {
        BedFixture fixture;
        const auto bed_prefix = fixture.create_bed_files(5, 4, 0.0).first;

        {
            std::fstream bed_file(
                bed_prefix.string() + ".bed",
                std::ios::in | std::ios::out | std::ios::binary);
            bed_file.seekp(0);
            bed_file.put(0x00);
            bed_file.put(0x00);
            bed_file.put(0x00);
        }

        REQUIRE_THROWS_AS(open_bed_source(bed_prefix.string()), GelexException);
    }

    SECTION("truncated payload")
    {
        BedFixture fixture;
        const auto bed_prefix = fixture.create_bed_files(5, 4, 0.0).first;
        const auto bed_path = bed_prefix.string() + ".bed";

        fs::resize_file(bed_path, fs::file_size(bed_path) - 1);

        REQUIRE_THROWS_AS(open_bed_source(bed_prefix.string()), GelexException);
    }

    SECTION("extra payload")
    {
        BedFixture fixture;
        const auto bed_prefix = fixture.create_bed_files(5, 4, 0.0).first;

        {
            std::ofstream bed_file(
                bed_prefix.string() + ".bed", std::ios::binary | std::ios::app);
            bed_file.put(0x00);
        }

        REQUIRE_THROWS_AS(open_bed_source(bed_prefix.string()), GelexException);
    }
}
