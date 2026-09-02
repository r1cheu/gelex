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
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

#include "file_fixture.h"

namespace
{

namespace fs = std::filesystem;
namespace test = gelex::test;

}  // namespace

TEST_CASE(
    "BinaryWriter supports interleaved payload appends",
    "[io][binary_writer]")
{
    test::FileFixture fixture;
    const auto container_path = fixture.get_test_dir() / "interleaved.samples";
    const Eigen::MatrixXd expected_double{
        {1.0, 4.0, 7.0}, {2.0, 5.0, 8.0}, {3.0, 6.0, 9.0}};
    const Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>
        expected_uint8{{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};

    {
        gelex::BinaryWriter writer(container_path.string());
        auto double_payload
            = writer.reserve<double>("double", gelex::BinaryShape{3, 3});
        auto uint8_payload
            = writer.reserve<std::uint8_t>("uint8", gelex::BinaryShape{3, 3});

        for (Eigen::Index column = 0; column < expected_double.cols(); ++column)
        {
            double_payload.append(expected_double.col(column));
            uint8_payload.append(expected_uint8.col(column));
        }
        writer.close();
    }

    gelex::BinaryReader reader(container_path.string());
    REQUIRE(reader.to_map<double>("double").isApprox(expected_double));
    REQUIRE(reader.to_map<std::uint8_t>("uint8").isApprox(expected_uint8));
}

TEST_CASE("BinaryWriter rejects duplicate identifiers", "[io][binary_writer]")
{
    test::FileFixture fixture;
    gelex::BinaryWriter writer(
        (fixture.get_test_dir() / "duplicate.samples").string());
    [[maybe_unused]] const auto reserved
        = writer.reserve<double>("value", gelex::BinaryShape{1, 1});

    REQUIRE_THROWS_AS(
        (writer.reserve<double>("value", gelex::BinaryShape{1, 1})),
        gelex::GelexException);
}

TEST_CASE("BinaryWriter rejects empty identifiers", "[io][binary_writer]")
{
    test::FileFixture fixture;
    gelex::BinaryWriter writer(
        (fixture.get_test_dir() / "empty_identifier.samples").string());

    REQUIRE_THROWS_AS(
        (writer.reserve<double>("", gelex::BinaryShape{1, 1})),
        gelex::GelexException);
}

TEST_CASE("BinaryWriter rejects payload overflow", "[io][binary_writer]")
{
    test::FileFixture fixture;
    gelex::BinaryWriter writer(
        (fixture.get_test_dir() / "overflow.samples").string());
    auto payload = writer.reserve<double>("value", gelex::BinaryShape{2, 1});
    const std::array<double, 2> values{1.0, 2.0};

    payload.append(std::span<const double>{values});
    REQUIRE_THROWS_AS(payload.append(values), gelex::GelexException);
}

TEST_CASE("PayloadWriter transfers ownership on move", "[io][binary_writer]")
{
    STATIC_REQUIRE(!std::is_copy_constructible_v<gelex::PayloadWriter<double>>);
    STATIC_REQUIRE(
        std::is_nothrow_move_constructible_v<gelex::PayloadWriter<double>>);

    test::FileFixture fixture;
    gelex::BinaryWriter writer(
        (fixture.get_test_dir() / "move.samples").string());
    auto payload = writer.reserve<double>("value", gelex::BinaryShape{1, 1});
    auto moved_payload = std::move(payload);

    REQUIRE_THROWS_AS(payload.append(1.0), gelex::GelexException);
    moved_payload.append(1.0);
    writer.close();
}

TEST_CASE("BinaryWriter commits incomplete payloads", "[io][binary_writer]")
{
    test::FileFixture fixture;
    const auto container_path = fixture.get_test_dir() / "incomplete.samples";
    std::vector<std::string> warnings;
    gelex::set_sink(
        [&warnings](gelex::Level level, std::string_view message)
        {
            if (level == gelex::Level::Warn)
            {
                warnings.emplace_back(message);
            }
        });
    {
        gelex::BinaryWriter writer(container_path.string());
        auto full = writer.reserve<double>("full", gelex::BinaryShape{1, 1});
        auto partial
            = writer.reserve<double>("partial", gelex::BinaryShape{2, 3});
        full.append(3.0);
        partial.append(Eigen::VectorXd{{1.0, 2.0}});
        writer.close();
    }
    gelex::set_sink({});

    REQUIRE(warnings.size() == 1);
    CHECK_THAT(
        warnings.front(),
        Catch::Matchers::ContainsSubstring("1 of 2 payload(s)"));

    const gelex::BinaryReader reader(container_path.string());
    REQUIRE(reader.to_map<double>("full").isApprox(Eigen::MatrixXd{{3.0}}));
    REQUIRE(reader.to_map<double>("partial").isApprox(
        Eigen::MatrixXd{{1.0}, {2.0}}));
}

TEST_CASE("BinaryWriter rejects operations after close", "[io][binary_writer]")
{
    test::FileFixture fixture;
    gelex::BinaryWriter writer(
        (fixture.get_test_dir() / "closed.samples").string());
    auto payload = writer.reserve<double>("value", gelex::BinaryShape{1, 1});
    payload.append(1.0);
    writer.close();

    REQUIRE_THROWS_AS(payload.append(2.0), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (writer.reserve<double>("other", gelex::BinaryShape{1, 1})),
        gelex::GelexException);
    REQUIRE_THROWS_AS(writer.close(), gelex::GelexException);
}

TEST_CASE("BinaryWriter without close discards output", "[io][binary_writer]")
{
    test::FileFixture fixture;
    const auto container_path = fixture.get_test_dir() / "discarded.samples";
    std::vector<std::string> warnings;
    gelex::set_sink(
        [&warnings](gelex::Level level, std::string_view message)
        {
            if (level == gelex::Level::Warn)
            {
                warnings.emplace_back(message);
            }
        });
    {
        gelex::BinaryWriter writer(container_path.string());
        writer.reserve<double>("value", gelex::BinaryShape{1, 1}).append(1.0);
    }
    gelex::set_sink({});

    REQUIRE_FALSE(fs::exists(container_path));
    REQUIRE(warnings.size() == 1);
    CHECK_THAT(
        warnings.front(),
        Catch::Matchers::ContainsSubstring("destroyed without close()"));
}
