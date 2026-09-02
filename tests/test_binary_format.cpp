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
#include <algorithm>
#include <array>
#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

#include "file_fixture.h"

namespace
{

namespace fs = std::filesystem;
using namespace gelex;

}  // namespace

TEST_CASE("Container round-trip with mixed dtypes", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    Eigen::MatrixXd coeffs(3, 5);
    coeffs << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0;

    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic> tracker(3, 5);
    tracker << 0, 1, 2, 0, 1, 1, 2, 0, 1, 2, 2, 0, 1, 2, 0;

    Eigen::MatrixXd scalars(3, 5);
    scalars << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3,
        1.4, 1.5;

    auto container_path = dir / "test.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Additive/coeff", coeffs);
        writer.write("Additive/group/assignment", tracker);
        writer.write("Residual/variance", scalars);
    }

    // Read back
    BinaryReader reader(container_path.string());
    REQUIRE(reader.n_sections() == 3);

    SECTION("double section")
    {
        auto mat = reader.to_map<double>("Additive/coeff");
        REQUIRE(mat.rows() == 3);
        REQUIRE(mat.cols() == 5);
        REQUIRE(mat.isApprox(coeffs));
    }

    SECTION("int8_t section")
    {
        auto mat = reader.to_map<int8_t>("Additive/group/assignment");
        REQUIRE(mat.rows() == 3);
        REQUIRE(mat.cols() == 5);
        REQUIRE(mat == tracker);
    }

    SECTION("scalar section")
    {
        auto mat = reader.to_mat<double>("Residual/variance");
        REQUIRE(mat.rows() == 3);
        REQUIRE(mat.cols() == 5);
        REQUIRE(mat.isApprox(scalars));
    }
}

TEST_CASE("Container contains", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    Eigen::MatrixXd data(2, 2);
    data << 1.0, 2.0, 3.0, 4.0;

    auto container_path = dir / "test.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Additive/coeff", data);
    }

    BinaryReader reader(container_path.string());
    REQUIRE(reader.contains("Additive/coeff"));
    REQUIRE_FALSE(reader.contains("Dominance/coeff"));
}

TEST_CASE("Container scalar round-trip", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    auto container_path = dir / "scalar.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Residual/variance", 0.75);
        writer.write("Run/count", uint32_t{1});
    }

    BinaryReader reader(container_path.string());
    REQUIRE(reader.n_sections() == 2);

    const auto variance = reader.to_mat<double>("Residual/variance");
    REQUIRE(variance.rows() == 1);
    REQUIRE(variance.cols() == 1);
    REQUIRE_THAT(variance(0, 0), Catch::Matchers::WithinRel(0.75));

    const auto count = reader.to_mat<uint32_t>("Run/count");
    REQUIRE(count.rows() == 1);
    REQUIRE(count.cols() == 1);
    REQUIRE(count(0, 0) == 1);
}

TEST_CASE("Container section_paths", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    Eigen::MatrixXd data(2, 2);
    data << 1.0, 2.0, 3.0, 4.0;

    auto container_path = dir / "test.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Additive/coeff", data);
        writer.write("Additive/variance", data);
        writer.write("Residual/variance", data);
    }

    BinaryReader reader(container_path.string());
    auto paths = reader.section_paths();

    REQUIRE(paths.size() == 3);

    std::vector<std::string> sorted(paths.begin(), paths.end());
    std::sort(sorted.begin(), sorted.end());

    REQUIRE(sorted[0] == "Additive/coeff");
    REQUIRE(sorted[1] == "Additive/variance");
    REQUIRE(sorted[2] == "Residual/variance");
}

TEST_CASE("Container section not found throws", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    Eigen::MatrixXd data(2, 2);
    data << 1.0, 2.0, 3.0, 4.0;

    auto container_path = dir / "test.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Additive/coeff", data);
    }

    BinaryReader reader(container_path.string());
    REQUIRE_THROWS_AS(reader.to_map<double>("Dominance/coeff"), GelexException);
}

TEST_CASE("Container dtype mismatch throws", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    Eigen::MatrixXd data(2, 2);
    data << 1.0, 2.0, 3.0, 4.0;

    auto container_path = dir / "test.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Additive/coeff", data);
    }

    BinaryReader reader(container_path.string());
    REQUIRE_THROWS_AS(reader.to_map<int8_t>("Additive/coeff"), GelexException);
}

TEST_CASE("Container invalid magic throws", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    auto bad_path = dir / "bad.samples";
    std::ofstream out(bad_path, std::ios::binary);
    std::array<char, 64> garbage{};
    out.write(garbage.data(), static_cast<std::streamsize>(garbage.size()));
    out.close();

    REQUIRE_THROWS_AS(BinaryReader(bad_path.string()), GelexException);
}

TEST_CASE("Container file too small throws", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    auto tiny_path = dir / "tiny.samples";
    std::ofstream out(tiny_path, std::ios::binary);
    std::array<char, 10> data{};
    out.write(data.data(), static_cast<std::streamsize>(data.size()));
    out.close();

    REQUIRE_THROWS_AS(BinaryReader(tiny_path.string()), GelexException);
}

TEST_CASE("Container duplicate section path throws", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    BinaryWriter writer((dir / "test.samples").string());
    writer.reserve<double>("Additive/coeff", 2, 2);
    REQUIRE_THROWS_AS(
        (writer.reserve<double>("Additive/coeff", 2, 2)), GelexException);
}

TEST_CASE("Container sections with different indices", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    Eigen::MatrixXd random0(2, 4);
    random0 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0;

    Eigen::MatrixXd random1(3, 4);
    random1 << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2;

    Eigen::MatrixXd var0(1, 4);
    var0 << 0.5, 0.6, 0.7, 0.8;

    Eigen::MatrixXd var1(1, 4);
    var1 << 1.5, 1.6, 1.7, 1.8;

    auto container_path = dir / "multi_index.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Random/coeff/0", random0);
        writer.write("Random/coeff/1", random1);
        writer.write("Random/variance/0", var0);
        writer.write("Random/variance/1", var1);
    }

    BinaryReader reader(container_path.string());
    REQUIRE(reader.n_sections() == 4);

    REQUIRE(reader.contains("Random/coeff/0"));
    REQUIRE(reader.contains("Random/coeff/1"));
    REQUIRE_FALSE(reader.contains("Random/coeff/2"));

    SECTION("index 0 coeffs")
    {
        auto mat = reader.to_map<double>("Random/coeff/0");
        REQUIRE(mat.rows() == 2);
        REQUIRE(mat.cols() == 4);
        REQUIRE(mat.isApprox(random0));
    }

    SECTION("index 1 coeffs")
    {
        auto mat = reader.to_map<double>("Random/coeff/1");
        REQUIRE(mat.rows() == 3);
        REQUIRE(mat.cols() == 4);
        REQUIRE(mat.isApprox(random1));
    }

    SECTION("index 0 variance")
    {
        auto mat = reader.to_map<double>("Random/variance/0");
        REQUIRE(mat.rows() == 1);
        REQUIRE(mat.cols() == 4);
        REQUIRE(mat.isApprox(var0));
    }

    SECTION("index 1 variance")
    {
        auto mat = reader.to_map<double>("Random/variance/1");
        REQUIRE(mat.rows() == 1);
        REQUIRE(mat.cols() == 4);
        REQUIRE(mat.isApprox(var1));
    }
}

TEST_CASE("Reserve round-trip with column-wise writes", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    constexpr int rows = 5;
    constexpr int cols = 10;

    Eigen::MatrixXd expected_double(rows, cols);
    for (int c = 0; c < cols; ++c)
    {
        for (int r = 0; r < rows; ++r)
        {
            expected_double(r, c) = c * 100.0 + r;
        }
    }

    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic> expected_int8(
        rows, cols);
    for (int c = 0; c < cols; ++c)
    {
        for (int r = 0; r < rows; ++r)
        {
            expected_int8(r, c) = static_cast<int8_t>((c + r) % 3);
        }
    }

    auto container_path = dir / "reserve.samples";
    {
        BinaryWriter writer(container_path.string());
        auto h_double = writer.reserve<double>("Additive/coeff", rows, cols);
        auto h_int8
            = writer.reserve<int8_t>("Additive/group/assignment", rows, cols);

        for (int c = 0; c < cols; ++c)
        {
            writer.write(h_double, expected_double.col(c));
            writer.write(h_int8, expected_int8.col(c));
        }
    }

    BinaryReader reader(container_path.string());
    REQUIRE(reader.n_sections() == 2);

    auto mat_d = reader.to_map<double>("Additive/coeff");
    REQUIRE(mat_d.rows() == rows);
    REQUIRE(mat_d.cols() == cols);
    REQUIRE(mat_d.isApprox(expected_double));

    auto mat_i = reader.to_map<int8_t>("Additive/group/assignment");
    REQUIRE(mat_i.rows() == rows);
    REQUIRE(mat_i.cols() == cols);
    REQUIRE(mat_i == expected_int8);
}

TEST_CASE("Reserve write overflow throws", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    auto container_path = dir / "overflow.samples";
    BinaryWriter writer(container_path.string());
    auto h = writer.reserve<double>("Additive/coeff", 2, 1);

    Eigen::Vector3d too_large;
    too_large << 1.0, 2.0, 3.0;
    REQUIRE_THROWS_AS((writer.write(h, too_large)), GelexException);
}

TEST_CASE("Container accepts long section paths", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    const std::string long_path
        = "state/genetic_0/single/A/prior_state/"
          "scaled_mixture_gaussian/mixture/proportion";
    const Eigen::MatrixXd expected{
        {0.1, 0.4}, {0.2, 0.3}, {0.3, 0.2}, {0.4, 0.1}};

    auto container_path = dir / "long_path.draws";
    {
        BinaryWriter writer(container_path.string());
        writer.write(long_path, expected);
    }

    BinaryReader reader(container_path.string());
    REQUIRE(reader.contains(long_path));
    REQUIRE(reader.to_map<double>(long_path).isApprox(expected));
}

TEST_CASE(
    "BinaryWriter close reports incomplete sections",
    "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    auto container_path = dir / "incomplete.draws";
    BinaryWriter writer(container_path.string());
    writer.reserve<double>("state/fixed/coeffs", 1, 2);

    REQUIRE_THROWS_AS(writer.close(), GelexException);
}

TEST_CASE("BinaryWriter close commits complete container", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    auto container_path = dir / "closed.draws";
    BinaryWriter writer(container_path.string());
    writer.write("state/fixed/coeffs", Eigen::MatrixXd{{1.0, 2.0}});
    writer.close();

    BinaryReader reader(container_path.string());
    REQUIRE(reader.to_map<double>("state/fixed/coeffs")
                .isApprox(Eigen::MatrixXd{{1.0, 2.0}}));
}

TEST_CASE("Container string section round-trip", "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    const std::vector<std::string_view> names
        = {"intercept", "sex", "age_group", "pop_stratum"};

    Eigen::MatrixXd coeff(2, 3);
    coeff << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    auto container_path = dir / "string_test.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write_strings("Fixed/names", names);
        writer.write("Fixed/coeff", coeff);
    }

    BinaryReader reader(container_path.string());
    REQUIRE(reader.n_sections() == 2);

    SECTION("string section round-trips correctly")
    {
        auto result = reader.to_strings("Fixed/names");
        REQUIRE(result.size() == names.size());
        for (size_t i = 0; i < names.size(); ++i)
        {
            REQUIRE(result[i] == names[i]);
        }
    }

    SECTION("numeric section still intact")
    {
        auto mat = reader.to_map<double>("Fixed/coeff");
        REQUIRE(mat.rows() == 2);
        REQUIRE(mat.cols() == 3);
        REQUIRE(mat.isApprox(coeff));
    }
}

TEST_CASE(
    "Container to_strings on non-string section throws",
    "[binary_container]")
{
    test::FileFixture fixture;
    const auto& dir = fixture.get_test_dir();

    Eigen::MatrixXd data(2, 2);
    data << 1.0, 2.0, 3.0, 4.0;

    auto container_path = dir / "wrong_dtype.samples";
    {
        BinaryWriter writer(container_path.string());
        writer.write("Additive/coeff", data);
    }

    BinaryReader reader(container_path.string());
    REQUIRE_THROWS_AS(reader.to_strings("Additive/coeff"), GelexException);
}
