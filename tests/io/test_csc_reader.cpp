// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <string>
#include <type_traits>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/detail/binary_wire.h"

#include "file_fixture.h"

namespace
{

namespace fs = std::filesystem;
namespace test = gelex::test;

auto read_file(const fs::path& path) -> std::vector<char>
{
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input), {}};
}

auto write_file(const fs::path& path, const std::vector<char>& bytes) -> void
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

template <typename T>
auto patch(std::vector<char>& bytes, std::size_t offset, T value) -> void
{
    std::memcpy(bytes.data() + offset, &value, sizeof(value));
}

// Writes one 4x3 matrix "x".
auto write_sample(const fs::path& path) -> void
{
    auto writer = gelex::open_csc_writer(path.string());
    auto stream = writer.reserve<double>("x", {4, 3});
    const Eigen::MatrixXd values{{0, 0, 3}, {1, 0, 0}, {0, 0, 4}, {-2.5, 0, 0}};
    for (Eigen::Index column = 0; column < values.cols(); ++column)
    {
        stream << values.col(column);
    }
    writer.close();
}

}  // namespace

static_assert(!std::is_copy_constructible_v<gelex::CscReader>);
static_assert(std::is_nothrow_move_constructible_v<gelex::CscReader>);

TEMPLATE_TEST_CASE(
    "CscReader reads supported dtypes",
    "[io][csc_reader]",
    double,
    float,
    std::uint8_t)
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "dtype.csc";
    const Eigen::MatrixX<TestType> expected{
        {TestType{0}, TestType{2}, TestType{0}},
        {TestType{4}, TestType{0}, TestType{0}}};
    {
        auto writer = gelex::open_csc_writer(path.string());
        auto stream = writer.reserve<TestType>("values", {2, 3});
        for (Eigen::Index column = 0; column < expected.cols(); ++column)
        {
            stream << expected.col(column);
        }
        writer.close();
    }

    const gelex::CscReader reader(path.string());
    REQUIRE(reader.size() == 1);
    const auto& header = reader.info("values");
    REQUIRE(header.type == gelex::detail::binary_type_for<TestType>);
    REQUIRE(header.shape == (gelex::BinaryShape{2, 3}));
    REQUIRE(reader.nnz("values") == 2);

    const auto map = reader.to_map<TestType>("values");
    REQUIRE(map.rows() == 2);
    REQUIRE(map.cols() == 3);
    REQUIRE(map.nonZeros() == 2);
    REQUIRE(Eigen::MatrixX<TestType>{map} == expected);
    const auto matrix = reader.to_mat<TestType>("values");
    REQUIRE(matrix.nonZeros() == 2);
    REQUIRE(Eigen::MatrixX<TestType>{matrix} == expected);
}

TEST_CASE("CscReader exposes matrix metadata", "[io][csc_reader]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "metadata.csc";
    {
        auto writer = gelex::open_csc_writer(path.string());
        auto zeta = writer.reserve<double>("zeta", {1, 1});
        auto alpha = writer.reserve<float>("alpha", {1, 1});
        auto beta = writer.reserve<std::uint8_t>("beta", {1, 1});
        zeta << Eigen::VectorXd{{3.0}};
        alpha << Eigen::VectorXf{{0.0F}};
        beta << Eigen::VectorX<std::uint8_t>{{2}};
        writer.close();
    }

    const gelex::CscReader reader(path.string());
    REQUIRE(reader.size() == 3);
    REQUIRE(reader.contains("alpha"));
    REQUIRE_FALSE(reader.contains("missing"));
    REQUIRE(reader.nnz("alpha") == 0);
    REQUIRE(reader.nnz("beta") == 1);

    const auto payloads = reader.payloads();
    REQUIRE(payloads.size() == 3);
    REQUIRE(payloads[0].identifier == "alpha");
    REQUIRE(payloads[1].identifier == "beta");
    REQUIRE(payloads[2].identifier == "zeta");
    REQUIRE(payloads[1].type == gelex::BinaryType::uint8);
    REQUIRE_THROWS_AS(reader.info("missing"), gelex::GelexException);
    REQUIRE_THROWS_AS(reader.to_map<double>("missing"), gelex::GelexException);
}

TEST_CASE("CscReader reads empty and degenerate matrices", "[io][csc_reader]")
{
    test::FileFixture fixture;
    const auto empty_path = fixture.get_test_dir() / "empty.csc";
    {
        auto writer = gelex::open_csc_writer(empty_path.string());
        writer.close();
    }
    REQUIRE(gelex::CscReader(empty_path.string()).size() == 0);

    const auto path = fixture.get_test_dir() / "degenerate.csc";
    {
        auto writer = gelex::open_csc_writer(path.string());
        [[maybe_unused]] auto no_columns
            = writer.reserve<double>("no_columns", {3, 0});
        auto no_rows = writer.reserve<double>("no_rows", {0, 2});
        no_rows << Eigen::VectorXd{} << Eigen::VectorXd{};
        auto zeros = writer.reserve<double>("zeros", {2, 2});
        zeros << Eigen::VectorXd::Zero(2) << Eigen::VectorXd::Zero(2);
        writer.close();
    }
    const gelex::CscReader reader(path.string());
    REQUIRE(reader.to_mat<double>("no_columns").rows() == 3);
    REQUIRE(reader.to_mat<double>("no_columns").cols() == 0);
    REQUIRE(reader.to_mat<double>("no_rows").rows() == 0);
    REQUIRE(reader.to_mat<double>("no_rows").cols() == 2);
    REQUIRE(reader.nnz("zeros") == 0);
    REQUIRE(reader.to_mat<double>("zeros").nonZeros() == 0);
}

TEST_CASE("CscReader rejects dtype mismatch", "[io][csc_reader]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "mismatch.csc";
    write_sample(path);
    const gelex::CscReader reader(path.string());
    REQUIRE_THROWS_WITH(
        reader.to_map<float>("x"),
        Catch::Matchers::ContainsSubstring("dtype mismatch"));
}

TEST_CASE("CscReader rejects malformed files", "[io][csc_reader]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "malformed.csc";
    write_sample(path);
    auto bytes = read_file(path);
    const auto footer = bytes.size() - gelex::detail::footer_size;
    std::uint64_t index_offset{};
    std::memcpy(&index_offset, bytes.data() + footer + 8, sizeof(index_offset));

    SECTION("missing file")
    {
        REQUIRE_THROWS_WITH(
            gelex::CscReader((fixture.get_test_dir() / "none.csc").string()),
            Catch::Matchers::ContainsSubstring("not found"));
    }
    SECTION("file smaller than the footer")
    {
        bytes.resize(10);
        write_file(path, bytes);
        REQUIRE_THROWS_WITH(
            gelex::CscReader(path.string()),
            Catch::Matchers::ContainsSubstring("too small"));
    }
    SECTION("bad magic")
    {
        bytes[footer] = 'X';
        write_file(path, bytes);
        REQUIRE_THROWS_WITH(
            gelex::CscReader(path.string()),
            Catch::Matchers::ContainsSubstring("magic"));
    }
    SECTION("index offset beyond the footer")
    {
        patch<std::uint64_t>(bytes, footer + 8, footer + 64);
        write_file(path, bytes);
        REQUIRE_THROWS_WITH(
            gelex::CscReader(path.string()),
            Catch::Matchers::ContainsSubstring("index offset"));
    }
    SECTION("matrix count larger than the index")
    {
        patch<std::uint64_t>(bytes, footer + 16, 1000);
        write_file(path, bytes);
        REQUIRE_THROWS_WITH(
            gelex::CscReader(path.string()),
            Catch::Matchers::ContainsSubstring("matrix count"));
    }
    SECTION("array offset outside the data region")
    {
        // Entry layout: u32 name size, name "x", u8 type, rows, cols, nnz,
        // then the values offset.
        const auto values_offset = index_offset + 4 + 1 + 1 + 24;
        patch<std::uint64_t>(bytes, values_offset, index_offset);
        write_file(path, bytes);
        REQUIRE_THROWS_WITH(
            gelex::CscReader(path.string()),
            Catch::Matchers::ContainsSubstring("outside the data region"));
    }
    SECTION("overlapping arrays")
    {
        const auto indices_offset = index_offset + 4 + 1 + 1 + 32;
        patch<std::uint64_t>(bytes, indices_offset, std::uint64_t{0});
        write_file(path, bytes);
        REQUIRE_THROWS_WITH(
            gelex::CscReader(path.string()),
            Catch::Matchers::ContainsSubstring("overlap"));
    }
    SECTION("truncated index")
    {
        bytes.erase(
            bytes.begin() + static_cast<std::ptrdiff_t>(index_offset) + 4,
            bytes.begin() + static_cast<std::ptrdiff_t>(index_offset) + 5);
        write_file(path, bytes);
        REQUIRE_THROWS_AS(
            gelex::CscReader(path.string()), gelex::GelexException);
    }
}
