// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <string_view>

#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/dense_writer.h"
#include "gelex/io/detail/binary_wire.h"

#include "file_fixture.h"

namespace
{

namespace test = gelex::test;

template <typename T>
consteval auto expected_binary_type() -> gelex::BinaryType
{
    if constexpr (std::same_as<T, double>)
    {
        return gelex::BinaryType::float64;
    }
    else if constexpr (std::same_as<T, float>)
    {
        return gelex::BinaryType::float32;
    }
    else
    {
        static_assert(std::same_as<T, std::uint8_t>);
        return gelex::BinaryType::uint8;
    }
}

template <typename T>
auto write_matrix(
    gelex::DenseWriter& writer,
    std::string_view identifier,
    const Eigen::MatrixX<T>& matrix) -> void
{
    writer.reserve<T>(
        identifier,
        gelex::BinaryShape{
            static_cast<std::uint64_t>(matrix.rows()),
            static_cast<std::uint64_t>(matrix.cols())})
        << matrix.reshaped();
}

}  // namespace

TEMPLATE_TEST_CASE(
    "BinaryReader reads supported dtypes",
    "[io][binary_reader]",
    double,
    float,
    std::uint8_t)
{
    test::FileFixture fixture;
    const auto container_path = fixture.get_test_dir() / "dtype.samples";
    const Eigen::MatrixX<TestType> expected{
        {TestType{1}, TestType{2}, TestType{3}},
        {TestType{4}, TestType{5}, TestType{6}}};
    {
        auto writer = gelex::open_dense_writer(container_path.string());
        write_matrix(writer, "values", expected);
        writer.close();
    }

    gelex::BinaryReader reader(container_path.string());
    const auto& info = reader.info("values");
    REQUIRE(info.type == expected_binary_type<TestType>());
    REQUIRE(info.shape == (gelex::BinaryShape{2, 3}));

    const auto view = reader.to_map<TestType>("values");
    REQUIRE(
        reinterpret_cast<std::uintptr_t>(view.data())
            % gelex::detail::payload_alignment
        == 0);
    REQUIRE(view.isApprox(expected));

    const auto matrix = reader.to_mat<TestType>("values");
    REQUIRE(matrix.isApprox(expected));
    REQUIRE(matrix.data() != view.data());
}

TEST_CASE("BinaryReader exposes payload metadata", "[io][binary_reader]")
{
    test::FileFixture fixture;
    const auto container_path = fixture.get_test_dir() / "metadata.samples";
    {
        auto writer = gelex::open_dense_writer(container_path.string());
        writer.reserve<double>("zeta", gelex::BinaryShape{1, 1}) << 3.0;
        writer.reserve<double>("alpha", gelex::BinaryShape{1, 1}) << 1.0;
        writer.reserve<double>("beta", gelex::BinaryShape{1, 1}) << 2.0;
        writer.close();
    }

    gelex::BinaryReader reader(container_path.string());
    REQUIRE(reader.size() == 3);
    REQUIRE(reader.contains("alpha"));
    REQUIRE_FALSE(reader.contains("missing"));

    const auto& info = reader.info("beta");
    REQUIRE(info.identifier == "beta");
    REQUIRE(info.type == gelex::BinaryType::float64);
    REQUIRE(info.shape == (gelex::BinaryShape{1, 1}));

    const auto payloads = reader.payloads();
    REQUIRE(payloads.size() == 3);
    REQUIRE(payloads[0].identifier == "alpha");
    REQUIRE(payloads[1].identifier == "beta");
    REQUIRE(payloads[2].identifier == "zeta");
    REQUIRE_THROWS_AS(reader.info("missing"), gelex::GelexException);
}

TEST_CASE("BinaryReader rejects dtype mismatch", "[io][binary_reader]")
{
    test::FileFixture fixture;
    const auto container_path
        = fixture.get_test_dir() / "dtype_mismatch.samples";
    {
        auto writer = gelex::open_dense_writer(container_path.string());
        writer.reserve<double>("value", gelex::BinaryShape{1, 1}) << 1.0;
        writer.close();
    }

    gelex::BinaryReader reader(container_path.string());
    REQUIRE_THROWS_AS(
        reader.to_map<std::uint8_t>("value"), gelex::GelexException);
}

TEST_CASE("BinaryReader rejects malformed footer", "[io][binary_reader]")
{
    test::FileFixture fixture;
    const auto& directory = fixture.get_test_dir();

    SECTION("file is smaller than the footer")
    {
        const auto path = directory / "tiny.samples";
        std::ofstream output(path, std::ios::binary);
        const std::array<char, 10> bytes{};
        output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        output.close();

        REQUIRE_THROWS_AS(
            gelex::BinaryReader(path.string()), gelex::GelexException);
    }

    SECTION("footer magic is invalid")
    {
        const auto path = directory / "invalid_magic.samples";
        std::ofstream output(path, std::ios::binary);
        const std::array<char, 64> bytes{};
        output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        output.close();

        REQUIRE_THROWS_AS(
            gelex::BinaryReader(path.string()), gelex::GelexException);
    }
}
