// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
#include "gelex/io/dense_writer.h"

#include "file_fixture.h"

namespace
{

namespace fs = std::filesystem;
namespace test = gelex::test;

template <typename Writer>
concept CanReserve = requires(Writer&& writer) {
    std::forward<Writer>(writer).template reserve<double>(
        "value", gelex::BinaryShape{1, 1});
};

auto temporary_path(const fs::path& path) -> fs::path
{
    return path.string() + ".tmp";
}

template <typename Derived>
auto as_span(const Eigen::DenseBase<Derived>& values)
    -> std::span<const typename Derived::Scalar>
{
    static_assert(
        Derived::IsVectorAtCompileTime
        || Derived::InnerStrideAtCompileTime == 1);
    return {values.derived().data(), static_cast<std::size_t>(values.size())};
}

}  // namespace

static_assert(CanReserve<gelex::DenseWriter&>);
static_assert(!CanReserve<gelex::DenseWriter>);
static_assert(!std::is_copy_constructible_v<gelex::DenseWriter>);
static_assert(!std::is_move_constructible_v<gelex::DenseWriter>);
static_assert(!std::is_copy_constructible_v<gelex::DenseStream<double>>);
static_assert(std::is_nothrow_move_constructible_v<gelex::DenseStream<double>>);
static_assert(std::is_nothrow_move_assignable_v<gelex::DenseStream<double>>);

TEST_CASE(
    "DenseWriter supports interleaved matrix appends",
    "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "interleaved.samples";
    const Eigen::MatrixXd expected_double{
        {1.0, 4.0, 7.0}, {2.0, 5.0, 8.0}, {3.0, 6.0, 9.0}};
    const Eigen::Matrix<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>
        expected_uint8{{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
    {
        auto writer = gelex::open_dense_writer(path.string());
        auto doubles
            = writer.reserve<double>("double", gelex::BinaryShape{3, 3});
        auto bytes
            = writer.reserve<std::uint8_t>("uint8", gelex::BinaryShape{3, 3});
        for (Eigen::Index column = 0; column < expected_double.cols(); ++column)
        {
            doubles << as_span(expected_double.col(column));
            bytes << as_span(expected_uint8.col(column));
        }
        writer.close();
    }
    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<double>("double").isApprox(expected_double));
    REQUIRE(reader.to_map<std::uint8_t>("uint8").isApprox(expected_uint8));
}

TEST_CASE("DenseWriter writes whole matrices at once", "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "whole.samples";
    const Eigen::MatrixXf expected{{1.0F, 3.0F}, {2.0F, 4.0F}};
    {
        auto writer = gelex::open_dense_writer(path.string());
        auto stream = writer.reserve<float>("x", gelex::BinaryShape{2, 2});
        stream << as_span(expected);
        REQUIRE_THROWS_AS(stream << as_span(expected), gelex::GelexException);
        writer.close();
    }
    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("x").isApprox(expected));
}

TEST_CASE("DenseWriter validates reservations", "[io][dense_writer]")
{
    test::FileFixture fixture;
    auto writer = gelex::open_dense_writer(
        (fixture.get_test_dir() / "reserve.samples").string());
    [[maybe_unused]] const auto reserved
        = writer.reserve<double>("value", gelex::BinaryShape{1, 1});
    REQUIRE_THROWS_AS(
        (writer.reserve<double>("value", gelex::BinaryShape{1, 1})),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (writer.reserve<double>("", gelex::BinaryShape{1, 1})),
        gelex::GelexException);
    REQUIRE(writer.is_open());
}

TEST_CASE("DenseStream validates column and matrix sizes", "[io][dense_writer]")
{
    test::FileFixture fixture;
    auto writer = gelex::open_dense_writer(
        (fixture.get_test_dir() / "sizes.samples").string());
    auto stream = writer.reserve<double>("value", gelex::BinaryShape{2, 2});
    const std::array<double, 3> three{1.0, 2.0, 3.0};
    const std::array<double, 2> two{1.0, 2.0};
    REQUIRE_THROWS_AS(stream << three, gelex::GelexException);
    const std::array<double, 4> whole{1.0, 2.0, 3.0, 4.0};
    REQUIRE_THROWS_AS(stream << 1.0, gelex::GelexException);
    stream << two;
    REQUIRE_THROWS_AS(stream << whole, gelex::GelexException);
    stream << two;
    REQUIRE_THROWS_AS(stream << two, gelex::GelexException);
    REQUIRE(writer.is_open());
}

TEST_CASE("DenseStream transfers ownership on move", "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "move.samples";
    {
        gelex::DenseWriter writer(path.string());
        auto stream = writer.reserve<double>(
            "fixed/coefficients", gelex::BinaryShape{1, 1});
        for (std::size_t index = 0; index < 64; ++index)
        {
            [[maybe_unused]] auto other = writer.reserve<double>(
                "other/" + std::to_string(index), gelex::BinaryShape{1, 0});
        }
        REQUIRE(stream.identifier() == "fixed/coefficients");
        auto moved = std::move(stream);
        REQUIRE_THROWS_AS(stream << 1.0, gelex::GelexException);
        REQUIRE_THROWS_AS(stream.identifier(), gelex::GelexException);
        REQUIRE(moved.identifier() == "fixed/coefficients");
        moved << 1.0;
        writer.close();
    }
    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<double>("fixed/coefficients")
                .isApprox(Eigen::MatrixXd{{1.0}}));
}

TEST_CASE(
    "DenseWriter refuses to publish incomplete matrices",
    "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "incomplete.samples";
    {
        auto writer = gelex::open_dense_writer(path.string());
        auto full = writer.reserve<double>("full", gelex::BinaryShape{1, 1});
        auto partial
            = writer.reserve<double>("partial", gelex::BinaryShape{2, 3});
        full << 3.0;
        partial << as_span(Eigen::VectorXd{{1.0, 2.0}});
        REQUIRE_THROWS_WITH(
            writer.close(),
            Catch::Matchers::ContainsSubstring("\"partial\" is incomplete"));
    }
    REQUIRE_FALSE(fs::exists(path));
    REQUIRE_FALSE(fs::exists(temporary_path(path)));
}

TEST_CASE("DenseWriter commits an empty container", "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "empty.samples";
    {
        auto writer = gelex::open_dense_writer(path.string());
        writer.close();
    }
    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.size() == 0);
}

TEST_CASE(
    "DenseWriter closes explicitly and idempotently",
    "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "explicit_close.samples";
    auto writer = gelex::open_dense_writer(path.string());
    auto stream = writer.reserve<double>("value", gelex::BinaryShape{1, 1});
    stream << 1.0;
    REQUIRE(writer.is_open());
    writer.close();
    REQUIRE_FALSE(writer.is_open());
    REQUIRE_THROWS_WITH(
        stream << 2.0, Catch::Matchers::ContainsSubstring("is closed"));
    REQUIRE_THROWS_WITH(
        writer.reserve<double>("other", gelex::BinaryShape{1, 1}),
        Catch::Matchers::ContainsSubstring("is closed"));
    REQUIRE_NOTHROW(writer.close());
    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<double>("value").isApprox(Eigen::MatrixXd{{1.0}}));
}

TEST_CASE(
    "DenseWriter destruction discards uncommitted output",
    "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "discarded.samples";
    std::vector<std::string> errors;
    gelex::set_sink(
        [&errors](gelex::Level level, std::string_view message)
        {
            if (level == gelex::Level::Error)
            {
                errors.emplace_back(message);
            }
        });
    {
        auto writer = gelex::open_dense_writer(path.string());
        writer.reserve<double>("value", gelex::BinaryShape{1, 1}) << 1.0;
        REQUIRE(fs::exists(temporary_path(path)));
    }
    REQUIRE_FALSE(fs::exists(path));
    REQUIRE_FALSE(fs::exists(temporary_path(path)));
    // Forgetting close() on a normal path is reported.
    REQUIRE(errors.size() == 1);
    CHECK_THAT(errors[0], Catch::Matchers::ContainsSubstring("unclosed"));

    REQUIRE_THROWS_AS(
        [&]
        {
            gelex::DenseWriter writer(path.string());
            writer.reserve<double>("value", gelex::BinaryShape{1, 1}) << 1.0;
            throw gelex::GelexException("stop writing");
        }(),
        gelex::GelexException);
    gelex::set_sink({});
    REQUIRE_FALSE(fs::exists(path));
    REQUIRE_FALSE(fs::exists(temporary_path(path)));
    // Unwinding discards silently.
    REQUIRE(errors.size() == 1);
}

TEST_CASE(
    "DenseWriter failed publication discards output and cannot be retried",
    "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "failed_close.samples";
    auto writer = gelex::open_dense_writer(path.string());
    auto stream = writer.reserve<double>("value", gelex::BinaryShape{1, 1});
    stream << 1.0;
    REQUIRE(fs::create_directory(path));
    REQUIRE_THROWS_AS(writer.close(), gelex::GelexException);
    REQUIRE(fs::is_directory(path));
    REQUIRE_FALSE(fs::exists(temporary_path(path)));
    REQUIRE_THROWS_AS(writer.close(), gelex::GelexException);
    REQUIRE_THROWS_AS(stream << 1.0, gelex::GelexException);
}

TEST_CASE(
    "DenseWriter rejects a concurrent writer to the same path",
    "[io][dense_writer]")
{
    test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "exclusive.samples";
    auto writer = gelex::open_dense_writer(path.string());
    auto stream = writer.reserve<double>("value", gelex::BinaryShape{1, 2});
    stream << 1.0;
    REQUIRE(fs::exists(temporary_path(path)));
    REQUIRE_THROWS_AS(
        (gelex::DenseWriter{path.string()}), gelex::GelexException);
    REQUIRE(fs::exists(temporary_path(path)));
    stream << 2.0;
    writer.close();
    const gelex::BinaryReader reader(path.string());
    REQUIRE(
        reader.to_map<double>("value").isApprox(Eigen::MatrixXd{{1.0, 2.0}}));
}
