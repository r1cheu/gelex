// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/mapped_file.h"

#include "file_fixture.h"

#if defined(__linux__)
#include <fcntl.h>
#include <unistd.h>
#endif

namespace
{

struct StoredMatrix
{
    std::string identifier;
    std::uint8_t type{};
    gelex::BinaryShape shape{};
    std::uint64_t nnz{};
    std::array<std::uint64_t, 3> offsets{};
};

template <typename T>
auto read_integer(const gelex::MappedFile& file, std::size_t& cursor) -> T
{
    REQUIRE(cursor <= file.size());
    REQUIRE(sizeof(T) <= file.size() - cursor);
    T value{};
    std::memcpy(&value, file.data() + cursor, sizeof(T));
    cursor += sizeof(T);
    return value;
}

auto read_index(const gelex::MappedFile& file) -> std::vector<StoredMatrix>
{
    REQUIRE(file.size() >= 24);
    const auto footer = file.size() - 24;
    REQUIRE(
        std::string_view{reinterpret_cast<const char*>(file.data() + footer), 8}
        == "GELEXSC1");
    auto cursor = footer + 8;
    const auto index_offset = read_integer<std::uint64_t>(file, cursor);
    const auto count = read_integer<std::uint64_t>(file, cursor);
    REQUIRE(index_offset % 64 == 0);
    REQUIRE(index_offset <= footer);
    cursor = static_cast<std::size_t>(index_offset);
    std::vector<StoredMatrix> matrices;
    for (std::uint64_t index = 0; index < count; ++index)
    {
        StoredMatrix matrix;
        const auto name_size = read_integer<std::uint32_t>(file, cursor);
        REQUIRE(name_size <= file.size() - cursor);
        matrix.identifier = std::string{
            reinterpret_cast<const char*>(file.data() + cursor), name_size};
        cursor += name_size;
        matrix.type = read_integer<std::uint8_t>(file, cursor);
        matrix.shape[0] = read_integer<std::uint64_t>(file, cursor);
        matrix.shape[1] = read_integer<std::uint64_t>(file, cursor);
        matrix.nnz = read_integer<std::uint64_t>(file, cursor);
        for (auto& offset : matrix.offsets)
        {
            offset = read_integer<std::uint64_t>(file, cursor);
            REQUIRE(offset % 64 == 0);
            REQUIRE(offset <= index_offset);
        }
        matrices.push_back(std::move(matrix));
    }
    REQUIRE(cursor == footer);
    return matrices;
}

auto open_mapped(const std::filesystem::path& path) -> gelex::MappedFile
{
    gelex::MappedFile file;
    std::error_code ec;
    file.map(path.string(), ec);
    REQUIRE_FALSE(ec);
    return file;
}

template <typename T>
auto array_view(
    const gelex::MappedFile& file,
    std::uint64_t offset,
    std::uint64_t count) -> Eigen::Map<const Eigen::VectorX<T>>
{
    REQUIRE(offset <= file.size());
    REQUIRE(count <= (file.size() - offset) / sizeof(T));
    return {
        reinterpret_cast<const T*>(file.data() + offset),
        static_cast<Eigen::Index>(count)};
}

template <typename T>
auto to_dense(const gelex::MappedFile& file, const StoredMatrix& matrix)
    -> Eigen::MatrixX<T>
{
    const auto values = array_view<T>(file, matrix.offsets[0], matrix.nnz);
    const auto indices
        = array_view<std::int64_t>(file, matrix.offsets[1], matrix.nnz);
    const auto indptr = array_view<std::int64_t>(
        file, matrix.offsets[2], matrix.shape[1] + 1);
    const Eigen::Map<
        const Eigen::SparseMatrix<T, Eigen::ColMajor, std::int64_t>>
        sparse{
            static_cast<Eigen::Index>(matrix.shape[0]),
            static_cast<Eigen::Index>(matrix.shape[1]),
            static_cast<Eigen::Index>(matrix.nnz),
            indptr.data(),
            indices.data(),
            values.data()};
    return Eigen::MatrixX<T>{sparse};
}

// Array `number` is written to ".N.tmp" and committed to ".N" until it has
// been merged.
auto committed_spool_path(const std::filesystem::path& path, int number)
    -> std::filesystem::path
{
    return path.string() + "." + std::to_string(number);
}

auto spool_path(const std::filesystem::path& path, int number)
    -> std::filesystem::path
{
    return committed_spool_path(path, number).string() + ".tmp";
}

auto check_spools_removed(const std::filesystem::path& path, int count = 3)
    -> void
{
    for (int number = 1; number <= count; ++number)
    {
        CHECK_FALSE(std::filesystem::exists(spool_path(path, number)));
        CHECK_FALSE(
            std::filesystem::exists(committed_spool_path(path, number)));
    }
}

template <typename Writer>
concept CanReserve = requires(Writer&& writer) {
    std::forward<Writer>(writer).template reserve<double>(
        "x", gelex::BinaryShape{1, 1});
};

static_assert(!std::is_copy_constructible_v<gelex::CscWriter>);
static_assert(!std::is_move_constructible_v<gelex::CscWriter>);
static_assert(!std::is_copy_constructible_v<gelex::CscStream<double>>);
static_assert(std::is_move_constructible_v<gelex::CscStream<double>>);
static_assert(std::is_move_assignable_v<gelex::CscStream<double>>);
static_assert(CanReserve<gelex::CscWriter&>);
static_assert(!CanReserve<gelex::CscWriter>);

}  // namespace

TEST_CASE(
    "CscWriter interleaves matrices and preserves double precision",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "markers.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto coefficients = writer.reserve<double>("coefficients", {4, 3});
    auto assignments = writer.reserve<std::uint8_t>("assignments", {4, 3});
    for (int number = 1; number <= 6; ++number)
    {
        CHECK(std::filesystem::is_regular_file(spool_path(path, number)));
    }
    const Eigen::MatrixXd expected{
        {0, 0, 3}, {1.0000000001, 0, 0}, {0, 0, 4}, {-2.5, 0, 0}};
    const Eigen::MatrixX<std::uint8_t> classes{
        {0, 0, 2}, {1, 0, 0}, {0, 0, 3}, {0, 0, 0}};
    for (Eigen::Index column = 0; column < expected.cols(); ++column)
    {
        auto& result = coefficients << expected.col(column);
        REQUIRE(&result == &coefficients);
        assignments << classes.col(column);
    }
    CHECK_FALSE(std::filesystem::exists(path));
    writer.close();
    check_spools_removed(path, 6);

    const auto file = open_mapped(path);
    const auto matrices = read_index(file);
    REQUIRE(matrices.size() == 2);
    CHECK(matrices[0].identifier == "coefficients");
    CHECK(matrices[0].nnz == 4);
    CHECK(matrices[1].nnz == 3);
    CHECK(matrices[0].shape == gelex::BinaryShape{4, 3});
    CHECK(matrices[0].type == std::to_underlying(gelex::BinaryType::float64));
    CHECK(matrices[1].type == std::to_underlying(gelex::BinaryType::uint8));
    CHECK(to_dense<double>(file, matrices[0]).isApprox(expected));
    CHECK(to_dense<std::uint8_t>(file, matrices[1]).isApprox(classes));
    CHECK(
        array_view<std::int64_t>(file, matrices[0].offsets[1], 4)
            .isApprox(Eigen::VectorX<std::int64_t>{{1, 3, 0, 2}}));
    CHECK(
        array_view<std::int64_t>(file, matrices[0].offsets[2], 4)
            .isApprox(Eigen::VectorX<std::int64_t>{{0, 2, 2, 4}}));
}

TEST_CASE(
    "CscStream writes consecutive columns with different nonzero counts",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "large.csc";
    auto writer = gelex::open_csc_writer(path.string());
    constexpr Eigen::Index rows = 20000;
    auto stream = writer.reserve<double>("x", {rows, 2});
    const Eigen::VectorXd first
        = Eigen::VectorXd::LinSpaced(rows, 1.0, 20000.0);
    Eigen::VectorXd second = Eigen::VectorXd::Zero(rows);
    for (Eigen::Index row = 0; row < rows; row += 100)
    {
        second(row) = -0.125;
    }
    stream << first << second;
    CHECK_FALSE(std::filesystem::exists(path));
    writer.close();
    const auto file = open_mapped(path);
    const auto matrices = read_index(file);
    REQUIRE(matrices.size() == 1);
    CHECK(matrices[0].nnz == 20200);
    const auto dense = to_dense<double>(file, matrices[0]);
    CHECK(dense.col(0).isApprox(first));
    CHECK(dense.col(1).isApprox(second));
}

TEST_CASE("CscStream supports float and byte matrices", "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "types.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto floats = writer.reserve<float>("float", {3, 1});
    auto bytes = writer.reserve<std::uint8_t>("byte", {3, 1});
    const Eigen::VectorXf x{{0.0F, -0.25F, 2.5F}};
    const Eigen::VectorX<std::uint8_t> y{{3, 0, 7}};
    floats << x;
    bytes << y;
    writer.close();
    const auto file = open_mapped(path);
    const auto matrices = read_index(file);
    REQUIRE(matrices.size() == 2);
    CHECK(matrices[0].type == std::to_underlying(gelex::BinaryType::float32));
    CHECK(matrices[1].type == std::to_underlying(gelex::BinaryType::uint8));
    CHECK(to_dense<float>(file, matrices[0]).isApprox(x));
    CHECK(to_dense<std::uint8_t>(file, matrices[1]) == y);
}

TEST_CASE("CscStream drops only exact zeros", "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "special.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto stream = writer.reserve<double>("x", {6, 1});
    const auto tiny = std::numeric_limits<double>::denorm_min();
    const Eigen::VectorXd x{
        {0.0,
         -0.0,
         tiny,
         -tiny,
         std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::quiet_NaN()}};
    stream << x;
    writer.close();
    const auto file = open_mapped(path);
    const auto matrices = read_index(file);
    REQUIRE(matrices.size() == 1);
    CHECK(matrices[0].nnz == 4);
    const auto values = array_view<double>(file, matrices[0].offsets[0], 4);
    CHECK(values(0) == tiny);
    CHECK(values(1) == -tiny);
    CHECK(std::isinf(values(2)));
    CHECK(std::isnan(values(3)));
}

TEST_CASE(
    "CscWriter represents empty containers and zero dimensions",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "empty.csc";
    auto writer = gelex::open_csc_writer(path.string());
    SECTION("no matrices")
    {
        writer.close();
        const auto file = open_mapped(path);
        CHECK(file.size() == 24);
        CHECK(read_index(file).empty());
    }
    SECTION("empty dimensions and all-zero columns")
    {
        auto zero_rows = writer.reserve<double>("rows", {0, 3});
        [[maybe_unused]] auto zero_columns
            = writer.reserve<double>("cols", {3, 0});
        [[maybe_unused]] auto zero_both
            = writer.reserve<double>("both", {0, 0});
        auto zero_values = writer.reserve<double>("values", {4, 2});
        const Eigen::VectorXd empty;
        const Eigen::VectorXd zeros = Eigen::VectorXd::Zero(4);
        zero_rows << empty << empty << empty;
        zero_values << zeros << zeros;
        writer.close();
        const auto file = open_mapped(path);
        const auto matrices = read_index(file);
        REQUIRE(matrices.size() == 4);
        for (const auto& matrix : matrices)
        {
            CHECK(matrix.nnz == 0);
            CHECK(
                array_view<std::int64_t>(
                    file, matrix.offsets[2], matrix.shape[1] + 1)
                    .isZero());
            const auto dense = to_dense<double>(file, matrix);
            CHECK(dense.rows() == matrix.shape[0]);
            CHECK(dense.cols() == matrix.shape[1]);
            CHECK(dense.isZero());
        }
    }
}

TEST_CASE(
    "CscWriter rejects invalid identifiers without poisoning the writer",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    auto writer = gelex::open_csc_writer(
        (fixture.get_test_dir() / "reserve.csc").string());
    [[maybe_unused]] auto stream = writer.reserve<double>("x", {1, 0});
    REQUIRE_THROWS_AS(
        writer.reserve<double>("", {1, 0}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        writer.reserve<float>("x", {2, 0}), gelex::GelexException);
    writer.close();
}

TEST_CASE(
    "CscStream validates rows and refuses extra columns before writing",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "bounds.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto stream = writer.reserve<double>("x", {2, 1});
    const Eigen::VectorXd wrong{{1, 2, 3}};
    const Eigen::VectorXd correct{{1, 0}};
    REQUIRE_THROWS_AS(stream << wrong, gelex::GelexException);
    stream << correct;
    REQUIRE_THROWS_AS(stream << correct, gelex::GelexException);
    writer.close();
    REQUIRE_NOTHROW(writer.close());
    REQUIRE_THROWS_AS(stream << correct, gelex::GelexException);
    REQUIRE_THROWS_AS(
        writer.reserve<double>("y", {2, 0}), gelex::GelexException);
    const auto file = open_mapped(path);
    CHECK(to_dense<double>(file, read_index(file).front()).isApprox(correct));
}

TEST_CASE(
    "CscStream handles survive new reservations and transfer on move",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "move.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto original = writer.reserve<double>("original", {2, 2});
    const Eigen::VectorXd x{{0, 1}};
    original << x;
    for (int index = 0; index < 40; ++index)
    {
        [[maybe_unused]] auto other
            = writer.reserve<float>(std::to_string(index), {1, 0});
    }
    auto moved = std::move(original);
    REQUIRE_THROWS_AS(original << x, gelex::GelexException);
    auto destination = writer.reserve<double>("unused", {0, 0});
    destination = std::move(moved);
    REQUIRE_THROWS_AS(moved << x, gelex::GelexException);
    destination << x;
    writer.close();
    const auto file = open_mapped(path);
    const auto matrices = read_index(file);
    REQUIRE_FALSE(matrices.empty());
    CHECK(matrices.front().identifier == "original");
    CHECK(
        to_dense<double>(file, matrices.front())
            .isApprox(Eigen::MatrixXd{{0, 0}, {1, 1}}));
}

TEST_CASE(
    "CscWriter publishes matrices after their completed streams are destroyed",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "completed.csc";
    auto writer = gelex::open_csc_writer(path.string());
    const Eigen::MatrixXd expected{{0, 2}, {1, 0}, {0, -3}};
    {
        auto stream = writer.reserve<double>("x", {3, 2});
        stream << expected.col(0) << expected.col(1);
    }
    writer.close();
    const auto file = open_mapped(path);
    const auto matrices = read_index(file);
    REQUIRE(matrices.size() == 1);
    CHECK(matrices[0].nnz == 3);
    CHECK(to_dense<double>(file, matrices[0]).isApprox(expected));
}

TEST_CASE(
    "CscWriter refuses to publish until every matrix is complete",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path
        = fixture.create_named_text_file("incomplete.csc", "original");
    auto writer = gelex::open_csc_writer(path.string());
    auto full = writer.reserve<double>("full", {1, 1});
    auto partial = writer.reserve<double>("partial", {1, 2});
    const Eigen::VectorXd x{{1}};
    full << x;
    partial << x;
    REQUIRE_THROWS_AS(writer.close(), gelex::GelexException);
    {
        const auto file = open_mapped(path);
        CHECK(
            std::string_view{
                reinterpret_cast<const char*>(file.data()), file.size()}
            == "original");
    }
    partial << x;
    writer.close();
    check_spools_removed(path, 6);
    const auto file = open_mapped(path);
    const auto entries = read_index(file);
    REQUIRE(entries.size() == 2);
    CHECK(to_dense<double>(file, entries[1]).isApprox(Eigen::MatrixXd{{1, 1}}));
}

TEST_CASE(
    "CscWriter destruction discards even complete uncommitted matrices",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "discard.csc";
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
        auto writer = gelex::open_csc_writer(path.string());
        auto stream = writer.reserve<double>("x", {1, 1});
        stream << Eigen::VectorXd{{2}};
        CHECK_FALSE(std::filesystem::exists(path));
    }
    gelex::set_sink({});
    CHECK_FALSE(std::filesystem::exists(path));
    CHECK_FALSE(std::filesystem::exists(path.string() + ".tmp"));
    check_spools_removed(path);
    // Forgetting close() on a normal path is reported.
    REQUIRE(errors.size() == 1);
    CHECK(errors[0].contains("unclosed"));
}

TEST_CASE(
    "CscWriter failed publication discards temporary output",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "publish.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto stream = writer.reserve<double>("x", {1, 1});
    stream << Eigen::VectorXd{{2}};
    REQUIRE(std::filesystem::create_directory(path));
    REQUIRE_THROWS_AS(writer.close(), gelex::GelexException);
    CHECK(std::filesystem::is_directory(path));
    CHECK_FALSE(std::filesystem::exists(path.string() + ".tmp"));
    check_spools_removed(path);
}

TEST_CASE(
    "CscWriter failed reservation discards its arrays and keeps the writer "
    "usable",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "failure.csc";
    const auto existing
        = fixture.create_named_text_file("failure.csc.5.tmp", "original");
    {
        auto writer = gelex::open_csc_writer(path.string());
        auto stream = writer.reserve<double>("x", {1, 1});
        REQUIRE_THROWS_AS(
            writer.reserve<double>("other", {1, 1}), gelex::GelexException);
        CHECK_FALSE(std::filesystem::exists(spool_path(path, 4)));
        CHECK_FALSE(std::filesystem::exists(spool_path(path, 6)));
        stream << Eigen::VectorXd{{1}};
        writer.close();
    }
    CHECK_FALSE(std::filesystem::exists(path.string() + ".tmp"));
    check_spools_removed(path, 4);
    const auto file = open_mapped(path);
    const auto entries = read_index(file);
    REQUIRE(entries.size() == 1);
    CHECK(entries.front().identifier == "x");
    const auto preserved = open_mapped(existing);
    CHECK(
        std::string_view{
            reinterpret_cast<const char*>(preserved.data()), preserved.size()}
        == "original");
}

TEST_CASE(
    "CscWriter refuses to reserve over an existing committed array path",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "sibling.csc";
    const auto existing
        = fixture.create_named_text_file("sibling.csc.2", "original");
    {
        auto writer = gelex::open_csc_writer(path.string());
        REQUIRE_THROWS_AS(
            writer.reserve<double>("x", {1, 1}), gelex::GelexException);
        CHECK_FALSE(std::filesystem::exists(spool_path(path, 1)));
    }
    CHECK_FALSE(std::filesystem::exists(path.string() + ".tmp"));
    const auto file = open_mapped(existing);
    CHECK(
        std::string_view{
            reinterpret_cast<const char*>(file.data()), file.size()}
        == "original");
}

TEST_CASE(
    "Opening a CscWriter preserves existing temporary resources",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    REQUIRE_THROWS_AS(gelex::open_csc_writer(""), gelex::GelexException);
    REQUIRE_THROWS_AS(
        gelex::open_csc_writer(
            (fixture.get_test_dir() / "missing" / "x").string()),
        gelex::GelexException);
    const auto directory = fixture.get_test_dir() / "directory";
    REQUIRE(std::filesystem::create_directory(directory));
    REQUIRE_THROWS_AS(
        gelex::open_csc_writer(directory.string()), gelex::GelexException);
    CHECK(std::filesystem::is_directory(directory));
    check_spools_removed(directory);
    const auto path = fixture.get_test_dir() / "exclusive.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto stream = writer.reserve<double>("x", {1, 1});
    REQUIRE_THROWS_AS(
        gelex::open_csc_writer(path.string()), gelex::GelexException);
    stream << Eigen::VectorXd{{3}};
    writer.close();
    const auto file = open_mapped(path);
    CHECK(
        to_dense<double>(file, read_index(file).front())
            .isApprox(Eigen::VectorXd{{3}}));
}

TEST_CASE(
    "CscStream move assignment leaves an abandoned matrix incomplete",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "abandoned.csc";
    {
        auto writer = gelex::open_csc_writer(path.string());
        auto destination = writer.reserve<double>("abandoned", {1, 2});
        destination << Eigen::VectorXd{{1}};
        auto source = writer.reserve<double>("kept", {1, 1});
        destination = std::move(source);
        destination << Eigen::VectorXd{{2}};
        REQUIRE_THROWS_AS(
            source << Eigen::VectorXd{{3}}, gelex::GelexException);
        REQUIRE_THROWS_AS(writer.close(), gelex::GelexException);
    }
    CHECK_FALSE(std::filesystem::exists(path));
    check_spools_removed(path, 6);
}

#if defined(__linux__)
TEST_CASE(
    "CscStream write and close failures prevent publication",
    "[io][csc_writer]")
{
    if (!std::filesystem::exists("/dev/full")
        || !std::filesystem::exists("/proc/self/fd"))
    {
        SKIP("requires /dev/full and /proc/self/fd");
    }
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "stream_failure.csc";
    Eigen::Index rows = 1;
    int file_number = 1;
    SECTION("values write fails")
    {
        rows = 20000;
    }
    SECTION("merged output write fails")
    {
        rows = 20000;
        file_number = 0;
    }
    SECTION("values close fails") {}
    SECTION("indices close fails")
    {
        file_number = 2;
    }
    SECTION("indptr close fails")
    {
        file_number = 3;
    }
    {
        auto writer = gelex::open_csc_writer(path.string());
        auto stream = writer.reserve<double>(
            "x", {static_cast<std::uint64_t>(rows), 1});
        auto other = writer.reserve<double>("other", {1, 1});
        const auto target_path
            = file_number == 0 ? std::filesystem::path{path.string() + ".tmp"}
                               : spool_path(path, file_number);
        int target = -1;
        for (const auto& entry :
             std::filesystem::directory_iterator("/proc/self/fd"))
        {
            std::error_code ec;
            if (std::filesystem::read_symlink(entry.path(), ec) == target_path
                && !ec)
            {
                target = std::stoi(entry.path().filename().string());
                break;
            }
        }
        REQUIRE(target >= 0);
        // Redirect only this stream's descriptor to force a real I/O failure.
        const auto full = ::open("/dev/full", O_WRONLY);
        REQUIRE(full >= 0);
        const auto redirected = ::dup2(full, target);
        ::close(full);
        REQUIRE(redirected == target);
        const Eigen::VectorXd column = Eigen::VectorXd::Ones(rows);
        REQUIRE_THROWS_AS(stream << column, gelex::GelexException);
        REQUIRE_THROWS_AS(writer.close(), gelex::GelexException);
    }
    CHECK_FALSE(std::filesystem::exists(path));
    CHECK_FALSE(std::filesystem::exists(path.string() + ".tmp"));
    check_spools_removed(path, 6);
}
#endif

TEST_CASE(
    "CscStream merges on completion while the index keeps reservation order",
    "[io][csc_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "completion_order.csc";
    auto writer = gelex::open_csc_writer(path.string());
    auto first = writer.reserve<double>("first", {2, 2});
    first << Eigen::VectorXd{{1, 0}};
    {
        auto second = writer.reserve<float>("second", {2, 1});
        second << Eigen::VectorXf{{0, 3}};
        for (int number = 4; number <= 6; ++number)
        {
            CHECK_FALSE(std::filesystem::exists(spool_path(path, number)));
        }
        CHECK(std::filesystem::exists(spool_path(path, 1)));
        CHECK_FALSE(std::filesystem::exists(path));
    }
    first << Eigen::VectorXd{{0, 2}};
    check_spools_removed(path, 6);
    CHECK_FALSE(std::filesystem::exists(path));
    writer.close();
    const auto file = open_mapped(path);
    const auto matrices = read_index(file);
    REQUIRE(matrices.size() == 2);
    CHECK(matrices[0].identifier == "first");
    CHECK(matrices[1].identifier == "second");
    CHECK(matrices[1].offsets[0] < matrices[0].offsets[0]);
    CHECK(
        to_dense<double>(file, matrices[0])
            .isApprox(Eigen::MatrixXd{{1, 0}, {0, 2}}));
    CHECK(to_dense<float>(file, matrices[1]).isApprox(Eigen::VectorXf{{0, 3}}));
}
