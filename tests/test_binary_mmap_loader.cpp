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
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <string>
#include <type_traits>
#include <vector>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using gelex::test::FileFixture;

#if __has_include(                      \
    "../src/data/binary_mmap_loader.h") \
    && __has_include("../src/data/binary_writer.h")

#include "../src/data/binary_mmap_loader.h"
#include "../src/data/binary_writer.h"

namespace gelex::detail::test
{

constexpr size_t kMetaSize = BinaryWriter<double>::kMetaSize;

auto read_all_bytes_for_mmap_loader(const fs::path& path)
    -> std::vector<std::byte>
{
    std::ifstream file(path, std::ios::binary);
    REQUIRE(file.is_open());

    file.seekg(0, std::ios::end);
    const auto size = static_cast<size_t>(file.tellg());
    file.seekg(0, std::ios::beg);

    std::vector<std::byte> bytes(size);
    if (size > 0)
    {
        file.read(
            reinterpret_cast<char*>(bytes.data()),
            static_cast<std::streamsize>(size));
        REQUIRE(file.good());
    }

    return bytes;
}

template <typename eT>
auto make_vector(std::initializer_list<eT> values)
    -> Eigen::Matrix<eT, Eigen::Dynamic, 1>
{
    Eigen::Matrix<eT, Eigen::Dynamic, 1> vec(static_cast<int>(values.size()));
    int index = 0;
    for (const auto value : values)
    {
        vec(index++) = value;
    }
    return vec;
}

template <typename eT>
auto make_expected_matrix(
    const std::vector<Eigen::Matrix<eT, Eigen::Dynamic, 1>>& columns)
    -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>
{
    if (columns.empty())
    {
        return Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>(0, 0);
    }

    const auto rows = columns.front().size();
    const auto cols = static_cast<int>(columns.size());
    Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
    for (int col = 0; col < cols; ++col)
    {
        mat.col(col) = columns[static_cast<size_t>(col)];
    }

    return mat;
}

template <typename eT>
auto require_matrix_equal(
    const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>& actual,
    const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>& expected) -> void
{
    REQUIRE(actual.rows() == expected.rows());
    REQUIRE(actual.cols() == expected.cols());

    if constexpr (std::is_floating_point_v<eT>)
    {
        REQUIRE(actual.isApprox(expected));
    }
    else
    {
        for (Eigen::Index col = 0; col < actual.cols(); ++col)
        {
            for (Eigen::Index row = 0; row < actual.rows(); ++row)
            {
                REQUIRE(actual(row, col) == expected(row, col));
            }
        }
    }
}

}  // namespace gelex::detail::test

template <typename eT>
using TestVector = Eigen::Matrix<eT, Eigen::Dynamic, 1>;

using namespace gelex::detail;        // NOLINT
using namespace gelex::detail::test;  // NOLINT

#define GELEX_BINARY_MMAP_LOADER_TYPES uint8_t, float, double

TEMPLATE_TEST_CASE(
    "BinaryMmapLoader - round trip from BinaryWriter",
    "[data][binary_mmap_loader]",
    GELEX_BINARY_MMAP_LOADER_TYPES)
{
    FileFixture files;
    const auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    std::vector<TestVector<TestType>> columns;
    columns.push_back(
        make_vector<TestType>(
            {static_cast<TestType>(1),
             static_cast<TestType>(2),
             static_cast<TestType>(3)}));
    columns.push_back(
        make_vector<TestType>(
            {static_cast<TestType>(10),
             static_cast<TestType>(20),
             static_cast<TestType>(30)}));
    columns.push_back(
        make_vector<TestType>(
            {static_cast<TestType>(5),
             static_cast<TestType>(6),
             static_cast<TestType>(7)}));

    {
        BinaryWriter<TestType> writer(file_path_str);
        for (const auto& column : columns)
        {
            writer.write(column);
        }
        writer.finish();
    }

    BinaryMmapLoader<TestType> loader(file_path_str);
    const auto& mapped = loader.matrix();
    const auto expected = make_expected_matrix<TestType>(columns);

    REQUIRE(mapped.rows() == expected.rows());
    REQUIRE(mapped.cols() == expected.cols());
    require_matrix_equal<TestType>(mapped, expected);

    const auto copied = loader.load_copy();
    require_matrix_equal<TestType>(copied, expected);
}

TEMPLATE_TEST_CASE(
    "BinaryMmapLoader - empty matrix for header-only file",
    "[data][binary_mmap_loader]",
    GELEX_BINARY_MMAP_LOADER_TYPES)
{
    FileFixture files;
    const auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    {
        BinaryWriter<TestType> writer(file_path_str);
        writer.finish();
    }

    BinaryMmapLoader<TestType> loader(file_path_str);
    const auto& mapped = loader.matrix();
    REQUIRE(mapped.rows() == 0);
    REQUIRE(mapped.cols() == 0);

    const auto copied = loader.load_copy();
    REQUIRE(copied.rows() == 0);
    REQUIRE(copied.cols() == 0);
}

TEST_CASE(
    "BinaryMmapLoader - dtype mismatch throws",
    "[data][binary_mmap_loader]")
{
    FileFixture files;
    const auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    {
        BinaryWriter<double> writer(file_path_str);
        const auto record = make_vector<double>({1.0, 2.0, 3.0});
        writer.write(record);
        writer.finish();
    }

    REQUIRE_THROWS_AS(
        BinaryMmapLoader<float>(file_path_str),
        gelex::ArgumentValidationException);
}

TEST_CASE(
    "BinaryMmapLoader - invalid magic throws",
    "[data][binary_mmap_loader]")
{
    FileFixture files;
    const auto source_path = files.generate_random_file_path(".bin");
    const std::string source_path_str = source_path.string();

    {
        BinaryWriter<double> writer(source_path_str);
        const auto record = make_vector<double>({1.0, 2.0});
        writer.write(record);
        writer.finish();
    }

    auto bytes = read_all_bytes_for_mmap_loader(source_path);
    REQUIRE(bytes.size() >= kMetaSize);
    bytes[0] = std::byte{'X'};

    const auto bad_path
        = files.create_named_binary_file("bad_magic.bin", bytes);

    REQUIRE_THROWS_AS(
        BinaryMmapLoader<double>(bad_path.string()),
        gelex::FileFormatException);
}

TEST_CASE(
    "BinaryMmapLoader - truncated header throws",
    "[data][binary_mmap_loader]")
{
    FileFixture files;
    std::vector<std::byte> truncated(7, std::byte{0});
    const auto bad_path
        = files.create_named_binary_file("truncated_header.bin", truncated);

    REQUIRE_THROWS_AS(
        BinaryMmapLoader<double>(bad_path.string()),
        gelex::FileFormatException);
}

TEST_CASE(
    "BinaryMmapLoader - payload size mismatch throws",
    "[data][binary_mmap_loader]")
{
    FileFixture files;
    const auto source_path = files.generate_random_file_path(".bin");
    const std::string source_path_str = source_path.string();

    {
        BinaryWriter<float> writer(source_path_str);
        const auto record = make_vector<float>({1.0F, 2.0F, 3.0F, 4.0F});
        writer.write(record);
        writer.finish();
    }

    auto bytes = read_all_bytes_for_mmap_loader(source_path);
    REQUIRE(bytes.size() > kMetaSize);
    bytes.pop_back();
    const auto bad_path
        = files.create_named_binary_file("payload_mismatch.bin", bytes);

    REQUIRE_THROWS_AS(
        BinaryMmapLoader<float>(bad_path.string()), gelex::FileFormatException);
}

TEST_CASE(
    "BinaryMmapLoader - load_copy survives loader lifetime",
    "[data][binary_mmap_loader]")
{
    FileFixture files;
    const auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    {
        BinaryWriter<double> writer(file_path_str);
        writer.write(make_vector<double>({1.0, 2.0, 3.0}));
        writer.write(make_vector<double>({4.0, 5.0, 6.0}));
        writer.finish();
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> copied;
    {
        BinaryMmapLoader<double> loader(file_path_str);
        copied = loader.load_copy();
    }

    REQUIRE(copied.rows() == 3);
    REQUIRE(copied.cols() == 2);
    REQUIRE(copied(0, 0) == 1.0);
    REQUIRE(copied(1, 0) == 2.0);
    REQUIRE(copied(2, 0) == 3.0);
    REQUIRE(copied(0, 1) == 4.0);
    REQUIRE(copied(1, 1) == 5.0);
    REQUIRE(copied(2, 1) == 6.0);
}

#else

TEST_CASE(
    "BinaryMmapLoader - tests are gated until implementation header exists",
    "[data][binary_mmap_loader]")
{
    SUCCEED(
        "binary_mmap_loader.h is not present yet; full tests are compiled once "
        "it exists");
}

#endif
