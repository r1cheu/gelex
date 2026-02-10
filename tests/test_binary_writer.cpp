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
#include <cstring>
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

#if __has_include("../src/data/binary_writer.h")

#include "../src/data/binary_writer.h"

namespace gelex::detail::test
{

constexpr std::array<std::byte, 8> kExpectedMagic
    = {std::byte{'G'},
       std::byte{'E'},
       std::byte{'L'},
       std::byte{'E'},
       std::byte{'X'},
       std::byte{'B'},
       std::byte{'W'},
       std::byte{'1'}};
constexpr uint32_t kExpectedVersion = 1;
constexpr size_t kMetaSize = 32;

struct MetaView
{
    std::array<std::byte, 8> magic{};
    uint32_t version = 0;
    uint64_t n_rows = 0;
    uint64_t n_cols = 0;
    uint8_t dtype = 0;
};

auto read_all_bytes(const fs::path& path) -> std::vector<std::byte>
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

auto decode_le_u32(const std::byte* data) -> uint32_t
{
    return static_cast<uint32_t>(static_cast<uint8_t>(data[0]))
           | (static_cast<uint32_t>(static_cast<uint8_t>(data[1])) << 8)
           | (static_cast<uint32_t>(static_cast<uint8_t>(data[2])) << 16)
           | (static_cast<uint32_t>(static_cast<uint8_t>(data[3])) << 24);
}

auto decode_le_u64(const std::byte* data) -> uint64_t
{
    uint64_t value = 0;
    for (size_t i = 0; i < 8; ++i)
    {
        value |= static_cast<uint64_t>(static_cast<uint8_t>(data[i]))
                 << (i * 8);
    }
    return value;
}

auto parse_header_meta(const std::vector<std::byte>& bytes) -> MetaView
{
    REQUIRE(bytes.size() >= kMetaSize);

    const std::byte* raw = bytes.data();

    MetaView meta;
    for (size_t i = 0; i < meta.magic.size(); ++i)
    {
        meta.magic[i] = raw[i];
    }
    meta.version = decode_le_u32(raw + 8);
    meta.n_rows = decode_le_u64(raw + 12);
    meta.n_cols = decode_le_u64(raw + 20);
    meta.dtype = static_cast<uint8_t>(raw[28]);
    return meta;
}

auto extract_payload_after_header(const std::vector<std::byte>& bytes)
    -> std::vector<std::byte>
{
    REQUIRE(bytes.size() >= kMetaSize);
    return std::vector<std::byte>(
        bytes.begin() + static_cast<std::ptrdiff_t>(kMetaSize), bytes.end());
}

template <typename eT>
auto expected_dtype_code() -> uint8_t
{
    if constexpr (std::is_same_v<eT, uint8_t>)
    {
        return 1;
    }
    if constexpr (std::is_same_v<eT, float>)
    {
        return 2;
    }
    return 3;
}

template <typename eT>
auto to_bytes(const Eigen::Matrix<eT, Eigen::Dynamic, 1>& vector)
    -> std::vector<std::byte>
{
    const size_t bytes_count = static_cast<size_t>(vector.size()) * sizeof(eT);
    std::vector<std::byte> out(bytes_count);
    if (bytes_count > 0)
    {
        std::memcpy(out.data(), vector.data(), bytes_count);
    }
    return out;
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

}  // namespace gelex::detail::test

template <typename eT>
using TestVector = Eigen::Matrix<eT, Eigen::Dynamic, 1>;

using namespace gelex::detail;        // NOLINT
using namespace gelex::detail::test;  // NOLINT

#define GELEX_BINARY_WRITER_TYPES uint8_t, float, double

TEMPLATE_TEST_CASE(
    "BinaryWriter - explicit finish writes header meta and payload",
    "[data][binary_writer]",
    GELEX_BINARY_WRITER_TYPES)
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    TestVector<TestType> record = make_vector<TestType>(
        {static_cast<TestType>(1),
         static_cast<TestType>(2),
         static_cast<TestType>(3)});

    REQUIRE_NOTHROW(
        [&]()
        {
            BinaryWriter<TestType> writer(file_path_str);
            writer.write(record);
            writer.finish();
        }());

    const auto bytes = read_all_bytes(file_path);
    REQUIRE(bytes.size() == kMetaSize + to_bytes(record).size());

    const std::vector<std::byte> payload = extract_payload_after_header(bytes);
    REQUIRE(payload == to_bytes(record));

    const MetaView meta = parse_header_meta(bytes);
    REQUIRE(meta.magic == kExpectedMagic);
    REQUIRE(meta.version == kExpectedVersion);
    REQUIRE(meta.n_rows == static_cast<uint64_t>(record.size()));
    REQUIRE(meta.n_cols == 1);
    REQUIRE(meta.dtype == expected_dtype_code<TestType>());
}

TEMPLATE_TEST_CASE(
    "BinaryWriter - multiple write accumulates payload",
    "[data][binary_writer]",
    GELEX_BINARY_WRITER_TYPES)
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    TestVector<TestType> r1 = make_vector<TestType>(
        {static_cast<TestType>(10), static_cast<TestType>(11)});
    TestVector<TestType> r2 = make_vector<TestType>(
        {static_cast<TestType>(20), static_cast<TestType>(21)});
    TestVector<TestType> r3 = make_vector<TestType>(
        {static_cast<TestType>(30), static_cast<TestType>(31)});

    REQUIRE_NOTHROW(
        [&]()
        {
            BinaryWriter<TestType> writer(file_path_str);
            writer.write(r1);
            writer.write(r2);
            writer.write(r3);
            writer.finish();
        }());

    std::vector<std::byte> expected_payload;
    const auto b1 = to_bytes(r1);
    const auto b2 = to_bytes(r2);
    const auto b3 = to_bytes(r3);
    expected_payload.insert(expected_payload.end(), b1.begin(), b1.end());
    expected_payload.insert(expected_payload.end(), b2.begin(), b2.end());
    expected_payload.insert(expected_payload.end(), b3.begin(), b3.end());

    const auto bytes = read_all_bytes(file_path);
    REQUIRE(bytes.size() == kMetaSize + expected_payload.size());
    REQUIRE(extract_payload_after_header(bytes) == expected_payload);

    const MetaView meta = parse_header_meta(bytes);
    REQUIRE(meta.n_rows == static_cast<uint64_t>(r1.size()));
    REQUIRE(meta.n_cols == 3);
    REQUIRE(meta.dtype == expected_dtype_code<TestType>());
}

TEST_CASE(
    "BinaryWriter - zero write then finish writes only header meta",
    "[data][binary_writer]")
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    REQUIRE_NOTHROW(
        [&]()
        {
            BinaryWriter<double> writer(file_path_str);
            writer.finish();
        }());

    const auto bytes = read_all_bytes(file_path);
    REQUIRE(bytes.size() == kMetaSize);

    const MetaView meta = parse_header_meta(bytes);
    REQUIRE(meta.magic == kExpectedMagic);
    REQUIRE(meta.version == kExpectedVersion);
    REQUIRE(meta.n_rows == 0);
    REQUIRE(meta.n_cols == 0);
    REQUIRE(meta.dtype == expected_dtype_code<double>());
}

TEMPLATE_TEST_CASE(
    "BinaryWriter - destructor auto finish",
    "[data][binary_writer]",
    GELEX_BINARY_WRITER_TYPES)
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    TestVector<TestType> record = make_vector<TestType>(
        {static_cast<TestType>(7), static_cast<TestType>(9)});

    REQUIRE_NOTHROW(
        [&]()
        {
            BinaryWriter<TestType> writer(file_path_str);
            writer.write(record);
        }());

    const auto bytes = read_all_bytes(file_path);
    REQUIRE(bytes.size() == kMetaSize + to_bytes(record).size());

    const MetaView meta = parse_header_meta(bytes);
    REQUIRE(meta.n_rows == static_cast<uint64_t>(record.size()));
    REQUIRE(meta.n_cols == 1);
}

TEST_CASE(
    "BinaryWriter - explicit finish does not duplicate header meta",
    "[data][binary_writer]")
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    {
        BinaryWriter<double> writer(file_path_str);
        const auto record = make_vector<double>({1.0, 2.0, 3.0});
        writer.write(record);
        writer.finish();
    }

    const auto bytes = read_all_bytes(file_path);
    const size_t expected_size = kMetaSize + 3 * sizeof(double);
    REQUIRE(bytes.size() == expected_size);

    const MetaView meta = parse_header_meta(bytes);
    REQUIRE(meta.n_rows == 3);
    REQUIRE(meta.n_cols == 1);
    REQUIRE(meta.dtype == expected_dtype_code<double>());
}

TEMPLATE_TEST_CASE(
    "BinaryWriter - inconsistent record length throws",
    "[data][binary_writer]",
    GELEX_BINARY_WRITER_TYPES)
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    TestVector<TestType> first = make_vector<TestType>(
        {static_cast<TestType>(1),
         static_cast<TestType>(2),
         static_cast<TestType>(3)});
    TestVector<TestType> second = make_vector<TestType>(
        {static_cast<TestType>(4), static_cast<TestType>(5)});

    BinaryWriter<TestType> writer(file_path_str);
    writer.write(first);

    REQUIRE_THROWS_AS(writer.write(second), gelex::ArgumentValidationException);
}

TEMPLATE_TEST_CASE(
    "BinaryWriter - write after finish throws",
    "[data][binary_writer]",
    GELEX_BINARY_WRITER_TYPES)
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".bin");
    const std::string file_path_str = file_path.string();

    TestVector<TestType> record = make_vector<TestType>(
        {static_cast<TestType>(1), static_cast<TestType>(2)});

    BinaryWriter<TestType> writer(file_path_str);
    writer.write(record);
    writer.finish();

    REQUIRE_THROWS_AS(writer.write(record), gelex::InvalidOperationException);
}

TEST_CASE("BinaryWriter - directory path should throw", "[data][binary_writer]")
{
    FileFixture files;
    const auto dir_path = files.get_test_dir();
    const std::string dir_path_str = dir_path.string();

    REQUIRE_THROWS_AS(
        BinaryWriter<double>(std::string_view(dir_path_str)),
        gelex::FileOpenException);
}

#else

TEST_CASE(
    "BinaryWriter - tests are gated until implementation header exists",
    "[data][binary_writer]")
{
    SUCCEED(
        "binary_writer.h is not present yet; full tests are compiled once it "
        "exists");
}

#endif
