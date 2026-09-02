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
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <string>
#include <type_traits>
#include <utility>

#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"
#include "gelex/io/detail/binary_wire.h"

#include "file_fixture.h"

namespace
{

namespace fs = std::filesystem;
namespace test = gelex::test;

template <typename T>
    requires std::is_trivially_copyable_v<T>
auto read_value(std::ifstream& input) -> T
{
    T value{};
    input.read(
        reinterpret_cast<char*>(&value),
        static_cast<std::streamsize>(sizeof(value)));
    return value;
}

template <std::unsigned_integral T, std::size_t size>
auto read_integers(std::ifstream& input) -> std::array<T, size>
{
    std::array<T, size> values{};
    input.read(
        reinterpret_cast<char*>(values.data()),
        static_cast<std::streamsize>(sizeof(values)));
    return values;
}

}  // namespace

TEST_CASE("GELEXBF2 has a stable wire layout", "[io][binary_format]")
{
    test::FileFixture fixture;
    const auto container_path = fixture.get_test_dir() / "layout.samples";
    {
        gelex::BinaryWriter writer(container_path.string());
        writer.reserve<double>("x", gelex::BinaryShape{1, 1}).append(1.0);
        writer.close();
    }

    constexpr auto directory_offset = gelex::detail::payload_alignment;
    constexpr auto directory_size = gelex::detail::payload_entry_fixed_size + 1;
    REQUIRE(
        fs::file_size(container_path)
        == directory_offset + directory_size + gelex::detail::footer_size);

    std::ifstream input(container_path, std::ios::binary);
    REQUIRE(read_value<double>(input) == 1.0);

    input.seekg(static_cast<std::streamoff>(directory_offset));
    REQUIRE(read_value<std::uint32_t>(input) == 1);

    char identifier{};
    input.read(&identifier, 1);
    REQUIRE(identifier == 'x');
    REQUIRE(
        read_value<std::uint8_t>(input)
        == std::to_underlying(gelex::BinaryType::float64));

    const auto descriptor = read_integers<std::uint64_t, 4>(input);
    REQUIRE(
        descriptor == (std::array<std::uint64_t, 4>{1, 1, 0, sizeof(double)}));

    std::array<char, 8> magic{};
    input.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    REQUIRE(
        magic == (std::array<char, 8>{'G', 'E', 'L', 'E', 'X', 'B', 'F', '2'}));

    const auto footer = read_integers<std::uint64_t, 2>(input);
    REQUIRE(footer == (std::array<std::uint64_t, 2>{directory_offset, 1}));
    REQUIRE(input.peek() == std::ifstream::traits_type::eof());
}

TEST_CASE("GELEXBF2 dtype codes are stable", "[io][binary_format]")
{
    STATIC_REQUIRE(std::to_underlying(gelex::BinaryType::float64) == 1);
    STATIC_REQUIRE(std::to_underlying(gelex::BinaryType::float32) == 2);
    STATIC_REQUIRE(std::to_underlying(gelex::BinaryType::int32) == 3);
    STATIC_REQUIRE(std::to_underlying(gelex::BinaryType::uint8) == 4);

    REQUIRE(gelex::detail::binary_type_size(gelex::BinaryType::float64) == 8);
    REQUIRE(gelex::detail::binary_type_size(gelex::BinaryType::float32) == 4);
    REQUIRE(gelex::detail::binary_type_size(gelex::BinaryType::int32) == 4);
    REQUIRE(gelex::detail::binary_type_size(gelex::BinaryType::uint8) == 1);
}

TEST_CASE("GELEXBF2 stores variable-length identifiers", "[io][binary_format]")
{
    test::FileFixture fixture;
    const auto container_path = fixture.get_test_dir() / "identifier.samples";
    const std::string identifier(300, 'a');
    {
        gelex::BinaryWriter writer(container_path.string());
        writer.reserve<std::uint8_t>(identifier, gelex::BinaryShape{1, 1})
            .append(std::uint8_t{7});
        writer.close();
    }

    std::ifstream input(container_path, std::ios::binary);
    input.seekg(static_cast<std::streamoff>(gelex::detail::payload_alignment));
    REQUIRE(read_value<std::uint32_t>(input) == identifier.size());

    std::string decoded_identifier(identifier.size(), '\0');
    input.read(
        decoded_identifier.data(),
        static_cast<std::streamsize>(decoded_identifier.size()));
    REQUIRE(decoded_identifier == identifier);
}
