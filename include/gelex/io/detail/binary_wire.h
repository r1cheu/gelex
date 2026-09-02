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

#ifndef GELEX_IO_DETAIL_BINARY_WIRE_H_
#define GELEX_IO_DETAIL_BINARY_WIRE_H_

#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"

namespace gelex::detail
{

static_assert(std::endian::native == std::endian::little);
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559);
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
static_assert(sizeof(std::int32_t) == 4);
static_assert(sizeof(std::uint8_t) == 1);

inline constexpr std::array<std::byte, 8> binary_format_magic
    = {std::byte{'G'},
       std::byte{'E'},
       std::byte{'L'},
       std::byte{'E'},
       std::byte{'X'},
       std::byte{'B'},
       std::byte{'F'},
       std::byte{'2'}};
inline constexpr std::size_t footer_size = 24;
inline constexpr std::size_t payload_entry_fixed_size
    = sizeof(std::uint32_t) + sizeof(std::uint8_t)
      + (4 * sizeof(std::uint64_t));
inline constexpr std::uint64_t payload_alignment = 64;

struct PayloadEntry
{
    PayloadInfo info;
    std::uint64_t offset{};
    std::uint64_t size{};
};

template <std::unsigned_integral T>
inline auto read_integer(const std::byte* data) -> T
{
    T value;
    std::memcpy(&value, data, sizeof(value));
    return value;
}

inline auto decode_binary_type(std::byte byte) -> BinaryType
{
    const auto type
        = static_cast<BinaryType>(std::to_integer<std::uint8_t>(byte));
    switch (type)
    {
        case BinaryType::float64:
        case BinaryType::float32:
        case BinaryType::int32:
        case BinaryType::uint8:
            return type;
    }
    throw GelexException("unknown binary payload type");
}

inline auto binary_type_size(BinaryType type) -> std::uint64_t
{
    switch (type)
    {
        case BinaryType::float64:
            return sizeof(double);
        case BinaryType::float32:
            return sizeof(float);
        case BinaryType::int32:
            return sizeof(std::int32_t);
        case BinaryType::uint8:
            return sizeof(std::uint8_t);
    }
    throw GelexException("unknown binary payload type");
}

inline auto checked_product(std::uint64_t lhs, std::uint64_t rhs)
    -> std::uint64_t
{
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs)
    {
        throw GelexException("binary shape element count overflow");
    }
    return lhs * rhs;
}

inline auto checked_add(std::uint64_t lhs, std::uint64_t rhs) -> std::uint64_t
{
    if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs)
    {
        throw GelexException("binary payload size overflow");
    }
    return lhs + rhs;
}

}  // namespace gelex::detail

#endif  // GELEX_IO_DETAIL_BINARY_WIRE_H_
