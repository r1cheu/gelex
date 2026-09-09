// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_DETAIL_BINARY_WIRE_H_
#define GELEX_IO_DETAIL_BINARY_WIRE_H_

#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fmt/format.h>
#include <limits>
#include <span>
#include <string>
#include <string_view>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"

namespace gelex::detail
{

static_assert(std::endian::native == std::endian::little);
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559);
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);
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

// CSC v1: values, indices and indptr arrays per matrix, then an index and
// the same 24-byte footer layout as GELEXBF2.
inline constexpr std::array<std::byte, 8> csc_format_magic
    = {std::byte{'G'},
       std::byte{'E'},
       std::byte{'L'},
       std::byte{'E'},
       std::byte{'X'},
       std::byte{'S'},
       std::byte{'C'},
       std::byte{'1'}};
inline constexpr std::size_t csc_array_count = 3;
inline constexpr std::size_t csc_index_entry_fixed_size
    = sizeof(std::uint32_t) + sizeof(std::uint8_t)
      + ((3 + csc_array_count) * sizeof(std::uint64_t));

struct MatrixEntry
{
    MatrixHeader header;
    std::uint64_t offset{};
    std::uint64_t size{};
};

struct CscEntry
{
    MatrixHeader header;
    std::uint64_t nnz{};
    // values, indices, indptr
    std::array<std::uint64_t, csc_array_count> offsets{};
};

template <std::unsigned_integral T>
inline auto read_integer(const std::byte* data) -> T
{
    T value;
    std::memcpy(&value, data, sizeof(value));
    return value;
}

// Consumes a byte span front to back; running out throws with the field
// name.
class ByteCursor
{
   public:
    ByteCursor(std::span<const std::byte> bytes, std::string_view path) noexcept
        : remaining_(bytes), path_(path)
    {
    }

    auto read_bytes(std::size_t size, std::string_view field)
        -> std::span<const std::byte>
    {
        if (size > remaining_.size())
        {
            throw GelexException(fmt::format("{}: truncated {}", path_, field));
        }
        const auto bytes = remaining_.first(size);
        remaining_ = remaining_.subspan(size);
        return bytes;
    }

    template <std::unsigned_integral T>
    auto read(std::string_view field) -> T
    {
        return read_integer<T>(read_bytes(sizeof(T), field).data());
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t
    {
        return remaining_.size();
    }

    [[nodiscard]] auto empty() const noexcept -> bool
    {
        return remaining_.empty();
    }

   private:
    std::span<const std::byte> remaining_;
    std::string_view path_;
};

inline auto decode_binary_type(std::byte byte) -> BinaryType
{
    const auto type
        = static_cast<BinaryType>(std::to_integer<std::uint8_t>(byte));
    switch (type)
    {
        case BinaryType::float64:
        case BinaryType::float32:
        case BinaryType::uint8:
            return type;
    }
    throw GelexException("unknown binary payload type");
}

// Both index formats start an entry with the identifier, dtype and shape.
inline auto read_matrix_header(ByteCursor& cursor) -> MatrixHeader
{
    MatrixHeader header;
    const auto identifier_size = cursor.read<std::uint32_t>("index entry");
    const auto identifier_bytes
        = cursor.read_bytes(identifier_size, "matrix identifier");
    header.identifier.assign(
        reinterpret_cast<const char*>(identifier_bytes.data()),
        identifier_bytes.size());
    header.type
        = decode_binary_type(cursor.read_bytes(1, "matrix dtype").front());
    header.shape[0] = cursor.read<std::uint64_t>("matrix shape");
    header.shape[1] = cursor.read<std::uint64_t>("matrix shape");
    return header;
}

inline auto binary_type_size(BinaryType type) -> std::uint64_t
{
    switch (type)
    {
        case BinaryType::float64:
            return sizeof(double);
        case BinaryType::float32:
            return sizeof(float);
        case BinaryType::uint8:
            return sizeof(std::uint8_t);
    }
    throw GelexException("unknown binary payload type");
}

// Byte sizes of the values, indices and indptr arrays.
inline auto csc_array_sizes(
    BinaryType type,
    std::uint64_t columns,
    std::uint64_t nnz) -> std::array<std::uint64_t, csc_array_count>
{
    return {
        nnz * binary_type_size(type),
        nnz * sizeof(std::int64_t),
        (columns + 1) * sizeof(std::int64_t)};
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

inline auto align_payload_offset(std::uint64_t offset) -> std::uint64_t
{
    const auto remainder = offset % payload_alignment;
    return remainder == 0 ? offset : offset + payload_alignment - remainder;
}

}  // namespace gelex::detail

#endif  // GELEX_IO_DETAIL_BINARY_WIRE_H_
