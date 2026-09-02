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

#ifndef GELEX_IO_DETAIL_BINARY_FORMAT_H_
#define GELEX_IO_DETAIL_BINARY_FORMAT_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <type_traits>

namespace gelex::detail
{

inline constexpr std::array<std::byte, 8> binary_format_magic
    = {std::byte{'G'},
       std::byte{'E'},
       std::byte{'L'},
       std::byte{'E'},
       std::byte{'X'},
       std::byte{'B'},
       std::byte{'F'},
       std::byte{'1'}};
inline constexpr size_t footer_size = 24;
inline constexpr size_t path_size = 256;
inline constexpr size_t dtype_offset = path_size;
inline constexpr size_t offset_offset = dtype_offset + 8;
inline constexpr size_t size_offset = offset_offset + 8;
inline constexpr size_t rows_offset = size_offset + 8;
inline constexpr size_t cols_offset = rows_offset + 8;
inline constexpr size_t toc_entry_size = cols_offset + 8;
inline constexpr size_t page_alignment = 4096;
inline constexpr size_t max_path_length = path_size - 1;

// String dtype: value 0x01 cannot conflict with arithmetic types since
// the smallest arithmetic dtype encodes as (sizeof(T) << 2) >= 4.
inline constexpr uint8_t type_string = 0x01;

template <typename eT>
    requires std::is_arithmetic_v<eT>
inline constexpr uint8_t type_byte = static_cast<uint8_t>(
    (sizeof(eT) << 2U)
    | (static_cast<uint32_t>(std::is_floating_point_v<eT>) << 1U)
    | static_cast<uint32_t>(std::is_signed_v<eT>));

template <typename T>
    requires std::is_trivially_copyable_v<T>
inline auto decode(const std::byte* data) -> T
{
    T value;
    std::memcpy(&value, data, sizeof(value));
    return value;
}

template <typename T>
    requires std::is_trivially_copyable_v<T>
inline auto encode(T value, std::byte* out) -> void
{
    std::memcpy(out, &value, sizeof(value));
}

inline auto path_as_view(const std::array<char, path_size>& p)
    -> std::string_view
{
    auto len = std::find(p.begin(), p.end(), '\0') - p.begin();
    return {p.data(), static_cast<size_t>(len)};
}

struct TocEntry
{
    std::array<char, path_size> path{};
    uint8_t dtype{};
    uint64_t offset{};
    uint64_t size{};
    uint64_t rows{};
    uint64_t cols{};

    auto to_bytes() const -> std::array<std::byte, toc_entry_size>
    {
        std::array<std::byte, toc_entry_size> buf{};
        std::memcpy(buf.data(), path.data(), path_size);
        buf[dtype_offset] = static_cast<std::byte>(dtype);
        encode(offset, &buf[offset_offset]);
        encode(size, &buf[size_offset]);
        encode(rows, &buf[rows_offset]);
        encode(cols, &buf[cols_offset]);
        return buf;
    }

    static auto from_bytes(const std::byte* buf) -> TocEntry
    {
        TocEntry e;
        std::memcpy(e.path.data(), buf, path_size);
        e.dtype = static_cast<uint8_t>(buf[dtype_offset]);
        e.offset = decode<uint64_t>(buf + offset_offset);
        e.size = decode<uint64_t>(buf + size_offset);
        e.rows = decode<uint64_t>(buf + rows_offset);
        e.cols = decode<uint64_t>(buf + cols_offset);
        return e;
    }
};

}  // namespace gelex::detail

#endif  // GELEX_IO_DETAIL_BINARY_FORMAT_H_
