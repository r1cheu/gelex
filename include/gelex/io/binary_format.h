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

#ifndef GELEX_IO_BINARY_FORMAT_H_
#define GELEX_IO_BINARY_FORMAT_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <type_traits>

namespace gelex::io::detail
{

inline constexpr std::array<std::byte, 8> kBinaryFormatMagic
    = {std::byte{'G'},
       std::byte{'E'},
       std::byte{'L'},
       std::byte{'E'},
       std::byte{'X'},
       std::byte{'B'},
       std::byte{'F'},
       std::byte{'2'}};
inline constexpr size_t kFooterSize = 24;
inline constexpr size_t kTocEntrySize = 104;
inline constexpr size_t kPageAlignment = 4096;
inline constexpr size_t kMaxPathLength = 63;

// String dtype: value 0x01 cannot conflict with arithmetic types since
// the smallest arithmetic dtype encodes as (sizeof(T) << 2) >= 4.
inline constexpr uint8_t kTypeString = 0x01;

template <typename eT>
    requires std::is_arithmetic_v<eT>
inline constexpr uint8_t kTypeByte = static_cast<uint8_t>(
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

inline auto path_as_view(const std::array<char, 64>& p) -> std::string_view
{
    auto len = std::find(p.begin(), p.end(), '\0') - p.begin();
    return {p.data(), static_cast<size_t>(len)};
}

struct TocEntry
{
    std::array<char, 64> path{};
    uint8_t dtype{};
    uint64_t offset{};
    uint64_t size{};
    uint64_t rows{};
    uint64_t cols{};

    auto to_bytes() const -> std::array<std::byte, kTocEntrySize>
    {
        std::array<std::byte, kTocEntrySize> buf{};
        std::memcpy(buf.data(), path.data(), 64);
        buf[64] = static_cast<std::byte>(dtype);
        // buf[65..71] padding — already zero from aggregate init
        encode(offset, &buf[72]);
        encode(size, &buf[80]);
        encode(rows, &buf[88]);
        encode(cols, &buf[96]);
        return buf;
    }

    static auto from_bytes(const std::byte* buf) -> TocEntry
    {
        TocEntry e;
        std::memcpy(e.path.data(), buf, 64);
        e.dtype = static_cast<uint8_t>(buf[64]);
        e.offset = decode<uint64_t>(buf + 72);
        e.size = decode<uint64_t>(buf + 80);
        e.rows = decode<uint64_t>(buf + 88);
        e.cols = decode<uint64_t>(buf + 96);
        return e;
    }
};

}  // namespace gelex::io::detail

#endif  // GELEX_IO_BINARY_FORMAT_H_
