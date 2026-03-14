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

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

#include "gelex/types/genetic_effect_type.h"

namespace gelex::detail::binary_format
{
inline constexpr std::array<std::byte, 8> kBinaryFormatMagic
    = {std::byte{'G'},
       std::byte{'E'},
       std::byte{'L'},
       std::byte{'E'},
       std::byte{'X'},
       std::byte{'B'},
       std::byte{'F'},
       std::byte{'1'}};
inline constexpr size_t kFooterSize = 24;
inline constexpr size_t kTocEntrySize = 40;
inline constexpr size_t kPageAlignment = 4096;

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

}  // namespace gelex::detail::binary_format

namespace gelex::detail
{

enum class DataKind : uint8_t
{
    Coeff = 0,
    MixtureTracker = 1,
    SignTracker = 2,
    Variance = 3,
    MixtureProportion = 4,
    SignProportion = 5,
    Genotype = 6,
    LociStats = 7,
    MonoIndices = 8,
    GenoMethod = 9,
};

struct SectionKey
{
    gelex::EffectType effect{};
    DataKind kind{};
    uint8_t index{0};

    auto operator==(const SectionKey&) const -> bool = default;
};

struct SectionKeyHash
{
    auto operator()(const SectionKey& k) const noexcept -> size_t
    {
        return (static_cast<size_t>(k.effect.to_byte()) << 16U)
               | (static_cast<size_t>(k.kind) << 8U)
               | static_cast<size_t>(k.index);
    }
};

struct TocEntry
{
    SectionKey key;
    uint8_t dtype{};
    uint64_t offset{};
    uint64_t size{};
    uint64_t rows{};
    uint64_t cols{};

    auto to_bytes() const -> std::array<std::byte, binary_format::kTocEntrySize>
    {
        std::array<std::byte, binary_format::kTocEntrySize> buf{};
        buf[0] = static_cast<std::byte>(key.effect.to_byte());
        buf[1] = static_cast<std::byte>(key.kind);
        buf[2] = static_cast<std::byte>(dtype);
        buf[3] = static_cast<std::byte>(key.index);
        binary_format::encode(offset, &buf[8]);
        binary_format::encode(size, &buf[16]);
        binary_format::encode(rows, &buf[24]);
        binary_format::encode(cols, &buf[32]);
        return buf;
    }

    static auto from_bytes(const std::byte* buf) -> TocEntry
    {
        return {
            .key
            = SectionKey{.effect = gelex::EffectType::from_byte(static_cast<uint8_t>(buf[0])), .kind = static_cast<DataKind>(buf[1]), .index = static_cast<uint8_t>(buf[3])},
            .dtype = static_cast<uint8_t>(buf[2]),
            .offset = binary_format::decode<uint64_t>(buf + 8),
            .size = binary_format::decode<uint64_t>(buf + 16),
            .rows = binary_format::decode<uint64_t>(buf + 24),
            .cols = binary_format::decode<uint64_t>(buf + 32)};
    }
};

}  // namespace gelex::detail

#endif  // GELEX_IO_BINARY_FORMAT_H_
