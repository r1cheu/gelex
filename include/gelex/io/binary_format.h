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
#include <format>
#include <limits>
#include <string_view>
#include <type_traits>

#include "gelex/exception.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::detail::binary_format
{

// Encodes type identity into a single byte: (sizeof << 2) | (float_bit << 1) |
// sign_bit Guarantees unique codes for all arithmetic types of different
// size/signedness/category.
template <typename eT>
    requires std::is_arithmetic_v<eT>
constexpr auto dtype_code() -> uint8_t
{
    return static_cast<uint8_t>(
        (sizeof(eT) << 2U)
        | (static_cast<uint8_t>(std::is_floating_point_v<eT>) << 1U)
        | static_cast<uint8_t>(std::is_signed_v<eT>));
}

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

// --- Shared validation utilities ---

template <size_t N>
inline auto validate_magic(
    const std::byte* data,
    const std::array<std::byte, N>& expected,
    std::string_view path,
    std::string_view label) -> void
{
    if (!std::equal(expected.begin(), expected.end(), data))
    {
        throw FileFormatException(
            std::format("{}: invalid {} magic", path, label));
    }
}

inline auto
checked_mul(size_t a, size_t b, std::string_view path, std::string_view ctx)
    -> size_t
{
    if (a == 0 || b == 0)
    {
        return 0;
    }
    if (a > (std::numeric_limits<size_t>::max() / b))
    {
        throw FileFormatException(
            std::format("{}: size overflow while computing {}", path, ctx));
    }
    return a * b;
}

inline auto
checked_add(size_t a, size_t b, std::string_view path, std::string_view ctx)
    -> size_t
{
    if (a > (std::numeric_limits<size_t>::max() - b))
    {
        throw FileFormatException(
            std::format("{}: size overflow while computing {}", path, ctx));
    }
    return a + b;
}

// --- Container format (GELEXBW2) constants ---

inline constexpr std::array<std::byte, 8> kBinaryFormatMagic
    = {std::byte{'G'},
       std::byte{'E'},
       std::byte{'L'},
       std::byte{'E'},
       std::byte{'X'},
       std::byte{'B'},
       std::byte{'W'},
       std::byte{'2'}};
inline constexpr uint32_t kBinaryFormatVersion = 2;
inline constexpr size_t kFooterSize = 32;
inline constexpr size_t kTocEntrySize = 40;
inline constexpr size_t kPageAlignment = 4096;

}  // namespace gelex::detail::binary_format

namespace gelex::detail
{

enum class EffectType : uint8_t
{
    Add = 0,
    Dom = 1,
    Fixed = 2,
    Random = 3,
    Residual = 4,
};

enum class DataKind : uint8_t
{
    Coeff = 0,
    Mixture = 1,
    Sign = 2,
    Variance = 3,
    Pi = 4,
    PositiveProb = 5,
    Genotype = 6,
    SnpStats = 7,
    MonoIndices = 8,
};

struct SectionKey
{
    EffectType effect{};
    DataKind kind{};
    uint8_t index{0};

    auto operator==(const SectionKey&) const -> bool = default;
};

struct SectionKeyHash
{
    auto operator()(const SectionKey& k) const noexcept -> size_t
    {
        return (static_cast<size_t>(k.effect) << 16U)
               | (static_cast<size_t>(k.kind) << 8U)
               | static_cast<size_t>(k.index);
    }
};

constexpr auto to_section_effect_type(GeneticEffectType type) -> EffectType
{
    switch (type)
    {
        case GeneticEffectType::Add:
            return EffectType::Add;
        case GeneticEffectType::Dom:
            return EffectType::Dom;
    }
    __builtin_unreachable();
}

struct TocEntry
{
    SectionKey key;
    uint8_t dtype{};
    uint64_t offset{};
    uint64_t size{};
    uint64_t rows{};
    uint64_t cols{};
};

}  // namespace gelex::detail

#endif  // GELEX_IO_BINARY_FORMAT_H_
