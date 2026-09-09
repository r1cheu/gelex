// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_BINARY_FORMAT_H_
#define GELEX_IO_BINARY_FORMAT_H_

#include <array>
#include <concepts>
#include <cstdint>
#include <string>
#include <type_traits>

namespace gelex
{

using BinaryShape = std::array<std::uint64_t, 2>;

enum class BinaryType : std::uint8_t
{
    float64 = 1,
    float32 = 2,
    uint8 = 4,
};

// One named matrix in a dense or CSC file: the index entry minus offsets.
struct MatrixHeader
{
    std::string identifier;
    BinaryType type{};
    BinaryShape shape{};

    auto operator==(const MatrixHeader&) const -> bool = default;
};

namespace detail
{

template <typename T>
concept SupportedDtype = std::same_as<std::remove_cv_t<T>, double>
                         || std::same_as<std::remove_cv_t<T>, float>
                         || std::same_as<std::remove_cv_t<T>, std::uint8_t>;

template <SupportedDtype T>
inline constexpr BinaryType binary_type_for = []
{
    using Value = std::remove_cv_t<T>;
    if constexpr (std::same_as<Value, double>)
    {
        return BinaryType::float64;
    }
    else if constexpr (std::same_as<Value, float>)
    {
        return BinaryType::float32;
    }
    else
    {
        return BinaryType::uint8;
    }
}();

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_IO_BINARY_FORMAT_H_
