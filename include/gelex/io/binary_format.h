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
    int32 = 3,
    uint8 = 4,
};

struct PayloadDescriptor
{
    BinaryType type{};
    BinaryShape shape{};

    auto operator==(const PayloadDescriptor&) const -> bool = default;
};

struct PayloadInfo
{
    std::string identifier;
    PayloadDescriptor descriptor;

    auto operator==(const PayloadInfo&) const -> bool = default;
};

namespace detail
{

template <typename T>
concept SupportedDtype = std::same_as<std::remove_cv_t<T>, double>
                         || std::same_as<std::remove_cv_t<T>, float>
                         || std::same_as<std::remove_cv_t<T>, std::int32_t>
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
    else if constexpr (std::same_as<Value, std::int32_t>)
    {
        return BinaryType::int32;
    }
    else
    {
        return BinaryType::uint8;
    }
}();

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_IO_BINARY_FORMAT_H_
