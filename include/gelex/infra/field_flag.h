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

#ifndef GELEX_INFRA_FIELD_FLAG_H_
#define GELEX_INFRA_FIELD_FLAG_H_

#include <cstdint>

namespace gelex
{

enum class FieldFlag : std::uint8_t
{
    none = 0,
    checkpoint = 1U << 0U,
    trace = 1U << 1U,
    report = 1U << 2U,
    summary = 1U << 3U,
};

constexpr auto operator|(FieldFlag lhs, FieldFlag rhs) -> FieldFlag
{
    return static_cast<FieldFlag>(
        static_cast<std::uint8_t>(lhs) | static_cast<std::uint8_t>(rhs));
}

constexpr auto operator|=(FieldFlag& lhs, FieldFlag rhs) -> FieldFlag&
{
    lhs = lhs | rhs;
    return lhs;
}

constexpr auto has(FieldFlag flags, FieldFlag flag) -> bool
{
    const auto bits = static_cast<std::uint8_t>(flags);
    const auto mask = static_cast<std::uint8_t>(flag);
    return mask != 0 && (bits & mask) == mask;
}

}  // namespace gelex

#endif  // GELEX_INFRA_FIELD_FLAG_H_
