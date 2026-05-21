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

#ifndef GELEX_MODEL_BAYES_LEGACY_BAYES_BASE_H_
#define GELEX_MODEL_BAYES_LEGACY_BAYES_BASE_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <string_view>
#include <utility>

#include <fmt/base.h>
#include <fmt/format.h>

namespace gelex
{

enum class BayesBase : uint8_t
{
    A,
    B,
    C,
    R,
    RR,
    kCount,
};

inline constexpr std::array<
    std::pair<BayesBase, std::string_view>,
    std::to_underlying(BayesBase::kCount)>
    kBayesBaseNames = {{
        {BayesBase::A, "A"},
        {BayesBase::B, "B"},
        {BayesBase::C, "C"},
        {BayesBase::R, "R"},
        {BayesBase::RR, "RR"},
    }};

inline auto get_bayes_base(std::string_view sv) -> std::optional<BayesBase>
{
    const auto* it = std::ranges::find_if(
        kBayesBaseNames, [sv](const auto& p) { return p.second == sv; });
    if (it != kBayesBaseNames.end())
    {
        return it->first;
    }
    return std::nullopt;
}

}  // namespace gelex

namespace fmt
{
template <>
struct formatter<gelex::BayesBase> : formatter<string_view>
{
    auto format(gelex::BayesBase b, format_context& ctx) const
        -> format_context::iterator
    {
        return formatter<string_view>::format(to_string_view(b), ctx);
    }

   private:
    static constexpr auto to_string_view(gelex::BayesBase b) -> string_view
    {
        for (const auto& [value, name] : gelex::kBayesBaseNames)
        {
            if (value == b)
            {
                return name;
            }
        }
        return "unknown";
    }
};

}  // namespace fmt

#endif  // GELEX_MODEL_BAYES_LEGACY_BAYES_BASE_H_
