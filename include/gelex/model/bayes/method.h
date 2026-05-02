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

#ifndef GELEX_MODEL_BAYES_METHOD_H_
#define GELEX_MODEL_BAYES_METHOD_H_
#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

#include <fmt/base.h>
#include <fmt/format.h>

#include "gelex/types/genetic_effect_type.h"

namespace gelex
{
enum class BayesBase : uint8_t
{
    A,
    B,
    C,
    R,
    RR,
};

struct BayesMethodConfig
{
    BayesBase base{};
    GeneticMode mode = GeneticMode::A;
    bool asymmetric = false;
    bool estimate_pi = false;

    constexpr auto operator==(const BayesMethodConfig&) const -> bool = default;
};

// {base, mode, asymmetric, estimate_pi}
inline constexpr auto kValidMethods = std::array<BayesMethodConfig, 22>{{
    {BayesBase::A, GeneticMode::A, false, false},
    {BayesBase::A, GeneticMode::D, false, false},
    {BayesBase::A, GeneticMode::AD, false, false},
    {BayesBase::B, GeneticMode::A, false, false},
    {BayesBase::B, GeneticMode::D, false, false},
    {BayesBase::B, GeneticMode::AD, false, false},
    {BayesBase::B, GeneticMode::A, false, true},
    {BayesBase::B, GeneticMode::D, false, true},
    {BayesBase::B, GeneticMode::AD, false, true},
    {BayesBase::C, GeneticMode::A, false, false},
    {BayesBase::C, GeneticMode::D, false, false},
    {BayesBase::C, GeneticMode::AD, false, false},
    {BayesBase::C, GeneticMode::A, false, true},
    {BayesBase::C, GeneticMode::D, false, true},
    {BayesBase::C, GeneticMode::AD, false, true},
    {BayesBase::R, GeneticMode::A, false, false},
    {BayesBase::R, GeneticMode::D, false, false},
    {BayesBase::R, GeneticMode::AD, false, false},
    {BayesBase::R, GeneticMode::AD, true, false},
    {BayesBase::RR, GeneticMode::A, false, false},
    {BayesBase::RR, GeneticMode::D, false, false},
    {BayesBase::RR, GeneticMode::AD, false, false},
}};

constexpr auto is_valid_method(const BayesMethodConfig& m) -> bool
{
    return std::ranges::find(kValidMethods, m) != kValidMethods.end();
}

inline auto get_bayes_base(std::string_view sv) -> std::optional<BayesBase>
{
    static const std::unordered_map<std::string_view, BayesBase> string_to_base
        = {
            {"A", BayesBase::A},
            {"B", BayesBase::B},
            {"C", BayesBase::C},
            {"R", BayesBase::R},
            {"RR", BayesBase::RR},
        };

    auto it = string_to_base.find(sv);
    if (it != string_to_base.end())
    {
        return it->second;
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
        string_view name = "unknown";
        switch (b)
        {
            case gelex::BayesBase::A:
                name = "A";
                break;
            case gelex::BayesBase::B:
                name = "B";
                break;
            case gelex::BayesBase::C:
                name = "C";
                break;
            case gelex::BayesBase::R:
                name = "R";
                break;
            case gelex::BayesBase::RR:
                name = "RR";
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

template <>
struct formatter<gelex::BayesMethodConfig> : formatter<string_view>
{
    static auto format(const gelex::BayesMethodConfig& c, format_context& ctx)
        -> format_context::iterator
    {
        auto name = fmt::format("Bayes{}", c.base);
        if (c.estimate_pi)
        {
            name += "pi";
        }
        if (c.asymmetric)
        {
            name += " + asymmetric dominance";
        }
        else if (c.mode == gelex::GeneticMode::D)
        {
            name += " (dominance only)";
        }
        else if (c.mode == gelex::GeneticMode::AD)
        {
            name += " + dominance";
        }
        return fmt::format_to(ctx.out(), "{}", name);
    }
};

}  // namespace fmt

#endif  // GELEX_MODEL_BAYES_METHOD_H_
