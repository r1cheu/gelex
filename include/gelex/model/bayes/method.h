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

#include <cstdint>
#include <optional>
#include <string_view>
#include <unordered_map>

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
            case gelex::BayesBase::kCount:
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

}  // namespace fmt

#endif  // GELEX_MODEL_BAYES_METHOD_H_
