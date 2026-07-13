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

#include "cli/completion/choice_descriptions.h"

#include <array>
#include <span>
#include <string_view>

namespace cli
{
namespace
{
constexpr std::array<ChoiceEntry, 10> CODING{
    {{"SH", "HWE standardize"},
     {"CH", "HWE center"},
     {"OSH", "HWE orthogonalized standardize"},
     {"OCH", "HWE orthogonalized center"},
     {"S", "Sample standardize"},
     {"C", "Sample center"},
     {"OS", "Sample orthogonalized standardize"},
     {"OC", "Sample orthogonalized center"},
     {"NS", "NOIA standardize"},
     {"NC", "NOIA center"}}};

constexpr std::array<ChoiceEntry, 6> MODEL{
    {{"RR", "BayesRR"},
     {"A", "BayesA"},
     {"B", "BayesB"},
     {"C", "BayesC"},
     {"R", "BayesR"},
     {"CD", "BayesCD"}}};

constexpr std::array<ChoiceEntry, 3> MODE{
    {{"A", "Additive"}, {"D", "Dominance"}, {"AD", "Additive + Dominance"}}};

constexpr std::array<ChoiceEntry, 3> TRANSFORM{
    {{"none", "No rank-INT"},
     {"dint", "Direct rank-INT"},
     {"iint", "Indirect rank-INT"}}};
}  // namespace

auto choice_group(std::string_view type_token) -> std::span<const ChoiceEntry>
{
    if (type_token == "CODING")
    {
        return CODING;
    }
    if (type_token == "MODEL")
    {
        return MODEL;
    }
    if (type_token == "MODE")
    {
        return MODE;
    }
    if (type_token == "TRANSFORM")
    {
        return TRANSFORM;
    }
    return {};
}
}  // namespace cli
