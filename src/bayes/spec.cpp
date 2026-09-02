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

#include "gelex/bayes/spec.h"

#include <cmath>
#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <string_view>

#include "gelex/exception.h"

namespace gelex
{

namespace
{

auto validate_open_probability(double value, std::string_view name) -> void
{
    if (!std::isfinite(value) || value <= 0.0 || value >= 1.0)
    {
        throw GelexException(
            fmt::format(
                "{} must lie in the open interval (0, 1), got {}",
                name,
                value));
    }
}

template <std::size_t Size>
auto validate_probability_simplex(
    const std::array<double, Size>& probabilities,
    std::string_view name) -> void
{
    auto total = 0.0;
    for (const auto [index, probability] :
         probabilities | std::views::enumerate)
    {
        if (!std::isfinite(probability) || probability <= 0.0)
        {
            throw GelexException(
                fmt::format(
                    "{}[{}] must be finite and positive, got {}",
                    name,
                    index,
                    probability));
        }
        total += probability;
    }

    constexpr double simplex_tolerance = 1e-9;
    if (!std::isfinite(total) || std::abs(total - 1.0) > simplex_tolerance)
    {
        throw GelexException(
            fmt::format("{} must sum to 1, got {}", name, total));
    }
}

}  // namespace

HalfNormal::HalfNormal() : HalfNormal{default_positive_probability} {}

HalfNormal::HalfNormal(double positive_probability)
    : positive_probability_{positive_probability}
{
    validate_open_probability(
        positive_probability_, "half-normal positive probability");
}

SpikeSlab::SpikeSlab() : SpikeSlab{default_probability} {}

SpikeSlab::SpikeSlab(double probability) : probability_{probability}
{
    validate_open_probability(probability_, "spike-slab inclusion probability");
}

ScaledMixture::ScaledMixture()
    : ScaledMixture{default_probabilities, default_scales}
{
}

ScaledMixture::ScaledMixture(std::array<double, 5> probabilities)
    : ScaledMixture{probabilities, default_scales}
{
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
ScaledMixture::ScaledMixture(
    std::array<double, 5> probabilities,
    std::array<double, 5> scales)
    : probabilities_{probabilities}, scales_{scales}
{
    validate_probability_simplex(
        probabilities_, "scaled-mixture probabilities");
    if (scales_.front() != 0.0)
    {
        throw GelexException(
            fmt::format(
                "scaled-mixture scales[0] must be zero, got {}",
                scales_.front()));
    }
    for (const auto [index, scale] :
         scales_ | std::views::drop(1) | std::views::enumerate)
    {
        if (!std::isfinite(scale) || scale <= 0.0)
        {
            throw GelexException(
                fmt::format(
                    "scaled-mixture scales[{}] must be finite and positive, "
                    "got {}",
                    index + 1,
                    scale));
        }
    }
}
// NOLINTEND(bugprone-easily-swappable-parameters)

JointSpikeSlab::JointSpikeSlab() : JointSpikeSlab{default_probabilities} {}

JointSpikeSlab::JointSpikeSlab(std::array<double, 4> probabilities)
    : probabilities_{probabilities}
{
    validate_probability_simplex(
        probabilities_, "joint spike-slab probabilities");
}

}  // namespace gelex
