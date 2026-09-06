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

#ifndef GELEX_BAYES_SPEC_H_
#define GELEX_BAYES_SPEC_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <string_view>

#include "gelex/bayes/genetic_policy.h"
#include "gelex/exception.h"

namespace gelex
{

namespace detail
{

inline auto validate_probability_simplex(
    std::span<const double> probabilities,
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

}  // namespace detail

template <VarianceLayout Kind = VarianceLayout::Pooled>
struct GaussianSpec
{
};

template <
    VarianceLayout Kind = VarianceLayout::Pooled,
    MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
class SpikeSlabSpec
{
    static constexpr double default_probability = 0.01;

   public:
    SpikeSlabSpec() : SpikeSlabSpec{default_probability} {}

    explicit SpikeSlabSpec(double probability) : probability_{probability}
    {
        if (!std::isfinite(probability_) || probability_ <= 0.0
            || probability_ >= 1.0)
        {
            throw GelexException(
                fmt::format(
                    "spike-slab inclusion probability must lie in the open "
                    "interval (0, 1), got {}",
                    probability_));
        }
    }

    [[nodiscard]] auto probability() const noexcept -> double
    {
        return probability_;
    }

   private:
    double probability_;
};

struct HalfNormalSpec
{
};

template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
class ScaledMixtureSpec
{
    static constexpr std::array default_probabilities{
        0.99,
        0.005,
        0.003,
        0.001,
        0.001};
    static constexpr std::array default_scales{0.0, 0.001, 0.01, 0.1, 1.0};

   public:
    static constexpr std::size_t class_count
        = 5;  // null, small, medium, large, xlarge

    ScaledMixtureSpec()
        : ScaledMixtureSpec{default_probabilities, default_scales}
    {
    }

    explicit ScaledMixtureSpec(std::array<double, class_count> probabilities)
        : ScaledMixtureSpec{probabilities, default_scales}
    {
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    ScaledMixtureSpec(
        std::array<double, class_count> probabilities,
        std::array<double, class_count> scales)
        : probabilities_{probabilities}, scales_{scales}
    {
        detail::validate_probability_simplex(
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
                        "scaled-mixture scales[{}] must be finite and "
                        "positive, "
                        "got {}",
                        index + 1,
                        scale));
            }
        }
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    [[nodiscard]] auto probabilities() const noexcept
        -> const std::array<double, class_count>&
    {
        return probabilities_;
    }

    [[nodiscard]] auto scales() const noexcept
        -> const std::array<double, class_count>&
    {
        return scales_;
    }

   private:
    std::array<double, class_count> probabilities_;
    std::array<double, class_count> scales_;
};

template <MixtureWeightUpdate WeightUpdate = MixtureWeightUpdate::Enabled>
class JointSpikeSlabSpec
{
    static constexpr std::array default_probabilities{
        0.99,
        1.0 / 300,
        1.0 / 300,
        1.0 / 300};

   public:
    static constexpr std::size_t class_count = 4;  // null, A, D, AD

    JointSpikeSlabSpec() : JointSpikeSlabSpec{default_probabilities} {}

    explicit JointSpikeSlabSpec(std::array<double, class_count> probabilities)
        : probabilities_{probabilities}
    {
        detail::validate_probability_simplex(
            probabilities_, "joint spike-slab probabilities");
    }

    [[nodiscard]] auto probabilities() const noexcept
        -> const std::array<double, class_count>&
    {
        return probabilities_;
    }

   private:
    std::array<double, class_count> probabilities_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_SPEC_H_
