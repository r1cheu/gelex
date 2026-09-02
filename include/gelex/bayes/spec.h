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
#include <cstddef>

namespace gelex
{

struct Gaussian
{
};

class SpikeSlab
{
    static constexpr double default_probability = 0.01;

   public:
    SpikeSlab();
    explicit SpikeSlab(double probability);

    [[nodiscard]] auto probability() const noexcept -> double
    {
        return probability_;
    }

   private:
    double probability_;
};

class HalfNormal
{
    static constexpr double default_positive_probability = 0.5;

   public:
    HalfNormal();
    explicit HalfNormal(double positive_probability);

    [[nodiscard]] auto positive_probability() const noexcept -> double
    {
        return positive_probability_;
    }

   private:
    double positive_probability_;
};

class ScaledMixture
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

    ScaledMixture();
    explicit ScaledMixture(std::array<double, class_count> probabilities);
    ScaledMixture(
        std::array<double, class_count> probabilities,
        std::array<double, class_count> scales);

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

class JointSpikeSlab
{
    static constexpr std::array default_probabilities{
        0.99,
        1.0 / 300,
        1.0 / 300,
        1.0 / 300};

   public:
    static constexpr std::size_t class_count = 4;  // null, A, D, AD

    JointSpikeSlab();
    explicit JointSpikeSlab(std::array<double, class_count> probabilities);

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
