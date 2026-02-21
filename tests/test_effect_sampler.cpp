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

#include <cmath>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/exception.h"

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;

TEST_CASE("EffectSampler - validation", "[effect_sampler]")
{
    SECTION("Valid config does not throw")
    {
        std::mt19937_64 rng(42);
        REQUIRE_NOTHROW(EffectSampler({{1.0, 1.0}}, {{1.0, 1.0}}, rng));
    }

    SECTION("Empty effect classes throws")
    {
        std::mt19937_64 rng(42);
        REQUIRE_THROWS_AS(
            EffectSampler({}, {{1.0, 1.0}}, rng), ArgumentValidationException);
    }

    SECTION("Proportions must sum to 1")
    {
        std::mt19937_64 rng(42);
        REQUIRE_THROWS_AS(
            EffectSampler({{0.3, 1.0}, {0.3, 1.0}}, {{1.0, 1.0}}, rng),
            ArgumentValidationException);
    }

    SECTION("Zero variance is allowed")
    {
        std::mt19937_64 rng(42);
        REQUIRE_NOTHROW(EffectSampler({{1.0, 0.0}}, {{1.0, 1.0}}, rng));
    }

    SECTION("Negative variance throws")
    {
        std::mt19937_64 rng(42);
        REQUIRE_THROWS_AS(
            EffectSampler({{1.0, -0.1}}, {{1.0, 1.0}}, rng),
            ArgumentValidationException);
    }

    SECTION("Negative proportion throws")
    {
        std::mt19937_64 rng(42);
        REQUIRE_THROWS_AS(
            EffectSampler({{-0.5, 1.0}, {1.5, 1.0}}, {{1.0, 1.0}}, rng),
            ArgumentValidationException);
    }

    SECTION("Dominance classes validated when provided")
    {
        std::mt19937_64 rng(42);
        REQUIRE_THROWS_AS(
            EffectSampler({{1.0, 1.0}}, {{0.3, 1.0}, {0.3, 1.0}}, rng),
            ArgumentValidationException);
    }

    SECTION("Dominance classes skipped when vector is empty")
    {
        std::mt19937_64 rng(42);
        REQUIRE_NOTHROW(EffectSampler({{1.0, 1.0}}, {}, rng));
    }
}

TEST_CASE("EffectSampler - sampling", "[effect_sampler]")
{
    SECTION("Single class produces all same class assignments")
    {
        std::mt19937_64 rng(42);
        EffectSampler sampler({{1.0, 1.0}}, {}, rng);

        auto effects = sampler.sample(100);

        REQUIRE(effects.size() == 100);
        for (Eigen::Index i = 0; i < effects.size(); ++i)
        {
            REQUIRE(effects.add_class(i) == 0);
        }
    }

    SECTION("Multi-class proportions are approximately correct")
    {
        std::mt19937_64 rng(42);
        EffectSampler sampler({{0.5, 0.001}, {0.3, 0.01}, {0.2, 1.0}}, {}, rng);

        constexpr Eigen::Index N_SNPS = 1000;
        auto effects = sampler.sample(N_SNPS);

        REQUIRE(effects.size() == N_SNPS);

        std::array<int, 3> class_counts = {0, 0, 0};
        for (Eigen::Index i = 0; i < effects.size(); ++i)
        {
            const int add_class = effects.add_class(i);
            REQUIRE(add_class >= 0);
            REQUIRE(add_class < 3);
            class_counts[add_class]++;
        }

        // Class proportions should be approximately 50%, 30%, 20%
        REQUIRE_THAT(
            static_cast<double>(class_counts[0]) / N_SNPS,
            WithinAbs(0.5, 0.05));
        REQUIRE_THAT(
            static_cast<double>(class_counts[1]) / N_SNPS,
            WithinAbs(0.3, 0.05));
        REQUIRE_THAT(
            static_cast<double>(class_counts[2]) / N_SNPS,
            WithinAbs(0.2, 0.05));
    }

    SECTION("Zero variance class produces zero effects")
    {
        std::mt19937_64 rng(42);
        EffectSampler sampler({{1.0, 0.0}}, {}, rng);

        auto effects = sampler.sample(100);

        for (Eigen::Index i = 0; i < effects.size(); ++i)
        {
            REQUIRE(effects.additive(i) == 0.0);
        }
    }

    SECTION("Dominance effects sampled when dominance classes exist")
    {
        std::mt19937_64 rng(42);
        EffectSampler sampler({{1.0, 1.0}}, {{1.0, 1.0}}, rng);

        auto effects = sampler.sample(100);

        bool has_nonzero_dominance = false;
        for (Eigen::Index i = 0; i < effects.size(); ++i)
        {
            if (std::abs(effects.dominance(i)) > 1e-10)
            {
                has_nonzero_dominance = true;
                break;
            }
        }

        REQUIRE(has_nonzero_dominance);
    }

    SECTION("Reproducibility with same seed")
    {
        std::mt19937_64 rng1(123);
        EffectSampler sampler1({{0.5, 1.0}, {0.5, 0.1}}, {}, rng1);
        auto effects1 = sampler1.sample(50);

        std::mt19937_64 rng2(123);
        EffectSampler sampler2({{0.5, 1.0}, {0.5, 0.1}}, {}, rng2);
        auto effects2 = sampler2.sample(50);

        for (Eigen::Index i = 0; i < 50; ++i)
        {
            REQUIRE(effects1.additive(i) == effects2.additive(i));
            REQUIRE(effects1.add_class(i) == effects2.add_class(i));
        }
    }
}
