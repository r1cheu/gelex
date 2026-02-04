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

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/effect_sampler.h"
#include "gelex/exception.h"

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;

TEST_CASE("EffectSampler - validation", "[effect_sampler]")
{
    SECTION("Valid config does not throw")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, 1.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        REQUIRE_NOTHROW(EffectSampler(config));
    }

    SECTION("Empty effect classes throws")
    {
        EffectSampler::Config config{
            .add_classes = {},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        REQUIRE_THROWS_AS(EffectSampler(config), ArgumentValidationException);
    }

    SECTION("Proportions must sum to 1")
    {
        EffectSampler::Config config{
            .add_classes = {{0.3, 1.0}, {0.3, 1.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        REQUIRE_THROWS_AS(EffectSampler(config), ArgumentValidationException);
    }

    SECTION("Zero variance is allowed")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, 0.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        REQUIRE_NOTHROW(EffectSampler(config));
    }

    SECTION("Negative variance throws")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, -0.1}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        REQUIRE_THROWS_AS(EffectSampler(config), ArgumentValidationException);
    }

    SECTION("Negative proportion throws")
    {
        EffectSampler::Config config{
            .add_classes = {{-0.5, 1.0}, {1.5, 1.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        REQUIRE_THROWS_AS(EffectSampler(config), ArgumentValidationException);
    }

    SECTION("Dominance classes validated when has_dominance is true")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, 1.0}},
            .dom_classes = {{0.3, 1.0}, {0.3, 1.0}},
            .has_dominance = true,
            .seed = 42,
        };
        REQUIRE_THROWS_AS(EffectSampler(config), ArgumentValidationException);
    }

    SECTION("Dominance classes not validated when has_dominance is false")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, 1.0}},
            .dom_classes = {{0.3, 1.0}, {0.3, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        REQUIRE_NOTHROW(EffectSampler(config));
    }
}

TEST_CASE("EffectSampler - sampling", "[effect_sampler]")
{
    SECTION("Single class produces all same class assignments")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, 1.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        EffectSampler sampler(config);

        auto effects = sampler.sample(100);

        REQUIRE(effects.size() == 100);
        for (const auto& [idx, effect] : effects)
        {
            REQUIRE(effect.add_class == 0);
        }
    }

    SECTION("Multi-class proportions are approximately correct")
    {
        EffectSampler::Config config{
            .add_classes = {{0.5, 0.001}, {0.3, 0.01}, {0.2, 1.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        EffectSampler sampler(config);

        constexpr Eigen::Index N_SNPS = 1000;
        auto effects = sampler.sample(N_SNPS);

        REQUIRE(effects.size() == N_SNPS);

        std::array<int, 3> class_counts = {0, 0, 0};
        for (const auto& [idx, effect] : effects)
        {
            REQUIRE(effect.add_class >= 0);
            REQUIRE(effect.add_class < 3);
            class_counts[effect.add_class]++;
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
        EffectSampler::Config config{
            .add_classes = {{1.0, 0.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        EffectSampler sampler(config);

        auto effects = sampler.sample(100);

        for (const auto& [idx, effect] : effects)
        {
            REQUIRE(effect.additive == 0.0);
        }
    }

    SECTION("Dominance effects sampled when has_dominance is true")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, 1.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = true,
            .seed = 42,
        };
        EffectSampler sampler(config);

        auto effects = sampler.sample(100);

        bool has_nonzero_dominance = std::ranges::any_of(
            effects, [](const auto& kv) {
                return std::abs(kv.second.dominance) > 1e-10;
            });

        REQUIRE(has_nonzero_dominance);
    }

    SECTION("No dominance effects when has_dominance is false")
    {
        EffectSampler::Config config{
            .add_classes = {{1.0, 1.0}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 42,
        };
        EffectSampler sampler(config);

        auto effects = sampler.sample(100);

        for (const auto& [idx, effect] : effects)
        {
            REQUIRE(effect.dominance == 0.0);
        }
    }

    SECTION("Reproducibility with same seed")
    {
        EffectSampler::Config config{
            .add_classes = {{0.5, 1.0}, {0.5, 0.1}},
            .dom_classes = {{1.0, 1.0}},
            .has_dominance = false,
            .seed = 123,
        };

        EffectSampler sampler1(config);
        auto effects1 = sampler1.sample(50);

        EffectSampler sampler2(config);
        auto effects2 = sampler2.sample(50);

        for (Eigen::Index i = 0; i < 50; ++i)
        {
            REQUIRE(effects1[i].additive == effects2[i].additive);
            REQUIRE(effects1[i].add_class == effects2[i].add_class);
        }
    }
}
