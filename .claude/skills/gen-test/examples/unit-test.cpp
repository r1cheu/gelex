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
}
