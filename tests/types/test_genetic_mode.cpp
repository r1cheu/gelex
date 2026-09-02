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

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <ranges>
#include <string_view>
#include <vector>

#include "gelex/types/genetic_mode.h"

using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::name_of;

namespace
{

template <GeneticModeSet>
inline constexpr bool ACCEPTS_GENETIC_MODE_SET_NTTP = true;

static_assert(ACCEPTS_GENETIC_MODE_SET_NTTP<GeneticMode::A | GeneticMode::D>);

static_assert(GeneticModeSet{}.size() == 0);
static_assert(!GeneticModeSet{}.contains(GeneticMode::A));
static_assert(!GeneticModeSet{}.contains(GeneticMode::D));

inline constexpr auto STRAY_BIT_MODES = []
{
    auto modes = GeneticModeSet{GeneticMode::D};
    modes.bits |= 0b1000'0000U;
    return modes;
}();

static_assert(STRAY_BIT_MODES.size() == 1);
static_assert(STRAY_BIT_MODES.contains(GeneticMode::D));
static_assert(!STRAY_BIT_MODES.contains(GeneticMode::A));
static_assert(STRAY_BIT_MODES.index_of(GeneticMode::D) == 0);

static_assert((GeneticMode::A | GeneticMode::D).index_of(GeneticMode::A) == 0);
static_assert((GeneticMode::A | GeneticMode::D).index_of(GeneticMode::D) == 1);
static_assert(GeneticModeSet{GeneticMode::D}.index_of(GeneticMode::D) == 0);

// An absent mode reports size(), mirroring std::find returning end().
static_assert(GeneticModeSet{GeneticMode::D}.index_of(GeneticMode::A) == 1);
static_assert(GeneticModeSet{}.index_of(GeneticMode::A) == 0);

static_assert(name_of(GeneticModeSet{GeneticMode::A}) == "A");
static_assert(name_of(GeneticModeSet{GeneticMode::D}) == "D");
static_assert(name_of(GeneticMode::A | GeneticMode::D) == "AD");

// Undefined bits are ignored here too, so a stray write cannot turn a labelled
// set into an unnamed one.
static_assert(name_of(STRAY_BIT_MODES) == "D");
static_assert(name_of(GeneticModeSet{}).empty());

}  // namespace

TEST_CASE("GeneticModeSet defaults to the empty set", "[types][genetic_mode]")
{
    const auto modes = GeneticModeSet{};

    REQUIRE(modes.size() == 0);
    REQUIRE_FALSE(modes.contains(GeneticMode::A));
    REQUIRE((modes.each() | std::ranges::to<std::vector>()).empty());
}

TEST_CASE("GeneticModeSet iterates in canonical order", "[types][genetic_mode]")
{
    const auto modes = GeneticMode::D | GeneticMode::A;
    const auto values = modes.each() | std::ranges::to<std::vector>();

    REQUIRE(values == std::vector{GeneticMode::A, GeneticMode::D});
}

TEST_CASE("GeneticModeSet::index_of inverts each()", "[types][genetic_mode]")
{
    const auto modes = GeneticMode::D | GeneticMode::A;

    for (const auto [index, mode] : modes.each() | std::views::enumerate)
    {
        REQUIRE(modes.index_of(mode) == static_cast<std::size_t>(index));
    }
}

TEST_CASE("GeneticModeSet reports membership and size", "[types][genetic_mode]")
{
    auto modes = GeneticModeSet{GeneticMode::A};
    modes |= GeneticModeSet{GeneticMode::D};

    REQUIRE(modes.contains(GeneticMode::A));
    REQUIRE(modes.contains(GeneticMode::D));
    REQUIRE(modes.size() == 2);
}
