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
#include <concepts>
#include <cstddef>
#include <ranges>
#include <string_view>
#include <vector>

#include "gelex/genetic_mode.h"

using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::name_of;

namespace
{

template <GeneticModeSet>
inline constexpr bool accepts_genetic_mode_set_nttp = true;

static_assert(accepts_genetic_mode_set_nttp<GeneticMode::A | GeneticMode::D>);

inline constexpr auto stray_bit_modes = []
{
    auto modes = GeneticModeSet{GeneticMode::D};
    modes.bits |= 0b1000'0000U;
    return modes;
}();

static_assert(stray_bit_modes.size() == 1);
static_assert(stray_bit_modes.contains(GeneticMode::D));
static_assert(!stray_bit_modes.contains(GeneticMode::A));
static_assert(stray_bit_modes.index_of(GeneticMode::D) == 0);

static_assert((GeneticMode::A | GeneticMode::D).index_of(GeneticMode::A) == 0);
static_assert((GeneticMode::A | GeneticMode::D).index_of(GeneticMode::D) == 1);
static_assert(GeneticModeSet{GeneticMode::D}.index_of(GeneticMode::D) == 0);

// An absent mode reports size(), mirroring std::find returning end().
static_assert(GeneticModeSet{GeneticMode::D}.index_of(GeneticMode::A) == 1);

static_assert(name_of(GeneticModeSet{GeneticMode::A}) == "A");
static_assert(name_of(GeneticModeSet{GeneticMode::D}) == "D");
static_assert(name_of(GeneticMode::A | GeneticMode::D) == "AD");

// Undefined bits are ignored here too, so a stray write cannot turn a labelled
// set into an unnamed one.
static_assert(name_of(stray_bit_modes) == "D");

// The type has no default constructor, so a mode set can never start out empty.
static_assert(!std::default_initializable<GeneticModeSet>);

}  // namespace

TEST_CASE("GeneticModeSet iterates in canonical order", "[genetic_mode]")
{
    const auto modes = GeneticMode::D | GeneticMode::A;
    const auto values = modes.each() | std::ranges::to<std::vector>();

    REQUIRE(values == std::vector{GeneticMode::A, GeneticMode::D});
}

TEST_CASE("GeneticModeSet::index_of inverts each()", "[genetic_mode]")
{
    const auto modes = GeneticMode::D | GeneticMode::A;

    for (const auto [index, mode] : modes.each() | std::views::enumerate)
    {
        REQUIRE(modes.index_of(mode) == static_cast<std::size_t>(index));
    }
}

TEST_CASE("GeneticModeSet reports membership and size", "[genetic_mode]")
{
    auto modes = GeneticModeSet{GeneticMode::A};
    modes |= GeneticModeSet{GeneticMode::D};

    REQUIRE(modes.contains(GeneticMode::A));
    REQUIRE(modes.contains(GeneticMode::D));
    REQUIRE(modes.size() == 2);
}
