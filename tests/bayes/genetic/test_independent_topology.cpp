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

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <memory>
#include <type_traits>
#include <vector>

#include "gelex/bayes/genetic/independent_topology.h"

using gelex::generate_mode_values;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::IndependentTopology;
using gelex::transform_mode_values;

namespace
{

struct TestValue
{
    int value;
};

constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using AdditiveDominanceTopology = IndependentTopology<mode_ad, TestValue>;
using DominanceAdditiveTopology
    = IndependentTopology<GeneticMode::D | GeneticMode::A, TestValue>;

static_assert(
    AdditiveDominanceTopology::modes == (GeneticMode::A | GeneticMode::D));
static_assert(
    std::same_as<AdditiveDominanceTopology, DominanceAdditiveTopology>);
static_assert(std::same_as<AdditiveDominanceTopology::value_type, TestValue>);
static_assert(
    std::constructible_from<AdditiveDominanceTopology, TestValue, TestValue>);
static_assert(std::constructible_from<
              AdditiveDominanceTopology,
              std::array<TestValue, 2>>);

static_assert(
    AdditiveDominanceTopology{TestValue{1}, TestValue{2}}
        .get<GeneticMode::A>()
        .value
    == 1);
static_assert(
    AdditiveDominanceTopology{TestValue{1}, TestValue{2}}
        .get<GeneticMode::D>()
        .value
    == 2);

constexpr auto generated_topology = generate_mode_values<mode_ad>(
    [](GeneticMode mode) { return TestValue{mode == GeneticMode::A ? 1 : 2}; });

static_assert(std::same_as<
              std::remove_cvref_t<decltype(generated_topology)>,
              AdditiveDominanceTopology>);
static_assert(generated_topology.get<GeneticMode::A>().value == 1);
static_assert(generated_topology.get<GeneticMode::D>().value == 2);

constexpr auto transformed_topology = transform_mode_values(
    generated_topology,
    [](GeneticMode mode, const TestValue& value)
    { return value.value + (mode == GeneticMode::A ? 10 : 20); });

static_assert(std::same_as<
              std::remove_cvref_t<decltype(transformed_topology)>,
              IndependentTopology<mode_ad, int>>);
static_assert(transformed_topology.get<GeneticMode::A>() == 11);
static_assert(transformed_topology.get<GeneticMode::D>() == 22);

}  // namespace

TEST_CASE(
    "IndependentTopology stores values in canonical mode order",
    "[bayes][genetic][independent_topology]")
{
    auto topology = IndependentTopology<
        GeneticMode::D | GeneticMode::A,
        std::unique_ptr<int>>{
        std::make_unique<int>(1), std::make_unique<int>(2)};

    std::vector<GeneticMode> visited_modes;
    std::vector<int> visited_values;
    for (const auto& [mode, value] : topology.each())
    {
        visited_modes.push_back(mode);
        visited_values.push_back(*value);
    }

    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(visited_values == std::vector{1, 2});
}

TEST_CASE(
    "IndependentTopology provides mode-indexed value access",
    "[bayes][genetic][independent_topology]")
{
    auto topology
        = IndependentTopology<GeneticModeSet{GeneticMode::D}, TestValue>{
            TestValue{2}};

    topology.get<GeneticMode::D>().value = 3;
    const auto& constant_topology = topology;

    REQUIRE(constant_topology.get<GeneticMode::D>().value == 3);
}

TEST_CASE(
    "IndependentTopology maps every mode to its constructor argument",
    "[bayes][genetic][independent_topology]")
{
    const auto topology = AdditiveDominanceTopology{TestValue{1}, TestValue{2}};

    REQUIRE(topology.get<GeneticMode::A>().value == 1);
    REQUIRE(topology.get<GeneticMode::D>().value == 2);

    std::vector<int> visited_values;
    for (const auto& [mode, value] : topology.each())
    {
        visited_values.push_back(value.value);
    }

    REQUIRE(visited_values == std::vector{1, 2});
}

TEST_CASE(
    "generate_mode_values invokes the factory once in canonical order",
    "[bayes][genetic][independent_topology]")
{
    std::vector<GeneticMode> visited_modes;

    const auto topology = generate_mode_values<mode_ad>(
        [&](GeneticMode mode)
        {
            visited_modes.push_back(mode);
            return std::make_unique<int>(mode == GeneticMode::A ? 1 : 2);
        });

    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(*topology.get<GeneticMode::A>() == 1);
    REQUIRE(*topology.get<GeneticMode::D>() == 2);
}

TEST_CASE(
    "transform_mode_values returns owning values without changing its source",
    "[bayes][genetic][independent_topology]")
{
    const auto topology = AdditiveDominanceTopology{TestValue{1}, TestValue{2}};
    std::vector<GeneticMode> visited_modes;

    const auto transformed = transform_mode_values(
        topology,
        [&](GeneticMode mode, const TestValue& value) -> const int&
        {
            visited_modes.push_back(mode);
            return value.value;
        });

    static_assert(std::same_as<
                  std::remove_cvref_t<decltype(transformed)>::value_type,
                  int>);
    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(transformed.get<GeneticMode::A>() == 1);
    REQUIRE(transformed.get<GeneticMode::D>() == 2);
    REQUIRE(topology.get<GeneticMode::A>().value == 1);
    REQUIRE(topology.get<GeneticMode::D>().value == 2);
}

TEST_CASE(
    "mode value algorithms preserve a singleton dominance mode",
    "[bayes][genetic][independent_topology]")
{
    constexpr auto mode_d = GeneticModeSet{GeneticMode::D};
    std::vector<GeneticMode> generated_modes;
    std::vector<GeneticMode> transformed_modes;

    const auto generated = generate_mode_values<mode_d>(
        [&](GeneticMode mode)
        {
            generated_modes.push_back(mode);
            return TestValue{2};
        });
    const auto transformed = transform_mode_values(
        generated,
        [&](GeneticMode mode, const TestValue& value)
        {
            transformed_modes.push_back(mode);
            return value.value * 2;
        });

    REQUIRE(generated_modes == std::vector{GeneticMode::D});
    REQUIRE(transformed_modes == std::vector{GeneticMode::D});
    REQUIRE(generated.get<GeneticMode::D>().value == 2);
    REQUIRE(transformed.get<GeneticMode::D>() == 4);
}
