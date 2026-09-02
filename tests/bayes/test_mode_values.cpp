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
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/mode_values.h"

using gelex::generate_mode_values;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::HomogeneousModeValues;
using gelex::JointModeValues;
using gelex::ModeValues;
using gelex::transform_mode_values;

namespace
{

struct TestValue
{
    int value;
};

struct NonDefaultValue
{
    NonDefaultValue() = delete;
    explicit NonDefaultValue(int value) : value{value} {}

    int value;
};

constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using AdditiveDominanceValues = ModeValues<mode_ad, TestValue, TestValue>;
using DominanceAdditiveValues
    = ModeValues<GeneticMode::D | GeneticMode::A, TestValue, TestValue>;
using HeterogeneousValues = ModeValues<mode_ad, int, std::string>;
using HomogeneousValues = HomogeneousModeValues<mode_ad, TestValue>;
using HomogeneousNonDefaultValues
    = HomogeneousModeValues<mode_ad, NonDefaultValue>;
using TestJointModeValues
    = JointModeValues<AdditiveDominanceValues, std::string>;

template <typename ModeValuesType, typename JointValue>
concept ValidJointModeValues
    = requires { typename JointModeValues<ModeValuesType, JointValue>; };

static_assert(
    AdditiveDominanceValues::modes == (GeneticMode::A | GeneticMode::D));
static_assert(std::same_as<AdditiveDominanceValues, DominanceAdditiveValues>);
static_assert(std::same_as<HomogeneousValues, AdditiveDominanceValues>);
static_assert(std::same_as<
              HomogeneousNonDefaultValues,
              ModeValues<mode_ad, NonDefaultValue, NonDefaultValue>>);
static_assert(std::same_as<
              AdditiveDominanceValues::mode_value_type<GeneticMode::A>,
              TestValue>);
static_assert(std::same_as<
              AdditiveDominanceValues::mode_value_type<GeneticMode::D>,
              TestValue>);
static_assert(
    std::constructible_from<AdditiveDominanceValues, TestValue, TestValue>);
static_assert(
    std::same_as<HeterogeneousValues::mode_value_type<GeneticMode::A>, int>);
static_assert(std::same_as<
              HeterogeneousValues::mode_value_type<GeneticMode::D>,
              std::string>);
static_assert(
    std::same_as<
        decltype(std::declval<HeterogeneousValues&>().get<GeneticMode::A>()),
        int&>);
static_assert(std::same_as<
              decltype(std::declval<const HeterogeneousValues&>()
                           .get<GeneticMode::D>()),
              const std::string&>);

static_assert(
    AdditiveDominanceValues{TestValue{1}, TestValue{2}}
        .get<GeneticMode::A>()
        .value
    == 1);
static_assert(
    AdditiveDominanceValues{TestValue{1}, TestValue{2}}
        .get<GeneticMode::D>()
        .value
    == 2);

constexpr auto generated_values = generate_mode_values<mode_ad>(
    []<GeneticMode Mode>()
    { return TestValue{Mode == GeneticMode::A ? 1 : 2}; });

static_assert(std::same_as<
              std::remove_cvref_t<decltype(generated_values)>,
              AdditiveDominanceValues>);
static_assert(generated_values.get<GeneticMode::A>().value == 1);
static_assert(generated_values.get<GeneticMode::D>().value == 2);

constexpr auto transformed_values = transform_mode_values(
    generated_values,
    []<GeneticMode Mode>(const TestValue& value)
    { return value.value + (Mode == GeneticMode::A ? 10 : 20); });

static_assert(std::same_as<
              std::remove_cvref_t<decltype(transformed_values)>,
              ModeValues<mode_ad, int, int>>);
static_assert(transformed_values.get<GeneticMode::A>() == 11);
static_assert(transformed_values.get<GeneticMode::D>() == 22);

constexpr auto heterogeneous_generated = generate_mode_values<mode_ad>(
    []<GeneticMode Mode>()
    {
        if constexpr (Mode == GeneticMode::A)
        {
            return 1;
        }
        else
        {
            return std::string_view{"dominance"};
        }
    });

static_assert(std::same_as<
              std::remove_cvref_t<decltype(heterogeneous_generated)>,
              ModeValues<mode_ad, int, std::string_view>>);

constexpr auto heterogeneous_transformed = transform_mode_values(
    heterogeneous_generated,
    []<GeneticMode Mode>(const auto& value)
    {
        if constexpr (Mode == GeneticMode::A)
        {
            return static_cast<double>(value);
        }
        else
        {
            return value.size();
        }
    });

static_assert(std::same_as<
              std::remove_cvref_t<decltype(heterogeneous_transformed)>,
              ModeValues<mode_ad, double, std::size_t>>);
static_assert(heterogeneous_transformed.get<GeneticMode::A>() == 1.0);
static_assert(heterogeneous_transformed.get<GeneticMode::D>() == 9);

static_assert(TestJointModeValues::modes == mode_ad);
static_assert(std::same_as<
              TestJointModeValues::mode_value_type<GeneticMode::A>,
              TestValue>);
static_assert(std::same_as<
              TestJointModeValues::mode_value_type<GeneticMode::D>,
              TestValue>);
static_assert(std::same_as<TestJointModeValues::joint_value_type, std::string>);
static_assert(std::same_as<
              TestJointModeValues::mode_values_type,
              AdditiveDominanceValues>);
static_assert(!ValidJointModeValues<TestJointModeValues, int>);
static_assert(!ValidJointModeValues<AdditiveDominanceValues, const int&>);
static_assert(std::constructible_from<
              TestJointModeValues,
              AdditiveDominanceValues,
              std::string>);
static_assert(
    !std::constructible_from<TestJointModeValues, AdditiveDominanceValues>);
static_assert(std::default_initializable<TestJointModeValues>);
static_assert(std::same_as<
              decltype(JointModeValues{
                  AdditiveDominanceValues{TestValue{1}, TestValue{2}},
                  3}),
              JointModeValues<AdditiveDominanceValues, int>>);
static_assert(std::same_as<
              decltype(JointModeValues{ModeValues<mode_ad, int, int>{1, 2}, 3}),
              JointModeValues<ModeValues<mode_ad, int, int>, int>>);
static_assert(
    JointModeValues{AdditiveDominanceValues{TestValue{1}, TestValue{2}}, 3}
        .mode_values()
        .get<GeneticMode::A>()
        .value
    == 1);
static_assert(
    JointModeValues{AdditiveDominanceValues{TestValue{1}, TestValue{2}}, 3}
        .mode_values()
        .get<GeneticMode::D>()
        .value
    == 2);
static_assert(
    JointModeValues{AdditiveDominanceValues{TestValue{1}, TestValue{2}}, 3}
        .joint()
    == 3);
static_assert(
    JointModeValues{ModeValues<mode_ad, int, int>{1, 2}, 3}
        .mode_values()
        .get<GeneticMode::D>()
    == 2);
static_assert(
    JointModeValues{ModeValues<mode_ad, int, int>{1, 2}, 3}.joint() == 3);

}  // namespace

TEST_CASE(
    "ModeValues stores values in canonical mode order",
    "[bayes][mode_values]")
{
    auto values = ModeValues<
        GeneticMode::D | GeneticMode::A,
        std::unique_ptr<int>,
        std::string>{std::make_unique<int>(1), "dominance"};

    std::vector<GeneticMode> visited_modes;
    std::vector<std::string> visited_values;
    values.for_each(
        [&]<GeneticMode Mode>(const auto& value)
        {
            visited_modes.push_back(Mode);
            if constexpr (Mode == GeneticMode::A)
            {
                visited_values.push_back(std::to_string(*value));
            }
            else
            {
                visited_values.push_back(value);
            }
        });

    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(visited_values == std::vector<std::string>{"1", "dominance"});
}

TEST_CASE(
    "ModeValues provides mode-indexed value access",
    "[bayes][mode_values]")
{
    auto values
        = ModeValues<GeneticModeSet{GeneticMode::D}, TestValue>{TestValue{2}};

    values.get<GeneticMode::D>().value = 3;
    const auto& constant_values = values;

    REQUIRE(constant_values.get<GeneticMode::D>().value == 3);
}

TEST_CASE(
    "ModeValues maps every mode to its constructor argument",
    "[bayes][mode_values]")
{
    const auto values = AdditiveDominanceValues{TestValue{1}, TestValue{2}};

    REQUIRE(values.get<GeneticMode::A>().value == 1);
    REQUIRE(values.get<GeneticMode::D>().value == 2);

    std::vector<int> visited_values;
    values.for_each([&]<GeneticMode /*Mode*/>(const TestValue& value)
                    { visited_values.push_back(value.value); });

    REQUIRE(visited_values == std::vector{1, 2});
}

TEST_CASE(
    "generate_mode_values invokes the factory once in canonical order",
    "[bayes][mode_values]")
{
    std::vector<GeneticMode> visited_modes;

    const auto values = generate_mode_values<mode_ad>(
        [&]<GeneticMode Mode>()
        {
            visited_modes.push_back(Mode);
            return std::make_unique<int>(Mode == GeneticMode::A ? 1 : 2);
        });

    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(*values.get<GeneticMode::A>() == 1);
    REQUIRE(*values.get<GeneticMode::D>() == 2);
}

TEST_CASE(
    "transform_mode_values returns owning values without changing its source",
    "[bayes][mode_values]")
{
    const auto values = AdditiveDominanceValues{TestValue{1}, TestValue{2}};
    std::vector<GeneticMode> visited_modes;

    const auto transformed = transform_mode_values(
        values,
        [&]<GeneticMode Mode>(const TestValue& value) -> const int&
        {
            visited_modes.push_back(Mode);
            return value.value;
        });

    static_assert(std::same_as<
                  std::remove_cvref_t<decltype(transformed)>::mode_value_type<
                      GeneticMode::A>,
                  int>);
    static_assert(std::same_as<
                  std::remove_cvref_t<decltype(transformed)>::mode_value_type<
                      GeneticMode::D>,
                  int>);
    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(transformed.get<GeneticMode::A>() == 1);
    REQUIRE(transformed.get<GeneticMode::D>() == 2);
    REQUIRE(values.get<GeneticMode::A>().value == 1);
    REQUIRE(values.get<GeneticMode::D>().value == 2);
}

TEST_CASE(
    "mode value algorithms preserve a singleton dominance mode",
    "[bayes][mode_values]")
{
    constexpr auto mode_d = GeneticModeSet{GeneticMode::D};
    std::vector<GeneticMode> generated_modes;
    std::vector<GeneticMode> transformed_modes;

    const auto generated = generate_mode_values<mode_d>(
        [&]<GeneticMode Mode>()
        {
            generated_modes.push_back(Mode);
            return TestValue{2};
        });
    const auto transformed = transform_mode_values(
        generated,
        [&]<GeneticMode Mode>(const TestValue& value)
        {
            transformed_modes.push_back(Mode);
            return value.value * 2;
        });

    REQUIRE(generated_modes == std::vector{GeneticMode::D});
    REQUIRE(transformed_modes == std::vector{GeneticMode::D});
    REQUIRE(generated.get<GeneticMode::D>().value == 2);
    REQUIRE(transformed.get<GeneticMode::D>() == 4);
}

TEST_CASE(
    "JointModeValues stores mode values in canonical order",
    "[bayes][mode_values]")
{
    auto values = JointModeValues{
        ModeValues<mode_ad, std::unique_ptr<int>, std::unique_ptr<int>>{
            std::make_unique<int>(1), std::make_unique<int>(2)},
        std::string{"joint"}};

    std::vector<GeneticMode> visited_modes;
    std::vector<int> visited_values;
    values.mode_values().for_each(
        [&]<GeneticMode Mode>(const auto& value)
        {
            visited_modes.push_back(Mode);
            visited_values.push_back(*value);
        });

    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(visited_values == std::vector{1, 2});
    REQUIRE(values.joint() == "joint");
}

TEST_CASE(
    "JointModeValues provides separate mode and joint access",
    "[bayes][mode_values]")
{
    auto values = JointModeValues{
        ModeValues<mode_ad, int, int>{1, 2}, std::string{"joint"}};

    values.mode_values().get<GeneticMode::A>() = 3;
    values.mode_values().get<GeneticMode::D>() = 4;
    values.joint() = "updated";

    const auto& constant_values = values;
    REQUIRE(constant_values.mode_values().get<GeneticMode::A>() == 3);
    REQUIRE(constant_values.mode_values().get<GeneticMode::D>() == 4);
    REQUIRE(constant_values.joint() == "updated");
}

TEST_CASE(
    "JointModeValues traversal of a const container yields const mode values",
    "[bayes][mode_values]")
{
    const auto values = JointModeValues{
        AdditiveDominanceValues{TestValue{1}, TestValue{2}},
        std::string{"joint"}};

    std::vector<int> visited_values;
    values.mode_values().for_each(
        [&]<GeneticMode /*Mode*/>(const auto& value)
        {
            static_assert(
                std::is_const_v<std::remove_reference_t<decltype(value)>>);
            visited_values.push_back(value.value);
        });

    REQUIRE(visited_values == std::vector{1, 2});
}
