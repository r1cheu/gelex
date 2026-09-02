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
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"

using gelex::GeneticMode;
using gelex::IndependentTopology;
using gelex::JointTopology;

namespace
{

struct TestValue
{
    int value;
};

constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using TestModeTopology = IndependentTopology<mode_ad, TestValue, TestValue>;
using TestJointTopology = JointTopology<TestModeTopology, std::string>;

template <typename ModeTopology, typename JointValue>
concept ValidJointTopology
    = requires { typename JointTopology<ModeTopology, JointValue>; };

static_assert(TestJointTopology::modes == (GeneticMode::A | GeneticMode::D));
static_assert(std::same_as<
              TestJointTopology::mode_value_type<GeneticMode::A>,
              TestValue>);
static_assert(std::same_as<
              TestJointTopology::mode_value_type<GeneticMode::D>,
              TestValue>);
static_assert(std::same_as<TestJointTopology::joint_value_type, std::string>);
static_assert(
    std::same_as<TestJointTopology::mode_topology_type, TestModeTopology>);
static_assert(!ValidJointTopology<TestJointTopology, int>);
static_assert(!ValidJointTopology<TestModeTopology, const int&>);

static_assert(
    std::constructible_from<TestJointTopology, TestModeTopology, std::string>);

// The value constructor still requires the joint payload, while the topology
// composes the defaults provided by its leaves.
static_assert(!std::constructible_from<TestJointTopology, TestModeTopology>);
static_assert(std::default_initializable<TestJointTopology>);

// Mode and joint types deduce independently, including when they coincide.
static_assert(std::same_as<
              decltype(JointTopology{
                  TestModeTopology{TestValue{1}, TestValue{2}},
                  3}),
              JointTopology<TestModeTopology, int>>);
static_assert(std::same_as<
              decltype(JointTopology{
                  IndependentTopology<mode_ad, int, int>{1, 2},
                  3}),
              JointTopology<IndependentTopology<mode_ad, int, int>, int>>);

static_assert(
    JointTopology{TestModeTopology{TestValue{1}, TestValue{2}}, 3}
        .mode_values()
        .get<GeneticMode::A>()
        .value
    == 1);
static_assert(
    JointTopology{TestModeTopology{TestValue{1}, TestValue{2}}, 3}
        .mode_values()
        .get<GeneticMode::D>()
        .value
    == 2);
static_assert(
    JointTopology{TestModeTopology{TestValue{1}, TestValue{2}}, 3}.joint()
    == 3);

// A coinciding joint type must not shift into the mode storage.
static_assert(
    JointTopology{IndependentTopology<mode_ad, int, int>{1, 2}, 3}
        .mode_values()
        .get<GeneticMode::D>()
    == 2);
static_assert(
    JointTopology{IndependentTopology<mode_ad, int, int>{1, 2}, 3}.joint()
    == 3);

}  // namespace

TEST_CASE(
    "JointTopology stores mode values in canonical order",
    "[bayes][genetic][joint_topology]")
{
    auto topology = JointTopology{
        IndependentTopology<
            mode_ad,
            std::unique_ptr<int>,
            std::unique_ptr<int>>{
            std::make_unique<int>(1), std::make_unique<int>(2)},
        std::string{"joint"}};

    std::vector<GeneticMode> visited_modes;
    std::vector<int> visited_values;
    topology.mode_values().for_each(
        [&]<GeneticMode Mode>(const auto& value)
        {
            visited_modes.push_back(Mode);
            visited_values.push_back(*value);
        });

    REQUIRE(visited_modes == std::vector{GeneticMode::A, GeneticMode::D});
    REQUIRE(visited_values == std::vector{1, 2});
    REQUIRE(topology.joint() == "joint");
}

TEST_CASE(
    "JointTopology provides separate mode and joint access",
    "[bayes][genetic][joint_topology]")
{
    auto topology = JointTopology{
        IndependentTopology<mode_ad, int, int>{1, 2}, std::string{"joint"}};

    topology.mode_values().get<GeneticMode::A>() = 3;
    topology.mode_values().get<GeneticMode::D>() = 4;
    topology.joint() = "updated";

    const auto& constant_topology = topology;
    REQUIRE(constant_topology.mode_values().get<GeneticMode::A>() == 3);
    REQUIRE(constant_topology.mode_values().get<GeneticMode::D>() == 4);
    REQUIRE(constant_topology.joint() == "updated");
}

TEST_CASE(
    "JointTopology traversal of a const topology yields const mode values",
    "[bayes][genetic][joint_topology]")
{
    const auto topology = JointTopology{
        TestModeTopology{TestValue{1}, TestValue{2}}, std::string{"joint"}};

    std::vector<int> visited_values;
    topology.mode_values().for_each(
        [&]<GeneticMode /*Mode*/>(const auto& value)
        {
            static_assert(
                std::is_const_v<std::remove_reference_t<decltype(value)>>);
            visited_values.push_back(value.value);
        });

    REQUIRE(visited_values == std::vector{1, 2});
}
