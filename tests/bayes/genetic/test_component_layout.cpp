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

#include "gelex/bayes/genetic/component_layout.h"
#include "gelex/types/genetic_mode.h"

using gelex::ComponentLayout;
using gelex::GeneticMode;
using gelex::JointZeroInflatedComponentLayout;
using gelex::SingleComponentLayout;
using gelex::ZeroInflatedComponentLayout;

static_assert(ComponentLayout<SingleComponentLayout>);
static_assert(ComponentLayout<ZeroInflatedComponentLayout<5>>);
static_assert(ComponentLayout<JointZeroInflatedComponentLayout>);

TEST_CASE(
    "single and zero-inflated component layouts map active classes",
    "[bayes][genetic][component_layout]")
{
    REQUIRE(SingleComponentLayout::component_count == 1);
    REQUIRE(SingleComponentLayout::component_index(GeneticMode::A, 0) == 0);
    REQUIRE(
        SingleComponentLayout::component_index(GeneticMode::D, 1)
        == SingleComponentLayout::no_component);

    using Layout = ZeroInflatedComponentLayout<5>;
    REQUIRE(Layout::component_count == 4);
    REQUIRE(Layout::component_index(GeneticMode::A, 0) == Layout::no_component);
    REQUIRE(Layout::component_index(GeneticMode::D, 1) == 0);
    REQUIRE(Layout::component_index(GeneticMode::D, 4) == 3);
}

TEST_CASE(
    "joint component layout maps classes independently for each mode",
    "[bayes][genetic][component_layout]")
{
    using Layout = JointZeroInflatedComponentLayout;

    REQUIRE(Layout::component_count == 2);
    REQUIRE(Layout::component_index(GeneticMode::A, 0) == Layout::no_component);
    REQUIRE(Layout::component_index(GeneticMode::A, 1) == 0);
    REQUIRE(Layout::component_index(GeneticMode::A, 2) == Layout::no_component);
    REQUIRE(Layout::component_index(GeneticMode::A, 3) == 1);

    REQUIRE(Layout::component_index(GeneticMode::D, 0) == Layout::no_component);
    REQUIRE(Layout::component_index(GeneticMode::D, 1) == Layout::no_component);
    REQUIRE(Layout::component_index(GeneticMode::D, 2) == 0);
    REQUIRE(Layout::component_index(GeneticMode::D, 3) == 1);
}
