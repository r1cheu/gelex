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
#include <initializer_list>
#include <span>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/model/bayes/legacy_algorithm_shape.h"
#include "gelex/model/bayes/legacy_bayes_policy.h"
#include "gelex/types/genetic_effect_type.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::bayes::AlgorithmShape;
using gelex::bayes::BayesPolicy;
using gelex::bayes::resolve_shape;
using gelex::bayes::to_file_suffix;
using gelex::bayes::to_heritability_label;
using gelex::bayes::to_variance_label;

namespace
{

inline constexpr std::array kIndepShapes{
    AlgorithmShape::a_only,
    AlgorithmShape::d_only,
    AlgorithmShape::ad_independent,
};

inline constexpr std::array kJointShapes{
    AlgorithmShape::ad_joint,
};

auto resolve(
    std::span<const AlgorithmShape> policy_shapes,
    std::initializer_list<GeneticMode> effects) -> AlgorithmShape
{
    const BayesPolicy policy{.shapes = policy_shapes};
    const std::vector<GeneticMode> effect_vec{effects};
    return resolve_shape(policy, effect_vec);
}

}  // namespace

TEST_CASE("resolve_shape returns expected shape", "[algorithm_shape]")
{
    using enum AlgorithmShape;
    using enum GeneticMode;

    SECTION("{A} indep -> a_only")
    {
        REQUIRE(resolve(kIndepShapes, {A}) == a_only);
    }
    SECTION("{D} indep -> d_only")
    {
        REQUIRE(resolve(kIndepShapes, {D}) == d_only);
    }
    SECTION("{A,D} indep -> ad_independent")
    {
        REQUIRE(resolve(kIndepShapes, {A, D}) == ad_independent);
    }
    SECTION("{D,A} order does not matter")
    {
        REQUIRE(resolve(kIndepShapes, {D, A}) == ad_independent);
    }
    SECTION("{A,D} joint -> ad_joint")
    {
        REQUIRE(resolve(kJointShapes, {A, D}) == ad_joint);
    }
}

TEST_CASE("resolve_shape throws on unsupported request", "[algorithm_shape]")
{
    using enum GeneticMode;

    SECTION("{A} with joint-only policy")
    {
        REQUIRE_THROWS_AS(resolve(kJointShapes, {A}), GelexException);
    }
    SECTION("{D} with joint-only policy")
    {
        REQUIRE_THROWS_AS(resolve(kJointShapes, {D}), GelexException);
    }
    SECTION("{A} with empty policy")
    {
        REQUIRE_THROWS_AS(resolve({}, {A}), GelexException);
    }
    SECTION("{D} with empty policy")
    {
        REQUIRE_THROWS_AS(resolve({}, {D}), GelexException);
    }
    SECTION("{A,D} with empty policy")
    {
        REQUIRE_THROWS_AS(resolve({}, {A, D}), GelexException);
    }
}

TEST_CASE("label helpers map each shape to fixed strings", "[algorithm_shape]")
{
    using enum AlgorithmShape;

    SECTION("a_only")
    {
        REQUIRE(to_variance_label(a_only) == "σ²_add");
        REQUIRE(to_heritability_label(a_only) == "h²");
        REQUIRE(to_file_suffix(a_only) == "add");
    }
    SECTION("d_only")
    {
        REQUIRE(to_variance_label(d_only) == "σ²_dom");
        REQUIRE(to_heritability_label(d_only) == "δ²");
        REQUIRE(to_file_suffix(d_only) == "dom");
    }
    SECTION("ad_independent")
    {
        REQUIRE(to_variance_label(ad_independent) == "σ²_g");
        REQUIRE(to_heritability_label(ad_independent) == "H²");
        REQUIRE(to_file_suffix(ad_independent) == "ad");
    }
    SECTION("ad_joint")
    {
        REQUIRE(to_variance_label(ad_joint) == "σ²_g");
        REQUIRE(to_heritability_label(ad_joint) == "H²");
        REQUIRE(to_file_suffix(ad_joint) == "ad");
    }
}
