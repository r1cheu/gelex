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
#include <catch2/matchers/catch_matchers_string.hpp>
#include <limits>
#include <string>

#include "gelex/bayes/detail/recipe_validation.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

using Catch::Matchers::ContainsSubstring;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::IndependentTopology;
using gelex::JointSpikeSlab;
using gelex::NoParameters;
using gelex::ScaledMixture;
using gelex::SpikeSlab;
using gelex::VarianceBudget;
using gelex::detail::validate_recipe_inputs;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();

using SpikeSlabA = IndependentTopology<mode_a, SpikeSlab>;
using SpikeSlabAD = IndependentTopology<mode_ad, SpikeSlab>;
using ScaledMixtureA = IndependentTopology<mode_a, ScaledMixture>;
using ScaledMixtureAD = IndependentTopology<mode_ad, ScaledMixture>;

const auto BUDGET_A = VarianceBudget{{.additive = 0.5}};
const auto BUDGET_AD = VarianceBudget{{.additive = 0.5, .dominance = 0.2}};

// Reports the message a failing validation produces, so that each test asserts
// on what the user actually reads.
auto message_of(auto&& validate) -> std::string
{
    try
    {
        validate();
    }
    catch (const GelexException& error)
    {
        return error.what();
    }
    return {};
}

}  // namespace

TEST_CASE(
    "validate_recipe_inputs stays silent on valid input",
    "[bayes][recipe_validation]")
{
    REQUIRE_NOTHROW(validate_recipe_inputs(NoParameters{}, BUDGET_A, mode_a));
    REQUIRE_NOTHROW(validate_recipe_inputs(
        ScaledMixtureAD{ScaledMixture{}, ScaledMixture{}}, BUDGET_AD, mode_ad));
    REQUIRE_NOTHROW(
        validate_recipe_inputs(JointSpikeSlab{}, BUDGET_AD, mode_ad));
}

TEST_CASE(
    "validate_recipe_inputs reports every issue at once",
    "[bayes][recipe_validation]")
{
    const auto message = message_of(
        []
        {
            validate_recipe_inputs(
                SpikeSlabAD{
                    SpikeSlab{.probability = 1.5},
                    SpikeSlab{.probability = 0.0},
                },
                VarianceBudget{{.additive = 0.7, .dominance = 0.4}},
                mode_ad);
        });

    REQUIRE_THAT(
        message, ContainsSubstring("A inclusion probability must lie"));
    REQUIRE_THAT(
        message, ContainsSubstring("D inclusion probability must lie"));
    REQUIRE_THAT(
        message, ContainsSubstring("variance shares must sum to less than 1"));
}

TEST_CASE(
    "validate_recipe_inputs labels an issue with its own mode",
    "[bayes][recipe_validation]")
{
    const auto message = message_of(
        []
        {
            validate_recipe_inputs(
                SpikeSlabAD{
                    SpikeSlab{},
                    SpikeSlab{.probability = 1.5},
                },
                BUDGET_AD,
                mode_ad);
        });

    REQUIRE_THAT(
        message, ContainsSubstring("D inclusion probability must lie"));
    REQUIRE_THAT(message, !ContainsSubstring("A inclusion probability"));
}

TEST_CASE(
    "validate_recipe_inputs rejects a closed unit interval",
    "[bayes][recipe_validation]")
{
    auto probability = double{};

    SECTION("zero is out")
    {
        probability = 0.0;
    }

    SECTION("one is out")
    {
        probability = 1.0;
    }

    SECTION("a non-number is out")
    {
        probability = not_a_number;
    }

    const auto message = message_of(
        [probability]
        {
            validate_recipe_inputs(
                SpikeSlabA{SpikeSlab{.probability = probability}},
                BUDGET_A,
                mode_a);
        });

    REQUIRE_THAT(
        message,
        ContainsSubstring(
            "A inclusion probability must lie in the open interval (0, 1)"));
}

TEST_CASE(
    "validate_recipe_inputs rejects non-simplex weights",
    "[bayes][recipe_validation]")
{
    auto weights = ScaledMixture{}.probabilities;

    SECTION("components must be strictly positive")
    {
        weights = {1.0, 0.0, 0.0, 0.0, 0.0};
    }

    SECTION("components must be finite")
    {
        weights = {0.5, not_a_number, 0.2, 0.2, 0.1};
    }

    const auto message = message_of(
        [weights]
        {
            validate_recipe_inputs(
                ScaledMixtureA{ScaledMixture{.probabilities = weights}},
                BUDGET_A,
                mode_a);
        });

    REQUIRE_THAT(
        message,
        ContainsSubstring("A mixture weights[1] must be a finite positive"));
}

TEST_CASE(
    "validate_recipe_inputs holds the simplex sum to a tolerance",
    "[bayes][recipe_validation]")
{
    const auto validate = [](double last)
    {
        return [last]
        {
            validate_recipe_inputs(
                ScaledMixtureA{
                    ScaledMixture{.probabilities = {0.2, 0.2, 0.2, 0.2, last}}},
                BUDGET_A,
                mode_a);
        };
    };

    REQUIRE_NOTHROW(validate(0.2 - 1e-12)());
    REQUIRE_THAT(
        message_of(validate(0.2 - 1e-6)),
        ContainsSubstring("A mixture weights must sum to 1"));
}

TEST_CASE(
    "validate_recipe_inputs rejects malformed multipliers",
    "[bayes][recipe_validation]")
{
    SECTION("the null scale must come first")
    {
        const auto message = message_of(
            []
            {
                validate_recipe_inputs(
                    ScaledMixtureA{
                        ScaledMixture{.scales = {0.001, 0.01, 0.1, 1.0, 10.0}}},
                    BUDGET_A,
                    mode_a);
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring("A variance multipliers[0] must be zero"));
    }

    SECTION("the remaining scales must be finite and positive")
    {
        const auto message = message_of(
            []
            {
                validate_recipe_inputs(
                    ScaledMixtureAD{
                        ScaledMixture{},
                        ScaledMixture{.scales = {0.0, -0.01, 0.1, 1.0, 10.0}}},
                    BUDGET_AD,
                    mode_ad);
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring(
                "D variance multipliers[1] must be a finite positive "
                "multiplier"));
    }
}

TEST_CASE(
    "validate_recipe_inputs checks the joint allocation",
    "[bayes][recipe_validation]")
{
    const auto message = message_of(
        []
        {
            validate_recipe_inputs(
                JointSpikeSlab{
                    .probabilities = {0.9, 0.05, 0.02, 0.02},
                    .positive_probability = not_a_number,
                },
                BUDGET_AD,
                mode_ad);
        });

    REQUIRE_THAT(
        message, ContainsSubstring("AD mixture weights must sum to 1"));
    REQUIRE_THAT(
        message,
        ContainsSubstring("AD dominance positive-sign probability must lie"));
}

TEST_CASE(
    "validate_recipe_inputs cross-checks the budget against modes",
    "[bayes][recipe_validation]")
{
    SECTION("a present mode needs a positive share")
    {
        const auto message = message_of(
            []
            {
                validate_recipe_inputs(
                    NoParameters{}, VarianceBudget{{.additive = 0.5}}, mode_ad);
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring("D variance share must be finite and positive"));
    }

    SECTION("an absent mode needs a zero share")
    {
        const auto message = message_of(
            [] { validate_recipe_inputs(NoParameters{}, BUDGET_AD, mode_a); });

        REQUIRE_THAT(
            message, ContainsSubstring("D variance share must be zero"));
    }

    SECTION("shares must leave room for the residual")
    {
        const auto message = message_of(
            []
            {
                validate_recipe_inputs(
                    NoParameters{},
                    VarianceBudget{{.additive = 0.6, .dominance = 0.4}},
                    mode_ad);
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring("variance shares must sum to less than 1"));
    }

    SECTION("the random share must be finite and non-negative")
    {
        const auto message = message_of(
            []
            {
                validate_recipe_inputs(
                    NoParameters{},
                    VarianceBudget{{.additive = 0.5, .random = -0.1}},
                    mode_a);
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring(
                "random variance share must be finite and non-negative"));
    }

    SECTION("an unlabelled issue carries no mode prefix")
    {
        const auto message = message_of(
            []
            {
                validate_recipe_inputs(
                    NoParameters{},
                    VarianceBudget{{.additive = 0.5, .random = -0.1}},
                    mode_a);
            });

        REQUIRE_THAT(message, !ContainsSubstring("A random variance share"));
    }
}
