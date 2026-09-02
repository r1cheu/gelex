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
#include <catch2/matchers/catch_matchers_string.hpp>
#include <concepts>
#include <limits>
#include <string>

#include "gelex/bayes/spec.h"
#include "gelex/exception.h"

using Catch::Matchers::ContainsSubstring;
using gelex::GelexException;
using gelex::HalfNormal;
using gelex::JointSpikeSlab;
using gelex::ScaledMixture;
using gelex::SpikeSlab;

namespace
{

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();

static_assert(ScaledMixture::class_count == 5);
static_assert(JointSpikeSlab::class_count == 4);
static_assert(std::default_initializable<HalfNormal>);

auto message_of(auto&& construct) -> std::string
{
    try
    {
        construct();
    }
    catch (const GelexException& error)
    {
        return error.what();
    }
    return {};
}

}  // namespace

TEST_CASE("Bayes structural specs provide defaults", "[bayes][spec]")
{
    const auto spike_slab = SpikeSlab{};
    REQUIRE(spike_slab.probability() == 0.01);

    const auto scaled_mixture = ScaledMixture{};
    REQUIRE(
        scaled_mixture.probabilities()
        == std::array{0.99, 0.005, 0.003, 0.001, 0.001});
    REQUIRE(scaled_mixture.scales() == std::array{0.0, 0.001, 0.01, 0.1, 1.0});

    const auto joint_spike_slab = JointSpikeSlab{};
    REQUIRE(
        joint_spike_slab.probabilities()
        == std::array{0.99, 1.0 / 300, 1.0 / 300, 1.0 / 300});
}

TEST_CASE("Bayes structural specs accept resolved values", "[bayes][spec]")
{
    const auto spike_slab = SpikeSlab{0.2};
    REQUIRE(spike_slab.probability() == 0.2);

    const auto scaled_mixture = ScaledMixture{
        {0.8, 0.05, 0.05, 0.05, 0.05}, {0.0, 0.01, 0.1, 1.0, 10.0}};
    REQUIRE(
        scaled_mixture.probabilities()
        == std::array{0.8, 0.05, 0.05, 0.05, 0.05});
    REQUIRE(scaled_mixture.scales() == std::array{0.0, 0.01, 0.1, 1.0, 10.0});

    const auto joint_spike_slab = JointSpikeSlab{{0.7, 0.1, 0.1, 0.1}};
    REQUIRE(joint_spike_slab.probabilities() == std::array{0.7, 0.1, 0.1, 0.1});
}

TEST_CASE("SpikeSlab rejects invalid probabilities", "[bayes][spec]")
{
    auto probability = 0.0;

    SECTION("zero")
    {
        probability = 0.0;
    }
    SECTION("one")
    {
        probability = 1.0;
    }
    SECTION("non-finite")
    {
        probability = not_a_number;
    }

    REQUIRE_THAT(
        message_of([probability] { return SpikeSlab{probability}; }),
        ContainsSubstring("must lie in the open interval (0, 1)"));
}

TEST_CASE("ScaledMixture rejects invalid probabilities", "[bayes][spec]")
{
    SECTION("a probability is not positive")
    {
        REQUIRE_THAT(
            message_of([] { return ScaledMixture{{1.0, 0.0, 0.0, 0.0, 0.0}}; }),
            ContainsSubstring("probabilities[1] must be finite and positive"));
    }

    SECTION("probabilities do not sum to one")
    {
        REQUIRE_THAT(
            message_of([] { return ScaledMixture{{0.2, 0.2, 0.2, 0.2, 0.1}}; }),
            ContainsSubstring("probabilities must sum to 1"));
    }
}

TEST_CASE("ScaledMixture rejects invalid scales", "[bayes][spec]")
{
    const auto defaults = ScaledMixture{};
    const auto probabilities = defaults.probabilities();

    SECTION("the null scale is not first")
    {
        REQUIRE_THAT(
            message_of(
                [probabilities]
                {
                    return ScaledMixture{
                        probabilities, {0.01, 0.1, 1.0, 10.0, 100.0}};
                }),
            ContainsSubstring("scales[0] must be zero"));
    }

    SECTION("an active scale is not positive")
    {
        REQUIRE_THAT(
            message_of(
                [probabilities]
                {
                    return ScaledMixture{
                        probabilities, {0.0, -0.01, 0.1, 1.0, 10.0}};
                }),
            ContainsSubstring("scales[1] must be finite and positive"));
    }
}

TEST_CASE(
    "Bayes structural specs report the first invalid field",
    "[bayes][spec]")
{
    SECTION("scaled mixture probabilities precede scales")
    {
        const auto message = message_of(
            []
            {
                return ScaledMixture{
                    {0.9, 0.05, 0.02, 0.02, 0.0}, {1.0, 0.1, 1.0, 10.0, 100.0}};
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring("probabilities[4] must be finite and positive"));
        REQUIRE_THAT(message, !ContainsSubstring("scales[0]"));
    }
}

TEST_CASE("JointSpikeSlab rejects invalid inputs", "[bayes][spec]")
{
    SECTION("allocation probabilities do not form a simplex")
    {
        REQUIRE_THAT(
            message_of([] { return JointSpikeSlab{{0.9, 0.05, 0.02, 0.02}}; }),
            ContainsSubstring("probabilities must sum to 1"));
    }
}
