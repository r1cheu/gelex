// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <concepts>
#include <limits>
#include <string>

#include "gelex/bayes/spec.h"
#include "gelex/exception.h"

using Catch::Matchers::ContainsSubstring;
using gelex::GelexException;
using gelex::HalfNormalSpec;
using gelex::JointSpikeSlabSpec;
using gelex::MixtureWeightUpdate;
using gelex::ScaledMixtureSpec;
using gelex::SpikeSlabSpec;
using gelex::VarianceLayout;

namespace
{

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();

static_assert(ScaledMixtureSpec<>::class_count == 5);
static_assert(JointSpikeSlabSpec<>::class_count == 4);
static_assert(std::default_initializable<HalfNormalSpec>);

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
    const auto spike_slab = SpikeSlabSpec<>{};
    REQUIRE(spike_slab.probability() == 0.01);

    const auto scaled_mixture = ScaledMixtureSpec<>{};
    REQUIRE(
        scaled_mixture.probabilities()
        == std::array{0.99, 0.005, 0.003, 0.001, 0.001});
    REQUIRE(scaled_mixture.scales() == std::array{0.0, 0.001, 0.01, 0.1, 1.0});

    const auto joint_spike_slab = JointSpikeSlabSpec<>{};
    REQUIRE(
        joint_spike_slab.probabilities()
        == std::array{0.99, 1.0 / 300, 1.0 / 300, 1.0 / 300});
}

TEST_CASE("Bayes structural specs accept resolved values", "[bayes][spec]")
{
    const auto spike_slab = SpikeSlabSpec<>{0.2};
    REQUIRE(spike_slab.probability() == 0.2);

    const auto scaled_mixture = ScaledMixtureSpec<>{
        {0.8, 0.05, 0.05, 0.05, 0.05}, {0.0, 0.01, 0.1, 1.0, 10.0}};
    REQUIRE(
        scaled_mixture.probabilities()
        == std::array{0.8, 0.05, 0.05, 0.05, 0.05});
    REQUIRE(scaled_mixture.scales() == std::array{0.0, 0.01, 0.1, 1.0, 10.0});

    const auto joint_spike_slab = JointSpikeSlabSpec<>{{0.7, 0.1, 0.1, 0.1}};
    REQUIRE(joint_spike_slab.probabilities() == std::array{0.7, 0.1, 0.1, 0.1});
}

TEMPLATE_TEST_CASE(
    "SpikeSlab rejects invalid probabilities",
    "[bayes][spec]",
    (SpikeSlabSpec<VarianceLayout::Pooled, MixtureWeightUpdate::Enabled>),
    (SpikeSlabSpec<VarianceLayout::Unpooled, MixtureWeightUpdate::Enabled>),
    (SpikeSlabSpec<VarianceLayout::Pooled, MixtureWeightUpdate::Disabled>),
    (SpikeSlabSpec<VarianceLayout::Unpooled, MixtureWeightUpdate::Disabled>))
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
        message_of([probability] { return TestType{probability}; }),
        ContainsSubstring("must lie in the open interval (0, 1)"));
}

TEMPLATE_TEST_CASE(
    "ScaledMixture rejects invalid probabilities",
    "[bayes][spec]",
    ScaledMixtureSpec<MixtureWeightUpdate::Enabled>,
    ScaledMixtureSpec<MixtureWeightUpdate::Disabled>)
{
    SECTION("a probability is not positive")
    {
        REQUIRE_THAT(
            message_of([] { return TestType{{1.0, 0.0, 0.0, 0.0, 0.0}}; }),
            ContainsSubstring("probabilities[1] must be finite and positive"));
    }

    SECTION("probabilities do not sum to one")
    {
        REQUIRE_THAT(
            message_of([] { return TestType{{0.2, 0.2, 0.2, 0.2, 0.1}}; }),
            ContainsSubstring("probabilities must sum to 1"));
    }
}

TEMPLATE_TEST_CASE(
    "ScaledMixture rejects invalid scales",
    "[bayes][spec]",
    ScaledMixtureSpec<MixtureWeightUpdate::Enabled>,
    ScaledMixtureSpec<MixtureWeightUpdate::Disabled>)
{
    const auto defaults = TestType{};
    const auto probabilities = defaults.probabilities();

    SECTION("the null scale is not first")
    {
        REQUIRE_THAT(
            message_of(
                [probabilities]
                {
                    return TestType{
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
                    return TestType{
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
                return ScaledMixtureSpec<>{
                    {0.9, 0.05, 0.02, 0.02, 0.0}, {1.0, 0.1, 1.0, 10.0, 100.0}};
            });

        REQUIRE_THAT(
            message,
            ContainsSubstring("probabilities[4] must be finite and positive"));
        REQUIRE_THAT(message, !ContainsSubstring("scales[0]"));
    }
}

TEMPLATE_TEST_CASE(
    "JointSpikeSlab rejects invalid inputs",
    "[bayes][spec]",
    JointSpikeSlabSpec<MixtureWeightUpdate::Enabled>,
    JointSpikeSlabSpec<MixtureWeightUpdate::Disabled>)
{
    SECTION("allocation probabilities do not form a simplex")
    {
        REQUIRE_THAT(
            message_of([] { return TestType{{0.9, 0.05, 0.02, 0.02}}; }),
            ContainsSubstring("probabilities must sum to 1"));
    }
}
