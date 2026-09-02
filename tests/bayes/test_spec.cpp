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

#include "gelex/bayes/spec.h"

using gelex::JointSpikeSlab;
using gelex::ScaledMixture;
using gelex::SpikeSlab;
TEST_CASE("Bayes structural specs provide defaults", "[bayes][spec]")
{
    const auto spike_slab = SpikeSlab{};
    REQUIRE(spike_slab.probability == 0.01);

    const auto scaled_mixture = ScaledMixture{};
    REQUIRE(
        scaled_mixture.probabilities
        == std::array{0.99, 0.005, 0.003, 0.001, 0.001});
    REQUIRE(scaled_mixture.scales == std::array{0.0, 0.001, 0.01, 0.1, 1.0});

    const auto joint_spike_slab = JointSpikeSlab{};
    REQUIRE(
        joint_spike_slab.probabilities
        == std::array{0.99, 1.0 / 300, 1.0 / 300, 1.0 / 300});
    REQUIRE(joint_spike_slab.positive_probability == 0.5);
}

TEST_CASE("Bayes structural specs accept resolved values", "[bayes][spec]")
{
    const auto spike_slab = SpikeSlab{.probability = 0.2};
    REQUIRE(spike_slab.probability == 0.2);

    const auto joint_spike_slab = JointSpikeSlab{
        .probabilities = {0.7, 0.1, 0.1, 0.1}, .positive_probability = 0.6};
    REQUIRE(joint_spike_slab.probabilities == std::array{0.7, 0.1, 0.1, 0.1});
    REQUIRE(joint_spike_slab.positive_probability == 0.6);
}
