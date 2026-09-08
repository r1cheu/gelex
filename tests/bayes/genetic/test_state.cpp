// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <limits>

#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/exception.h"

namespace
{
template <typename State>
concept MutableProbability
    = requires(State& state) { state.set_probability(0.5); };

template <typename State>
concept MutableProbabilities = requires(State& state) {
    state.set_probabilities(state.probabilities());
};

using gelex::MixtureWeightUpdate;
using gelex::VarianceLayout;

static_assert(
    MutableProbability<gelex::SpikeSlabState<VarianceLayout::Pooled>>);
static_assert(!MutableProbability<gelex::SpikeSlabState<
                  VarianceLayout::Pooled,
                  MixtureWeightUpdate::Disabled>>);
static_assert(MutableProbabilities<gelex::ScaledMixtureState<>>);
static_assert(!MutableProbabilities<
              gelex::ScaledMixtureState<MixtureWeightUpdate::Disabled>>);
static_assert(MutableProbabilities<gelex::JointSpikeSlabState<>>);
static_assert(!MutableProbabilities<
              gelex::JointSpikeSlabState<MixtureWeightUpdate::Disabled>>);
static_assert(requires(gelex::SpikeSlabState<VarianceLayout::Pooled>& state) {
    { state.probability() } -> std::same_as<double>;
});
static_assert(requires(gelex::ScaledMixtureState<>& state) {
    { state.probabilities() } -> std::same_as<const std::array<double, 5>&>;
});
static_assert(requires(gelex::JointSpikeSlabState<>& state) {
    { state.probabilities() } -> std::same_as<const std::array<double, 4>&>;
});
}  // namespace

TEST_CASE(
    "Spike slab state validates probability changes",
    "[bayes][genetic][state]")
{
    const gelex::GeneticDimensions dimensions{.individual = 3, .marker = 2};
    gelex::SpikeSlabState<VarianceLayout::Pooled> state{1.0, 0.25, dimensions};
    state.set_probability(0.75);
    REQUIRE(state.probability() == 0.75);
    for (double invalid :
         {-0.1,
          0.0,
          1.0,
          1.1,
          std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::quiet_NaN()})
    {
        REQUIRE_THROWS_AS(
            state.set_probability(invalid), gelex::GelexException);
        REQUIRE(state.probability() == 0.75);
        REQUIRE_THROWS_AS(
            (gelex::SpikeSlabState<VarianceLayout::Pooled>{
                1.0, invalid, dimensions}),
            gelex::GelexException);
        REQUIRE_THROWS_AS(
            (gelex::SpikeSlabState<
                VarianceLayout::Pooled,
                MixtureWeightUpdate::Disabled>{1.0, invalid, dimensions}),
            gelex::GelexException);
    }
}

TEST_CASE("Mixture states validate simplex changes", "[bayes][genetic][state]")
{
    const gelex::GeneticDimensions dimensions{.individual = 3, .marker = 2};
    gelex::ScaledMixtureState<> scaled{
        1.0, {0.6, 0.1, 0.1, 0.1, 0.1}, dimensions};
    gelex::JointSpikeSlabState<> joint{{0.7, 0.1, 0.1, 0.1}, dimensions};
    const auto check = [](auto& state)
    {
        auto probabilities = state.probabilities();
        probabilities[0] -= 0.1;
        probabilities[1] += 0.1;
        state.set_probabilities(probabilities);
        REQUIRE(state.probabilities() == probabilities);
        for (double invalid :
             {-0.1,
              0.0,
              2.0,
              std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::quiet_NaN()})
        {
            auto invalid_probabilities = probabilities;
            invalid_probabilities[0] = invalid;
            REQUIRE_THROWS_AS(
                state.set_probabilities(invalid_probabilities),
                gelex::GelexException);
            REQUIRE(state.probabilities() == probabilities);
        }
    };
    check(scaled);
    check(joint);
    REQUIRE_THROWS_AS(
        (gelex::ScaledMixtureState<>{1.0, {}, dimensions}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::ScaledMixtureState<MixtureWeightUpdate::Disabled>{
            1.0, {}, dimensions}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::JointSpikeSlabState<>{{}, dimensions}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::JointSpikeSlabState<MixtureWeightUpdate::Disabled>{
            {}, dimensions}),
        gelex::GelexException);
}
