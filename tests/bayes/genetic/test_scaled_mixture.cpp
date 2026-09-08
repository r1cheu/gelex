// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <random>

#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/scaled_mixture_kernel.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"

TEMPLATE_TEST_CASE_SIG(
    "Scaled mixture maintains assignments counts and kernel residual",
    "[bayes][genetic][scaled_mixture]",
    ((gelex::MixtureWeightUpdate Update), Update),
    gelex::MixtureWeightUpdate::Enabled,
    gelex::MixtureWeightUpdate::Disabled)
{
    using state_type = gelex::ScaledMixtureState<Update>;
    const Eigen::MatrixXd genotypes{
        {0.0, 2.0}, {1.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}};
    const auto design = gelex::test::make_genetic_design(genotypes);
    const Eigen::MatrixXd x = genotypes.rowwise() - genotypes.colwise().mean();
    const std::array<double, 5> probabilities{0.6, 0.1, 0.1, 0.1, 0.1};
    const gelex::ScaledMixturePrior<Update> prior{
        .variance = {1.0, gelex::make_scaled_inv_chi2_prior(4.0, 1.0)},
        .probabilities = gelex::detail::make_parameter<Update>(
            probabilities, gelex::make_uniform_dirichlet_prior<5>()),
        .scales = {0.0, 0.01, 0.1, 1.0, 10.0}};
    const gelex::GeneticDimensions dimensions{.individual = 4, .marker = 2};
    const Eigen::VectorXd response{{3.0, -1.0, 2.0, 0.5}};
    const auto check = [](const state_type& state)
    {
        std::array<std::size_t, state_type::class_count> counts{};
        for (Eigen::Index marker = 0; marker < 2; ++marker)
        {
            const auto assignment = state.assignments()(marker);
            ++counts.at(assignment);
            if (assignment == 0)
            {
                REQUIRE(state.coefficients()(marker) == 0.0);
            }
        }

        REQUIRE(state.class_counts() == counts);
    };

    SECTION(
        "All class transitions including equal coefficients and zero effects")
    {
        for (std::uint8_t old_class = 0; old_class < state_type::class_count;
             ++old_class)
        {
            for (std::uint8_t new_class = 0;
                 new_class < state_type::class_count;
                 ++new_class)
            {
                for (double new_value : {0.75, -1.25, 0.0})
                {
                    CAPTURE(old_class, new_class, new_value);
                    auto state = gelex::make_state(prior, dimensions);
                    state.transition(0, 0.75, old_class);
                    state.transition(1, 0.5, 4);
                    check(state);
                    state.transition(0, new_value, new_class);
                    REQUIRE(state.assignments()(0) == new_class);
                    REQUIRE(
                        state.coefficients()(0)
                        == (new_class == 0 ? 0.0 : new_value));
                    check(state);
                }
            }
        }
    }
    SECTION("Kernel maintains residual and assignment invariants across sweeps")
    {
        auto state = gelex::make_state(prior, dimensions);
        gelex::ResidualState residual{response, 1.0};
        auto kernel = gelex::make_kernel(prior);
        std::mt19937_64 rng{42};
        for (int sweep = 0; sweep < 5; ++sweep)
        {
            kernel.template step<gelex::GeneticMode::A>(
                design, state, residual, rng);
            check(state);
            REQUIRE(residual.adjusted_response.isApprox(
                response - (x * state.coefficients())));
            if constexpr (Update == gelex::MixtureWeightUpdate::Disabled)
            {
                REQUIRE(state.probabilities() == probabilities);
            }
        }
    }
}
