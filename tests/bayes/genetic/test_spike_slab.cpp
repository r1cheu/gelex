// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <random>

#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/spike_slab_kernel.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"

TEMPLATE_TEST_CASE_SIG(
    "Spike slab maintains assignments, counts and kernel residual",
    "[bayes][genetic][spike_slab]",
    ((gelex::VarianceLayout Kind, gelex::MixtureWeightUpdate Update),
     Kind,
     Update),
    (gelex::VarianceLayout::Pooled, gelex::MixtureWeightUpdate::Enabled),
    (gelex::VarianceLayout::Unpooled, gelex::MixtureWeightUpdate::Enabled),
    (gelex::VarianceLayout::Pooled, gelex::MixtureWeightUpdate::Disabled),
    (gelex::VarianceLayout::Unpooled, gelex::MixtureWeightUpdate::Disabled))
{
    const Eigen::MatrixXd genotypes{
        {0.0, 2.0}, {1.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}};
    const auto design = gelex::test::make_genetic_design(genotypes);
    const Eigen::MatrixXd x = genotypes.rowwise() - genotypes.colwise().mean();
    const gelex::SpikeSlabPrior<Kind, Update> prior{
        .variance = {1.0, gelex::make_scaled_inv_chi2_prior(4.0, 1.0)},
        .probability = gelex::detail::make_parameter<Update>(
            0.5, gelex::make_beta_prior(1.0, 1.0))};
    auto state = gelex::make_state(
        prior, gelex::GeneticDimensions{.individual = 4, .marker = 2});
    const Eigen::VectorXd response{{3.0, -1.0, 2.0, 0.5}};
    state.transition(0, 0.25, true);

    SECTION("Activation deactivation and repeated assignment")
    {
        state.transition(0, 1.5, true);
        REQUIRE(state.class_counts() == std::array<std::size_t, 2>{1, 1});
        REQUIRE(state.coefficients().isApprox(Eigen::VectorXd{{1.5, 0.0}}));

        state.transition(1, -0.75, true);
        REQUIRE(state.class_counts() == std::array<std::size_t, 2>{0, 2});

        state.transition(0, 9.0, false);
        REQUIRE(state.coefficients().isApprox(Eigen::VectorXd{{0.0, -0.75}}));
        REQUIRE(state.class_counts() == std::array<std::size_t, 2>{1, 1});

        state.transition(0, 5.0, false);
        REQUIRE(state.coefficients().isApprox(Eigen::VectorXd{{0.0, -0.75}}));
        REQUIRE(state.class_counts() == std::array<std::size_t, 2>{1, 1});
        REQUIRE(state.assignments().template cast<double>().isApprox(
            Eigen::VectorXd{{0.0, 1.0}}));
    }
    SECTION("Kernel preserves state invariants across sweeps")
    {
        auto kernel = gelex::make_kernel(prior);
        gelex::ResidualState residual_state{
            response - (x * state.coefficients()), 1.0};
        std::mt19937_64 rng{42};
        for (int sweep = 0; sweep < 5; ++sweep)
        {
            kernel.template step<gelex::GeneticMode::A>(
                design, state, residual_state, rng);

            REQUIRE(residual_state.adjusted_response.isApprox(
                response - (x * state.coefficients())));
            REQUIRE(
                state.class_counts()[1]
                == state.assignments().template cast<std::size_t>().sum());
            REQUIRE(state.class_counts()[0] + state.class_counts()[1] == 2);
            for (Eigen::Index marker = 0; marker < 2; ++marker)
            {
                if (state.assignments()(marker) == 0)
                {
                    REQUIRE(state.coefficients()(marker) == 0.0);
                }
            }
            if constexpr (Update == gelex::MixtureWeightUpdate::Disabled)
            {
                REQUIRE(state.probability() == 0.5);
            }
        }
    }
}
