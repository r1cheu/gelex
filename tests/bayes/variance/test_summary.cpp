// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <limits>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/bayes/variance/summary.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "bayes/bayes_model_fixture.h"

TEST_CASE(
    "Variance summary includes covariance between genetic modes",
    "[bayes][variance_summary]")
{
    constexpr auto modes = gelex::GeneticMode::A | gelex::GeneticMode::D;
    const auto model = gelex::test::make_random_effect_model(modes);
    const auto prior = gelex::make_prior(
        gelex::BuiltinBayesRecipe<modes, gelex::BayesMethod::RR>{
            gelex::VarianceBudget{
                {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    state.genetic().get<gelex::GeneticMode::A>().transition(
        Eigen::VectorXd{{0, 1, 2}});
    state.genetic().get<gelex::GeneticMode::D>().transition(
        Eigen::VectorXd{{0, 1, 2}});
    state.residual().variance = 2.0;
    const auto summary = gelex::make_variance_summary(state);
    REQUIRE(
        summary.genetic<gelex::GeneticMode::A>() == Catch::Approx(2.0 / 3.0));
    REQUIRE(summary.genetic_total() == Catch::Approx(8.0 / 3.0));
    REQUIRE(summary.total_heritability() == Catch::Approx(4.0 / 7.0));
}

TEST_CASE(
    "Variance summary rejects invalid variance components",
    "[bayes][variance_summary]")
{
    constexpr auto modes = gelex::GeneticModeSet{gelex::GeneticMode::A};
    using Values = gelex::HomogeneousModeValues<modes, double>;
    using Summary = gelex::VarianceSummary<modes>;
    REQUIRE_THROWS_AS((Summary{Values{0.0}, 0.0, 0.0}), gelex::GelexException);
    REQUIRE_THROWS_AS((Summary{Values{1.0}, 1.0, -1.0}), gelex::GelexException);
    REQUIRE_THROWS_AS((Summary{Values{-1.0}, 1.0, 1.0}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (Summary{Values{1.0}, std::numeric_limits<double>::infinity(), 1.0}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (Summary{Values{std::numeric_limits<double>::quiet_NaN()}, 1.0, 1.0}),
        gelex::GelexException);
}
