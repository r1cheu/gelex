// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/draws.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"

#include "bayes_model_fixture.h"
#include "file_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;

}  // namespace

TEST_CASE(
    "Bayes draws serialize state and commit a short run",
    "[bayes][draws]")
{
    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "short.draws").string();
    const auto model = gelex::test::make_random_effect_model(mode_ad);
    const auto prior = gelex::make_prior(
        gelex::BuiltinBayesRecipe<mode_ad, gelex::BayesMethod::RR>{
            gelex::VarianceBudget{
                {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    state.fixed().coefficients = Eigen::VectorXd{{1.25}};
    state.random()[0].coefficients = Eigen::VectorXd{{2.0, -3.0}};
    state.random()[0].variance = 4.0;
    state.residual().variance = 5.0;
    state.genetic().get<gelex::GeneticMode::A>().transition(0, 0.5);
    state.genetic().get<gelex::GeneticMode::D>().transition(1, -0.25);
    {
        auto draws = gelex::BayesDraws{state, model, path, 3};
        draws.append(state);
        REQUIRE_FALSE(std::filesystem::exists(path));
    }
    const gelex::BinaryReader reader{path};
    REQUIRE(reader.to_map<double>("fixed/coefficients")
                .isApprox(Eigen::MatrixXd{{1.25}}));
    REQUIRE(reader.to_map<float>("random/batch/coefficients")
                .isApprox(Eigen::VectorXf{{2.0F, -3.0F}}));
    REQUIRE(reader.to_map<double>("random/batch/variance")
                .isApprox(Eigen::MatrixXd{{4.0}}));
    REQUIRE(reader.to_map<double>("residual/variance")
                .isApprox(Eigen::MatrixXd{{5.0}}));
    REQUIRE(reader.to_map<float>("genetic/A/coefficients")
                .isApprox(Eigen::VectorXf{{0.5F, 0.0F}}));
    REQUIRE(reader.to_map<float>("genetic/D/coefficients")
                .isApprox(Eigen::VectorXf{{0.0F, -0.25F}}));
}

TEST_CASE(
    "Bayes draws reject samples beyond the reserved count",
    "[bayes][draws]")
{
    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "bounded.draws").string();
    const auto model = gelex::test::make_random_effect_model(mode_a);
    const auto prior = gelex::make_prior(
        gelex::BuiltinBayesRecipe<mode_a, gelex::BayesMethod::RR>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.1}}},
        model);
    const auto state = gelex::make_state(prior, model);
    {
        auto draws = gelex::BayesDraws{state, model, path, 1};
        draws.append(state);
        REQUIRE_THROWS_AS(draws.append(state), gelex::GelexException);
    }
    const gelex::BinaryReader reader{path};
    REQUIRE(reader.to_map<float>("genetic/A/coefficients").cols() == 1);
    REQUIRE(reader.to_map<double>("residual/variance").cols() == 1);
}
