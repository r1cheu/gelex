// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <cstdint>
#include <filesystem>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/operations.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/random_design.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/data/fixed_design.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"

#include "bayes_model_fixture.h"
#include "compact_genotype_fixture.h"
#include "file_fixture.h"
#include "random_design_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;

using JointSpikeSlabSpecAD = gelex::JointModeValues<
    gelex::ModeValues<mode_ad, gelex::GaussianSpec<>, gelex::HalfNormalSpec>,
    gelex::JointSpikeSlabSpec<>>;
using ScaledMixtureSpecAD
    = gelex::HomogeneousModeValues<mode_ad, gelex::ScaledMixtureSpec<>>;

}  // namespace

// Scaled mixtures keep a per-class decomposition, so summing the modes folds
// two lazy row-sum expressions; a dangling operand there would surface as
// garbage rather than a compile error.
TEST_CASE("BayesDraws decomposes per-class genetic values", "[bayes][draws]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "mixture_decomposition.draws";
    const auto model = gelex::test::make_random_effect_model(mode_ad);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, ScaledMixtureSpecAD>{gelex::VarianceBudget{
            {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);

    {
        auto draws = gelex::BayesDraws{state, model, path.string(), 1};

        // Each mode's class columns sum row-wise to {0, 3, 0}.
        state.genetic().for_each(
            []<gelex::GeneticMode Mode>(auto& genetic)
            {
                auto first = std::get<gelex::bayes::AxpyTarget>(
                    genetic.transition(0, 1.0, 1));
                Eigen::Map<Eigen::VectorXd>(
                    first.target.data(),
                    static_cast<Eigen::Index>(first.target.size()))
                    += first.scale * Eigen::VectorXd{{0.0, 1.0, 0.0}};
                auto second = std::get<gelex::bayes::AxpyTarget>(
                    genetic.transition(1, 2.0, 2));
                Eigen::Map<Eigen::VectorXd>(
                    second.target.data(),
                    static_cast<Eigen::Index>(second.target.size()))
                    += second.scale * Eigen::VectorXd{{0.0, 1.0, 0.0}};
            });
        state.residual().variance = 2.0;
        draws.append(state);
    }

    const gelex::BinaryReader reader(path.string());
    // The identical modes have total variance 8; residual variance adds 2.
    REQUIRE(reader.to_map<double>("genetic/A/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0}}));
    REQUIRE(reader.to_map<double>("genetic/total/explained_variance")
                .isApprox(Eigen::MatrixXd{{8.0}}));
    REQUIRE(reader.to_map<double>("genetic/total/heritability")
                .isApprox(Eigen::MatrixXd{{0.8}}));
}

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
