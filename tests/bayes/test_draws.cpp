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

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <filesystem>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance_budget.h"
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

TEST_CASE("BayesDraws records every state component", "[bayes][draws]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "bayes.draws";
    const auto model = gelex::test::make_random_effect_model(mode_a);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);

    {
        gelex::BayesDraws draws(prior, model, path.string(), 2);
        static_assert(!std::movable<decltype(draws)>);
        auto& additive = state.genetic().get<gelex::GeneticMode::A>();

        state.fixed().coefficients = Eigen::VectorXd{{1.0}};
        state.random()[0].coefficients = Eigen::VectorXd{{2.0, 3.0}};
        state.random()[0].variance = 4.0;
        additive.coefficients = Eigen::VectorXd{{5.0, 6.0}};
        additive.family_state.variance = 0.5;
        state.residual().variance = 7.0;
        draws.append(state);

        state.fixed().coefficients = Eigen::VectorXd{{3.0}};
        state.random()[0].coefficients = Eigen::VectorXd{{4.0, 5.0}};
        state.random()[0].variance = 6.0;
        additive.coefficients = Eigen::VectorXd{{7.0, 8.0}};
        additive.family_state.variance = 1.5;
        state.residual().variance = 9.0;
        draws.append(state);

        REQUIRE_FALSE(std::filesystem::exists(path));

        REQUIRE(draws.fixed().result().mean.isApprox(Eigen::VectorXd{{2.0}}));
        REQUIRE(draws.random()[0].coefficients().result().mean.isApprox(
            Eigen::VectorXd{{3.0, 4.0}}));
        REQUIRE(
            draws.random()[0].variance().result().mean == Catch::Approx(5.0));
        REQUIRE(draws.genetic()
                    .get<gelex::GeneticMode::A>()
                    .coefficients.result()
                    .mean.isApprox(Eigen::VectorXd{{6.0, 7.0}}));
        REQUIRE(draws.residual().result().mean == Catch::Approx(8.0));
    }

    REQUIRE(std::filesystem::exists(path));
    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("fixed/coefficients")
                .isApprox(Eigen::MatrixXf{{1.0, 3.0}}));
    REQUIRE(reader.to_map<float>("random/batch/coefficients")
                .isApprox(Eigen::MatrixXf{{2.0, 4.0}, {3.0, 5.0}}));
    REQUIRE(reader.to_map<double>("random/batch/variance")
                .isApprox(Eigen::MatrixXd{{4.0, 6.0}}));
    REQUIRE(reader.to_map<float>("genetic/A/coefficients")
                .isApprox(Eigen::MatrixXf{{5.0, 7.0}, {6.0, 8.0}}));
    REQUIRE(reader.to_map<double>("genetic/A/variance")
                .isApprox(Eigen::MatrixXd{{0.5, 1.5}}));
    REQUIRE(reader.to_map<double>("residual/variance")
                .isApprox(Eigen::MatrixXd{{7.0, 9.0}}));
}

TEST_CASE("BayesDraws bounds the number of appended draws", "[bayes][draws]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "exceeded.draws";
    const auto model = gelex::test::make_random_effect_model(mode_a);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);

    gelex::BayesDraws draws(prior, model, path.string(), 1);
    draws.append(state);

    REQUIRE_THROWS_AS(draws.append(state), gelex::GelexException);
}

TEST_CASE("BayesDraws records the variance decomposition", "[bayes][draws]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "decomposition.draws";
    const auto model = gelex::test::make_random_effect_model(mode_a);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);

    {
        gelex::BayesDraws draws(prior, model, path.string(), 1);

        state.genetic().get<gelex::GeneticMode::A>().family_state.fitted_values
            = Eigen::VectorXd{{0.0, 3.0, 0.0}};
        state.random()[0].fitted_values = Eigen::VectorXd{{0.0, 3.0, 0.0}};
        state.residual().variance = 4.0;
        draws.append(state);

        const auto& summary = draws.variance_summary();
        REQUIRE(
            summary.explained_variance<gelex::GeneticMode::A>().result().mean
            == Catch::Approx(2.0));
        REQUIRE(
            summary.total_heritability().result().mean == Catch::Approx(0.25));
        REQUIRE(
            draws.random()[0].explained_variance().result().mean
            == Catch::Approx(2.0));
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<double>("genetic/A/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0}}));
    REQUIRE(reader.to_map<double>("genetic/A/heritability")
                .isApprox(Eigen::MatrixXd{{0.25}}));
    // A single-mode model still writes the totals, so readers see one layout.
    REQUIRE(reader.to_map<double>("genetic/total/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0}}));
    REQUIRE(reader.to_map<double>("genetic/total/heritability")
                .isApprox(Eigen::MatrixXd{{0.25}}));
    REQUIRE(reader.to_map<double>("random/batch/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0}}));
}

// Joint families reach the mode states through JointModeValues, a path no
// other draws test covers.
TEST_CASE("BayesDraws decomposes a joint spike-slab state", "[bayes][draws]")
{
    using Family = gelex::JointSpikeSlabFamily<>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "joint_decomposition.draws";
    const auto model = gelex::test::make_random_effect_model(mode_ad);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, Family>{gelex::VarianceBudget{
            {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);

    {
        gelex::BayesDraws draws(prior, model, path.string(), 1);

        state.genetic().get<gelex::GeneticMode::A>().family_state.fitted_values
            = Eigen::VectorXd{{1.0, 2.0, 3.0}};
        state.genetic().get<gelex::GeneticMode::D>().family_state.fitted_values
            = Eigen::VectorXd{{0.0, 1.0, 2.0}};
        state.random()[0].fitted_values = Eigen::VectorXd{{0.0, 3.0, 0.0}};
        state.residual().variance = 10.0 / 3.0;
        draws.append(state);
    }

    const gelex::BinaryReader reader(path.string());
    // Each mode alone varies by 2/3, but their sum {1, 3, 5} varies by 8/3:
    // the total carries the A-D covariance that the parts cannot.
    REQUIRE(reader.to_map<double>("genetic/A/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0 / 3.0}}));
    REQUIRE(reader.to_map<double>("genetic/D/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0 / 3.0}}));
    REQUIRE(reader.to_map<double>("genetic/total/explained_variance")
                .isApprox(Eigen::MatrixXd{{8.0 / 3.0}}));
    // Phenotypic variance is 8/3 + 2 + 10/3 = 8.
    REQUIRE(reader.to_map<double>("genetic/A/heritability")
                .isApprox(Eigen::MatrixXd{{1.0 / 12.0}}));
    REQUIRE(reader.to_map<double>("genetic/total/heritability")
                .isApprox(Eigen::MatrixXd{{1.0 / 3.0}}));
    REQUIRE(reader.to_map<double>("random/batch/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0}}));
}

// Scaled mixtures keep a per-class decomposition, so summing the modes folds
// two lazy row-sum expressions; a dangling operand there would surface as
// garbage rather than a compile error.
TEST_CASE("BayesDraws decomposes per-class genetic values", "[bayes][draws]")
{
    using Family = gelex::ScaledMixtureFamily<>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "mixture_decomposition.draws";
    const auto model = gelex::test::make_random_effect_model(mode_ad);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, Family>{gelex::VarianceBudget{
            {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);

    {
        gelex::BayesDraws draws(prior, model, path.string(), 1);

        // Each mode's class columns sum row-wise to {0, 3, 0}.
        const Eigen::MatrixXd classes{
            {0.0, 0.0, 0.0, 0.0}, {1.0, 1.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
        state.genetic().get<gelex::GeneticMode::A>().family_state.fitted_values
            = classes;
        state.genetic().get<gelex::GeneticMode::D>().family_state.fitted_values
            = classes;
        state.random()[0].fitted_values = Eigen::VectorXd{{0.0, 3.0, 0.0}};
        state.residual().variance = 2.0;
        draws.append(state);
    }

    const gelex::BinaryReader reader(path.string());
    // Each mode varies by 2, their sum {0, 6, 0} by 8; phenotypic is 8+2+2.
    REQUIRE(reader.to_map<double>("genetic/A/explained_variance")
                .isApprox(Eigen::MatrixXd{{2.0}}));
    REQUIRE(reader.to_map<double>("genetic/total/explained_variance")
                .isApprox(Eigen::MatrixXd{{8.0}}));
    REQUIRE(reader.to_map<double>("genetic/total/heritability")
                .isApprox(Eigen::MatrixXd{{2.0 / 3.0}}));
}

TEST_CASE("BayesDraws commits a short run", "[bayes][draws]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "short.draws";
    const auto model = gelex::test::make_random_effect_model(mode_a);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);
    state.residual().variance = 5.0;

    {
        gelex::BayesDraws draws(prior, model, path.string(), 3);
        draws.append(state);
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<double>("residual/variance")
                .isApprox(Eigen::MatrixXd{{5.0}}));
}
