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
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};

auto make_model() -> gelex::BayesModel
{
    auto genetic = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, mode_a);
    std::vector<gelex::bayes::RandomDesign> random;
    random.emplace_back(
        "batch",
        std::vector<std::string>{"batch_1", "batch_2"},
        Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}});
    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        std::move(random),
        std::move(genetic)};
}

}  // namespace

TEST_CASE("make_draws records every state component", "[bayes][draws]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "bayes.draws";
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);

    {
        auto draws = gelex::make_draws(prior, model, path.string(), 2);
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
        draws.close();
        REQUIRE(std::filesystem::exists(path));

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

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<float>("fixed/coefficients")
                .isApprox(Eigen::MatrixXf{{1.0, 3.0}}));
    REQUIRE(reader.to_map<float>("random/0/coefficients")
                .isApprox(Eigen::MatrixXf{{2.0, 4.0}, {3.0, 5.0}}));
    REQUIRE(reader.to_map<double>("random/0/variance")
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
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);

    auto draws = gelex::make_draws(prior, model, path.string(), 1);
    draws.append(state);

    REQUIRE_THROWS_AS(draws.append(state), gelex::GelexException);
    draws.close();
}

TEST_CASE("BayesDraws commits a short run", "[bayes][draws]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "short.draws";
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);
    state.residual().variance = 5.0;

    {
        auto draws = gelex::make_draws(prior, model, path.string(), 3);
        draws.append(state);
        draws.close();
    }

    const gelex::BinaryReader reader(path.string());
    REQUIRE(reader.to_map<double>("residual/variance")
                .isApprox(Eigen::MatrixXd{{5.0}}));
}
