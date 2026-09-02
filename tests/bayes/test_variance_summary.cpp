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
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/bayes/variance_summary.h"
#include "gelex/exception.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"
#include "gelex/types/mode_values.h"

#include "compact_genotype_fixture.h"

using Catch::Approx;

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;
using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

template <gelex::GeneticModeSet Modes, typename... Priors>
using GaussianStateFor
    = gelex::detail::genetic_state_t<gelex::ModeValues<Modes, Priors...>>;
using AdditiveState = GaussianStateFor<
    mode_a,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>>;
using AdditiveDominanceState = GaussianStateFor<
    mode_ad,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>>;

// The summed genetic value must stay a view onto the state: materializing it
// would allocate a length-n vector on every draw.
static_assert(std::is_reference_v<decltype(gelex::total_genetic_value<mode_a>(
                  std::declval<const AdditiveState&>()))>);
static_assert(!std::is_same_v<
              std::remove_cvref_t<decltype(gelex::total_genetic_value<mode_ad>(
                  std::declval<const AdditiveDominanceState&>()))>,
              Eigen::VectorXd>);

auto make_ad_model_with_random() -> gelex::BayesModel
{
    auto genetic = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        mode_ad);
    std::vector<gelex::bayes::RandomDesign> random;
    random.emplace_back(
        "batch",
        std::vector<std::string>{"batch"},
        Eigen::MatrixXd{{0.0}, {1.0}, {0.0}, {1.0}});
    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        gelex::FixedDesign::make(4),
        std::move(random),
        std::move(genetic)};
}

}  // namespace

TEST_CASE(
    "genetic_value reduces a per-class decomposition to the mode total",
    "[bayes][variance_summary]")
{
    gelex::ScaledMixtureState<3> state;
    state.fitted_values = Eigen::MatrixXd{{1.0, 2.0}, {3.0, 4.0}};

    REQUIRE(gelex::genetic_value(state).isApprox(Eigen::VectorXd{{3.0, 7.0}}));
}

TEST_CASE(
    "make_variance_summary keeps the mode covariance out of the parts",
    "[bayes][variance_summary]")
{
    const auto model = make_ad_model_with_random();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, Family>{gelex::VarianceBudget{
            {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);

    state.genetic().get<gelex::GeneticMode::A>().family_state.fitted_values
        = Eigen::VectorXd{{1.0, 2.0, 3.0, 4.0}};
    state.genetic().get<gelex::GeneticMode::D>().family_state.fitted_values
        = Eigen::VectorXd{{0.0, 1.0, 0.0, 1.0}};
    state.random()[0].fitted_values = Eigen::VectorXd{{2.0, 0.0, 2.0, 0.0}};
    state.residual().variance = 3.0;

    const auto summary = gelex::make_variance_summary(state);

    REQUIRE(summary.genetic<gelex::GeneticMode::A>() == Approx(1.25));
    REQUIRE(summary.genetic<gelex::GeneticMode::D>() == Approx(0.25));
    // The modes are positively correlated, so the total {1, 3, 3, 5} exceeds
    // the sum of the parts.
    REQUIRE(summary.genetic_total() == Approx(2.0));
    REQUIRE(summary.random_total() == Approx(1.0));
    REQUIRE(summary.residual() == Approx(3.0));
    REQUIRE(summary.phenotypic() == Approx(6.0));
    REQUIRE(
        summary.heritability<gelex::GeneticMode::A>() == Approx(1.25 / 6.0));
    REQUIRE(
        summary.heritability<gelex::GeneticMode::D>() == Approx(0.25 / 6.0));
    REQUIRE(summary.total_heritability() == Approx(2.0 / 6.0));
}

TEST_CASE(
    "make_variance_summary rejects a draw its ratios are undefined on",
    "[bayes][variance_summary]")
{
    const auto model = make_ad_model_with_random();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, Family>{gelex::VarianceBudget{
            {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);

    const auto summarize = [&] { return gelex::make_variance_summary(state); };

    SECTION("no phenotypic variance to divide by")
    {
        // Every component is zero, so heritability would be 0/0.
        state.residual().variance = 0.0;
        REQUIRE_THROWS_AS(summarize(), gelex::GelexException);
    }

    SECTION("a diverged genetic value")
    {
        state.genetic().get<gelex::GeneticMode::A>().family_state.fitted_values
            = Eigen::VectorXd{
                {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0}};
        REQUIRE_THROWS_AS(summarize(), gelex::GelexException);
    }

    SECTION("a negative residual variance")
    {
        state.residual().variance = -1.0;
        REQUIRE_THROWS_AS(summarize(), gelex::GelexException);
    }
}
