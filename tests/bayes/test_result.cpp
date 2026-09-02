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
#include <algorithm>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/result_factory.h"
#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/result.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/fixed_design.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"

#include "bayes_model_fixture.h"
#include "compact_genotype_fixture.h"
#include "file_fixture.h"
#include "random_design_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;

template <gelex::GeneticModeSet Modes, typename Family, typename Configure>
auto collect_result(const std::filesystem::path& path, Configure configure)
{
    auto model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        Modes);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<Modes, Family>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    configure(state);

    auto draws = gelex::make_draws(prior, model, path.string(), 1);
    draws.append(state);
    return gelex::make_result(model, draws);
}

template <typename T>
concept EmptyVariance = requires(T value) {
    requires std::same_as<
        std::remove_cvref_t<decltype(value.variance)>,
        gelex::EmptyPosteriorResult>;
};

template <typename T>
concept EmptyProbability = requires(T value) {
    requires std::same_as<
        std::remove_cvref_t<decltype(value.probability)>,
        gelex::EmptyPosteriorResult>;
};

template <typename T>
concept EmptyProbabilities = requires(T value) {
    requires std::same_as<
        std::remove_cvref_t<decltype(value.probabilities)>,
        gelex::EmptyPosteriorResult>;
};

}  // namespace

TEST_CASE(
    "Posterior result constructors validate stored invariants",
    "[bayes][result]")
{
    REQUIRE_THROWS_AS(
        (gelex::ScalarPosteriorResult{"", {}}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::VectorPosteriorResult{"", {}}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::VectorPosteriorResult{
            "posterior/vector",
            gelex::VectorRunningStatsResult{
                .mean = Eigen::VectorXd{{1.0, 2.0}},
                .stddev = Eigen::VectorXd{{0.1}}}}),
        gelex::GelexException);

    auto posterior = gelex::VectorPosteriorResult{
        "fixed/coefficients",
        gelex::VectorRunningStatsResult{
            .mean = Eigen::VectorXd{{1.0, 2.0}},
            .stddev = Eigen::VectorXd{{0.1, 0.2}}}};
    REQUIRE_THROWS_AS(
        (gelex::CoefficientPosteriorResult{
            std::move(posterior), {std::string{gelex::intercept_name}}}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::MarkerPipResult{Eigen::VectorXd{}}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::MarkerPipResult{Eigen::VectorXd{{-0.1, 0.5}}}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::MarkerPipResult{Eigen::VectorXd{{0.5, 1.1}}}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::MarkerPveResult{Eigen::VectorXd{}}), gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::MarkerPveResult{Eigen::VectorXd{{-0.1, 0.5}}}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        (gelex::MarkerPveResult{
            Eigen::VectorXd{{0.5, std::numeric_limits<double>::quiet_NaN()}}}),
        gelex::GelexException);
    REQUIRE_NOTHROW(gelex::MarkerPveResult(Eigen::VectorXd{{0.5, 1.1}}));
}

TEST_CASE(
    "BayesResult owns named fixed random and variance summaries",
    "[bayes][result]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "owned_result.draws";
    const auto result = [&]
    {
        auto model = gelex::test::make_random_effect_model(mode_a);
        const auto prior = gelex::make_prior(
            gelex::BayesRecipe<mode_a, Family>{
                gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
            model);
        auto state = gelex::make_state(prior, model);
        state.fixed().coefficients = Eigen::VectorXd{{2.0}};
        state.random()[0].coefficients = Eigen::VectorXd{{3.0, 4.0}};
        state.random()[0].variance = 5.0;
        state.random()[0].fitted_values = Eigen::VectorXd{{0.0, 3.0, 0.0}};
        state.genetic().get<gelex::GeneticMode::A>().coefficients
            = Eigen::VectorXd{{6.0, 7.0}};
        state.genetic().get<gelex::GeneticMode::A>().family_state.variance
            = 0.5;
        state.residual().variance = 4.0;

        auto draws = gelex::make_draws(prior, model, path.string(), 1);
        draws.append(state);
        static_assert(std::same_as<
                      decltype(gelex::make_result(model, draws)),
                      gelex::BayesResult<
                          typename decltype(prior)::genetic_prior_type>>);
        return gelex::make_result(model, draws);
    }();

    REQUIRE(result.fixed().identifier() == "fixed/coefficients");
    REQUIRE(
        std::ranges::equal(
            result.fixed().column_names(),
            std::array{std::string{gelex::intercept_name}}));
    REQUIRE(result.fixed().statistics().mean.isApprox(Eigen::VectorXd{{2.0}}));

    REQUIRE(result.random().size() == 1);
    const auto& random = result.random()[0];
    REQUIRE(random.coefficients().identifier() == "random/batch/coefficients");
    REQUIRE(
        std::ranges::equal(
            random.coefficients().column_names(),
            std::array{std::string{"batch_1"}, std::string{"batch_2"}}));
    REQUIRE(random.coefficients().statistics().mean.isApprox(
        Eigen::VectorXd{{3.0, 4.0}}));
    REQUIRE(random.variance().identifier() == "random/batch/variance");
    REQUIRE(
        random.explained_variance().identifier()
        == "random/batch/explained_variance");

    const auto& marker_effect
        = result.marker_effects().get<gelex::GeneticMode::A>();
    REQUIRE(
        marker_effect.coefficients().identifier() == "genetic/A/coefficients");
    REQUIRE(marker_effect.coefficients().statistics().mean.isApprox(
        Eigen::VectorXd{{6.0, 7.0}}));

    const auto& genetic
        = result.genetic_parameters().get<gelex::GeneticMode::A>();
    REQUIRE(genetic.variance.identifier() == "genetic/A/variance");
    REQUIRE(genetic.variance.statistics().mean == Catch::Approx(0.5));
    STATIC_REQUIRE(
        std::same_as<
            std::remove_cvref_t<decltype(marker_effect.pip())>,
            gelex::EmptyPosteriorResult>);
    REQUIRE(marker_effect.pve().values().isApprox(
        Eigen::VectorXd{{36.0, 49.0 / 3.0}}));
    REQUIRE(result.residual().identifier() == "residual/variance");
    REQUIRE(result.residual().statistics().mean == Catch::Approx(4.0));

    const auto& summary = result.variance_summary();
    REQUIRE(
        summary.explained_variance<gelex::GeneticMode::A>().identifier()
        == "genetic/A/explained_variance");
    REQUIRE(
        summary.heritability<gelex::GeneticMode::A>().identifier()
        == "genetic/A/heritability");
    REQUIRE(
        summary.total_explained_variance().identifier()
        == "genetic/total/explained_variance");
    REQUIRE(
        summary.total_heritability().identifier()
        == "genetic/total/heritability");
}

TEST_CASE(
    "BayesResult excludes unpooled marker variance without materializing its "
    "statistics",
    "[bayes][result]")
{
    gelex::test::FileFixture fixture;
    gelex::BinaryWriter writer(
        (fixture.get_test_dir() / "unmaterialized.draws").string());
    gelex::GaussianDraws<gelex::VarianceLayout::Unpooled> draws{
        .variance = gelex::VectorDraw{writer.reserve<float>(
            "genetic/A/variance", gelex::BinaryShape{2, 0})}};

    const auto result = gelex::detail::make_result(draws);

    STATIC_REQUIRE(EmptyVariance<decltype(result)>);
}

TEST_CASE(
    "BayesResult projects pooled unpooled sampled and fixed spike-slab fields",
    "[bayes][result]")
{
    gelex::test::FileFixture fixture;

    SECTION("sampled unpooled")
    {
        using Family = gelex::SpikeSlabFamily<
            gelex::VarianceLayout::Unpooled,
            gelex::MixtureWeightUpdate::Enabled>;
        const auto result = collect_result<mode_a, Family>(
            fixture.get_test_dir() / "sampled_unpooled.draws",
            [](auto& state)
            {
                auto& family = state.genetic()
                                   .template get<gelex::GeneticMode::A>()
                                   .family_state;
                family.variance = Eigen::VectorXd{{0.1, 0.2}};
                family.assignment = Eigen::VectorX<std::uint8_t>{{0, 1}};
                family.probability = 0.35;
            });
        const auto& family
            = result.genetic_parameters().template get<gelex::GeneticMode::A>();

        STATIC_REQUIRE(EmptyVariance<std::remove_cvref_t<decltype(family)>>);
        REQUIRE(family.probability.identifier() == "genetic/A/probability");
        REQUIRE(family.probability.statistics().mean == Catch::Approx(0.35));
        REQUIRE(result.marker_effects()
                    .template get<gelex::GeneticMode::A>()
                    .pip()
                    .probabilities()
                    .isApprox(Eigen::VectorXd{{0.0, 1.0}}));
    }

    SECTION("fixed pooled")
    {
        using Family = gelex::SpikeSlabFamily<
            gelex::VarianceLayout::Pooled,
            gelex::MixtureWeightUpdate::Disabled>;
        const auto result = collect_result<mode_a, Family>(
            fixture.get_test_dir() / "fixed_pooled.draws",
            [](auto& state)
            {
                state.genetic()
                    .template get<gelex::GeneticMode::A>()
                    .family_state.variance = 0.75;
                state.genetic()
                    .template get<gelex::GeneticMode::A>()
                    .family_state.assignment
                    = Eigen::VectorX<std::uint8_t>{{1, 0}};
            });
        const auto& family
            = result.genetic_parameters().template get<gelex::GeneticMode::A>();

        STATIC_REQUIRE(EmptyProbability<std::remove_cvref_t<decltype(family)>>);
        REQUIRE(family.variance.identifier() == "genetic/A/variance");
        REQUIRE(family.variance.statistics().mean == Catch::Approx(0.75));
        REQUIRE(result.marker_effects()
                    .template get<gelex::GeneticMode::A>()
                    .pip()
                    .probabilities()
                    .isApprox(Eigen::VectorXd{{1.0, 0.0}}));
    }
}

TEST_CASE(
    "BayesResult projects scaled-mixture sampled and fixed fields",
    "[bayes][result]")
{
    gelex::test::FileFixture fixture;

    SECTION("sampled")
    {
        using Family = gelex::ScaledMixtureFamily<>;
        const auto result = collect_result<mode_a, Family>(
            fixture.get_test_dir() / "sampled_mixture.draws",
            [](auto& state)
            {
                auto& family = state.genetic()
                                   .template get<gelex::GeneticMode::A>()
                                   .family_state;
                family.variance = 0.8;
                family.assignment = Eigen::VectorX<std::uint8_t>{{0, 4}};
                family.probabilities = {0.5, 0.2, 0.15, 0.1, 0.05};
                family.fitted_values = Eigen::MatrixXd{
                    {0.0, 1.0, 2.0, 0.0},
                    {3.0, 1.0, 0.0, 0.0},
                    {0.0, 1.0, 4.0, 0.0}};
            });
        const auto& family
            = result.genetic_parameters().template get<gelex::GeneticMode::A>();

        REQUIRE(family.variance.identifier() == "genetic/A/variance");
        REQUIRE(family.probabilities.identifier() == "genetic/A/probabilities");
        REQUIRE(family.probabilities.statistics().mean.isApprox(
            Eigen::VectorXd{{0.5, 0.2, 0.15, 0.1, 0.05}}));
        REQUIRE(
            family.component_explained_variance.identifier()
            == "genetic/A/component_explained_variance");
        REQUIRE(
            family.component_explained_variance.statistics().mean.size() == 4);
        REQUIRE(result.marker_effects()
                    .template get<gelex::GeneticMode::A>()
                    .pip()
                    .probabilities()
                    .isApprox(Eigen::VectorXd{{0.0, 1.0}}));
    }

    SECTION("fixed")
    {
        using Family
            = gelex::ScaledMixtureFamily<gelex::MixtureWeightUpdate::Disabled>;
        const auto result = collect_result<mode_a, Family>(
            fixture.get_test_dir() / "fixed_mixture.draws",
            [](auto& state)
            {
                state.genetic()
                    .template get<gelex::GeneticMode::A>()
                    .family_state.fitted_values.setZero();
            });
        const auto& family
            = result.genetic_parameters().template get<gelex::GeneticMode::A>();

        STATIC_REQUIRE(
            EmptyProbabilities<std::remove_cvref_t<decltype(family)>>);
        REQUIRE(
            family.component_explained_variance.identifier()
            == "genetic/A/component_explained_variance");
    }
}

TEST_CASE("BayesResult projects half-normal joint fields", "[bayes][result]")
{
    gelex::test::FileFixture fixture;

    SECTION("count")
    {
        using Family = gelex::JointSpikeSlabFamily<>;
        const auto result = collect_result<mode_ad, Family>(
            fixture.get_test_dir() / "count_joint.draws",
            [](auto& state)
            {
                auto& dominance = state.genetic()
                                      .template get<gelex::GeneticMode::D>()
                                      .family_state;
                dominance.variance = 0.25;
                dominance.probit_coefficients = Eigen::Vector2d{{0.7, -0.2}};
                auto& joint = state.genetic().joint();
                joint.assignment = Eigen::VectorX<std::uint8_t>{{1, 3}};
                joint.probabilities = {0.7, 0.1, 0.1, 0.1};
                joint.fitted_values = Eigen::MatrixXd{
                    {1.0, 0.0, 2.0, 0.0},
                    {2.0, 0.0, 2.0, 4.0},
                    {3.0, 3.0, 2.0, 2.0}};
            });
        const auto& dominance
            = result.genetic_parameters().template get<gelex::GeneticMode::D>();

        REQUIRE(dominance.variance.identifier() == "genetic/D/variance");
        REQUIRE(
            dominance.probit_coefficients.identifier()
            == "genetic/D/probit_coefficients");
        REQUIRE(dominance.probit_coefficients.statistics().mean.isApprox(
            Eigen::Vector2d{{0.7, -0.2}}));
        REQUIRE(
            result.genetic_parameters().joint().probabilities.identifier()
            == "genetic/joint/probabilities");
        REQUIRE(
            result.genetic_parameters()
                .joint()
                .component_explained_variance.statistics()
                .mean.size()
            == 4);
        REQUIRE(result.marker_effects()
                    .template get<gelex::GeneticMode::A>()
                    .pip()
                    .probabilities()
                    .isApprox(Eigen::VectorXd{{1.0, 1.0}}));
        REQUIRE(result.marker_effects()
                    .template get<gelex::GeneticMode::D>()
                    .pip()
                    .probabilities()
                    .isApprox(Eigen::VectorXd{{0.0, 1.0}}));
        REQUIRE(result.marker_effects().joint().pip().probabilities().isApprox(
            Eigen::VectorXd{{1.0, 1.0}}));
        REQUIRE(result.marker_effects().joint().pve().values().isApprox(
            Eigen::VectorXd::Zero(2)));
    }

    SECTION("fixed joint probabilities")
    {
        using Family
            = gelex::JointSpikeSlabFamily<gelex::MixtureWeightUpdate::Disabled>;
        const auto result = collect_result<mode_ad, Family>(
            fixture.get_test_dir() / "fixed_joint.draws",
            [](auto& state)
            { state.genetic().joint().fitted_values.setZero(); });
        const auto& joint = result.genetic_parameters().joint();

        STATIC_REQUIRE(
            EmptyProbabilities<std::remove_cvref_t<decltype(joint)>>);
        REQUIRE(
            joint.component_explained_variance.identifier()
            == "genetic/joint/component_explained_variance");
    }
}

TEST_CASE(
    "BayesResult validates model and coefficient draw alignment",
    "[bayes][result]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    auto source_model = gelex::test::make_random_effect_model(mode_a);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        source_model);
    auto state = gelex::make_state(prior, source_model);
    auto draws = gelex::make_draws(
        prior,
        source_model,
        (fixture.get_test_dir() / "mismatched.draws").string(),
        1);
    draws.append(state);

    auto no_random_model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        mode_a);
    REQUIRE_THROWS_AS(
        gelex::make_result(no_random_model, draws), gelex::GelexException);

    auto fixed_design = gelex::FixedDesign::make(
        gelex::QuantitativeCovariate{
            .names = {"age"}, .X = Eigen::MatrixXd{{1.0}, {2.0}, {3.0}}});
    auto fixed_mismatch = gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        std::move(fixed_design),
        {},
        gelex::test::make_genetic_design(
            Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, mode_a)};
    REQUIRE_THROWS_AS(
        gelex::make_result(fixed_mismatch, draws), gelex::GelexException);

    std::vector<gelex::bayes::RandomDesign> one_column_random;
    one_column_random.push_back(
        gelex::test::make_random_design(
            "batch",
            std::array{std::string{"batch"}},
            Eigen::MatrixXd{{1.0}, {0.0}, {1.0}}));
    auto random_mismatch = gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        std::move(one_column_random),
        gelex::test::make_genetic_design(
            Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, mode_a)};
    REQUIRE_THROWS_AS(
        gelex::make_result(random_mismatch, draws), gelex::GelexException);
}

TEST_CASE(
    "BayesResult rejects draws without any recorded sample",
    "[bayes][result]")
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    gelex::test::FileFixture fixture;
    auto model = gelex::test::make_random_effect_model(mode_a);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto draws = gelex::make_draws(
        prior, model, (fixture.get_test_dir() / "empty.draws").string(), 1);

    REQUIRE_THROWS_AS(gelex::make_result(model, draws), gelex::GelexException);
}
