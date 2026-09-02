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
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <ranges>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior_compilation.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"

using Catch::Approx;
using gelex::bayes_prior_t;
using gelex::BayesModel;
using gelex::BayesPrior;
using gelex::BayesRecipe;
using gelex::compile;
using gelex::GaussianMethod;
using gelex::GaussianPrior;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::IndependentTopology;
using gelex::JointSpikeSlabMethod;
using gelex::JointSpikeSlabPrior;
using gelex::JointTopology;
using gelex::ScaledMixture;
using gelex::ScaledMixtureMethod;
using gelex::ScaledMixturePrior;
using gelex::SpikeSlab;
using gelex::SpikeSlabMethod;
using gelex::SpikeSlabPrior;
using gelex::UpdatePolicy;
using gelex::VarianceBudget;
using gelex::VarianceLayout;

namespace
{

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using SpikeSlabAD = gelex::IndependentTopology<mode_ad, SpikeSlab>;
using ScaledMixtureAD = gelex::IndependentTopology<mode_ad, ScaledMixture>;
using PooledGaussianMethod = GaussianMethod<VarianceLayout::Pooled>;
using UnpooledGaussianMethod = GaussianMethod<VarianceLayout::Unpooled>;
using PooledSpikeSlabMethod = SpikeSlabMethod<VarianceLayout::Pooled>;
using UnpooledSpikeSlabMethod = SpikeSlabMethod<VarianceLayout::Unpooled>;
using FixedUnpooledSpikeSlabMethod
    = SpikeSlabMethod<VarianceLayout::Unpooled, UpdatePolicy::Fixed>;
using DefaultScaledMixtureMethod = ScaledMixtureMethod<>;
using FixedScaledMixtureMethod = ScaledMixtureMethod<UpdatePolicy::Fixed>;
using DefaultJointSpikeSlabMethod = JointSpikeSlabMethod<>;
using MixedJointSpikeSlabMethod
    = JointSpikeSlabMethod<UpdatePolicy::Fixed, UpdatePolicy::Sampled>;

// Each recipe type admits exactly one prior type, and the five independent
// families differ only in their leaf.
static_assert(std::same_as<
              bayes_prior_t<BayesRecipe<PooledGaussianMethod, mode_ad>>,
              BayesPrior<IndependentTopology<
                  mode_ad,
                  GaussianPrior<VarianceLayout::Pooled>>>>);
static_assert(std::same_as<
              bayes_prior_t<BayesRecipe<UnpooledGaussianMethod, mode_a>>,
              BayesPrior<IndependentTopology<
                  mode_a,
                  GaussianPrior<VarianceLayout::Unpooled>>>>);
static_assert(std::same_as<
              bayes_prior_t<BayesRecipe<UnpooledSpikeSlabMethod, mode_ad>>,
              BayesPrior<IndependentTopology<
                  mode_ad,
                  SpikeSlabPrior<VarianceLayout::Unpooled>>>>);
static_assert(std::same_as<
              bayes_prior_t<BayesRecipe<PooledSpikeSlabMethod, mode_ad>>,
              BayesPrior<IndependentTopology<
                  mode_ad,
                  SpikeSlabPrior<VarianceLayout::Pooled>>>>);
static_assert(std::same_as<
              bayes_prior_t<BayesRecipe<DefaultScaledMixtureMethod, mode_ad>>,
              BayesPrior<IndependentTopology<mode_ad, ScaledMixturePrior<>>>>);
static_assert(std::same_as<
              bayes_prior_t<BayesRecipe<DefaultJointSpikeSlabMethod, mode_ad>>,
              BayesPrior<JointTopology<
                  GaussianPrior<VarianceLayout::Pooled>,
                  JointSpikeSlabPrior<>>>>);
static_assert(
    std::same_as<
        bayes_prior_t<BayesRecipe<FixedUnpooledSpikeSlabMethod, mode_ad>>,
        BayesPrior<IndependentTopology<
            mode_ad,
            SpikeSlabPrior<VarianceLayout::Unpooled, UpdatePolicy::Fixed>>>>);

template <typename Recipe>
concept Compilable = requires(const Recipe& recipe, const BayesModel& model) {
    compile(recipe, model);
};

static_assert(Compilable<BayesRecipe<PooledGaussianMethod, mode_a>>);
static_assert(Compilable<BayesRecipe<DefaultJointSpikeSlabMethod, mode_ad>>);

template <typename T>
concept HasHyperprior = requires(const T& parameter) { parameter.hyperprior; };

auto make_model(GeneticModeSet modes) -> BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        modes);
}

auto make_model_with_random(std::vector<gelex::bayes::RandomDesign> random)
    -> BayesModel
{
    auto genetic = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, mode_a);
    return BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        std::move(random),
        std::move(genetic)};
}

auto projection_variance(const Eigen::MatrixXd& design) -> double
{
    double variance = 0.0;
    for (Eigen::Index column = 0; column < design.cols(); ++column)
    {
        const auto values = design.col(column);
        variance += (values.array() - values.mean()).square().mean();
    }
    return variance;
}

// The calibration target restated from ADR 0002, so the test states the
// contract rather than the arithmetic the implementation happens to use.
auto expected_variance(
    const BayesModel& model,
    GeneticMode mode,
    double share,
    double activity) -> double
{
    return model.phenotype_variance() * share
           / (activity * model.genetic().projection(mode).col_var().sum());
}

}  // namespace

TEST_CASE(
    "compile calibrates each mode against its own share",
    "[bayes][prior_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<PooledGaussianMethod, mode_ad>{
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}}};

    const auto prior = compile(recipe, model);

    REQUIRE(
        prior.genetic().get<GeneticMode::A>().variance.initial_value()
        == Approx(expected_variance(model, GeneticMode::A, 0.4, 1.0)));
    REQUIRE(
        prior.genetic().get<GeneticMode::D>().variance.initial_value()
        == Approx(expected_variance(model, GeneticMode::D, 0.1, 1.0)));
}

TEST_CASE(
    "compile calibrates the variance prior mean",
    "[bayes][prior_compilation]")
{
    const auto model = make_model(mode_a);
    const auto prior = compile(
        BayesRecipe<UnpooledGaussianMethod, mode_a>::defaults(), model);

    const auto& variance = prior.genetic().get<GeneticMode::A>().variance;

    REQUIRE(variance.prior().degrees_of_freedom() == 4.0);
    REQUIRE(variance.prior().scale() == Approx(0.5 * variance.initial_value()));
}

TEST_CASE(
    "compile spreads the share over the included markers only",
    "[bayes][prior_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<PooledSpikeSlabMethod, mode_ad>{
        SpikeSlabAD{
            SpikeSlab{.probability = 0.05},
            SpikeSlab{.probability = 0.2},
        },
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
    };

    const auto prior = compile(recipe, model);

    // A smaller inclusion probability concentrates the same share on fewer
    // markers, so each included marker carries more variance.
    REQUIRE(
        prior.genetic().get<GeneticMode::A>().variance.initial_value()
        == Approx(expected_variance(model, GeneticMode::A, 0.4, 0.05)));
    REQUIRE(
        prior.genetic().get<GeneticMode::D>().variance.initial_value()
        == Approx(expected_variance(model, GeneticMode::D, 0.1, 0.2)));

    REQUIRE(prior.genetic().get<GeneticMode::A>().probability.initial == 0.05);
    REQUIRE(prior.genetic().get<GeneticMode::D>().probability.initial == 0.2);
}

TEST_CASE(
    "compile gives fixed and sampled methods distinct parameter types",
    "[bayes][prior_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto fixed_recipe
        = BayesRecipe<FixedUnpooledSpikeSlabMethod, mode_ad>{
            SpikeSlabAD{
                SpikeSlab{.probability = 0.01},
                SpikeSlab{.probability = 0.02},
            },
            VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
        };
    const auto sampled_recipe = BayesRecipe<UnpooledSpikeSlabMethod, mode_ad>{
        SpikeSlabAD{
            SpikeSlab{.probability = 0.01},
            SpikeSlab{.probability = 0.02},
        },
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
    };

    const auto fixed_prior = compile(fixed_recipe, model);
    const auto sampled_prior = compile(sampled_recipe, model);
    const auto& fixed_probability
        = fixed_prior.genetic().get<GeneticMode::A>().probability;
    const auto& sampled_probability
        = sampled_prior.genetic().get<GeneticMode::A>().probability;

    static_assert(!HasHyperprior<decltype(fixed_probability)>);
    static_assert(HasHyperprior<decltype(sampled_probability)>);
    REQUIRE(fixed_probability.initial == 0.01);
    REQUIRE(sampled_probability.hyperprior.alpha == 1.0);
    REQUIRE(sampled_probability.hyperprior.beta == 1.0);
}

TEST_CASE(
    "compile uses class probabilities to weight activity",
    "[bayes][prior_compilation]")
{
    const auto model = make_model(mode_a);
    const auto recipe = BayesRecipe<DefaultScaledMixtureMethod, mode_a>{
        gelex::IndependentTopology<mode_a, ScaledMixture>{ScaledMixture{}},
        VarianceBudget{{.additive = 0.5}},
    };

    const auto prior = compile(recipe, model);
    const auto fixed_prior = compile(
        BayesRecipe<FixedScaledMixtureMethod, mode_a>{
            gelex::IndependentTopology<mode_a, ScaledMixture>{ScaledMixture{}},
            VarianceBudget{{.additive = 0.5}}},
        model);
    const auto& leaf = prior.genetic().get<GeneticMode::A>();
    const auto& fixed_probabilities
        = fixed_prior.genetic().get<GeneticMode::A>().probabilities;

    // Default probabilities and scales: 0.005 * 0.001 + 0.003 * 0.01
    // + 0.001 * 0.1 + 0.001 * 1, with the null class contributing nothing.
    constexpr double activity = 0.001135;

    REQUIRE(
        leaf.variance.initial_value()
        == Approx(expected_variance(model, GeneticMode::A, 0.5, activity)));
    REQUIRE(leaf.scales == ScaledMixture{}.scales);
    REQUIRE(leaf.probabilities.initial == ScaledMixture{}.probabilities);
    REQUIRE(
        leaf.probabilities.hyperprior.concentration
        == std::array<double, 5>{1.0, 1.0, 1.0, 1.0, 1.0});
    static_assert(!HasHyperprior<decltype(fixed_probabilities)>);
    REQUIRE(fixed_probabilities.initial == ScaledMixture{}.probabilities);
}

TEST_CASE(
    "compile derives joint marginal activity",
    "[bayes][prior_compilation]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<MixedJointSpikeSlabMethod, mode_ad>{
        gelex::JointSpikeSlab{
            .probabilities = {0.8, 0.1, 0.05, 0.05},
            .positive_probability = 0.6,
        },
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
    };

    const auto prior = compile(recipe, model);

    REQUIRE(
        prior.genetic()
            .mode_values()
            .get<GeneticMode::A>()
            .variance.initial_value()
        == Approx(expected_variance(model, GeneticMode::A, 0.4, 0.15)));
    REQUIRE(
        prior.genetic()
            .mode_values()
            .get<GeneticMode::D>()
            .variance.initial_value()
        == Approx(expected_variance(model, GeneticMode::D, 0.1, 0.1)));
    static_assert(
        !HasHyperprior<decltype(prior.genetic().joint().probabilities)>);
    REQUIRE(prior.genetic().joint().positive_probability.initial == 0.6);
    REQUIRE(
        prior.genetic().joint().positive_probability.hyperprior.alpha == 1.0);
}

TEST_CASE(
    "compile owns the residual part of the variance budget",
    "[bayes][prior_compilation]")
{
    const auto model = make_model(mode_a);
    const auto recipe = BayesRecipe<PooledGaussianMethod, mode_a>{
        VarianceBudget{{.additive = 0.4}}};

    const auto prior = compile(recipe, model);

    REQUIRE(prior.random().empty());
    REQUIRE(
        prior.residual().initial_value()
        == Approx(model.phenotype_variance() * 0.6));
    REQUIRE(prior.residual().prior().degrees_of_freedom() == 4.0);
    REQUIRE(
        prior.residual().prior().scale()
        == Approx(0.5 * prior.residual().initial_value()));
}

TEST_CASE(
    "compile splits the total random share equally across random blocks",
    "[bayes][prior_compilation]")
{
    auto model = make_model_with_random(
        {gelex::bayes::RandomDesign{
             "first", {"first"}, Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}},
         gelex::bayes::RandomDesign{
             "second", {"second"}, Eigen::MatrixXd{{0.0}, {2.0}, {4.0}}}});
    const auto recipe = BayesRecipe<PooledGaussianMethod, mode_a>{
        VarianceBudget{{.additive = 0.4, .random = 0.2}}};

    const auto prior = compile(recipe, model);
    const double block_target = model.phenotype_variance() * 0.2 / 2.0;

    REQUIRE(prior.random().size() == model.random().size());
    for (const auto& [parameter, design] :
         std::views::zip(prior.random(), model.random()))
    {
        REQUIRE(
            parameter.initial_value() * projection_variance(design.X)
            == Approx(block_target));
        REQUIRE(parameter.prior().degrees_of_freedom() == 4.0);
        REQUIRE(
            parameter.prior().scale()
            == Approx(0.5 * parameter.initial_value()));
    }
}

TEST_CASE(
    "compile requires the random share and random designs to agree",
    "[bayes][prior_compilation]")
{
    SECTION("a share without a design")
    {
        const auto model = make_model(mode_a);
        const auto recipe = BayesRecipe<PooledGaussianMethod, mode_a>{
            VarianceBudget{{.additive = 0.4, .random = 0.1}}};

        REQUIRE_THROWS_AS(compile(recipe, model), gelex::GelexException);
    }

    SECTION("a design without a share")
    {
        const auto model = make_model_with_random({gelex::bayes::RandomDesign{
            "random", {"random"}, Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}}});

        REQUIRE_THROWS_AS(
            compile(
                BayesRecipe<PooledGaussianMethod, mode_a>::defaults(), model),
            gelex::GelexException);
    }
}

TEST_CASE(
    "compile rejects a random block with zero projection variance",
    "[bayes][prior_compilation]")
{
    const auto model = make_model_with_random({gelex::bayes::RandomDesign{
        "constant", {"constant"}, Eigen::MatrixXd{{1.0}, {1.0}, {1.0}}}});
    const auto recipe = BayesRecipe<PooledGaussianMethod, mode_a>{
        VarianceBudget{{.additive = 0.4, .random = 0.1}}};

    REQUIRE_THROWS_AS(compile(recipe, model), gelex::GelexException);
}
