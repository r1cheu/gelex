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

#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/data/fixed_design.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "random_design_fixture.h"

using Catch::Approx;
using gelex::BayesModel;
using gelex::BayesPrior;
using gelex::BayesRecipe;
using gelex::Gaussian;
using gelex::GaussianFamily;
using gelex::GaussianPrior;
using gelex::GeneticMode;
using gelex::GeneticModeSet;
using gelex::HalfNormal;
using gelex::HalfNormalPrior;
using gelex::JointModeValues;
using gelex::JointSpikeSlab;
using gelex::JointSpikeSlabFamily;
using gelex::JointSpikeSlabPrior;
using gelex::make_prior;
using gelex::MixtureWeightUpdate;
using gelex::ModeValues;
using gelex::ScaledMixture;
using gelex::ScaledMixtureFamily;
using gelex::ScaledMixturePrior;
using gelex::SpikeSlab;
using gelex::SpikeSlabFamily;
using gelex::SpikeSlabPrior;
using gelex::VarianceBudget;
using gelex::VarianceLayout;

namespace
{

template <typename Recipe>
using prior_result_t = decltype(make_prior(
    std::declval<const Recipe&>(),
    std::declval<const BayesModel&>()));

constexpr auto mode_a = GeneticModeSet{GeneticMode::A};
constexpr auto mode_ad = GeneticMode::A | GeneticMode::D;

using SpikeSlabAD = gelex::ModeValues<mode_ad, SpikeSlab, SpikeSlab>;
using ScaledMixtureAD
    = gelex::ModeValues<mode_ad, ScaledMixture, ScaledMixture>;
using JointModeSpecs = ModeValues<mode_ad, Gaussian, HalfNormal>;
using JointSpikeSlabAD = JointModeValues<JointModeSpecs, JointSpikeSlab>;
using PooledGaussianFamily = GaussianFamily<VarianceLayout::Pooled>;
using UnpooledGaussianFamily = GaussianFamily<VarianceLayout::Unpooled>;
using PooledSpikeSlabFamily = SpikeSlabFamily<VarianceLayout::Pooled>;
using UnpooledSpikeSlabFamily = SpikeSlabFamily<VarianceLayout::Unpooled>;
using FixedUnpooledSpikeSlabFamily
    = SpikeSlabFamily<VarianceLayout::Unpooled, MixtureWeightUpdate::Disabled>;
using DefaultScaledMixtureFamily = ScaledMixtureFamily<>;
using FixedScaledMixtureFamily
    = ScaledMixtureFamily<MixtureWeightUpdate::Disabled>;
using DefaultJointSpikeSlabFamily = JointSpikeSlabFamily<>;
using FixedJointSpikeSlabFamily
    = JointSpikeSlabFamily<MixtureWeightUpdate::Disabled>;

// Each recipe type admits exactly one prior type, and the five independent
// families differ only in their leaf.
static_assert(std::same_as<
              prior_result_t<BayesRecipe<mode_ad, PooledGaussianFamily>>,
              BayesPrior<ModeValues<
                  mode_ad,
                  GaussianPrior<VarianceLayout::Pooled>,
                  GaussianPrior<VarianceLayout::Pooled>>>>);
static_assert(
    std::same_as<
        prior_result_t<BayesRecipe<mode_a, UnpooledGaussianFamily>>,
        BayesPrior<
            ModeValues<mode_a, GaussianPrior<VarianceLayout::Unpooled>>>>);
static_assert(std::same_as<
              prior_result_t<BayesRecipe<mode_ad, UnpooledSpikeSlabFamily>>,
              BayesPrior<ModeValues<
                  mode_ad,
                  SpikeSlabPrior<VarianceLayout::Unpooled>,
                  SpikeSlabPrior<VarianceLayout::Unpooled>>>>);
static_assert(std::same_as<
              prior_result_t<BayesRecipe<mode_ad, PooledSpikeSlabFamily>>,
              BayesPrior<ModeValues<
                  mode_ad,
                  SpikeSlabPrior<VarianceLayout::Pooled>,
                  SpikeSlabPrior<VarianceLayout::Pooled>>>>);
static_assert(std::same_as<
              prior_result_t<BayesRecipe<mode_ad, DefaultScaledMixtureFamily>>,
              BayesPrior<ModeValues<
                  mode_ad,
                  ScaledMixturePrior<ScaledMixture::class_count>,
                  ScaledMixturePrior<ScaledMixture::class_count>>>>);
static_assert(std::same_as<
              prior_result_t<BayesRecipe<mode_ad, DefaultJointSpikeSlabFamily>>,
              BayesPrior<JointModeValues<
                  ModeValues<
                      mode_ad,
                      GaussianPrior<VarianceLayout::Pooled>,
                      HalfNormalPrior>,
                  JointSpikeSlabPrior<JointSpikeSlab::class_count>>>>);
static_assert(
    std::same_as<
        prior_result_t<BayesRecipe<mode_ad, FixedUnpooledSpikeSlabFamily>>,
        BayesPrior<ModeValues<
            mode_ad,
            SpikeSlabPrior<
                VarianceLayout::Unpooled,
                MixtureWeightUpdate::Disabled>,
            SpikeSlabPrior<
                VarianceLayout::Unpooled,
                MixtureWeightUpdate::Disabled>>>>);

template <typename Recipe>
concept CanMakePrior = requires(const Recipe& recipe, const BayesModel& model) {
    make_prior(recipe, model);
};

static_assert(CanMakePrior<BayesRecipe<mode_a, PooledGaussianFamily>>);
static_assert(CanMakePrior<BayesRecipe<mode_ad, DefaultJointSpikeSlabFamily>>);

template <typename T>
concept HasPrior = requires(const T& parameter) { parameter.prior; };

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
    "make_prior calibrates each mode against its own share",
    "[bayes][prior]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<mode_ad, PooledGaussianFamily>{
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}}};

    const auto prior = make_prior(recipe, model);

    REQUIRE(
        prior.genetic().get<GeneticMode::A>().variance.initial
        == Approx(expected_variance(model, GeneticMode::A, 0.4, 1.0)));
    REQUIRE(
        prior.genetic().get<GeneticMode::D>().variance.initial
        == Approx(expected_variance(model, GeneticMode::D, 0.1, 1.0)));
}

TEST_CASE("make_prior calibrates the variance prior mean", "[bayes][prior]")
{
    const auto model = make_model(mode_a);
    const auto prior = make_prior(
        BayesRecipe<mode_a, UnpooledGaussianFamily>::defaults(), model);

    const auto& variance = prior.genetic().get<GeneticMode::A>().variance;

    const auto parameters = variance.prior.scaled_inv_chi2_parameters();
    REQUIRE(parameters.degrees_of_freedom() == 4.0);
    REQUIRE(parameters.scale() == Approx(0.5 * variance.initial));
}

TEST_CASE(
    "make_prior spreads the share over the included markers only",
    "[bayes][prior]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<mode_ad, PooledSpikeSlabFamily>{
        SpikeSlabAD{
            SpikeSlab{0.05},
            SpikeSlab{0.2},
        },
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
    };

    const auto prior = make_prior(recipe, model);

    // A smaller inclusion probability concentrates the same share on fewer
    // markers, so each included marker carries more variance.
    REQUIRE(
        prior.genetic().get<GeneticMode::A>().variance.initial
        == Approx(expected_variance(model, GeneticMode::A, 0.4, 0.05)));
    REQUIRE(
        prior.genetic().get<GeneticMode::D>().variance.initial
        == Approx(expected_variance(model, GeneticMode::D, 0.1, 0.2)));

    REQUIRE(prior.genetic().get<GeneticMode::A>().probability.initial == 0.05);
    REQUIRE(prior.genetic().get<GeneticMode::D>().probability.initial == 0.2);
}

TEST_CASE(
    "make_prior gives fixed and sampled families distinct parameter types",
    "[bayes][prior]")
{
    const auto model = make_model(mode_ad);
    const auto fixed_recipe
        = BayesRecipe<mode_ad, FixedUnpooledSpikeSlabFamily>{
            SpikeSlabAD{
                SpikeSlab{0.01},
                SpikeSlab{0.02},
            },
            VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
        };
    const auto sampled_recipe = BayesRecipe<mode_ad, UnpooledSpikeSlabFamily>{
        SpikeSlabAD{
            SpikeSlab{0.01},
            SpikeSlab{0.02},
        },
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
    };

    const auto fixed_prior = make_prior(fixed_recipe, model);
    const auto sampled_prior = make_prior(sampled_recipe, model);
    const auto& fixed_probability
        = fixed_prior.genetic().get<GeneticMode::A>().probability;
    const auto& sampled_probability
        = sampled_prior.genetic().get<GeneticMode::A>().probability;

    static_assert(!HasPrior<decltype(fixed_probability)>);
    static_assert(HasPrior<decltype(sampled_probability)>);
    REQUIRE(fixed_probability.initial == 0.01);
    REQUIRE(
        sampled_probability.prior.dirichlet_parameters().concentrations()
        == std::array<double, 2>{1.0, 1.0});
}

TEST_CASE(
    "make_prior uses class probabilities to weight activity",
    "[bayes][prior]")
{
    const auto model = make_model(mode_a);
    const auto recipe = BayesRecipe<mode_a, DefaultScaledMixtureFamily>{
        gelex::ModeValues<mode_a, ScaledMixture>{ScaledMixture{}},
        VarianceBudget{{.additive = 0.5}},
    };

    const auto prior = make_prior(recipe, model);
    const auto fixed_prior = make_prior(
        BayesRecipe<mode_a, FixedScaledMixtureFamily>{
            gelex::ModeValues<mode_a, ScaledMixture>{ScaledMixture{}},
            VarianceBudget{{.additive = 0.5}}},
        model);
    const auto& leaf = prior.genetic().get<GeneticMode::A>();
    const auto& fixed_probabilities
        = fixed_prior.genetic().get<GeneticMode::A>().probabilities;
    const auto defaults = ScaledMixture{};

    // Default probabilities and scales: 0.005 * 0.001 + 0.003 * 0.01
    // + 0.001 * 0.1 + 0.001 * 1, with the null class contributing nothing.
    constexpr double activity = 0.001135;

    REQUIRE(
        leaf.variance.initial
        == Approx(expected_variance(model, GeneticMode::A, 0.5, activity)));
    REQUIRE(leaf.scales == defaults.scales());
    REQUIRE(leaf.probabilities.initial == defaults.probabilities());
    REQUIRE(
        leaf.probabilities.prior.dirichlet_parameters().concentrations()
        == std::array<double, 5>{1.0, 1.0, 1.0, 1.0, 1.0});
    static_assert(!HasPrior<decltype(fixed_probabilities)>);
    REQUIRE(fixed_probabilities.initial == defaults.probabilities());
}

TEST_CASE("make_prior derives joint marginal activity", "[bayes][prior]")
{
    const auto model = make_model(mode_ad);
    const auto recipe = BayesRecipe<mode_ad, FixedJointSpikeSlabFamily>{
        JointSpikeSlabAD{
            JointModeSpecs{Gaussian{}, HalfNormal{}},
            JointSpikeSlab{{0.8, 0.1, 0.05, 0.05}}},
        VarianceBudget{{.additive = 0.4, .dominance = 0.1}},
    };

    const auto prior = make_prior(recipe, model);

    REQUIRE(
        prior.genetic().mode_values().get<GeneticMode::A>().variance.initial
        == Approx(expected_variance(model, GeneticMode::A, 0.4, 0.15)));
    REQUIRE(
        prior.genetic().mode_values().get<GeneticMode::D>().variance.initial
        == Approx(expected_variance(model, GeneticMode::D, 0.1, 0.1)));
    static_assert(!HasPrior<decltype(prior.genetic().joint().probabilities)>);
}

TEST_CASE(
    "make_prior owns the residual part of the variance budget",
    "[bayes][prior]")
{
    const auto model = make_model(mode_a);
    const auto recipe = BayesRecipe<mode_a, PooledGaussianFamily>{
        VarianceBudget{{.additive = 0.4}}};

    const auto prior = make_prior(recipe, model);

    REQUIRE(prior.random().empty());
    REQUIRE(
        prior.residual().initial == Approx(model.phenotype_variance() * 0.6));
    const auto residual_parameters
        = prior.residual().prior.scaled_inv_chi2_parameters();
    REQUIRE(residual_parameters.degrees_of_freedom() == 4.0);
    REQUIRE(
        residual_parameters.scale() == Approx(0.5 * prior.residual().initial));
}

TEST_CASE(
    "make_prior splits the total random share equally across random blocks",
    "[bayes][prior]")
{
    auto model = make_model_with_random(
        {gelex::test::make_random_design(
             "first",
             std::vector<std::string>{"first"},
             Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}),
         gelex::test::make_random_design(
             "second",
             std::vector<std::string>{"second"},
             Eigen::MatrixXd{{0.0}, {2.0}, {4.0}})});
    const auto recipe = BayesRecipe<mode_a, PooledGaussianFamily>{
        VarianceBudget{{.additive = 0.4, .random = 0.2}}};

    const auto prior = make_prior(recipe, model);
    const double block_target = model.phenotype_variance() * 0.2 / 2.0;

    REQUIRE(prior.random().size() == model.random().size());
    for (const auto& [parameter, design] :
         std::views::zip(prior.random(), model.random()))
    {
        REQUIRE(
            parameter.initial * projection_variance(design.X())
            == Approx(block_target));
        const auto parameters = parameter.prior.scaled_inv_chi2_parameters();
        REQUIRE(parameters.degrees_of_freedom() == 4.0);
        REQUIRE(parameters.scale() == Approx(0.5 * parameter.initial));
    }
}

TEST_CASE(
    "make_prior requires the random share and random designs to agree",
    "[bayes][prior]")
{
    SECTION("a share without a design")
    {
        const auto model = make_model(mode_a);
        const auto recipe = BayesRecipe<mode_a, PooledGaussianFamily>{
            VarianceBudget{{.additive = 0.4, .random = 0.1}}};

        REQUIRE_THROWS_AS(make_prior(recipe, model), gelex::GelexException);
    }

    SECTION("a design without a share")
    {
        const auto model
            = make_model_with_random({gelex::test::make_random_design(
                "random",
                std::vector<std::string>{"random"},
                Eigen::MatrixXd{{0.0}, {1.0}, {2.0}})});

        REQUIRE_THROWS_AS(
            make_prior(
                BayesRecipe<mode_a, PooledGaussianFamily>::defaults(), model),
            gelex::GelexException);
    }
}

TEST_CASE(
    "make_prior rejects a random block with zero projection variance",
    "[bayes][prior]")
{
    const auto model = make_model_with_random({gelex::test::make_random_design(
        "constant",
        std::vector<std::string>{"constant"},
        Eigen::MatrixXd{{1.0}, {1.0}, {1.0}})});
    const auto recipe = BayesRecipe<mode_a, PooledGaussianFamily>{
        VarianceBudget{{.additive = 0.4, .random = 0.1}}};

    REQUIRE_THROWS_AS(make_prior(recipe, model), gelex::GelexException);
}
