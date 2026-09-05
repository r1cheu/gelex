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
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/kernel.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/data/fixed_design.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "random_design_fixture.h"

using Catch::Approx;

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;
using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;
using UnpooledFamily = gelex::GaussianFamily<gelex::VarianceLayout::Unpooled>;
using PooledSpikeSlabFamily
    = gelex::SpikeSlabFamily<gelex::VarianceLayout::Pooled>;
using UnpooledSpikeSlabFamily
    = gelex::SpikeSlabFamily<gelex::VarianceLayout::Unpooled>;
using FixedUnpooledSpikeSlabFamily = gelex::SpikeSlabFamily<
    gelex::VarianceLayout::Unpooled,
    gelex::MixtureWeightUpdate::Disabled>;
using SampledScaledMixtureFamily = gelex::ScaledMixtureFamily<>;
using FixedScaledMixtureFamily
    = gelex::ScaledMixtureFamily<gelex::MixtureWeightUpdate::Disabled>;
using AdditiveGeneticPrior = gelex::
    ModeValues<mode_a, gelex::GaussianPrior<gelex::VarianceLayout::Pooled>>;
using AdditiveDominanceGeneticPrior = gelex::ModeValues<
    mode_ad,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>>;
using UnpooledAdditiveGeneticPrior = gelex::
    ModeValues<mode_a, gelex::GaussianPrior<gelex::VarianceLayout::Unpooled>>;
using PooledSpikeSlabGeneticPrior = gelex::
    ModeValues<mode_a, gelex::SpikeSlabPrior<gelex::VarianceLayout::Pooled>>;
using PooledSpikeSlabADGeneticPrior = gelex::ModeValues<
    mode_ad,
    gelex::SpikeSlabPrior<gelex::VarianceLayout::Pooled>,
    gelex::SpikeSlabPrior<gelex::VarianceLayout::Pooled>>;
using UnpooledSpikeSlabGeneticPrior = gelex::
    ModeValues<mode_a, gelex::SpikeSlabPrior<gelex::VarianceLayout::Unpooled>>;
using SampledScaledMixtureGeneticPrior = gelex::ModeValues<
    mode_a,
    gelex::ScaledMixturePrior<gelex::ScaledMixture::class_count>>;
using FixedScaledMixtureGeneticPrior = gelex::ModeValues<
    mode_a,
    gelex::ScaledMixturePrior<
        gelex::ScaledMixture::class_count,
        gelex::MixtureWeightUpdate::Disabled>>;
using HeterogeneousGeneticPrior = gelex::ModeValues<
    mode_ad,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::SpikeSlabPrior<
        gelex::VarianceLayout::Unpooled,
        gelex::MixtureWeightUpdate::Disabled>>;

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<const gelex::BayesPrior<AdditiveGeneticPrior>&>())),
        gelex::BayesKernel<AdditiveGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<AdditiveDominanceGeneticPrior>&>())),
        gelex::BayesKernel<AdditiveDominanceGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<UnpooledAdditiveGeneticPrior>&>())),
        gelex::BayesKernel<UnpooledAdditiveGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<PooledSpikeSlabGeneticPrior>&>())),
        gelex::BayesKernel<PooledSpikeSlabGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<PooledSpikeSlabADGeneticPrior>&>())),
        gelex::BayesKernel<PooledSpikeSlabADGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<UnpooledSpikeSlabGeneticPrior>&>())),
        gelex::BayesKernel<UnpooledSpikeSlabGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<SampledScaledMixtureGeneticPrior>&>())),
        gelex::BayesKernel<SampledScaledMixtureGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<FixedScaledMixtureGeneticPrior>&>())),
        gelex::BayesKernel<FixedScaledMixtureGeneticPrior>>);

static_assert(std::same_as<
              decltype(gelex::make_kernel(
                  std::declval<
                      const gelex::BayesPrior<HeterogeneousGeneticPrior>&>())),
              gelex::BayesKernel<HeterogeneousGeneticPrior>>);
static_assert(
    std::is_destructible_v<gelex::BayesKernel<HeterogeneousGeneticPrior>>);

auto make_model() -> gelex::BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        mode_a);
}

auto make_ad_model() -> gelex::BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        mode_ad);
}

auto make_model_with_random() -> gelex::BayesModel
{
    auto genetic = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        mode_a);
    std::vector<gelex::bayes::RandomDesign> random;
    random.push_back(
        gelex::test::make_random_design(
            "batch",
            std::vector<std::string>{"batch"},
            Eigen::MatrixXd{{0.0}, {1.0}, {0.0}, {1.0}}));
    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        gelex::FixedDesign::make(4),
        std::move(random),
        std::move(genetic)};
}

auto reconstruct_genetic_fitted(
    const gelex::bayes::GeneticDesign& design,
    gelex::GeneticMode mode,
    const Eigen::VectorXd& coefficients) -> Eigen::VectorXd
{
    Eigen::VectorXd fitted = Eigen::VectorXd::Zero(design.rows());
    const auto& projection = design.projection(mode);
    for (Eigen::Index marker = 0; marker < coefficients.size(); ++marker)
    {
        projection.axpy(marker, coefficients(marker), fitted);
    }
    return fitted;
}

auto reconstruct_scaled_mixture_fitted(
    const gelex::bayes::GeneticDesign& design,
    gelex::GeneticMode mode,
    const Eigen::VectorXd& coefficients,
    const Eigen::VectorX<std::uint8_t>& assignment) -> Eigen::MatrixXd
{
    Eigen::MatrixXd fitted = Eigen::MatrixXd::Zero(
        design.rows(),
        static_cast<Eigen::Index>(gelex::ScaledMixture::class_count - 1));
    const auto& projection = design.projection(mode);
    for (Eigen::Index marker = 0; marker < coefficients.size(); ++marker)
    {
        const auto class_index = static_cast<std::size_t>(assignment(marker));
        if (class_index != 0)
        {
            projection.axpy(
                marker,
                coefficients(marker),
                fitted.col(static_cast<Eigen::Index>(class_index - 1)));
        }
    }
    return fitted;
}

}  // namespace

TEST_CASE(
    "pooled Gaussian kernel preserves adjusted-response and fitted-cache "
    "invariants",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE(genetic.family_state.fitted_values.isApprox(reconstructed));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(genetic.coefficients(1) == 0.0);
    REQUIRE(state.residual().variance > 0.0);
    REQUIRE(genetic.family_state.variance > 0.0);
}

TEST_CASE(
    "independent AD pooled Gaussian kernel preserves mode caches and the "
    "adjusted-response invariant",
    "[bayes][kernel]")
{
    const auto model = make_ad_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, Family>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& additive = state.genetic().get<gelex::GeneticMode::A>();
    const auto& dominance = state.genetic().get<gelex::GeneticMode::D>();
    const auto additive_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients);
    const auto dominance_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE(additive.family_state.fitted_values.isApprox(additive_fitted));
    REQUIRE(dominance.family_state.fitted_values.isApprox(dominance_fitted));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + additive_fitted
             + dominance_fitted)
                .isApprox(model.phenotype()));
    REQUIRE(additive.coefficients(1) == 0.0);
    REQUIRE(dominance.coefficients(1) == 0.0);
    REQUIRE(additive.family_state.variance > 0.0);
    REQUIRE(dominance.family_state.variance > 0.0);
}

TEST_CASE(
    "unpooled Gaussian kernel preserves adjusted-response and fitted-cache "
    "invariants",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, UnpooledFamily>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    const double invalid_marker_variance
        = state.genetic().get<gelex::GeneticMode::A>().family_state.variance(1);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE(genetic.family_state.fitted_values.isApprox(reconstructed));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(genetic.coefficients(1) == 0.0);
    REQUIRE(genetic.family_state.variance.size() == model.genetic().cols());
    REQUIRE((genetic.family_state.variance.array() > 0.0).all());
    REQUIRE(genetic.family_state.variance(1) == invalid_marker_variance);
}

TEST_CASE(
    "pooled Gaussian kernel preserves random-effect and residual invariants",
    "[bayes][kernel]")
{
    const auto model = make_model_with_random();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    // The fitted cache is rebuilt each sweep, so a missing reset only shows up
    // from the second one onwards.
    for (int sweep = 0; sweep < 3; ++sweep)
    {
        kernel.step(model, state, rng);
    }

    const auto& random = state.random().front();
    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;
    const Eigen::VectorXd random_fitted
        = model.random().front().X() * random.coefficients;
    const auto genetic_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);

    REQUIRE(random.coefficients.allFinite());
    REQUIRE(random.variance > 0.0);
    REQUIRE(random.fitted_values.isApprox(random_fitted));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + random_fitted
             + genetic_fitted)
                .isApprox(model.phenotype()));
}

TEST_CASE(
    "pooled spike-slab kernel preserves collapsed-state invariants",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, PooledSpikeSlabFamily>{
            gelex::ModeValues<mode_a, gelex::SpikeSlab>{gelex::SpikeSlab{0.5}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE(genetic.family_state.fitted_values.isApprox(reconstructed));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(family.assignment.size() == model.genetic().cols());
    for (Eigen::Index marker = 0; marker < family.assignment.size(); ++marker)
    {
        REQUIRE(family.assignment(marker) <= 1);
        if (family.assignment(marker) == 0)
        {
            REQUIRE(genetic.coefficients(marker) == 0.0);
        }
    }
    REQUIRE(family.assignment(1) == 0);
    REQUIRE(genetic.coefficients(1) == 0.0);
    REQUIRE(family.variance > 0.0);
    REQUIRE(family.probability > 0.0);
    REQUIRE(family.probability < 1.0);
    REQUIRE(family.probability != 0.5);
}

TEST_CASE(
    "independent AD pooled spike-slab kernel preserves mode caches",
    "[bayes][kernel]")
{
    const auto model = make_ad_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, PooledSpikeSlabFamily>{
            gelex::ModeValues<mode_ad, gelex::SpikeSlab, gelex::SpikeSlab>{
                gelex::SpikeSlab{0.5}, gelex::SpikeSlab{0.5}},
            gelex::VarianceBudget{{.additive = 0.4, .dominance = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& additive = state.genetic().get<gelex::GeneticMode::A>();
    const auto& dominance = state.genetic().get<gelex::GeneticMode::D>();
    const auto additive_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients);
    const auto dominance_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE(additive.family_state.fitted_values.isApprox(additive_fitted));
    REQUIRE(dominance.family_state.fitted_values.isApprox(dominance_fitted));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + additive_fitted
             + dominance_fitted)
                .isApprox(model.phenotype()));
    for (Eigen::Index marker = 0; marker < model.genetic().cols(); ++marker)
    {
        if (additive.family_state.assignment(marker) == 0)
        {
            REQUIRE(additive.coefficients(marker) == 0.0);
        }
        if (dominance.family_state.assignment(marker) == 0)
        {
            REQUIRE(dominance.coefficients(marker) == 0.0);
        }
    }
}

TEST_CASE(
    "unpooled spike-slab kernel samples inactive variances and keeps a fixed "
    "probability",
    "[bayes][kernel]")
{
    const auto model = make_model();
    constexpr double fixed_probability = 1.0e-12;
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, FixedUnpooledSpikeSlabFamily>{
            gelex::ModeValues<mode_a, gelex::SpikeSlab>{
                gelex::SpikeSlab{fixed_probability}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    auto& initial_family
        = state.genetic().get<gelex::GeneticMode::A>().family_state;
    const Eigen::VectorXd initial_variance = initial_family.variance;
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE(genetic.coefficients.isZero());
    REQUIRE(family.assignment.isZero());
    REQUIRE(genetic.family_state.fitted_values.isZero());
    REQUIRE((state.residual().adjusted_response + fixed_fitted)
                .isApprox(model.phenotype()));
    REQUIRE((family.variance.array() > 0.0).all());
    REQUIRE(family.variance(0) != initial_variance(0));
    REQUIRE(family.variance(1) == initial_variance(1));
    REQUIRE(family.probability == fixed_probability);
}

TEST_CASE(
    "scaled-mixture kernel samples classes and probabilities",
    "[bayes][kernel]")
{
    const auto model = make_model();
    constexpr std::array probabilities{0.05, 0.1, 0.15, 0.2, 0.5};
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, SampledScaledMixtureFamily>{
            gelex::ModeValues<mode_a, gelex::ScaledMixture>{
                gelex::ScaledMixture{probabilities}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const auto fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const auto component_fitted = reconstruct_scaled_mixture_fitted(
        model.genetic(),
        gelex::GeneticMode::A,
        genetic.coefficients,
        family.assignment);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;
    REQUIRE((state.residual().adjusted_response + fixed_fitted + fitted)
                .isApprox(model.phenotype()));
    REQUIRE(family.fitted_values.isApprox(component_fitted));
    REQUIRE(family.assignment(1) == 0);
    REQUIRE(genetic.coefficients(1) == 0.0);
    for (Eigen::Index marker = 0; marker < family.assignment.size(); ++marker)
    {
        REQUIRE(family.assignment(marker) < gelex::ScaledMixture::class_count);
        if (family.assignment(marker) == 0)
        {
            REQUIRE(genetic.coefficients(marker) == 0.0);
        }
    }
    double probability_sum = 0.0;
    for (const double probability : family.probabilities)
    {
        REQUIRE(probability > 0.0);
        probability_sum += probability;
    }
    REQUIRE(probability_sum == Approx(1.0));
    REQUIRE(family.probabilities != probabilities);
    REQUIRE(family.variance > 0.0);
}

TEST_CASE(
    "fixed scaled-mixture kernel omits probability updates",
    "[bayes][kernel]")
{
    const auto model = make_model();
    constexpr std::array probabilities{
        1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0 - 4.0e-12};
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, FixedScaledMixtureFamily>{
            gelex::ModeValues<mode_a, gelex::ScaledMixture>{
                gelex::ScaledMixture{probabilities}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const auto component_fitted = reconstruct_scaled_mixture_fitted(
        model.genetic(),
        gelex::GeneticMode::A,
        genetic.coefficients,
        family.assignment);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE(family.probabilities == probabilities);
    REQUIRE(family.fitted_values.isApprox(component_fitted));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(family.variance > 0.0);
}
